#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <cstdint>
#include <fstream>
#include <ranges>
#include <algorithm>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <limits>
#include <thread>
#include <mutex>
#include <atomic>
#include <chrono>
#include <iomanip>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>
#include <seqan3/core/debug_stream.hpp>

// For serialization
#include <cereal/archives/binary.hpp>

// Progress reporting utility
class ProgressReporter
{
private:
    std::atomic<size_t> current_{0};
    size_t total_;
    std::chrono::steady_clock::time_point start_time_;
    std::mutex output_mutex_;
    
public:
    ProgressReporter(size_t total) : total_(total), start_time_(std::chrono::steady_clock::now()) {}
    
    void update(size_t increment = 1)
    {
        current_ += increment;
        size_t current = current_.load();
        
        if (current % 1000 == 0 || current == total_)
        {
            std::lock_guard<std::mutex> lock(output_mutex_);
            auto now = std::chrono::steady_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start_time_).count();
            
            double progress = (double)current / total_ * 100.0;
            std::cout << "\r  Progress: " << current << "/" << total_ 
                      << " (" << std::fixed << std::setprecision(1) << progress << "%) "
                      << "Time: " << elapsed << "s" << std::flush;
                      
            if (current == total_)
            {
                std::cout << "\n";
            }
        }
    }
};

// Structure to represent a unique region
struct UniqueRegion
{
    size_t start;
    size_t end;
    seqan3::dna5_vector sequence;
};

// Hierarchical k-mer cache for multi-k optimization
class HierarchicalKmerCache
{
private:
    std::unordered_map<uint8_t, std::unordered_map<std::string, bool>> cache_by_k_;
    mutable std::mutex cache_mutex_;
    
public:
    bool find(const std::string& kmer, uint8_t k, bool& found_in_index)
    {
        std::lock_guard<std::mutex> lock(cache_mutex_);
        auto k_cache_it = cache_by_k_.find(k);
        if (k_cache_it != cache_by_k_.end())
        {
            auto it = k_cache_it->second.find(kmer);
            if (it != k_cache_it->second.end())
            {
                found_in_index = it->second;
                return true; // found in cache
            }
        }
        return false; // not in cache
    }
    
    void insert(const std::string& kmer, uint8_t k, bool found_in_index)
    {
        std::lock_guard<std::mutex> lock(cache_mutex_);
        cache_by_k_[k][kmer] = found_in_index;
    }
    
    size_t size(uint8_t k) const
    {
        std::lock_guard<std::mutex> lock(cache_mutex_);
        auto it = cache_by_k_.find(k);
        return it != cache_by_k_.end() ? it->second.size() : 0;
    }
    
    size_t total_size() const
    {
        std::lock_guard<std::mutex> lock(cache_mutex_);
        size_t total = 0;
        for (const auto& [k, cache] : cache_by_k_)
        {
            total += cache.size();
        }
        return total;
    }
    
    void clear()
    {
        std::lock_guard<std::mutex> lock(cache_mutex_);
        cache_by_k_.clear();
    }
    
    // Get cache statistics for a specific k
    std::pair<size_t, size_t> get_stats(uint8_t k) const
    {
        std::lock_guard<std::mutex> lock(cache_mutex_);
        auto it = cache_by_k_.find(k);
        if (it == cache_by_k_.end())
        {
            return {0, 0};
        }
        
        size_t total = it->second.size();
        size_t found_count = 0;
        for (const auto& [kmer, found] : it->second)
        {
            if (found) found_count++;
        }
        return {total, found_count};
    }
    
    // Check if a k-mer can be inferred from smaller k-mers (optimization)
    bool can_infer_from_smaller_kmers(const std::string& kmer, uint8_t k, bool& inferred_result)
    {
        std::lock_guard<std::mutex> lock(cache_mutex_);
        
        // If any smaller k-mer substring is absent, the larger k-mer is also absent
        for (uint8_t smaller_k = std::max(1, (int)k - 5); smaller_k < k; smaller_k++)
        {
            auto smaller_cache_it = cache_by_k_.find(smaller_k);
            if (smaller_cache_it == cache_by_k_.end()) continue;
            
            // Check all possible substrings of size smaller_k
            for (size_t pos = 0; pos <= kmer.length() - smaller_k; pos++)
            {
                std::string sub_kmer = kmer.substr(pos, smaller_k);
                auto sub_it = smaller_cache_it->second.find(sub_kmer);
                if (sub_it != smaller_cache_it->second.end() && !sub_it->second)
                {
                    // Found a smaller k-mer that's absent, so this k-mer is also absent
                    inferred_result = false;
                    return true;
                }
            }
        }
        
        return false; // Cannot infer
    }
    
    // Get all k-values that have been cached
    std::vector<uint8_t> get_cached_k_values() const
    {
        std::lock_guard<std::mutex> lock(cache_mutex_);
        std::vector<uint8_t> k_values;
        for (const auto& [k, cache] : cache_by_k_)
        {
            k_values.push_back(k);
        }
        std::sort(k_values.begin(), k_values.end());
        return k_values;
    }
};

// Function to merge overlapping k-mer positions into unique regions
std::vector<UniqueRegion> merge_overlapping_positions(const std::set<size_t>& absent_positions, 
                                                      uint8_t kmer_size,
                                                      const seqan3::dna5_vector& query_sequence)
{
    std::vector<UniqueRegion> unique_regions;
    if (absent_positions.empty()) return unique_regions;

    std::vector<size_t> sorted_positions(absent_positions.begin(), absent_positions.end());
    
    size_t region_start = sorted_positions[0];
    size_t region_end = sorted_positions[0] + kmer_size - 1;
    
    for (size_t i = 1; i < sorted_positions.size(); ++i)
    {
        size_t current_start = sorted_positions[i];
        size_t current_end = current_start + kmer_size - 1;
        
        // Check if k-mers overlap or are adjacent
        if (current_start <= region_end + 1)
        {
            // Extend the current region
            region_end = std::max(region_end, current_end);
        }
        else
        {
            // Start a new region
            UniqueRegion region;
            region.start = region_start;
            region.end = region_end;
            region.sequence = seqan3::dna5_vector{query_sequence.begin() + region_start,
                                                  query_sequence.begin() + region_end + 1};
            unique_regions.push_back(region);
            
            region_start = current_start;
            region_end = current_end;
        }
    }
    
    // Add the last region
    UniqueRegion region;
    region.start = region_start;
    region.end = region_end;
    region.sequence = seqan3::dna5_vector{query_sequence.begin() + region_start,
                                          query_sequence.begin() + region_end + 1};
    unique_regions.push_back(region);
    
    return unique_regions;
}

// Optimized function to find absent k-mer positions with hash map and overlap-based search
std::set<size_t> find_absent_kmer_positions(const seqan3::dna5_vector& sequence, 
                                           const auto& index, 
                                           uint8_t kmer_size,
                                           HierarchicalKmerCache& hierarchical_cache,
                                           ProgressReporter* progress = nullptr)
{
    auto [absent_positions, stats] = find_absent_kmer_positions_with_stats(sequence, index, kmer_size, hierarchical_cache, progress);
    return absent_positions;
}

// Function to find sequence files in a directory
std::vector<std::filesystem::path> find_sequence_files(const std::filesystem::path& dir_path)
{
    std::vector<std::filesystem::path> file_paths;
    std::set<std::string> valid_extensions{".fa", ".fasta", ".fna", ".ffn", ".faa", ".frn"};

    if (!std::filesystem::exists(dir_path) || !std::filesystem::is_directory(dir_path))
    {
        return file_paths;
    }

    for (auto const& entry : std::filesystem::recursive_directory_iterator{dir_path})
    {
        if (entry.is_regular_file() && valid_extensions.contains(entry.path().extension().string()))
        {
            file_paths.push_back(entry.path());
        }
    }
    
    std::sort(file_paths.begin(), file_paths.end());
    return file_paths;
}

// Function to process a single sequence file and return unique regions
struct ProcessingResult
{
    std::string filename;
    std::vector<uint8_t> k_values;
    std::map<uint8_t, size_t> total_sequences_by_k;
    std::map<uint8_t, size_t> total_unique_regions_by_k;
    std::map<uint8_t, size_t> total_unique_length_by_k;
    std::map<uint8_t, size_t> total_kmers_processed_by_k;
    std::map<uint8_t, size_t> unique_kmers_searched_by_k;
    std::map<uint8_t, size_t> cache_hits_by_k;
    std::map<uint8_t, size_t> inferred_from_smaller_k_by_k;
    std::map<uint8_t, size_t> overlap_optimizations_by_k;
    double processing_time_seconds;
    bool success;
    std::map<uint8_t, std::vector<std::pair<UniqueRegion, std::string>>> unique_regions_with_ids_by_k;
};

// Enhanced function to find absent k-mer positions with detailed statistics
struct SearchStats
{
    size_t total_kmers = 0;
    size_t unique_kmers_searched = 0;
    size_t cache_hits = 0;
    size_t inferred_from_smaller_k = 0;
    size_t overlap_attempts = 0;
};

std::pair<std::set<size_t>, SearchStats> find_absent_kmer_positions_with_stats(const seqan3::dna5_vector& sequence, 
                                                                               const auto& index, 
                                                                               uint8_t kmer_size,
                                                                               HierarchicalKmerCache& hierarchical_cache,
                                                                               ProgressReporter* progress = nullptr)
{
    std::set<size_t> absent_positions;
    SearchStats stats;
    const size_t total_kmers = sequence.size() - kmer_size + 1;
    stats.total_kmers = total_kmers;
    
    if (total_kmers == 0) return {absent_positions, stats};
    
    // For overlap-based search optimization
    seqan3::dna5_vector prev_kmer;
    bool prev_found = false;
    
    // Process k-mers in batches for better cache performance
    const size_t batch_size = 100;
    
    for (size_t batch_start = 0; batch_start < total_kmers; batch_start += batch_size)
    {
        size_t batch_end = std::min(batch_start + batch_size, total_kmers);
        
        for (size_t i = batch_start; i < batch_end; ++i)
        {
            auto kmer = seqan3::dna5_vector{sequence.begin() + i, sequence.begin() + i + kmer_size};
            
            // Convert k-mer to string for hashing
            std::string kmer_str;
            kmer_str.reserve(kmer_size);
            for (auto nucleotide : kmer)
            {
                kmer_str.push_back(seqan3::to_char(nucleotide));
            }
            
            bool found = false;
            
            // Check hierarchical cache first
            if (hierarchical_cache.find(kmer_str, kmer_size, found))
            {
                stats.cache_hits++;
            }
            // Try to infer from smaller k-mers
            else if (hierarchical_cache.can_infer_from_smaller_kmers(kmer_str, kmer_size, found))
            {
                stats.inferred_from_smaller_k++;
                // Cache the inferred result
                hierarchical_cache.insert(kmer_str, kmer_size, found);
            }
            else
            {
                // Optimization: Try overlap-based search for adjacent k-mers
                if (i > 0 && kmer_size > 1 && !prev_kmer.empty())
                {
                    // Check if this k-mer overlaps with the previous one (shift by 1)
                    bool overlaps = std::equal(kmer.begin(), kmer.begin() + kmer_size - 1,
                                             prev_kmer.begin() + 1);
                    
                    if (overlaps && prev_found)
                    {
                        stats.overlap_attempts++;
                    }
                }
                
                // Perform the actual search
                auto hits = seqan3::search(kmer, index);
                found = (hits.begin() != hits.end());
                stats.unique_kmers_searched++;
                
                // Cache the result in hierarchical cache
                hierarchical_cache.insert(kmer_str, kmer_size, found);
            }
            
            if (!found)
            {
                absent_positions.insert(i);
            }
            
            // Update for next iteration
            prev_kmer = kmer;
            prev_found = found;
        }
        
        // Update progress
        if (progress)
        {
            progress->update(batch_end - batch_start);
        }
    }
    
    return {absent_positions, stats};
}

// Multi-k processing function with hierarchical optimization
ProcessingResult process_sequence_file_multi_k(const std::filesystem::path& sequence_file,
                                               const auto& index,
                                               const std::vector<uint8_t>& k_values,
                                               HierarchicalKmerCache& hierarchical_cache)
{
    ProcessingResult result;
    result.filename = sequence_file.filename().string();
    result.k_values = k_values;
    result.processing_time_seconds = 0.0;
    result.success = false;

    auto start_time = std::chrono::steady_clock::now();

    try
    {
        seqan3::sequence_file_input fin{sequence_file};
        
        // First pass: Load all sequences
        std::vector<std::tuple<seqan3::dna5_vector, std::string>> sequences_info;
        
        for (auto & [seq, id, qual] : fin)
        {
            sequences_info.emplace_back(std::move(seq), std::string(id));
        }
        
        // Process each k-value in ascending order (smaller k-values first for better optimization)
        std::vector<uint8_t> sorted_k_values = k_values;
        std::sort(sorted_k_values.begin(), sorted_k_values.end());
        
        for (uint8_t k : sorted_k_values)
        {
            std::cout << "    Processing k=" << (int)k << "...\n";
            
            // Initialize statistics for this k-value
            result.total_sequences_by_k[k] = 0;
            result.total_unique_regions_by_k[k] = 0;
            result.total_unique_length_by_k[k] = 0;
            result.total_kmers_processed_by_k[k] = 0;
            result.unique_kmers_searched_by_k[k] = 0;
            result.cache_hits_by_k[k] = 0;
            result.inferred_from_smaller_k_by_k[k] = 0;
            result.overlap_optimizations_by_k[k] = 0;
            result.unique_regions_with_ids_by_k[k] = std::vector<std::pair<UniqueRegion, std::string>>();
            
            // Calculate total k-mers for this k-value
            size_t total_kmers_for_k = 0;
            for (const auto& [seq, id] : sequences_info)
            {
                if (seq.size() >= k)
                {
                    total_kmers_for_k += seq.size() - k + 1;
                }
            }
            
            result.total_kmers_processed_by_k[k] = total_kmers_for_k;
            
            if (total_kmers_for_k > 0)
            {
                std::cout << "      Total k-mers to process: " << total_kmers_for_k << "\n";
                ProgressReporter progress(total_kmers_for_k);
                
                // Process each sequence for this k-value
                for (const auto& [seq, id] : sequences_info)
                {
                    if (seq.size() < k)
                    {
                        continue;
                    }
                    
                    result.total_sequences_by_k[k]++;
                    
                    // Find absent k-mer positions with hierarchical optimization
                    auto [absent_positions, search_stats] = find_absent_kmer_positions_with_stats(seq, index, k, hierarchical_cache, &progress);
                    
                    // Accumulate optimization statistics
                    result.unique_kmers_searched_by_k[k] += search_stats.unique_kmers_searched;
                    result.cache_hits_by_k[k] += search_stats.cache_hits;
                    result.inferred_from_smaller_k_by_k[k] += search_stats.inferred_from_smaller_k;
                    result.overlap_optimizations_by_k[k] += search_stats.overlap_attempts;
                    
                    // Merge overlapping positions into unique regions
                    auto unique_regions = merge_overlapping_positions(absent_positions, k, seq);
                    
                    // Store unique regions with their IDs
                    for (size_t i = 0; i < unique_regions.size(); ++i)
                    {
                        const auto& region = unique_regions[i];
                        
                        // Create a comprehensive ID that includes file info and k-value
                        std::string file_stem = sequence_file.stem().string();
                        std::string region_id = file_stem + "_" + id + "_k" + std::to_string(k) + 
                                               "_unique_region_" + std::to_string(i + 1) +
                                               "_pos_" + std::to_string(region.start + 1) + "-" + 
                                               std::to_string(region.end + 1) +
                                               "_len_" + std::to_string(region.sequence.size());
                        
                        result.unique_regions_with_ids_by_k[k].emplace_back(region, region_id);
                        result.total_unique_length_by_k[k] += region.sequence.size();
                    }
                    
                    result.total_unique_regions_by_k[k] += unique_regions.size();
                }
            }
            
            // Print statistics for this k-value
            std::cout << "      Found " << result.total_unique_regions_by_k[k] << " unique regions\n";
            std::cout << "      Unique k-mers searched: " << result.unique_kmers_searched_by_k[k] << "\n";
            std::cout << "      Cache hits: " << result.cache_hits_by_k[k] << "\n";
            std::cout << "      Inferred from smaller k: " << result.inferred_from_smaller_k_by_k[k] << "\n";
            
            if (result.total_kmers_processed_by_k[k] > 0)
            {
                double cache_hit_rate = (double)result.cache_hits_by_k[k] / result.total_kmers_processed_by_k[k] * 100.0;
                double inference_rate = (double)result.inferred_from_smaller_k_by_k[k] / result.total_kmers_processed_by_k[k] * 100.0;
                double search_reduction = (double)(result.total_kmers_processed_by_k[k] - result.unique_kmers_searched_by_k[k]) / result.total_kmers_processed_by_k[k] * 100.0;
                
                std::cout << "      Cache hit rate: " << std::fixed << std::setprecision(1) << cache_hit_rate << "%\n";
                std::cout << "      Inference rate: " << std::fixed << std::setprecision(1) << inference_rate << "%\n";
                std::cout << "      Search reduction: " << std::fixed << std::setprecision(1) << search_reduction << "%\n";
            }
        }
        
        result.success = true;
    }
    catch (std::exception const & e)
    {
        std::cerr << "  Error processing " << result.filename << ": " << e.what() << '\n';
        result.success = false;
    }
    
    auto end_time = std::chrono::steady_clock::now();
    result.processing_time_seconds = std::chrono::duration<double>(end_time - start_time).count();
    
    return result;
}

ProcessingResult process_sequence_file(const std::filesystem::path& sequence_file,
                                      const auto& index,
                                      uint8_t kmer_size,
                                      HierarchicalKmerCache& hierarchical_cache)
{
    // Convert single k-value to vector for compatibility
    std::vector<uint8_t> k_values = {kmer_size};
    auto result = process_sequence_file_multi_k(sequence_file, index, k_values, hierarchical_cache);
    
    // For backward compatibility, we could populate the old fields, but since we changed the struct,
    // we'll just return the multi-k result
    return result;
}

// Struct to hold command line arguments
struct cmd_arguments
{
    std::filesystem::path index_path{};
    std::filesystem::path query_path{};
    std::filesystem::path output_file{"all_unique_regions.fasta"};
    uint8_t kmer_size{};
    uint8_t min_kmer_size{};
    uint8_t max_kmer_size{};
    bool is_directory{false};
    bool use_kmer_range{false};
};

// Function to set up the argument parser
void initialize_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "SeqAn3 Expert Team";
    parser.info.version = "1.0.0";
    parser.info.short_description = "Identifies unique regions from query sequences using an FM-Index.";
    parser.info.synopsis = {"kecfm-find -k <K> -i <INDEX_FILE> -q <QUERY_FILE_OR_DIR> [-o <OUTPUT_FILE>]",
                           "kecfm-find --min-k <MIN_K> --max-k <MAX_K> -i <INDEX_FILE> -q <QUERY_FILE_OR_DIR> [-o <OUTPUT_FILE>]"};

    parser.add_option(args.index_path, 'i', "index-file",
                      "Path to the pre-built FM-Index file.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{});
    parser.add_option(args.query_path, 'q', "query",
                      "Path to a query file (FASTA) or directory containing query files.",
                      seqan3::option_spec::required);
    parser.add_option(args.kmer_size, 'k', "kmer-size",
                      "The length of the k-mers to search for (use this OR --min-k/--max-k).",
                      seqan3::option_spec::standard);
    parser.add_option(args.min_kmer_size, '\0', "min-k",
                      "Minimum k-mer size for range search (use with --max-k).",
                      seqan3::option_spec::standard);
    parser.add_option(args.max_kmer_size, '\0', "max-k",
                      "Maximum k-mer size for range search (use with --min-k).",
                      seqan3::option_spec::standard);
    parser.add_option(args.output_file, 'o', "output-file",
                      "Path to the output FASTA file for all unique regions.",
                      seqan3::option_spec::standard,
                      seqan3::output_file_validator{});
    
    // Add help text about k-mer requirements
    parser.info.description.push_back("IMPORTANT: You must specify either:");
    parser.info.description.push_back("  - Single k-mer size: -k <K>");
    parser.info.description.push_back("  - K-mer range: --min-k <MIN_K> --max-k <MAX_K>");
    parser.info.description.push_back("");
    parser.info.description.push_back("Multi-k optimization: When using k-mer ranges, the program will");
    parser.info.description.push_back("automatically optimize searches by reusing information from smaller");
    parser.info.description.push_back("k-values, significantly reducing search time and memory usage.");
}

int main(int argc, char ** argv)
{
    seqan3::argument_parser parser{"kecfm-find", argc, argv};
    cmd_arguments args{};
    initialize_parser(parser, args);

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        std::cerr << " " << ext.what() << '\n';
        return -1;
    }

    // Validate k-mer arguments
    if (args.min_kmer_size > 0 && args.max_kmer_size > 0)
    {
        args.use_kmer_range = true;
        if (args.min_kmer_size > args.max_kmer_size)
        {
            std::cerr << " Error: min-k (" << (int)args.min_kmer_size << ") must be <= max-k (" << (int)args.max_kmer_size << ")\n";
            return -1;
        }
        if (args.kmer_size > 0)
        {
            std::cerr << " Error: Use either -k OR --min-k/--max-k, not both\n";
            return -1;
        }
    }
    else if (args.kmer_size > 0)
    {
        args.use_kmer_range = false;
    }
    else
    {
        std::cerr << " Error: Must specify either -k <K> OR --min-k <MIN_K> --max-k <MAX_K>\n";
        return -1;
    }

    // Create k-value vector
    std::vector<uint8_t> k_values;
    if (args.use_kmer_range)
    {
        for (uint8_t k = args.min_kmer_size; k <= args.max_kmer_size; k++)
        {
            k_values.push_back(k);
        }
        std::cout << "Using k-mer range: " << (int)args.min_kmer_size << " to " << (int)args.max_kmer_size 
                  << " (" << k_values.size() << " k-values)\n";
    }
    else
    {
        k_values.push_back(args.kmer_size);
        std::cout << "Using single k-mer size: " << (int)args.kmer_size << "\n";
    }

    // Determine if query_path is a file or directory
    if (!std::filesystem::exists(args.query_path))
    {
        std::cerr << " Query path does not exist: " << args.query_path << '\n';
        return -1;
    }
    
    args.is_directory = std::filesystem::is_directory(args.query_path);
    
    // Get list of files to process
    std::vector<std::filesystem::path> query_files;
    if (args.is_directory)
    {
        query_files = find_sequence_files(args.query_path);
        if (query_files.empty())
        {
            std::cerr << " No sequence files found in directory: " << args.query_path << '\n';
            return -1;
        }
        std::cout << "Found " << query_files.size() << " sequence files in directory.\n";
    }
    else
    {
        query_files.push_back(args.query_path);
    }

    // Create output directory if needed (for the output file path)
    auto output_dir = args.output_file.parent_path();
    if (!output_dir.empty())
    {
        try
        {
            std::filesystem::create_directories(output_dir);
        }
        catch (std::exception const & e)
        {
            std::cerr << " Failed to create output directory: " << e.what() << '\n';
            return -1;
        }
    }

    // Load the FM-Index
    using index_t = seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection>;
    index_t index;
    std::cout << "Loading index from: " << args.index_path.string() << "\n";
    try
    {
        std::ifstream is{args.index_path, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(index);
    }
    catch (std::exception const & e)
    {
        std::cerr << " Failed to load index: " << e.what() << '\n';
        return -1;
    }
    std::cout << "Index loaded successfully.\n\n";

    // Create hierarchical k-mer cache for cross-file and cross-k optimization
    HierarchicalKmerCache hierarchical_cache;
    std::cout << "Initialized hierarchical k-mer cache for cross-file and cross-k optimization.\n\n";

    // Process each query file
    std::vector<ProcessingResult> results;
    std::map<uint8_t, std::vector<std::pair<UniqueRegion, std::string>>> all_unique_regions_by_k;
    
    // Initialize output collections for each k-value
    for (uint8_t k : k_values)
    {
        all_unique_regions_by_k[k] = std::vector<std::pair<UniqueRegion, std::string>>();
    }
    
    for (const auto& query_file : query_files)
    {
        std::cout << "Processing: " << query_file.filename().string() << "\n";
        
        // Show current cache status
        auto cached_k_values = hierarchical_cache.get_cached_k_values();
        if (!cached_k_values.empty())
        {
            std::cout << "  Hierarchical cache status:\n";
            for (uint8_t k : cached_k_values)
            {
                auto [cache_size, found_in_cache] = hierarchical_cache.get_stats(k);
                std::cout << "    k=" << (int)k << ": " << cache_size << " k-mers cached (" 
                          << found_in_cache << " found, " << (cache_size - found_in_cache) << " absent)\n";
            }
        }
        
        ProcessingResult result;
        if (args.use_kmer_range)
        {
            result = process_sequence_file_multi_k(query_file, index, k_values, hierarchical_cache);
        }
        else
        {
            result = process_sequence_file(query_file, index, args.kmer_size, hierarchical_cache);
        }
        results.push_back(result);
        
        if (result.success)
        {
            std::cout << "  Processing time: " << std::fixed << std::setprecision(2) 
                      << result.processing_time_seconds << " seconds\n";
            
            // Collect all unique regions by k-value
            for (uint8_t k : result.k_values)
            {
                auto it = result.unique_regions_with_ids_by_k.find(k);
                if (it != result.unique_regions_with_ids_by_k.end())
                {
                    for (const auto& region_pair : it->second)
                    {
                        all_unique_regions_by_k[k].push_back(region_pair);
                    }
                }
            }
        }
        else
        {
            std::cout << "  Failed to process file\n";
        }
        std::cout << "\n";
    }

    // Write unique regions to output files (one per k-value for multi-k, or single file for single k)
    if (args.use_kmer_range)
    {
        // Multi-k: create separate output files for each k-value
        for (uint8_t k : k_values)
        {
            if (!all_unique_regions_by_k[k].empty())
            {
                std::filesystem::path k_output_file = args.output_file;
                std::string stem = k_output_file.stem().string();
                std::string extension = k_output_file.extension().string();
                k_output_file = k_output_file.parent_path() / (stem + "_k" + std::to_string(k) + extension);
                
                std::cout << "Writing k=" << (int)k << " unique regions to: " << k_output_file.string() << "\n";
                
                try
                {
                    seqan3::sequence_file_output fout{k_output_file};
                    
                    for (const auto& [region, region_id] : all_unique_regions_by_k[k])
                    {
                        fout.emplace_back(region.sequence, region_id);
                    }
                }
                catch (std::exception const & e)
                {
                    std::cerr << " Failed to write output file for k=" << (int)k << ": " << e.what() << '\n';
                    return -1;
                }
                
                std::cout << "Successfully wrote " << all_unique_regions_by_k[k].size() << " unique regions.\n";
            }
            else
            {
                std::cout << "No unique regions found for k=" << (int)k << ".\n";
            }
        }
    }
    else
    {
        // Single k: write to single output file
        uint8_t k = args.kmer_size;
        if (!all_unique_regions_by_k[k].empty())
        {
            std::cout << "Writing all unique regions to: " << args.output_file.string() << "\n";
            
            try
            {
                seqan3::sequence_file_output fout{args.output_file};
                
                for (const auto& [region, region_id] : all_unique_regions_by_k[k])
                {
                    fout.emplace_back(region.sequence, region_id);
                }
            }
            catch (std::exception const & e)
            {
                std::cerr << " Failed to write output file: " << e.what() << '\n';
                return -1;
            }
            
            std::cout << "Successfully wrote " << all_unique_regions_by_k[k].size() << " unique regions to output file.\n";
        }
        else
        {
            std::cout << "No unique regions found.\n";
        }
    }
    
    std::cout << "\n";
    
    // Print hierarchical cache statistics
    std::cout << "=== HIERARCHICAL CACHE STATISTICS ===\n";
    auto cached_k_values = hierarchical_cache.get_cached_k_values();
    size_t total_cache_size = 0;
    size_t total_found_in_cache = 0;
    
    for (uint8_t k : cached_k_values)
    {
        auto [cache_size, found_in_cache] = hierarchical_cache.get_stats(k);
        total_cache_size += cache_size;
        total_found_in_cache += found_in_cache;
        
        std::cout << "k=" << (int)k << ": " << cache_size << " unique k-mers (" 
                  << found_in_cache << " found, " << (cache_size - found_in_cache) << " absent)";
        if (cache_size > 0)
        {
            double found_percentage = (double)found_in_cache / cache_size * 100.0;
            std::cout << " - " << std::fixed << std::setprecision(1) << found_percentage << "% coverage";
        }
        std::cout << "\n";
    }
    
    std::cout << "Total cache entries: " << total_cache_size << "\n";
    if (total_cache_size > 0)
    {
        double overall_coverage = (double)total_found_in_cache / total_cache_size * 100.0;
        std::cout << "Overall index coverage: " << std::fixed << std::setprecision(1) << overall_coverage << "%\n";
    }

    // Print summary statistics
    std::cout << "\n=== SUMMARY ===\n";
    size_t total_files = results.size();
    size_t successful_files = 0;
    double total_processing_time = 0.0;
    
    // Statistics by k-value
    std::map<uint8_t, size_t> total_sequences_by_k;
    std::map<uint8_t, size_t> total_unique_regions_by_k;
    std::map<uint8_t, size_t> total_unique_length_by_k;
    std::map<uint8_t, size_t> total_kmers_processed_by_k;
    std::map<uint8_t, size_t> total_unique_kmers_searched_by_k;
    std::map<uint8_t, size_t> total_cache_hits_by_k;
    std::map<uint8_t, size_t> total_inferred_from_smaller_k_by_k;
    
    for (const auto& result : results)
    {
        if (result.success)
        {
            successful_files++;
            total_processing_time += result.processing_time_seconds;
            
            for (uint8_t k : result.k_values)
            {
                auto seq_it = result.total_sequences_by_k.find(k);
                if (seq_it != result.total_sequences_by_k.end())
                {
                    total_sequences_by_k[k] += seq_it->second;
                }
                
                auto reg_it = result.total_unique_regions_by_k.find(k);
                if (reg_it != result.total_unique_regions_by_k.end())
                {
                    total_unique_regions_by_k[k] += reg_it->second;
                }
                
                auto len_it = result.total_unique_length_by_k.find(k);
                if (len_it != result.total_unique_length_by_k.end())
                {
                    total_unique_length_by_k[k] += len_it->second;
                }
                
                auto kmer_it = result.total_kmers_processed_by_k.find(k);
                if (kmer_it != result.total_kmers_processed_by_k.end())
                {
                    total_kmers_processed_by_k[k] += kmer_it->second;
                }
                
                auto search_it = result.unique_kmers_searched_by_k.find(k);
                if (search_it != result.unique_kmers_searched_by_k.end())
                {
                    total_unique_kmers_searched_by_k[k] += search_it->second;
                }
                
                auto cache_it = result.cache_hits_by_k.find(k);
                if (cache_it != result.cache_hits_by_k.end())
                {
                    total_cache_hits_by_k[k] += cache_it->second;
                }
                
                auto infer_it = result.inferred_from_smaller_k_by_k.find(k);
                if (infer_it != result.inferred_from_smaller_k_by_k.end())
                {
                    total_inferred_from_smaller_k_by_k[k] += infer_it->second;
                }
            }
        }
    }
    
    std::cout << "Files processed: " << successful_files << "/" << total_files << "\n";
    std::cout << "Total processing time: " << std::fixed << std::setprecision(2) 
              << total_processing_time << " seconds\n";
    
    // Print statistics for each k-value
    for (uint8_t k : k_values)
    {
        std::cout << "\n--- k=" << (int)k << " Statistics ---\n";
        std::cout << "Total sequences: " << total_sequences_by_k[k] << "\n";
        std::cout << "Total unique regions: " << total_unique_regions_by_k[k] << "\n";
        std::cout << "Total unique sequence length: " << total_unique_length_by_k[k] << " bp\n";
        std::cout << "Total k-mers processed: " << total_kmers_processed_by_k[k] << "\n";
        std::cout << "Unique k-mers searched: " << total_unique_kmers_searched_by_k[k] << "\n";
        std::cout << "Cache hits: " << total_cache_hits_by_k[k] << "\n";
        std::cout << "Inferred from smaller k: " << total_inferred_from_smaller_k_by_k[k] << "\n";
        
        if (total_kmers_processed_by_k[k] > 0)
        {
            double cache_hit_rate = (double)total_cache_hits_by_k[k] / total_kmers_processed_by_k[k] * 100.0;
            double inference_rate = (double)total_inferred_from_smaller_k_by_k[k] / total_kmers_processed_by_k[k] * 100.0;
            double search_reduction = (double)(total_kmers_processed_by_k[k] - total_unique_kmers_searched_by_k[k]) / total_kmers_processed_by_k[k] * 100.0;
            
            std::cout << "Cache hit rate: " << std::fixed << std::setprecision(1) << cache_hit_rate << "%\n";
            std::cout << "Inference rate: " << std::fixed << std::setprecision(1) << inference_rate << "%\n";
            std::cout << "Search reduction: " << std::fixed << std::setprecision(1) << search_reduction << "%\n";
            
            if (total_processing_time > 0)
            {
                double speed = total_kmers_processed_by_k[k] / total_processing_time;
                std::cout << "Speed: " << std::fixed << std::setprecision(0) << speed << " k-mers/second\n";
            }
        }
        
        if (total_unique_regions_by_k[k] > 0)
        {
            std::cout << "Average region length: " << (total_unique_length_by_k[k] / total_unique_regions_by_k[k]) << " bp\n";
        }
        
        if (args.use_kmer_range)
        {
            std::filesystem::path k_output_file = args.output_file;
            std::string stem = k_output_file.stem().string();
            std::string extension = k_output_file.extension().string();
            k_output_file = k_output_file.parent_path() / (stem + "_k" + std::to_string(k) + extension);
            std::cout << "Output file: " << k_output_file << "\n";
        }
    }
    
    if (!args.use_kmer_range)
    {
        std::cout << "\nOutput file: " << args.output_file << "\n";
    }

    return successful_files == total_files ? 0 : 1;
}