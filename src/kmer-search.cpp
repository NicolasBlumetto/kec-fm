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

// Global k-mer cache for cross-file optimization
class GlobalKmerCache
{
private:
    std::unordered_map<std::string, bool> cache_;
    mutable std::mutex cache_mutex_;
    
public:
    bool find(const std::string& kmer, bool& found_in_index)
    {
        std::lock_guard<std::mutex> lock(cache_mutex_);
        auto it = cache_.find(kmer);
        if (it != cache_.end())
        {
            found_in_index = it->second;
            return true; // found in cache
        }
        return false; // not in cache
    }
    
    void insert(const std::string& kmer, bool found_in_index)
    {
        std::lock_guard<std::mutex> lock(cache_mutex_);
        cache_[kmer] = found_in_index;
    }
    
    size_t size() const
    {
        std::lock_guard<std::mutex> lock(cache_mutex_);
        return cache_.size();
    }
    
    void clear()
    {
        std::lock_guard<std::mutex> lock(cache_mutex_);
        cache_.clear();
    }
    
    // Get cache statistics
    std::pair<size_t, size_t> get_stats() const
    {
        std::lock_guard<std::mutex> lock(cache_mutex_);
        size_t total = cache_.size();
        size_t found_count = 0;
        for (const auto& [kmer, found] : cache_)
        {
            if (found) found_count++;
        }
        return {total, found_count};
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
                                           GlobalKmerCache& global_cache,
                                           ProgressReporter* progress = nullptr)
{
    auto [absent_positions, stats] = find_absent_kmer_positions_with_stats(sequence, index, kmer_size, global_cache, progress);
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
    size_t total_sequences;
    size_t total_unique_regions;
    size_t total_unique_length;
    size_t total_kmers_processed;
    size_t unique_kmers_searched;  // Number of unique k-mers actually searched
    size_t cache_hits;            // Number of k-mers found in cache
    size_t overlap_optimizations; // Number of overlap-based optimizations attempted
    double processing_time_seconds;
    bool success;
    std::vector<std::pair<UniqueRegion, std::string>> unique_regions_with_ids; // region and full sequence ID
};

// Enhanced function to find absent k-mer positions with detailed statistics
struct SearchStats
{
    size_t total_kmers = 0;
    size_t unique_kmers_searched = 0;
    size_t cache_hits = 0;
    size_t overlap_attempts = 0;
};

std::pair<std::set<size_t>, SearchStats> find_absent_kmer_positions_with_stats(const seqan3::dna5_vector& sequence, 
                                                                               const auto& index, 
                                                                               uint8_t kmer_size,
                                                                               GlobalKmerCache& global_cache,
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
            
            // Check global cache first
            if (global_cache.find(kmer_str, found))
            {
                stats.cache_hits++;
            }
            else
            {
                // Optimization 2: Try overlap-based search for adjacent k-mers
                if (i > 0 && kmer_size > 1 && !prev_kmer.empty())
                {
                    // Check if this k-mer overlaps with the previous one (shift by 1)
                    bool overlaps = std::equal(kmer.begin(), kmer.begin() + kmer_size - 1,
                                             prev_kmer.begin() + 1);
                    
                    if (overlaps && prev_found)
                    {
                        // Count overlap optimization attempts
                        stats.overlap_attempts++;
                    }
                }
                
                // Perform the actual search
                auto hits = seqan3::search(kmer, index);
                found = (hits.begin() != hits.end());
                stats.unique_kmers_searched++;
                
                // Cache the result in global cache
                global_cache.insert(kmer_str, found);
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

ProcessingResult process_sequence_file(const std::filesystem::path& sequence_file,
                                      const auto& index,
                                      uint8_t kmer_size,
                                      GlobalKmerCache& global_cache)
{
    ProcessingResult result;
    result.filename = sequence_file.filename().string();
    result.total_sequences = 0;
    result.total_unique_regions = 0;
    result.total_unique_length = 0;
    result.total_kmers_processed = 0;
    result.unique_kmers_searched = 0;
    result.cache_hits = 0;
    result.overlap_optimizations = 0;
    result.processing_time_seconds = 0.0;
    result.success = false;

    auto start_time = std::chrono::steady_clock::now();

    try
    {
        seqan3::sequence_file_input fin{sequence_file};
        
        // First pass: count sequences and total k-mers for progress reporting
        std::vector<std::tuple<seqan3::dna5_vector, std::string, size_t>> sequences_info; // seq, id, kmer_count
        size_t total_kmers_in_file = 0;
        
        for (auto & [seq, id, qual] : fin)
        {
            if (seq.size() >= kmer_size)
            {
                size_t kmers_in_seq = seq.size() - kmer_size + 1;
                sequences_info.emplace_back(std::move(seq), std::string(id), kmers_in_seq);
                total_kmers_in_file += kmers_in_seq;
            }
            else
            {
                sequences_info.emplace_back(std::move(seq), std::string(id), 0);
            }
        }
        
        result.total_sequences = sequences_info.size();
        result.total_kmers_processed = total_kmers_in_file;
        
        if (total_kmers_in_file > 0)
        {
            std::cout << "  Total k-mers to process: " << total_kmers_in_file << "\n";
            ProgressReporter progress(total_kmers_in_file);
            
            // Process each sequence
            for (const auto& [seq, id, kmer_count] : sequences_info)
            {
                if (seq.size() < kmer_size)
                {
                    std::cout << "  Warning: Sequence '" << id << "' is shorter than k-mer size, skipping.\n";
                    continue;
                }
                
                // Find absent k-mer positions with optimized search and collect stats
                auto [absent_positions, search_stats] = find_absent_kmer_positions_with_stats(seq, index, kmer_size, global_cache, &progress);
                
                // Accumulate optimization statistics
                result.unique_kmers_searched += search_stats.unique_kmers_searched;
                result.cache_hits += search_stats.cache_hits;
                result.overlap_optimizations += search_stats.overlap_attempts;
                
                // Merge overlapping positions into unique regions
                auto unique_regions = merge_overlapping_positions(absent_positions, kmer_size, seq);
                
                // Store unique regions with their IDs for later output
                for (size_t i = 0; i < unique_regions.size(); ++i)
                {
                    const auto& region = unique_regions[i];
                    
                    // Create a comprehensive ID that includes file info
                    std::string file_stem = sequence_file.stem().string();
                    std::string region_id = file_stem + "_" + id + "_unique_region_" + std::to_string(i + 1) +
                                           "_pos_" + std::to_string(region.start + 1) + "-" + 
                                           std::to_string(region.end + 1) +
                                           "_len_" + std::to_string(region.sequence.size());
                    
                    result.unique_regions_with_ids.emplace_back(region, region_id);
                    result.total_unique_length += region.sequence.size();
                }
                
                result.total_unique_regions += unique_regions.size();
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

// Struct to hold command line arguments
struct cmd_arguments
{
    std::filesystem::path index_path{};
    std::filesystem::path query_path{};
    std::filesystem::path output_file{"all_unique_regions.fasta"};
    uint8_t kmer_size{};
    bool is_directory{false};
};

// Function to set up the argument parser
void initialize_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "SeqAn3 Expert Team";
    parser.info.version = "1.0.0";
    parser.info.short_description = "Identifies unique regions from query sequences using an FM-Index.";
    parser.info.synopsis = {"kecfm-find -k <K> -i <INDEX_FILE> -q <QUERY_FILE_OR_DIR> [-o <OUTPUT_FILE>]"};

    parser.add_option(args.index_path, 'i', "index-file",
                      "Path to the pre-built FM-Index file.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{});
    parser.add_option(args.query_path, 'q', "query",
                      "Path to a query file (FASTA) or directory containing query files.",
                      seqan3::option_spec::required);
    parser.add_option(args.kmer_size, 'k', "kmer-size",
                      "The length of the k-mers to search for.",
                      seqan3::option_spec::required);
    parser.add_option(args.output_file, 'o', "output-file",
                      "Path to the output FASTA file for all unique regions.",
                      seqan3::option_spec::standard,
                      seqan3::output_file_validator{});
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

    // Create global k-mer cache for cross-file optimization
    GlobalKmerCache global_cache;
    std::cout << "Initialized global k-mer cache for cross-file optimization.\n\n";

    // Process each query file
    std::vector<ProcessingResult> results;
    std::vector<std::pair<UniqueRegion, std::string>> all_unique_regions;
    
    for (const auto& query_file : query_files)
    {
        std::cout << "Processing: " << query_file.filename().string() << "\n";
        
        // Show current cache status
        auto [cache_size, found_in_cache] = global_cache.get_stats();
        if (cache_size > 0)
        {
            std::cout << "  Global cache: " << cache_size << " k-mers cached (" 
                      << found_in_cache << " found, " << (cache_size - found_in_cache) << " absent)\n";
        }
        
        auto result = process_sequence_file(query_file, index, args.kmer_size, global_cache);
        results.push_back(result);
        
        if (result.success)
        {
            std::cout << "  Processed " << result.total_sequences << " sequences\n";
            std::cout << "  Found " << result.total_unique_regions << " unique regions\n";
            std::cout << "  Total unique sequence length: " << result.total_unique_length << " bp\n";
            std::cout << "  K-mers processed: " << result.total_kmers_processed << "\n";
            std::cout << "  Unique k-mers searched: " << result.unique_kmers_searched << "\n";
            std::cout << "  Cache hits: " << result.cache_hits << "\n";
            std::cout << "  Overlap optimizations: " << result.overlap_optimizations << "\n";
            if (result.total_kmers_processed > 0)
            {
                double cache_hit_rate = (double)result.cache_hits / result.total_kmers_processed * 100.0;
                double search_reduction = (double)(result.total_kmers_processed - result.unique_kmers_searched) / result.total_kmers_processed * 100.0;
                std::cout << "  Cache hit rate: " << std::fixed << std::setprecision(1) << cache_hit_rate << "%\n";
                std::cout << "  Search reduction: " << std::fixed << std::setprecision(1) << search_reduction << "%\n";
            }
            std::cout << "  Processing time: " << std::fixed << std::setprecision(2) 
                      << result.processing_time_seconds << " seconds\n";
            if (result.processing_time_seconds > 0)
            {
                double kmers_per_second = result.total_kmers_processed / result.processing_time_seconds;
                std::cout << "  Speed: " << std::fixed << std::setprecision(0) 
                          << kmers_per_second << " k-mers/second\n";
            }
            
            // Collect all unique regions
            for (const auto& region_pair : result.unique_regions_with_ids)
            {
                all_unique_regions.push_back(region_pair);
            }
        }
        else
        {
            std::cout << "  Failed to process file\n";
        }
        std::cout << "\n";
    }

    // Write all unique regions to a single output file
    if (!all_unique_regions.empty())
    {
        std::cout << "Writing all unique regions to: " << args.output_file.string() << "\n";
        
        try
        {
            seqan3::sequence_file_output fout{args.output_file};
            
            for (const auto& [region, region_id] : all_unique_regions)
            {
                fout.emplace_back(region.sequence, region_id);
            }
        }
        catch (std::exception const & e)
        {
            std::cerr << " Failed to write output file: " << e.what() << '\n';
            return -1;
        }
        
        std::cout << "Successfully wrote " << all_unique_regions.size() << " unique regions to output file.\n\n";
    }
    else
    {
        std::cout << "No unique regions found across all files.\n\n";
    }
    
    // Final global cache statistics
    auto [final_cache_size, final_found_in_cache] = global_cache.get_stats();
    std::cout << "\n=== GLOBAL CACHE STATISTICS ===\n";
    std::cout << "Total unique k-mers encountered: " << final_cache_size << "\n";
    std::cout << "K-mers found in index: " << final_found_in_cache << "\n";
    std::cout << "K-mers absent from index: " << (final_cache_size - final_found_in_cache) << "\n";
    if (final_cache_size > 0)
    {
        double found_percentage = (double)final_found_in_cache / final_cache_size * 100.0;
        std::cout << "Index coverage: " << std::fixed << std::setprecision(1) << found_percentage << "%\n";
    }

    // Print summary statistics
    std::cout << "\n=== SUMMARY ===\n";
    size_t total_files = results.size();
    size_t successful_files = 0;
    size_t total_sequences = 0;
    size_t total_unique_regions = 0;
    size_t total_unique_length = 0;
    size_t total_kmers_processed = 0;
    size_t total_unique_kmers_searched = 0;
    size_t total_cache_hits = 0;
    size_t total_overlap_optimizations = 0;
    double total_processing_time = 0.0;
    
    for (const auto& result : results)
    {
        if (result.success)
        {
            successful_files++;
            total_sequences += result.total_sequences;
            total_unique_regions += result.total_unique_regions;
            total_unique_length += result.total_unique_length;
            total_kmers_processed += result.total_kmers_processed;
            total_unique_kmers_searched += result.unique_kmers_searched;
            total_cache_hits += result.cache_hits;
            total_overlap_optimizations += result.overlap_optimizations;
            total_processing_time += result.processing_time_seconds;
        }
    }
    
    std::cout << "Files processed: " << successful_files << "/" << total_files << "\n";
    std::cout << "Total sequences: " << total_sequences << "\n";
    std::cout << "Total unique regions: " << total_unique_regions << "\n";
    std::cout << "Total unique sequence length: " << total_unique_length << " bp\n";
    std::cout << "Total k-mers processed: " << total_kmers_processed << "\n";
    std::cout << "Total unique k-mers searched: " << total_unique_kmers_searched << "\n";
    std::cout << "Total cache hits: " << total_cache_hits << "\n";
    std::cout << "Total overlap optimizations: " << total_overlap_optimizations << "\n";
    
    if (total_kmers_processed > 0)
    {
        double overall_cache_hit_rate = (double)total_cache_hits / total_kmers_processed * 100.0;
        double overall_search_reduction = (double)(total_kmers_processed - total_unique_kmers_searched) / total_kmers_processed * 100.0;
        std::cout << "Overall cache hit rate: " << std::fixed << std::setprecision(1) << overall_cache_hit_rate << "%\n";
        std::cout << "Overall search reduction: " << std::fixed << std::setprecision(1) << overall_search_reduction << "%\n";
        
        // Calculate global cache effectiveness
        if (final_cache_size > 0)
        {
            double cache_utilization = (double)total_cache_hits / (total_kmers_processed - total_unique_kmers_searched + total_cache_hits) * 100.0;
            double deduplication_ratio = (double)total_kmers_processed / final_cache_size;
            std::cout << "Global cache utilization: " << std::fixed << std::setprecision(1) << cache_utilization << "%\n";
            std::cout << "K-mer deduplication ratio: " << std::fixed << std::setprecision(1) << deduplication_ratio << ":1\n";
        }
    }
    
    std::cout << "Total processing time: " << std::fixed << std::setprecision(2) 
              << total_processing_time << " seconds\n";
    
    if (total_unique_regions > 0)
    {
        std::cout << "Average region length: " << (total_unique_length / total_unique_regions) << " bp\n";
    }
    
    if (total_processing_time > 0)
    {
        double overall_speed = total_kmers_processed / total_processing_time;
        std::cout << "Overall speed: " << std::fixed << std::setprecision(0) 
                  << overall_speed << " k-mers/second\n";
    }
    
    std::cout << "Output file: " << args.output_file << "\n";

    return successful_files == total_files ? 0 : 1;
}