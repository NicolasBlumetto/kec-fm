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
#include <numeric> // For std::iota

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
        
        // Update only periodically or at the end to avoid excessive console output
        if (total_ > 0 && (current % (total_ / 1000 + 1) == 0 || current == total_)) // Update at least 1000 times or for small totals
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
        for (const auto& pair : cache_) // Use const auto& pair
        {
            if (pair.second) found_count++; // pair.second is the boolean value
        }
        return {total, found_count};
    }
};

// Enhanced function to find absent k-mer positions with detailed statistics
struct SearchStats
{
    size_t total_kmers = 0;
    size_t unique_kmers_searched = 0; // K-mers that went to FM-index or cache
    size_t cache_hits = 0;
    size_t mem_skipped_kmers = 0; // K-mers skipped because they were covered by MEMs
};

// Simple function to find basic exact matches (disabled for now)
// We'll focus on optimizing the k-mer search instead of complex MEM finding

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

std::pair<std::set<size_t>, SearchStats> find_absent_kmer_positions_with_stats(const seqan3::dna5_vector& sequence, 
                                                                               const auto& index, 
                                                                               uint8_t kmer_size,
                                                                               GlobalKmerCache& global_cache,
                                                                               const std::vector<bool>& mem_covered_positions, // New: MEM coverage map
                                                                               ProgressReporter* progress = nullptr)
{
    std::set<size_t> absent_positions;
    SearchStats stats;
    const size_t total_kmers = sequence.size() >= kmer_size ? sequence.size() - kmer_size + 1 : 0;
    stats.total_kmers = total_kmers;
    
    if (total_kmers == 0) return {absent_positions, stats};
    
    for (size_t i = 0; i < total_kmers; ++i)
    {
        // Check if the current k-mer is fully covered by an MEM
        bool fully_mem_covered = true;
        for (size_t j = 0; j < kmer_size; ++j)
        {
            if (!mem_covered_positions[i + j])
            {
                fully_mem_covered = false;
                break;
            }
        }

        if (fully_mem_covered)
        {
            stats.mem_skipped_kmers++;
            // This k-mer is covered by an MEM, so it's considered "found"
            // We don't need to search the cache or index for it.
            if (progress) progress->update();
            continue; 
        }

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
            // Perform the actual search
            auto hits = seqan3::search(kmer, index);
            found = (hits.begin() != hits.end());
            stats.unique_kmers_searched++; // This counts k-mers that went to FM-index
            
            // Cache the result in global cache
            global_cache.insert(kmer_str, found);
        }
        
        if (!found)
        {
            absent_positions.insert(i);
        }
        
        // Update progress
        if (progress)
        {
            progress->update();
        }
    }
    
    return {absent_positions, stats};
}

// Advanced k-mer search using non-overlapping sampling and extension
std::pair<std::set<size_t>, SearchStats> find_absent_kmer_positions_with_extension(const seqan3::dna5_vector& sequence, 
                                                                                   const auto& index, 
                                                                                   uint8_t kmer_size,
                                                                                   GlobalKmerCache& global_cache,
                                                                                   ProgressReporter* progress = nullptr)
{
    std::set<size_t> absent_positions;
    SearchStats stats;
    const size_t total_kmers = sequence.size() >= kmer_size ? sequence.size() - kmer_size + 1 : 0;
    stats.total_kmers = total_kmers;
    
    if (total_kmers == 0) return {absent_positions, stats};
    
    // Track which positions we've already determined (found/not found)
    std::vector<int> position_status(total_kmers, -1); // -1 = unknown, 0 = absent, 1 = present
    
    // Step 1: Sample non-overlapping k-mers
    std::vector<size_t> sample_positions;
    for (size_t i = 0; i < total_kmers; i += kmer_size)
    {
        sample_positions.push_back(i);
    }
    
    std::cout << "    Sampling " << sample_positions.size() << " non-overlapping k-mers (vs " << total_kmers << " total)...\n";
    
    // Step 2: Search sampled k-mers and extend matches
    for (size_t pos : sample_positions)
    {
        if (position_status[pos] != -1) continue; // Already processed
        
        auto kmer = seqan3::dna5_vector{sequence.begin() + pos, sequence.begin() + pos + kmer_size};
        
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
            // Perform the actual search
            auto hits = seqan3::search(kmer, index);
            found = (hits.begin() != hits.end());
            stats.unique_kmers_searched++;
            
            // Cache the result
            global_cache.insert(kmer_str, found);
        }
        
        if (found)
        {
            // This k-mer is found in the index, so try to extend the match
            // to find the full region that exists in the index
            
            // Find the extent of the match region by extending carefully
            size_t left_start = pos;
            size_t right_end = pos + kmer_size - 1;
            
            // Extend to the right with binary search for efficiency
            size_t right_bound = std::min(pos + 100, total_kmers); // Limit extension to avoid excessive work
            for (size_t ext_pos = pos + 1; ext_pos < right_bound; ++ext_pos)
            {
                if (position_status[ext_pos] == 1) {
                    right_end = ext_pos + kmer_size - 1;
                    continue;
                } else if (position_status[ext_pos] == 0) {
                    break;
                }
                
                auto extended_kmer = seqan3::dna5_vector{sequence.begin() + ext_pos, sequence.begin() + ext_pos + kmer_size};
                
                std::string ext_kmer_str;
                ext_kmer_str.reserve(kmer_size);
                for (auto nucleotide : extended_kmer)
                {
                    ext_kmer_str.push_back(seqan3::to_char(nucleotide));
                }
                
                bool ext_found = false;
                if (global_cache.find(ext_kmer_str, ext_found))
                {
                    stats.cache_hits++;
                }
                else
                {
                    auto ext_hits = seqan3::search(extended_kmer, index);
                    ext_found = (ext_hits.begin() != ext_hits.end());
                    stats.unique_kmers_searched++;
                    global_cache.insert(ext_kmer_str, ext_found);
                }
                
                if (ext_found)
                {
                    right_end = ext_pos + kmer_size - 1;
                }
                else
                {
                    break;
                }
            }
            
            // Extend to the left with limits
            size_t left_bound = (pos >= 100) ? pos - 100 : 0;
            if (pos > 0)
            {
                for (size_t ext_pos = pos - 1; ext_pos >= left_bound && ext_pos != SIZE_MAX; --ext_pos)
                {
                    if (position_status[ext_pos] == 1) {
                        left_start = ext_pos;
                        continue;
                    } else if (position_status[ext_pos] == 0) {
                        break;
                    }
                    
                    auto extended_kmer = seqan3::dna5_vector{sequence.begin() + ext_pos, sequence.begin() + ext_pos + kmer_size};
                    
                    std::string ext_kmer_str;
                    ext_kmer_str.reserve(kmer_size);
                    for (auto nucleotide : extended_kmer)
                    {
                        ext_kmer_str.push_back(seqan3::to_char(nucleotide));
                    }
                    
                    bool ext_found = false;
                    if (global_cache.find(ext_kmer_str, ext_found))
                    {
                        stats.cache_hits++;
                    }
                    else
                    {
                        auto ext_hits = seqan3::search(extended_kmer, index);
                        ext_found = (ext_hits.begin() != ext_hits.end());
                        stats.unique_kmers_searched++;
                        global_cache.insert(ext_kmer_str, ext_found);
                    }
                    
                    if (ext_found)
                    {
                        left_start = ext_pos;
                    }
                    else
                    {
                        break;
                    }
                }
            }
            
            // Mark all positions in this extended region as found
            for (size_t mark_pos = left_start; mark_pos <= right_end - kmer_size + 1 && mark_pos < total_kmers; ++mark_pos)
            {
                position_status[mark_pos] = 1; // present
            }
        }
        else
        {
            // This k-mer is not found
            position_status[pos] = 0; // absent
            absent_positions.insert(pos);
        }
    }
    
    // Step 3: Fill in gaps - check positions that are still unknown
    size_t gaps_filled = 0;
    for (size_t i = 0; i < total_kmers; ++i)
    {
        if (position_status[i] == -1)
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
                // Perform the actual search
                auto hits = seqan3::search(kmer, index);
                found = (hits.begin() != hits.end());
                stats.unique_kmers_searched++;
                
                // Cache the result
                global_cache.insert(kmer_str, found);
            }
            
            if (!found)
            {
                absent_positions.insert(i);
            }
            
            gaps_filled++;
        }
    }
    
    if (gaps_filled > 0)
    {
        std::cout << "    Filled " << gaps_filled << " gap positions.\n";
    }
    
    // Update progress appropriately
    if (progress)
    {
        progress->update(total_kmers);
    }
    
    return {absent_positions, stats};
}

// Simple overlapping k-mer search - the original fast approach
std::pair<std::set<size_t>, SearchStats> find_absent_kmer_positions_simple(const seqan3::dna5_vector& sequence, 
                                                                           const auto& index, 
                                                                           uint8_t kmer_size,
                                                                           GlobalKmerCache& global_cache,
                                                                           ProgressReporter* progress = nullptr)
{
    std::set<size_t> absent_positions;
    SearchStats stats;
    const size_t total_kmers = sequence.size() >= kmer_size ? sequence.size() - kmer_size + 1 : 0;
    stats.total_kmers = total_kmers;
    
    if (total_kmers == 0) return {absent_positions, stats};
    
    for (size_t i = 0; i < total_kmers; ++i)
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
            // Perform the actual search
            auto hits = seqan3::search(kmer, index);
            found = (hits.begin() != hits.end());
            stats.unique_kmers_searched++;
            
            // Cache the result
            global_cache.insert(kmer_str, found);
        }
        
        if (!found)
        {
            absent_positions.insert(i);
        }
        
        // Update progress if reporter is provided
        if (progress)
        {
            progress->update();
        }
    }
    
    return {absent_positions, stats};
}

// Structure to hold statistics for MEM search
struct MemSearchStats
{
    size_t total_query_bases = 0;
    size_t total_mem_hits = 0;
    size_t total_mem_length = 0;
    size_t bases_covered_by_mems = 0;
};

// Function to process a single sequence file and return unique regions
struct ProcessingResult
{
    std::string filename;
    size_t total_sequences;
    size_t total_unique_regions;
    size_t total_unique_length;
    size_t total_kmers_processed;
    size_t unique_kmers_searched;  // Number of unique k-mers actually searched (index/cache)
    size_t cache_hits;            // Number of k-mers found in cache
    size_t mem_skipped_kmers;     // Number of k-mers skipped due to MEM coverage
    MemSearchStats mem_stats;     // Statistics for MEM finding
    double processing_time_seconds;
    bool success;
    std::vector<std::pair<UniqueRegion, std::string>> unique_regions_with_ids; // region and full sequence ID
};

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
    result.mem_skipped_kmers = 0;
    result.processing_time_seconds = 0.0;
    result.success = false;

    auto start_time = std::chrono::steady_clock::now();

    try
    {
        seqan3::sequence_file_input fin{sequence_file};
        
        std::vector<std::tuple<seqan3::dna5_vector, std::string, size_t>> sequences_info; // seq, id, kmer_count
        size_t total_kmers_in_file = 0;
        size_t total_bases_in_file = 0;
        
        for (auto & [seq, id, qual] : fin)
        {
            total_bases_in_file += seq.size();
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
        result.mem_stats.total_query_bases = total_bases_in_file;

        if (total_kmers_in_file == 0 && total_bases_in_file == 0)
        {
            std::cout << "  No valid sequences or k-mers found in file: " << result.filename << "\n";
            result.success = true; // Still a successful processing, just nothing to do
            auto end_time = std::chrono::steady_clock::now();
            result.processing_time_seconds = std::chrono::duration<double>(end_time - start_time).count();
            return result;
        }

        std::cout << "  Processing k-mers with cursor-based sampling and extension...\n";

        if (total_kmers_in_file > 0)
        {
            ProgressReporter progress(total_kmers_in_file);
            
            // Process each sequence for absent k-mers
            size_t current_seq_idx = 0;
            for (const auto& [seq, id, kmer_count] : sequences_info)
            {
                if (seq.size() < kmer_size)
                {
                    // Update progress for skipped k-mers if applicable
                    if (kmer_count > 0) progress.update(kmer_count);
                    std::cout << "  Warning: Sequence '" << id << "' is shorter than k-mer size, skipping k-mer search.\n";
                    current_seq_idx++;
                    continue;
                }
                
                // Find absent k-mer positions using cursor-based approach and collect stats
                auto [absent_positions, search_stats] = find_absent_kmer_positions_cursor(
                    seq, index, kmer_size, global_cache, &progress
                );
                
                // Accumulate optimization statistics
                result.unique_kmers_searched += search_stats.unique_kmers_searched;
                result.cache_hits += search_stats.cache_hits;
                result.mem_skipped_kmers += search_stats.mem_skipped_kmers;
                
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
                current_seq_idx++;
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
            std::cout << "  Total k-mers processed: " << result.total_kmers_processed << "\n";
            std::cout << "  Unique k-mers searched (FM-index): " << result.unique_kmers_searched << "\n";
            std::cout << "  Cache hits: " << result.cache_hits << "\n";
            
            if (result.total_kmers_processed > 0)
            {
                double cache_hit_rate = (result.total_kmers_processed > 0) ? (double)result.cache_hits / result.total_kmers_processed * 100.0 : 0.0;
                double search_reduction = (double)result.cache_hits / result.total_kmers_processed * 100.0;
                
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
            total_processing_time += result.processing_time_seconds;
        }
    }
    
    std::cout << "Files processed: " << successful_files << "/" << total_files << "\n";
    std::cout << "Total sequences: " << total_sequences << "\n";
    std::cout << "Total unique regions: " << total_unique_regions << "\n";
    std::cout << "Total unique sequence length: " << total_unique_length << " bp\n";
    std::cout << "Total k-mers processed: " << total_kmers_processed << "\n";
    std::cout << "Total unique k-mers searched (FM-index): " << total_unique_kmers_searched << "\n";
    std::cout << "Total cache hits: " << total_cache_hits << "\n";
    
    if (total_kmers_processed > 0)
    {
        double overall_cache_hit_rate = (double)total_cache_hits / total_kmers_processed * 100.0;
        double overall_search_reduction = (double)total_cache_hits / total_kmers_processed * 100.0;
        
        std::cout << "Overall cache hit rate: " << std::fixed << std::setprecision(1) << overall_cache_hit_rate << "%\n";
        std::cout << "Overall search reduction: " << std::fixed << std::setprecision(1) << overall_search_reduction << "%\n";
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

// Cursor-based k-mer search with non-overlapping sampling and efficient extension
std::pair<std::set<size_t>, SearchStats> find_absent_kmer_positions_cursor(const seqan3::dna5_vector& sequence, 
                                                                           const auto& index, 
                                                                           uint8_t kmer_size,
                                                                           GlobalKmerCache& global_cache,
                                                                           ProgressReporter* progress = nullptr)
{
    std::set<size_t> absent_positions;
    SearchStats stats;
    const size_t total_kmers = sequence.size() >= kmer_size ? sequence.size() - kmer_size + 1 : 0;
    stats.total_kmers = total_kmers;
    
    if (total_kmers == 0) return {absent_positions, stats};
    
    // Track which positions we've already determined (found/not found)
    std::vector<int> position_status(total_kmers, -1); // -1 = unknown, 0 = absent, 1 = present
    
    // Step 1: Sample non-overlapping k-mers
    std::vector<size_t> sample_positions;
    for (size_t i = 0; i < total_kmers; i += kmer_size)
    {
        sample_positions.push_back(i);
    }
    
    std::cout << "    Sampling " << sample_positions.size() << " non-overlapping k-mers (vs " << total_kmers << " total)...\n";
    
    // Step 2: Search sampled k-mers and use cursors to extend matches efficiently
    for (size_t pos : sample_positions)
    {
        if (position_status[pos] != -1) continue; // Already processed
        
        // Extract the k-mer at this position
        auto kmer = seqan3::dna5_vector{sequence.begin() + pos, sequence.begin() + pos + kmer_size};
        
        // Convert k-mer to string for caching
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
            // Use cursor to search for the k-mer
            auto cursor = index.cursor();
            bool cursor_success = cursor.extend_right(kmer);
            found = cursor_success && cursor.count() > 0;
            stats.unique_kmers_searched++;
            
            // Cache the result
            global_cache.insert(kmer_str, found);
            
            if (found)
            {
                // Use cursor to efficiently extend to the right to find the full matching region
                size_t extend_pos = pos + kmer_size;
                while (extend_pos < sequence.size() && extend_pos < pos + 100) // Limit extension
                {
                    // Try to extend the cursor by one more character
                    auto next_char = sequence[extend_pos];
                    if (cursor.extend_right(next_char) && cursor.count() > 0)
                    {
                        // Successfully extended - mark this k-mer position as found
                        if (extend_pos - kmer_size + 1 < total_kmers)
                        {
                            position_status[extend_pos - kmer_size + 1] = 1; // present
                        }
                        extend_pos++;
                    }
                    else
                    {
                        // Can't extend further
                        break;
                    }
                }
                
                // Mark all positions covered by the initial k-mer as found
                position_status[pos] = 1; // present
            }
        }
        
        // Mark the current position based on the result
        if (!found)
        {
            position_status[pos] = 0; // absent
            absent_positions.insert(pos);
        }
    }
    
    // Step 3: Fill in gaps - check positions that are still unknown
    size_t gaps_filled = 0;
    for (size_t i = 0; i < total_kmers; ++i)
    {
        if (position_status[i] == -1)
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
                // Use cursor for gap filling
                auto cursor = index.cursor();
                bool cursor_success = cursor.extend_right(kmer);
                found = cursor_success && cursor.count() > 0;
                stats.unique_kmers_searched++;
                
                // Cache the result
                global_cache.insert(kmer_str, found);
            }
            
            if (!found)
            {
                absent_positions.insert(i);
            }
            
            gaps_filled++;
        }
    }
    
    if (gaps_filled > 0)
    {
        std::cout << "    Filled " << gaps_filled << " gap positions using cursors.\n";
    }
    
    // Update progress appropriately
    if (progress)
    {
        progress->update(total_kmers);
    }
    
    return {absent_positions, stats};
}