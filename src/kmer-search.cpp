#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <cstdint>
#include <fstream>
#include <ranges>
#include <algorithm>
#include <set>
#include <limits>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>
#include <seqan3/core/debug_stream.hpp>

// For serialization
#include <cereal/archives/binary.hpp>

// Structure to represent a unique region
struct UniqueRegion
{
    size_t start;
    size_t end;
    seqan3::dna5_vector sequence;
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
    bool success;
    std::vector<std::pair<UniqueRegion, std::string>> unique_regions_with_ids; // region and full sequence ID
};

ProcessingResult process_sequence_file(const std::filesystem::path& sequence_file,
                                      const auto& index,
                                      uint8_t kmer_size)
{
    ProcessingResult result;
    result.filename = sequence_file.filename().string();
    result.total_sequences = 0;
    result.total_unique_regions = 0;
    result.total_unique_length = 0;
    result.success = false;

    try
    {
        seqan3::sequence_file_input fin{sequence_file};
        
        for (auto & [seq, id, qual] : fin)
        {
            result.total_sequences++;
            
            if (seq.size() < kmer_size)
            {
                std::cout << "  Warning: Sequence '" << id << "' is shorter than k-mer size, skipping.\n";
                continue;
            }
            
            // Find absent k-mer positions
            std::set<size_t> absent_positions;
            for (size_t i = 0; i <= seq.size() - kmer_size; ++i)
            {
                auto kmer = seqan3::dna5_vector{seq.begin() + i, seq.begin() + i + kmer_size};
                auto hits = seqan3::search(kmer, index);
                if (hits.begin() == hits.end())
                {
                    absent_positions.insert(i);
                }
            }
            
            // Merge overlapping positions into unique regions
            auto unique_regions = merge_overlapping_positions(absent_positions, kmer_size, seq);
            
            // Store unique regions with their IDs for later output
            for (size_t i = 0; i < unique_regions.size(); ++i)
            {
                const auto& region = unique_regions[i];
                
                // Create a comprehensive ID that includes file info
                std::string file_stem = sequence_file.stem().string();
                std::string region_id = file_stem + "_" + std::string(id) + "_unique_region_" + std::to_string(i + 1) +
                                       "_pos_" + std::to_string(region.start + 1) + "-" + 
                                       std::to_string(region.end + 1) +
                                       "_len_" + std::to_string(region.sequence.size());
                
                result.unique_regions_with_ids.emplace_back(region, region_id);
                result.total_unique_length += region.sequence.size();
            }
            
            result.total_unique_regions += unique_regions.size();
        }
        
        result.success = true;
    }
    catch (std::exception const & e)
    {
        std::cerr << "  Error processing " << result.filename << ": " << e.what() << '\n';
        result.success = false;
    }
    
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
    parser.info.synopsis = {"kmer-search -k <K> -i <INDEX_FILE> -q <QUERY_FILE_OR_DIR> [-o <OUTPUT_FILE>]"};

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
    seqan3::argument_parser parser{"kmer-search", argc, argv};
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

    // Process each query file
    std::vector<ProcessingResult> results;
    std::vector<std::pair<UniqueRegion, std::string>> all_unique_regions;
    
    for (const auto& query_file : query_files)
    {
        std::cout << "Processing: " << query_file.filename().string() << "\n";
        
        auto result = process_sequence_file(query_file, index, args.kmer_size);
        results.push_back(result);
        
        if (result.success)
        {
            std::cout << "  Processed " << result.total_sequences << " sequences\n";
            std::cout << "  Found " << result.total_unique_regions << " unique regions\n";
            std::cout << "  Total unique sequence length: " << result.total_unique_length << " bp\n";
            
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

    // Print summary statistics
    std::cout << "=== SUMMARY ===\n";
    size_t total_files = results.size();
    size_t successful_files = 0;
    size_t total_sequences = 0;
    size_t total_unique_regions = 0;
    size_t total_unique_length = 0;
    
    for (const auto& result : results)
    {
        if (result.success)
        {
            successful_files++;
            total_sequences += result.total_sequences;
            total_unique_regions += result.total_unique_regions;
            total_unique_length += result.total_unique_length;
        }
    }
    
    std::cout << "Files processed: " << successful_files << "/" << total_files << "\n";
    std::cout << "Total sequences: " << total_sequences << "\n";
    std::cout << "Total unique regions: " << total_unique_regions << "\n";
    std::cout << "Total unique sequence length: " << total_unique_length << " bp\n";
    
    if (total_unique_regions > 0)
    {
        std::cout << "Average region length: " << (total_unique_length / total_unique_regions) << " bp\n";
    }
    
    std::cout << "Output file: " << args.output_file << "\n";

    return successful_files == total_files ? 0 : 1;
}