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

// Struct to hold command line arguments
struct cmd_arguments
{
    std::filesystem::path index_path{};
    std::filesystem::path query_path{};
    std::filesystem::path output_path{"unique_regions.fasta"};
    uint8_t kmer_size{};
};

// Function to set up the argument parser
void initialize_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "SeqAn3 Expert Team";
    parser.info.version = "1.0.0";
    parser.info.short_description = "Identifies unique regions from a query file using an FM-Index.";
    parser.info.synopsis = {"kmer-search -k <K> -i <INDEX_FILE> -q <QUERY_FILE> [-o <OUTPUT_FILE>]"};

    parser.add_option(args.index_path, 'i', "index-file",
                      "Path to the pre-built FM-Index file.",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{});
    parser.add_option(args.query_path, 'q', "query-file",
                      "Path to the query file (e.g., FASTA).",
                      seqan3::option_spec::required,
                      seqan3::input_file_validator{});
    parser.add_option(args.kmer_size, 'k', "kmer-size",
                      "The length of the k-mers to search for.",
                      seqan3::option_spec::required);
    parser.add_option(args.output_path, 'o', "output-file",
                      "Path to the output FASTA file for unique regions.",
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

    // 1. Load (deserialize) the FM-Index
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
    std::cout << "Index loaded successfully.\n";

    // 2. Read the query sequence
    seqan3::dna5_vector query_sequence;
    try
    {
        seqan3::sequence_file_input fin{args.query_path};
        // We only process the first sequence in the query file
        if (fin.begin() == fin.end())
        {
             std::cerr << " Query file is empty.\n";
             return -1;
        }
        auto & [seq, id, qual] = *fin.begin();
        query_sequence = std::move(seq);
    }
    catch (std::exception const & e)
    {
        std::cerr << " Failed to read query file: " << e.what() << '\n';
        return -1;
    }

    if (query_sequence.size() < args.kmer_size)
    {
        std::cerr << " Query sequence is shorter than k-mer size.\n";
        return -1;
    }

    // 3. Iterate through k-mers and collect positions of absent ones
    std::set<size_t> absent_positions;
    std::cout << "Searching for " << query_sequence.size() - args.kmer_size + 1
              << " k-mers of size " << static_cast<int>(args.kmer_size) << "...\n";

    for (size_t i = 0; i <= query_sequence.size() - args.kmer_size; ++i)
    {
        // Use iterator-based construction instead of substr
        auto kmer = seqan3::dna5_vector{query_sequence.begin() + i, 
                                        query_sequence.begin() + i + args.kmer_size};
        auto hits = seqan3::search(kmer, index);
        if (hits.begin() == hits.end())
        {
            absent_positions.insert(i);
        }
    }

    // 4. Merge overlapping positions into unique regions
    auto unique_regions = merge_overlapping_positions(absent_positions, args.kmer_size, query_sequence);

    // 5. Output results
    if (unique_regions.empty())
    {
        std::cout << "\nAll k-mers from the query were found in the index. No unique regions identified.\n";
    }
    else
    {
        std::cout << "\nFound " << absent_positions.size() << " absent k-mers merged into " 
                  << unique_regions.size() << " unique regions.\n";
        std::cout << "Writing unique regions to: " << args.output_path.string() << "\n";

        // Read the original sequence ID for output
        std::string sequence_id = "query";
        try
        {
            seqan3::sequence_file_input fin{args.query_path};
            auto & [seq, id, qual] = *fin.begin();
            if (!id.empty())
            {
                sequence_id = std::string(id);
            }
        }
        catch (...)
        {
            // Use default ID if reading fails
        }

        // Write unique regions to FASTA file
        try
        {
            seqan3::sequence_file_output fout{args.output_path};
            
            for (size_t i = 0; i < unique_regions.size(); ++i)
            {
                const auto& region = unique_regions[i];
                std::string region_id = sequence_id + "_unique_region_" + std::to_string(i + 1) +
                                       "_pos_" + std::to_string(region.start + 1) + "-" + 
                                       std::to_string(region.end + 1) +
                                       "_len_" + std::to_string(region.sequence.size());
                
                fout.emplace_back(region.sequence, region_id);
            }
        }
        catch (std::exception const & e)
        {
            std::cerr << " Failed to write output file: " << e.what() << '\n';
            return -1;
        }
        
        std::cout << "Successfully wrote " << unique_regions.size() << " unique regions to output file.\n";
        
        // Print summary statistics
        size_t total_unique_length = 0;
        size_t min_length = std::numeric_limits<size_t>::max();
        size_t max_length = 0;
        
        for (const auto& region : unique_regions)
        {
            size_t length = region.sequence.size();
            total_unique_length += length;
            min_length = std::min(min_length, length);
            max_length = std::max(max_length, length);
        }
        
        std::cout << "\nSummary statistics:\n";
        std::cout << "  Total unique sequence length: " << total_unique_length << " bp\n";
        std::cout << "  Average region length: " << (total_unique_length / unique_regions.size()) << " bp\n";
        std::cout << "  Minimum region length: " << min_length << " bp\n";
        std::cout << "  Maximum region length: " << max_length << " bp\n";
    }

    return 0;
}