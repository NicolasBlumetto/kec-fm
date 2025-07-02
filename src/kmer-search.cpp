#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <cstdint>
#include <fstream>
#include <ranges>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>
#include <seqan3/core/debug_stream.hpp>

// For serialization
#include <cereal/archives/binary.hpp>

// Struct to hold command line arguments
struct cmd_arguments
{
    std::filesystem::path index_path{};
    std::filesystem::path query_path{};
    uint8_t kmer_size{};
};

// Function to set up the argument parser
void initialize_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "SeqAn3 Expert Team";
    parser.info.version = "1.0.0";
    parser.info.short_description = "Searches for k-mers from a query file in an FM-Index.";
    parser.info.synopsis = {"kmer-search -k <K> -i <INDEX_FILE> -q <QUERY_FILE>"};

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

    // 3. Iterate through k-mers, search, and collect absent ones
    std::vector<seqan3::dna5_vector> absent_kmers;
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
            absent_kmers.push_back(kmer);
        }
    }

    // 4. Report the results
    if (absent_kmers.empty())
    {
        std::cout << "\nAll k-mers from the query were found in the index.\n";
    }
    else
    {
        std::cout << "\nFound " << absent_kmers.size() << " k-mers that are absent from the index:\n";
        for (auto const & kmer : absent_kmers)
        {
            seqan3::debug_stream << kmer << '\n';
        }
    }

    return 0;
}