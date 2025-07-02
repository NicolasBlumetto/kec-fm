#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <set>
#include <fstream>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>

// For serialization
#include <cereal/archives/binary.hpp>

// Struct to hold command line arguments
struct cmd_arguments
{
    std::filesystem::path input_dir{};
    std::filesystem::path index_path{"out.fm_index"};
};

// Function to set up the argument parser
void initialize_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "SeqAn3 Expert Team";
    parser.info.version = "1.0.0";
    parser.info.short_description = "Creates an FM-Index from a directory of DNA files.";
    parser.info.synopsis = {"kecfm-index <INPUT_DIR>"};

    parser.add_positional_option(args.input_dir,
                                 "Path to a directory containing reference DNA files (e.g., FASTA).",
                                 seqan3::input_directory_validator{});
    parser.add_option(args.index_path, 'o', "output-file",
                      "Path to save the generated FM-Index.",
                      seqan3::option_spec::required,
                      seqan3::output_file_validator{});
}

// Function to find sequence files in a directory
std::vector<std::filesystem::path> find_sequence_files(std::filesystem::path const & dir_path)
{
    std::vector<std::filesystem::path> file_paths;
    std::set<std::string> valid_extensions{".fa", ".fasta", ".fna", ".ffn", ".faa", ".frn"};

    for (auto const & entry : std::filesystem::recursive_directory_iterator{dir_path})
    {
        if (entry.is_regular_file() && valid_extensions.contains(entry.path().extension().string()))
        {
            file_paths.push_back(entry.path());
        }
    }
    return file_paths;
}

int main(int argc, char ** argv)
{
    seqan3::argument_parser parser{"kecfm-index", argc, argv};
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

    // 1. Find all sequence files in the input directory
    std::vector<std::filesystem::path> file_paths;
    try
    {
        file_paths = find_sequence_files(args.input_dir);
        if (file_paths.empty())
        {
            std::cerr << " No sequence files found in " << args.input_dir << '\n';
            return 0;
        }
    }
    catch (std::runtime_error const & e)
    {
        std::cerr << " " << e.what() << '\n';
        return -1;
    }

    // 2. Read sequences from files and aggregate them
    std::vector<seqan3::dna5_vector> reference_sequences;
    std::cout << "Reading " << file_paths.size() << " reference files...\n";
    for (auto const & path : file_paths)
    {
        seqan3::sequence_file_input fin{path};
        for (auto & [seq, id, qual] : fin)
        {
            reference_sequences.push_back(std::move(seq));
        }
    }
    std::cout << "Finished reading. Total sequences loaded: " << reference_sequences.size() << "\n";

    // 3. Construct the FM-Index
    std::cout << "Constructing FM-Index...\n";
    seqan3::fm_index index{reference_sequences};
    std::cout << "Index construction complete.\n";

    // 4. Serialize the index to disk
    std::cout << "Serializing index to: " << args.index_path.string() << "\n";
    try
    {
        std::ofstream os{args.index_path, std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(index);
    }
    catch (std::exception const & e)
    {
        std::cerr << " Could not write index to file: " << e.what() << '\n';
        return -1;
    }
    std::cout << "Index successfully saved.\n";

    return 0;
}