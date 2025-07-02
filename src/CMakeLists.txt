cmake_minimum_required(VERSION 3.4)
project(kmer-tools LANGUAGES CXX)

# add seqan3 to search path
list (APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../seqan3/build_system")

find_package(seqan3 3.3.0 REQUIRED)

# Add compiler flags to suppress warnings
add_compile_definitions(SEQAN3_DISABLE_NEWER_COMPILER_DIAGNOSTIC)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-interference-size")

# Define the first executable for creating the index.
add_executable(kecfm-index create-index.cpp)

# Define the second executable for searching the index.
add_executable(kecfm-find kmer-search.cpp)

# Link both executables to the SeqAn3 library.
# seqan3::seqan3 is an INTERFACE target that provides the necessary
# include directories, dependencies (like cereal), and compiler flags.
target_link_libraries(kecfm-index PRIVATE seqan3::seqan3)
target_link_libraries(kecfm-find PRIVATE seqan3::seqan3)
