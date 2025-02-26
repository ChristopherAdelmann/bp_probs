cmake_minimum_required(VERSION 3.24)

project(base_pairing_probabilities VERSION 0.0.1)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RELEASE)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-psabi -Wno-interference-size")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
# set(CMAKE_CXX_FLAGS_RELEASE "-O3")


# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------
# # # CMake Modules
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# # # ViennaRNA
find_package(ViennaRNA )

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable (bp_probs  ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)

# # # ViennaRNA
message(STATUS "Finding ViennaRNA status: ${VIENNA_RNA_FOUND}")
if (NOT VIENNA_RNA_FOUND)
    message(STATUS "Building static ViennaRNA from source")
    include(viennaRNA)
    add_dependencies(bp_probs viennaRNA)
else()
    message(STATUS "Using ViennaRNA from system")
endif()

# # # Dependency Include Directories
target_include_directories (bp_probs SYSTEM PRIVATE   ${VIENNA_RNA_INCLUDE_DIR} ${VIENNA_RNA_INCLUDE_DIR}/ViennaRNA)

# # # Link Libraries
target_link_libraries (bp_probs PRIVATE   ${VIENNA_RNA_LIBRARY})
