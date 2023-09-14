cmake_minimum_required(VERSION 3.10)

# Set the project name
project(ExternalArchive)

# Add C++ 20 support
set(CMAKE_CXX_STANDARD 20)

# Specify the directories where the source and header files are located
set(SOURCE_DIR "${CMAKE_SOURCE_DIR}/src")
set(INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")

# Include directories (add the path where external_archive.h is located)
include_directories(${INCLUDE_DIR})

# Add the executable with your source files
add_executable(ExternalArchiveTest ${SOURCE_DIR}/external_archive.cpp ${SOURCE_DIR}/external_archive_test.cpp)

# Add the executable for pso test with its dependencies
add_executable(PSOTest 
    ${SOURCE_DIR}/pso_test.cpp 
    ${SOURCE_DIR}/pso.cpp 
    ${SOURCE_DIR}/particle.cpp 
    ${SOURCE_DIR}/external_archive.cpp
)

# Link necessary libraries if any
