cmake_minimum_required(VERSION 3.10)

# Set the project name
project(ExternalArchive)

# Add C++ 20 support
set(CMAKE_CXX_STANDARD 20)

find_package(spdlog REQUIRED)
find_package(crossguid CONFIG REQUIRED)
find_package(Arrow REQUIRED)
find_package(Parquet REQUIRED)

# Specify the directories where the source and header files are located
set(SOURCE_DIR "${CMAKE_SOURCE_DIR}/src")
set(INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")

set(LIB_REPORT_LOADS_HEADERS
        include/json.hpp
        include/csv.hpp
        )

# Include directories (add the path where external_archive.h is located)
include_directories(${INCLUDE_DIR})

# Add the executable with your source files
add_executable(ExternalArchiveTest ${SOURCE_DIR}/external_archive.cpp ${SOURCE_DIR}/external_archive_test.cpp ${SOURCE_DIR}/particle.cpp)
add_executable(Execute ${SOURCE_DIR}/exec_test.cpp ${SOURCE_DIR}/execute.cpp ${SOURCE_DIR}/misc_utilities.cpp)


# Add the executable for pso test with its dependencies
add_executable(PSOTest 
    ${SOURCE_DIR}/pso_test.cpp 
    ${SOURCE_DIR}/pso.cpp 
    ${SOURCE_DIR}/particle.cpp 
    ${SOURCE_DIR}/scenario.cpp 
    ${SOURCE_DIR}/misc_utilities.cpp
    ${SOURCE_DIR}/external_archive.cpp
    ${SOURCE_DIR}/amqp.cpp
    ${SOURCE_DIR}/execute.cpp
)

add_executable(ipopt_json_test
    ${SOURCE_DIR}/ipopt_json.cpp
    ${SOURCE_DIR}/pso.cpp 
    ${SOURCE_DIR}/particle.cpp 
    ${SOURCE_DIR}/scenario.cpp 
    ${SOURCE_DIR}/misc_utilities.cpp
    ${SOURCE_DIR}/external_archive.cpp
    ${SOURCE_DIR}/amqp.cpp
    ${SOURCE_DIR}/execute.cpp
)


add_executable(PSOCast
    ${SOURCE_DIR}/pso_cast.cpp 
    ${SOURCE_DIR}/scenario.cpp 
    ${SOURCE_DIR}/pso.cpp 
    ${SOURCE_DIR}/particle.cpp 
    ${SOURCE_DIR}/misc_utilities.cpp
    ${SOURCE_DIR}/external_archive.cpp
    ${SOURCE_DIR}/amqp.cpp
    ${SOURCE_DIR}/execute.cpp
)
add_executable(ScenarioTest
    ${SOURCE_DIR}/scenario_test.cpp 
    ${SOURCE_DIR}/scenario.cpp 
    ${SOURCE_DIR}/pso.cpp 
    ${SOURCE_DIR}/particle.cpp 
    ${SOURCE_DIR}/misc_utilities.cpp
    ${SOURCE_DIR}/external_archive.cpp
    ${SOURCE_DIR}/amqp.cpp
    ${SOURCE_DIR}/execute.cpp
)
target_link_libraries(ipopt_json_test PRIVATE arrow parquet fmt pthread crossguid hiredis redis++ SimpleAmqpClient)
target_link_libraries(PSOTest PRIVATE arrow parquet fmt pthread crossguid hiredis redis++ SimpleAmqpClient)
target_link_libraries(PSOCast PRIVATE arrow parquet fmt pthread crossguid hiredis redis++ SimpleAmqpClient)
target_link_libraries(ScenarioTest PRIVATE arrow parquet fmt pthread crossguid hiredis redis++ SimpleAmqpClient) 

target_link_libraries(ExternalArchiveTest PRIVATE fmt crossguid pthread) 
target_link_libraries(Execute PRIVATE fmt crossguid pthread) 
# Link necessary libraries if any
