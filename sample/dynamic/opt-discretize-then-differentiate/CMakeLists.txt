cmake_minimum_required(VERSION 3.22)
project(sample)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(icarat REQUIRED PATHS ../../../lib/ICARAT/lib/${CACHE_COMPILER}/${CACHE_TYPE})
get_filename_component(BIN_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
# dont change above

##### add source file(s) #####
set(SOURCES opt_discretize_then_differentiate.cpp)
##############################

# dont change below
add_executable(${BIN_NAME}.out ${SOURCES})
target_link_libraries(${BIN_NAME}.out icarat::icarat)
