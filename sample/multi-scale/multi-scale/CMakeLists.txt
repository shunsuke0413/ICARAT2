cmake_minimum_required(VERSION 3.22)
project(sample)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../)
find_package(icarat REQUIRED PATHS ../../../lib/ICARAT/lib/${CACHE_COMPILER}/${CACHE_TYPE})
get_filename_component(BIN_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
# dont change above

##### add source file(s) #####
set(SOURCES multiscale.cpp)
##############################

# dont change below
add_executable(${BIN_NAME}.out ${SOURCES})
target_link_libraries(${BIN_NAME}.out icarat::icarat)