# this test won't works well....
cmake_minimum_required(VERSION 3.22)
project(sample)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../)
find_package(icarat REQUIRED PATHS ../../lib/ICARAT/lib/${CACHE_COMPILER}/${CACHE_TYPE})
find_package(GTest REQUIRED PATHS ../../lib/googletest/lib/cmake)
get_filename_component(BIN_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)

# dont change above

# #### add source file(s) #####
set(SOURCES
    test-basis/make_element.cpp
    test-basis/test_bmatrix.cpp
    test-basis/test_nmatrix_surface.cpp
    test-io/test_mesh_bc.cpp
    test-optimization/test_filter.cpp
    test-optimization/test_opt_solver.cpp
)

# #############################

# dont change below
add_executable(${BIN_NAME}.out ${SOURCES})
target_link_libraries(${BIN_NAME}.out icarat::icarat GTest::gtest_main)