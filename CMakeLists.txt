# @file install_library.cpp
# @auhtor Takumi Sugiura
# @date Dec 3, 2022
# Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
# This software is released under tthickness MIT License.
# http:opensource.org/licenses/mit-license.php

# minimim version is set as 3.20 because intel oneapi supports from cmake 3.22
cmake_minimum_required(VERSION 3.22)

# set project name
project(ICARAT)

set(CMAKE_CXX_STANDARD 20) # specity using C++20
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# set default install prefix
set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/lib/ICARAT)

# include functions for this cmake
include(${PROJECT_SOURCE_DIR}/cmake/core/functions.cmake)

# source, include and header file(s)
set(SOURCES "") # source to make library
set(HEADERS "") # source to include this library

# set subdirectories for this project
add_subdirectory(src/analysis)
add_subdirectory(src/basis)
add_subdirectory(src/io)
add_subdirectory(src/misc)
add_subdirectory(src/optimization)
add_subdirectory(src/problem)

# set source files for this library
add_library(icarat SHARED ${SOURCES})

# add_library(icarat STATIC ${SOURCES})

# prepare for dependency libraries
include(${PROJECT_SOURCE_DIR}/cmake/core/install_library.cmake)

if(CACHE_MPI)
    include(${PROJECT_SOURCE_DIR}/cmake/core/install_mpi_library.cmake)
endif()

# add compiler options
target_compile_features(icarat PUBLIC cxx_std_14)
target_compile_options(icarat PUBLIC ${CACHE_FLAGS})
target_link_options(icarat PUBLIC ${CACHE_LDFLAGS})

# include directories to build this target
target_include_directories(icarat PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_include_directories(icarat INTERFACE ${CMAKE_INSTALL_PREFIX}/../headeronly)

# install step1: copy this library and dependent libraries
install(TARGETS icarat
    EXPORT icarat # export name
    LIBRARY DESTINATION lib/${CACHE_COMPILER}/${CACHE_TYPE} # destination of library
    INCLUDES DESTINATION include # destination for include files
)

# install step2: create config file for find_library()
install(EXPORT icarat
    FILE icarat-config.cmake
    NAMESPACE icarat::
    DESTINATION lib/${CACHE_COMPILER}/${CACHE_TYPE}
)

# install step3: copy header file of this project
set(SRC_DIR src/analysis src/basis src/io src/misc src/optimization src/problem)
install(DIRECTORY ${SRC_DIR}
    DESTINATION ./include
    PATTERN *.cpp EXCLUDE
    PATTERN CMakeLists.txt EXCLUDE
    PATTERN Makefile EXCLUDE
)
