# @file install_library.cpp
# @auhtor Takumi Sugiura
# @date Dec 2, 2022
# Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
# This software is released under tthickness MIT License.
# http:opensource.org/licenses/mit-license.php
#
# cmake cache file for MPI Compiler use
# how to use:
# cmake -B <build> -S <src> -C <this_file>
#

# compiler
set(CACHE_COMPILER "MPI" CACHE STRING "name" FORCE)
set(CACHE_TYPE "debug" CACHE STRING "type" FORCE)

# set compiler
set(CMAKE_C_COMPILER "mpicc" CACHE STRING "mpi compiler" FORCE)
set(CMAKE_CXX_COMPILER "mpicxx" CACHE STRING "mpi compiler" FORCE)

# set install library as mpi mode
set(CACHE_MPI ON CACHE BOOL "mpi flag" FORCE)

# set compiler flags
set(CACHE_FLAGS -DICARAT_MPI -Wall -Wextra -g CACHE STRING "compile debug flags" FORCE)

# check for main CMakeLists.txt
set(CACHE_IS_INCLUDED YES CACHE STRING "include flag" FORCE)