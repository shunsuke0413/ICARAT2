# @file install_library.cpp
# @auhtor Takumi Sugiura
# @date Dec 2, 2022
# Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
# This software is released under tthickness MIT License.
# http:opensource.org/licenses/mit-license.php
#
# cmake cache file for Intel Compiler use
# how to use:
# cmake -B <build> -S <src> -C <this_file>
#

# compiler
set(CACHE_COMPILER "Intel" CACHE STRING "name" FORCE)
set(CACHE_TYPE "debug" CACHE STRING "type" FORCE)

# set compiler
set(CMAKE_C_COMPILER "icc" CACHE STRING "intel compiler" FORCE)
set(CMAKE_CXX_COMPILER "icpc" CACHE STRING "intel compiler" FORCE)

# set install library as serial mode
set(CACHE_MPI OFF CACHE BOOL "mpi flag" FORCE)

# set compiler flags
set(CACHE_FLAGS -g -qmkl -DEIGEN_USE_MKL_ALL -DEIGEN_INITIALIZE_MATRICES_BY_ZERO CACHE STRING "flag" FORCE)

# set link flags
set(CACHE_LDFLAGS -qmkl CACHE STRING "link option" FORCE)

# check for main CMakeLists.txt
set(CACHE_IS_INCLUDED YES CACHE STRING "include flag" FORCE)