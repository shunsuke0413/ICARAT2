# @file install_library.cpp
# @auhtor Takumi Sugiura
# @date Dec 3, 2022
# Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
# This software is released under tthickness MIT License.
# http:opensource.org/licenses/mit-license.php
#
# cmake subdirectory file
#

project(ICARAT)

list(APPEND LOCAL_CXX
solver/CG.cpp
solver/solvers.cpp
)

list(APPEND LOCAL_INC
solver/solvers.hpp
analysis_utility.hpp
arc_length.hpp
explicit_dyna.hpp
general_eigen.hpp
homogenization.hpp
implicit_dyna.hpp
linear.hpp
)



convert_filenames_to_full_paths(LOCAL_CXX)
convert_filenames_to_full_paths(LOCAL_INCLUDE)

set(SOURCES ${SOURCES} ${LOCAL_CXX} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${LOCAL_INC} PARENT_SCOPE)