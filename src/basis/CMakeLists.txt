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
bmatrix.cpp
bmatrix_shell.cpp
bmatrix_total.cpp
bmatrix_updated.cpp
gauss_integration.cpp
nmatrix_surface.cpp
shape_function.cpp
)

list(APPEND LOCAL_INC
bmatrix.hpp
bmatrix_shell.hpp
bmatrix_total.hpp
bmatrix_updated.hpp
gauss_integration.hpp
newton_raphson.hpp
nmatrix_surface.hpp
shape_function.hpp

)



convert_filenames_to_full_paths(LOCAL_CXX)
convert_filenames_to_full_paths(LOCAL_INCLUDE)

set(SOURCES ${SOURCES} ${LOCAL_CXX} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${LOCAL_INC} PARENT_SCOPE)