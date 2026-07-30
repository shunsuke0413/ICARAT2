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
boundary_condition_base.cpp
boundary_condition_gid.cpp
)

list(APPEND LOCAL_INC
boundary_condition_base.hpp
boundary_condition_gid.hpp
boundary_condition_icarat.hpp
export_graph.hpp
export_vtk.hpp
file_reader.hpp
mesh_base.hpp
mesh_icarat.hpp
mesh_mdpa.hpp
mesh_vtk.hpp
)



convert_filenames_to_full_paths(LOCAL_CXX)
convert_filenames_to_full_paths(LOCAL_INCLUDE)

set(SOURCES ${SOURCES} ${LOCAL_CXX} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${LOCAL_INC} PARENT_SCOPE)