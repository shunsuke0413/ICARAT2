
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
elasticity.cpp
elasto_plast_iso.cpp
heat.cpp
fluid.cpp
mooney.cpp
multipli_elasto_plast_iso.cpp
neo_hooke.cpp
saint_venant.cpp
)

list(APPEND LOCAL_INC
base.hpp
elasticity.hpp
elasto_plast_iso.hpp
heat.hpp
fluid.hpp
generalized_maxwell.hpp
mooney.hpp
multipli_elasto_plast_iso.hpp
neo_hooke.hpp
saint_venant.hpp
)



convert_filenames_to_full_paths(LOCAL_CXX)
convert_filenames_to_full_paths(LOCAL_INCLUDE)

set(SOURCES ${SOURCES} ${LOCAL_CXX} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${LOCAL_INC} PARENT_SCOPE)