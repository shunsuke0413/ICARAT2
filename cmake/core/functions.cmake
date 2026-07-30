# @file install_library.cpp
# @auhtor Takumi Sugiura
# @date Dec 3, 2022
# Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
# This software is released under tthickness MIT License.
# http:opensource.org/licenses/mit-license.php
#
# a function which converts relative path to absolute path
#

function(convert_filenames_to_full_paths NAMES)
  unset(tmp_names)

  foreach(name ${${NAMES}})
    list(APPEND tmp_names ${CMAKE_CURRENT_SOURCE_DIR}/${name})
  endforeach()

  set(${NAMES} ${tmp_names} PARENT_SCOPE)
endfunction()