# @file install_library.cpp
# @auhtor Takumi Sugiura
# @date Dec 3, 2022
# Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
# This software is released under tthickness MIT License.
# http:opensource.org/licenses/mit-license.php
#
# cmake file to install MPI dependency libraries
#

include(ExternalProject)

# boost
# TODO: check build runs or not
set(BOOST_BUILD_DIR ${CMAKE_SOURCE_DIR}/lib/.cache/boost/)
set(BOOST_INSTALL_DIR ${PROJECT_SOURCE_DIR}/lib/boost)
set(BOOST_INCLUDE_DIR ${BOOST_INSTALL_DIR})

if(NOT OFFLINE)
        ExternalProject_Add(
                boost
                URL https://boostorg.jfrog.io/artifactory/main/release/1.76.0/source/boost_1_76_0.tar.bz2
                URL_HASH MD5=33334dd7f862e8ac9fe1cc7c6584fb6d
                PREFIX ${BOOST_BUILD_DIR}
                CONFIGURE_COMMAND ""
                UPDATE_COMMAND ""
                BUILD_COMMAND ""
                INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory ${BOOST_BUILD_DIR}/src/boost ${BOOST_INSTALL_DIR}/
                TEST_COMMAND ""
        )
endif()

if(NOT OFFLINE)
        # amgcl
        set(AMGCL_BUILD_DIR ${CMAKE_SOURCE_DIR}/lib/.cache/amgcl/)
        set(AMGCL_INSTALL_DIR ${PROJECT_SOURCE_DIR}/lib/amgcl)
        set(AMGCL_INCLUDE_DIR ${AMGCL_INSTALL_DIR})
        ExternalProject_Add(
                amgcl
                GIT_REPOSITORY https://github.com/ddemidov/amgcl.git
                GIT_TAG 1.4.3
                PREFIX ${AMGCL_BUILD_DIR}
                CONFIGURE_COMMAND ""
                BUILD_COMMAND ""
                UPDATE_COMMAND ""
                INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory ${AMGCL_BUILD_DIR}/src/amgcl/ ${AMGCL_INSTALL_DIR}/
                TEST_COMMAND ""
        )
endif()

if(NOT OFFLINE)
        # metis
        set(METIS_BUILD_DIR ${CMAKE_SOURCE_DIR}/lib/.cache/metis/)
        set(METIS_INSTALL_DIR ${PROJECT_SOURCE_DIR}/lib/metis)
        set(METIS_INCLUDE_DIR ${METIS_INSTALL_DIR}/metis)
        ExternalProject_Add(
                metis
                URL http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
                PREFIX ${METIS_BUILD_DIR}
                CONFIGURE_COMMAND ${CMAKE_COMMAND} -B ${METIS_BUILD_DIR}/src/metis-build -S ${METIS_BUILD_DIR}/src/metis -DGKLIB_PATH=${METIS_BUILD_DIR}/src/metis/GKlib
                BUILD_COMMAND ${CMAKE_COMMAND} --build ${METIS_BUILD_DIR}/src/metis-build -j
                UPDATE_COMMAND ""
                INSTALL_COMMAND ${CMAKE_COMMAND} --install ${METIS_BUILD_DIR}/src/metis-build --prefix ${METIS_INSTALL_DIR}
                TEST_COMMAND ""
        )
endif()

if(NOT OFFLINE)
        # zlib
        set(ZLIB_BUILD_DIR ${CMAKE_SOURCE_DIR}/lib/.cache/zlib/)
        set(ZLIB_INSTALL_DIR ${PROJECT_SOURCE_DIR}/lib/zlib)
        set(ZLIB_INCLUDE_DIR ${ZLIB_INSTALL_DIR}/zlib)
        ExternalProject_Add(
                zlib
                GIT_REPOSITORY https://github.com/madler/zlib.git
                GIT_TAG v1.2.13
                PREFIX ${ZLIB_BUILD_DIR}
                CONFIGURE_COMMAND ${CMAKE_COMMAND} -B ${ZLIB_BUILD_DIR}/src/zlib-build -S ${ZLIB_BUILD_DIR}/src/zlib -DCMAKE_INSTALL_PREFIX=${ZLIB_INSTALL_DIR}
                BUILD_COMMAND ${CMAKE_COMMAND} --build ${ZLIB_BUILD_DIR}/src/zlib-build -j
                UPDATE_COMMAND ""
                INSTALL_COMMAND ${CMAKE_COMMAND} --install ${ZLIB_BUILD_DIR}/src/zlib-build --prefix ${ZLIB_INSTALL_DIR}
                TEST_COMMAND ""
        )
endif()

if(NOT OFFLINE)
        # hypre
        set(HYPRE_BUILD_DIR ${CMAKE_SOURCE_DIR}/lib/.cache/hypre)
        set(HYPRE_INSTALL_DIR ${PROJECT_SOURCE_DIR}/lib/hypre/)
        set(HYPRE_INCLUDE_DIR ${HYPRE_INSTALL_DIR})
        ExternalProject_Add(
                hypre
                GIT_REPOSITORY https://github.com/hypre-space/hypre.git
                GIT_TAG v2.26.0
                PREFIX ${HYPRE_BUILD_DIR}
                CONFIGURE_COMMAND ${CMAKE_COMMAND} -B ${HYPRE_BUILD_DIR}/src/hypre-build/ -S ${HYPRE_BUILD_DIR}/src/hypre/src -DCMAKE_INSTALL_PREFIX=${HYPRE_INSTALL_DIR}/hypre
                UPDATE_COMMAND ""
                INSTALL_COMMAND ${CMAKE_COMMAND} --install ${HYPRE_BUILD_DIR}/src/hypre-build/
        )
endif()

target_include_directories(icarat PRIVATE ${BOOST_INCLUDE_DIR})
target_include_directories(icarat PRIVATE ${AMGCL_INCLUDE_DIR})
target_include_directories(icarat PRIVATE ${METIS_INCLUDE_DIR})
target_include_directories(icarat PRIVATE ${ZLIB_INCLUDE_DIR})
target_include_directories(icarat PRIVATE ${HYPRE_INCLUDE_DIR})

if(NOT OFFLINE)
        add_dependencies(icarat boost amgcl metis zlib hypre)
endif()
