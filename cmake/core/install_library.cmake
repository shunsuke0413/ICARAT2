# @file install_library.cpp
# @auhtor Takumi Sugiura
# @date Dec 3, 2022
# Copyright (c) 2021 SAODLab-Nagoya. All Rights reserved.
# This software is released under tthickness MIT License.
# http:opensource.org/licenses/mit-license.php
#
# cmake file to install dependency libraries
#

# # install header only libraries
include(ExternalProject)

# eigen
set(EIGEN_BUILD_DIR ${PROJECT_SOURCE_DIR}/lib/.cache/eigen)
set(EIGEN_INSTALL_DIR ${PROJECT_SOURCE_DIR}/lib/headeronly)
set(EIGEN_INCLUDE_DIR ${EIGEN_INSTALL_DIR})

if(NOT OFFLINE)
        ExternalProject_Add(
                eigen
		GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
		GIT_TAG 3.4.0
                PREFIX ${EIGEN_BUILD_DIR}
                CONFIGURE_COMMAND ""
                BUILD_COMMAND ""
                UPDATE_COMMAND ""
                INSTALL_COMMAND
                ${CMAKE_COMMAND} -E copy_directory ${EIGEN_BUILD_DIR}/src/eigen/Eigen ${EIGEN_INCLUDE_DIR}/Eigen
                && ${CMAKE_COMMAND} -E copy_directory ${EIGEN_BUILD_DIR}/src/eigen/unsupported ${EIGEN_INCLUDE_DIR}/unsupported
                TEST_COMMAND ""
        )
endif()

# toml11
# download toml11 library : https://github.com/ToruNiina/toml11.git
set(TOML11_BUILD_DIR ${CMAKE_SOURCE_DIR}/lib/.cache/toml11)
set(TOML11_INSTALL_DIR ${PROJECT_SOURCE_DIR}/lib/headeronly/toml)
set(TOML11_INCLUDE_DIR ${TOML11_INSTALL_DIR})

if(NOT OFFLINE)
        ExternalProject_Add(
                toml11
                GIT_REPOSITORY https://github.com/ToruNiina/toml11.git
                PREFIX ${TOML11_BUILD_DIR}
                CONFIGURE_COMMAND ""
                BUILD_COMMAND ""
                UPDATE_COMMAND ""
                INSTALL_COMMAND
                ${CMAKE_COMMAND} -E copy_directory ${TOML11_BUILD_DIR}/src/toml11/toml ${TOML11_INCLUDE_DIR}/toml
                && ${CMAKE_COMMAND} -E copy ${TOML11_BUILD_DIR}/src/toml11/toml.hpp ${TOML11_INCLUDE_DIR}
                TEST_COMMAND ""
        )
endif()

# spectra
set(SPECTRA_BUILD_DIR ${PROJECT_SOURCE_DIR}/lib/.cache/spectra)
set(SPECTRA_INSTALL_DIR ${PROJECT_SOURCE_DIR}/lib/headeronly/)
set(SPECTRA_INCLUDE_DIR ${SPECTRA_INSTALL_DIR})

if(NOT OFFLINE)
        ExternalProject_Add(
                spectra
                GIT_REPOSITORY https://github.com/yixuan/spectra.git
                GIT_TAG v1.0.1
                PREFIX ${SPECTRA_BUILD_DIR}
                CONFIGURE_COMMAND ""
                BUILD_COMMAND ""
                UPDATE_COMMAND ""
                INSTALL_COMMAND
                ${CMAKE_COMMAND} -E copy_directory ${SPECTRA_BUILD_DIR}/src/spectra/include/ ${SPECTRA_INCLUDE_DIR}
                TEST_COMMAND ""
        )
endif()

# googletest
set(GTEST_BUILD_DIR ${CMAKE_SOURCE_DIR}/lib/.cache/googletest/)
set(GTEST_INSTALL_DIR ${PROJECT_SOURCE_DIR}/lib/googletest)
set(GTEST_INCLUDE_DIR ${GTEST_INSTALL_DIR}/include)

if(NOT OFFLINE)
        ExternalProject_Add(
                googletest
                GIT_REPOSITORY https://github.com/google/googletest.git
                GIT_TAG release-1.12.1
                PREFIX ${GTEST_BUILD_DIR}
                CONFIGURE_COMMAND ${CMAKE_COMMAND} -B ${GTEST_BUILD_DIR}/src/googletest-build -S ${GTEST_BUILD_DIR}/src/googletest
                BUILD_COMMAND ${CMAKE_COMMAND} --build ${GTEST_BUILD_DIR}/src/googletest-build --parallel
                UPDATE_COMMAND ""
                INSTALL_COMMAND ${CMAKE_COMMAND} --install ${GTEST_BUILD_DIR}/src/googletest-build --prefix ${GTEST_INSTALL_DIR}
                TEST_COMMAND ""
        )
endif()

target_include_directories(icarat PRIVATE ${EIGEN_INCLUDE_DIR})
target_include_directories(icarat PRIVATE ${TOML11_INCLUDE_DIR})

if(NOT OFFLINE)
        add_dependencies(icarat eigen toml11 googletest spectra)
endif()
