# Copyright (c) 2015, Robert Lukierski.
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
# 
# Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
# 
# Neither the name of [project] nor the names of its
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# ------------------------------------------------------------------------------

cmake_minimum_required(VERSION 2.8)

project(CameraModels)

# ------------------------------------------------------------------------------
# Where to look for find-packages
# ------------------------------------------------------------------------------
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/")

# ------------------------------------------------------------------------------
# C++11
# ------------------------------------------------------------------------------
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++11" COMPILER_SUPPORTS_CXX11)
check_cxx_compiler_flag("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

# ----------------------------------------
# Dependencies
# ----------------------------------------
# Eigen
find_package(Eigen3 REQUIRED QUIET)
include_directories(${EIGEN3_INCLUDE_DIR})

# Sophus
find_package(Sophus REQUIRED QUIET)
include_directories(${Sophus_INCLUDE_DIR})

# ---------------------------------------------
# Our headers
# ---------------------------------------------
include_directories(.)

set(HEADERS
CameraModelHelpers.hpp
CameraModels.hpp
CameraPyramid.hpp
FisheyeCameraModel.hpp
FullGenericCameraModel.hpp
IdealFisheyeCameraModel.hpp
IdealGenericCameraModel.hpp
OpenGLProjectionMatrix.hpp
PinholeCameraModel.hpp
PinholeDisparityBrownConrady.hpp
PinholeDisparityCameraModel.hpp
PinholeDisparityDistortedCameraModel.hpp
PinholeDistortedCameraModel.hpp
SphericalCameraModel.hpp
SphericalPovRayCameraModel.hpp
)

# ----------------------------------------------------------------------
# Sources
# ----------------------------------------------------------------------
set(SOURCES
sources/CameraModels.cpp
)

# ---------------------------------------------
# Build
# ---------------------------------------------
add_library(camera_models SHARED ${SOURCES} ${HEADERS})

# --------------------------------------------------------------------------------
# Install
# --------------------------------------------------------------------------------
install(FILES 
${HEADERS}
DESTINATION include
)

install(TARGETS camera_models
RUNTIME DESTINATION bin
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib
)

set(CAMERA_MODELS_FOUND TRUE)
set(CAMERA_MODELS_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} ${EIGEN3_INCLUDE_DIR} )

configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/CameraModelsConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/CameraModelsConfig.cmake @ONLY IMMEDIATE)
export( PACKAGE CameraModels )

# ---------------------------------------------
# Unit Tests
# ---------------------------------------------
if(DEFINED BUILD_TESTS)
    if(BUILD_TESTS)
        add_subdirectory(tests)
    endif()
else()
    add_subdirectory(tests)
endif()
