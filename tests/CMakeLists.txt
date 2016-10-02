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

# ------------------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------------------
# Google Log
find_package(Glog REQUIRED QUIET)
include_directories(${GLOG_INCLUDE_DIRS})

# ------------------------------------------------------------------------------
# Internal headers for tests
# ------------------------------------------------------------------------------
include_directories(.)

# ----------------------------------------------------------------------
# UT sources
# ----------------------------------------------------------------------
set(TEST_SOURCES
tests_main.cpp
UT_CameraModels.cpp
UT_PolymorphicCameraModels.cpp
UT_CameraPyramid.cpp
)

# --------------------------------------------------
# Build UT executable
# --------------------------------------------------
add_executable(UT_CameraModels ${TEST_SOURCES})

# --------------------------------------------------
# Link external dependencies
# --------------------------------------------------
target_link_libraries(UT_CameraModels ${GLOG_LIBRARIES} ${GTEST_LIBRARIES} camera_models)

# --------------------------------------------------
# Add the test so it can be run with make test
# --------------------------------------------------
add_test(UT_CameraModels UT_CameraModels --gtest_output=xml:UT_CameraModels.xml)
