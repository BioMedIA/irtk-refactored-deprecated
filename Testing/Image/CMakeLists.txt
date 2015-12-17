# The Image Registration Toolkit (IRTK)
#
# Copyright 2008-2015 Imperial College London
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# ------------------------------------------------------------------------------
# Keep test executables separate from actual programs
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/Testing/bin")
set (EXECUTABLE_OUTPUT_PATH         "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
set (INPUT_DIR                      "${CMAKE_CURRENT_SOURCE_DIR}")

# ------------------------------------------------------------------------------
# Add any test source and required input arguments below

# Test names
set(TESTS
  irtkUnaryVoxelFunctionTest
  irtkConvolutionFunctionTest
  irtkDownsamplingTest
)

# Test arguments
#
# For each test which requires command-line arguments, set a CMake variable here
# named <test_source_name>_ARGS to the list of arguments that should be passed
# on to the test.
set(irtkUnaryVoxelFunctionTest_ARGS)

# ------------------------------------------------------------------------------
# Usually nothing has to be edited below to add a new test
if(BUILD_GTESTS)
  foreach(test IN LISTS TESTS)
    get_filename_component(name "${test}" NAME_WE)
    get_filename_component(ext  "${test}" EXT)
    if("${ext}" STREQUAL "")
      set(ext ".cc")
    endif()
    add_executable(${name} ${name}${ext})
    target_link_libraries(${name} ${GTEST_LIBRARIES})
    add_test(${name} "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${name}" ${${name}_ARGS})
  endforeach()
endif()

# ------------------------------------------------------------------------------
# CMake test driver (i.e., not based on GTest - DEPRECATED)
create_test_sourcelist(TEST_DRIVER_SRCS
  irtkImageTestDriver.cc
    irtkInterpolateExtrapolateImageFunctionTest.cc
    irtkDisplacementToVelocityFieldTest.cc
)

irtk_add_executable(irtkImageTestDriver ${TEST_DRIVER_SRCS})
macro(add_deprecated_test name)
  add_test(${name} "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/irtkImageTestDriver" ${name} ${ARGN})
endmacro()

add_deprecated_test(irtkInterpolateExtrapolateImageFunctionTest)
add_deprecated_test(irtkDisplacementToVelocityFieldTest)
