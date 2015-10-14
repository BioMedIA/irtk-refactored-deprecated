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

set(IRTK_VERSION_MAJOR "2")
set(IRTK_VERSION_MINOR "0")
set(IRTK_VERSION_PATCH "0")
set(IRTK_VERSION "${IRTK_VERSION_MAJOR}.${IRTK_VERSION_MINOR}")
set(IRTK_VERSION_FULL "${IRTK_VERSION}.${IRTK_VERSION_PATCH}")

# SOVERSION = 0 means no stable ABI.
set(IRTK_SOVERSION 0)

#execute_process(
#  COMMAND git log -1 --format=%h
#  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
#  OUTPUT_VARIABLE IRTK_REVISION
#  OUTPUT_STRIP_TRAILING_WHITESPACE
#)

if(NOT IRTK_REVISION)
  set(IRTK_REVISION "NULL")
endif()