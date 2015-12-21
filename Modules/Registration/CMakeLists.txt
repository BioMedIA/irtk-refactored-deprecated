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

set(IRTK_MODULE_NAME "Registration")
set(IRTK_MODULE_INSTALL_INCLUDEDIR
  "${IRTK_INSTALL_INCLUDEDIR}/irtk/${IRTK_MODULE_NAME}")

file(GLOB IRTK_MODULE_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h")

install(FILES ${IRTK_MODULE_HEADERS}
        DESTINATION "${IRTK_MODULE_INSTALL_INCLUDEDIR}"
        COMPONENT dev)

include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/include")

add_subdirectory(src)