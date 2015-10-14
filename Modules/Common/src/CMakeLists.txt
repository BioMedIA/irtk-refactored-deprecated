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

set(IRTK_MODULE_OUTPUT_NAME "irtk${IRTK_MODULE_NAME}")

set(IRTK_MODULE_SOURCES irtkBSpline.cc
                        irtkCommon.cc
                        irtkCifstream.cc
                        irtkCofstream.cc
                        irtkCxxLib.cc
                        irtkException.cc
                        irtkObserver.cc
                        irtkOptions.cc
                        irtkParallel.cc
                        irtkPath.cc
                        irtkProfiling.cc
                        irtkSinc.cc
                        irtkTerminal.cc
                        irtkZLib.cc)

if(WITH_MATLAB)
  list(APPEND IRTK_MODULE_SOURCES irtkMatlab.cc)
endif()

add_library(${IRTK_MODULE_NAME} ${IRTK_MODULE_SOURCES})

set_target_properties(${IRTK_MODULE_NAME} PROPERTIES
                      VERSION ${IRTK_VERSION}
                      SOVERSION ${IRTK_SOVERSION}
                      OUTPUT_NAME ${IRTK_MODULE_OUTPUT_NAME})

install(TARGETS ${IRTK_MODULE_NAME}
        DESTINATION "${IRTK_INSTALL_LIBDIR}"
        INCLUDES DESTINATION "${IRTK_MODULE_INSTALL_INCLUDEDIR}"
        COMPONENT lib)