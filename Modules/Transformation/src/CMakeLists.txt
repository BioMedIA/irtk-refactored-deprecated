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

include_directories("${CMAKE_SOURCE_DIR}/Modules/Common/include")
include_directories("${CMAKE_SOURCE_DIR}/Modules/Geometry/include")
include_directories("${CMAKE_SOURCE_DIR}/Modules/Image/include")
include_directories("${CMAKE_SOURCE_DIR}/Modules/Registration/include")

set(IRTK_MODULE_SOURCES irtkAffineTransformation.cc
                        irtkBSplineFreeFormTransformation3D.cc
                        irtkBSplineFreeFormTransformation4D.cc
                        irtkBSplineFreeFormTransformationStatistical.cc
                        irtkBSplineFreeFormTransformationSV.cc
                        irtkBSplineFreeFormTransformationTD.cc
                        irtkEigenFreeFormTransformation.cc
                        irtkFluidFreeFormTransformation.cc
                        irtkFreeFormTransformation3D.cc
                        irtkFreeFormTransformation4D.cc
                        irtkFreeFormTransformation.cc
                        irtkFreeFormTransformationRungeKutta.cc
                        irtkHomogeneousTransformation.cc
                        irtkImageFastFourierTransform.cc
                        irtkImageHomogeneousTransformation.cc
                        irtkImageTransformation2.cc
                        irtkImageTransformation.cc
                        irtkInverseAffineTransformation.cc
                        irtkLinearFreeFormTransformation3D.cc
                        irtkLinearFreeFormTransformation4D.cc
                        irtkLinearFreeFormTransformationTD.cc
                        irtkMultiLevelFreeFormTransformation.cc
                        irtkMultiLevelStationaryVelocityTransformation.cc
                        irtkMultiLevelTransformation.cc
                        irtkPartialAffineTransformation.cc
                        irtkPartialBSplineFreeFormTransformationSV.cc
                        irtkPartialMultiLevelStationaryVelocityTransformation.cc
                        irtkRigidTransformation.cc
                        irtkSimilarityTransformation.cc
                        irtkTransformationApproximationError.cc
                        irtkTransformation.cc
                        irtkTransformationInverse.cc)

add_library(${IRTK_MODULE_NAME} ${IRTK_MODULE_SOURCES})

target_link_libraries(${IRTK_MODULE_NAME} Common Geometry Image Registration)

set_target_properties(${IRTK_MODULE_NAME} PROPERTIES
                      VERSION ${IRTK_VERSION}
                      SOVERSION ${IRTK_SOVERSION}
                      OUTPUT_NAME ${IRTK_MODULE_OUTPUT_NAME})

install(TARGETS ${IRTK_MODULE_NAME}
        DESTINATION "${IRTK_INSTALL_LIBDIR}"
        INCLUDES DESTINATION "${IRTK_MODULE_INSTALL_INCLUDEDIR}"
        COMPONENT lib)