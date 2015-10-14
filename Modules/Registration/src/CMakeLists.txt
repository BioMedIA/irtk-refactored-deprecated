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
set(IRTK_MODULE_OUTPUT_NAME "irtk${IRTK_MODULE_NAME}")

include_directories("${CMAKE_SOURCE_DIR}/Modules/Common/include")
include_directories("${CMAKE_SOURCE_DIR}/Modules/Geometry/include")
include_directories("${CMAKE_SOURCE_DIR}/Modules/Image/include")
include_directories("${CMAKE_SOURCE_DIR}/Modules/Transformation/include")

set(IRTK_MODULE_SOURCES irtkAdaptiveLineSearch.cc
                        irtkConjugateGradientDescent.cc
                        irtkCosineOfNormalizedGradientField.cc
                        irtkDataFidelity.cc
                        irtkEnergyTerm.cc
                        irtkGenericRegistrationDebugger.cc
                        irtkGenericRegistrationFilter.cc
                        irtkGenericRegistrationFilter_registration2++.cc
                        irtkGenericRegistrationLogger.cc
                        irtkGradientDescent.cc
                        irtkGradientFieldSimilarity.cc
                        irtkImageSimilarity.cc
                        irtkInexactLineSearch.cc
                        irtkIntensityCrossCorrelation.cc
                        irtkJacobianConstraint.cc
                        irtkLineSearch.cc
                        irtkLocalOptimizer.cc
                        irtkLogJacobianConstraint.cc
                        irtkMaxStepLineSearch.cc
                        irtkMinJacobianConstraint.cc
                        irtkMutualImageInformation.cc
                        irtkNormalizedGradientFieldSimilarity.cc
                        irtkNormalizedIntensityCrossCorrelation.cc
                        irtkNormalizedMutualImageInformation.cc
                        irtkProbabilisticImageSimilarity.cc
                        irtkRegisteredImage.cc
                        irtkRegistrationEnergy.cc
                        irtkRegistrationUtils.cc
                        irtkSmoothnessConstraint.cc
                        irtkSparsityConstraint.cc
                        irtkSumOfSquaredIntensityDifferences.cc
                        irtkTopologyPreservationConstraint.cc
                        irtkTransformationConstraint.cc)

if(WITH_LBFGS)
  list(APPEND IRTK_MODULE_SOURCES irtkLimitedMemoryBFGSDescent.cc)
endif()

if(WITH_VTK)
  list(APPEND IRTK_MODULE_SOURCES irtkBalloonForce.cc
                                  irtkClosestCell.cc
                                  irtkClosestPoint.cc
                                  irtkClosestPointLabel.cc
                                  irtkCurrentsDistance.cc
                                  irtkCurvatureConstraint.cc
                                  irtkDeformableSurfaceDebugger.cc
                                  irtkDeformableSurfaceLogger.cc
                                  irtkDeformableSurfaceModel.cc
                                  irtkEdgeLengthConstraint.cc
                                  irtkExternalForce.cc
                                  irtkEulerMethod.cc
                                  irtkEulerMethodWithMomentum.cc
                                  irtkFiducialMatch.cc
                                  irtkFiducialRegistrationError.cc
                                  irtkFuzzyCorrespondence.cc
                                  irtkFuzzyCorrespondenceUtils.cc
                                  irtkImageEdgeForce.cc
                                  irtkInflationForce.cc
                                  irtkLaplacianConstraint.cc
                                  irtkNonSelfIntersectionConstraint.cc
                                  irtkPointCorrespondence.cc
                                  irtkPointCorrespondenceDistance.cc
                                  irtkPointLocator.cc
                                  irtkPointSetConstraint.cc
                                  irtkPointSetDistance.cc
                                  irtkPointSetForce.cc
                                  irtkRadialErrorFunction.cc
                                  irtkRegisteredPointSet.cc
                                  irtkRegisteredSurface.cc
                                  irtkRobustClosestPoint.cc
                                  irtkRobustPointMatch.cc
                                  irtkSpectralDecomposition.cc
                                  irtkSpectralMatch.cc
                                  irtkSurfaceCollisions.cc
                                  irtkSurfaceConstraint.cc
                                  irtkSurfaceDistance.cc
                                  irtkSurfaceForce.cc)
endif()

add_library(${IRTK_MODULE_NAME} ${IRTK_MODULE_SOURCES})

if(WITH_VTK)
  include_directories("${CMAKE_SOURCE_DIR}/Modules/PolyData/include")
  target_link_libraries(${IRTK_MODULE_NAME} PolyData ${VTK_LIBRARIES})
endif()

set_target_properties(${IRTK_MODULE_NAME} PROPERTIES
                      VERSION ${IRTK_VERSION}
                      SOVERSION ${IRTK_SOVERSION}
                      OUTPUT_NAME ${IRTK_MODULE_OUTPUT_NAME})

install(TARGETS ${IRTK_MODULE_NAME}
        DESTINATION "${IRTK_INSTALL_LIBDIR}"
        INCLUDES DESTINATION "${IRTK_MODULE_INSTALL_INCLUDEDIR}"
        COMPONENT lib)