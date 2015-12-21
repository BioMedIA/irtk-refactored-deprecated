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

include_directories(BEFORE "${CMAKE_SOURCE_DIR}/Modules/Common/include")
include_directories(BEFORE "${CMAKE_SOURCE_DIR}/Modules/Geometry/include")

set(IRTK_MODULE_SOURCES irtkAnisoDiffusion.cc
                        irtkBaseImage.cc
                        irtkConvolution_1D.cc
                        irtkConvolution_2D.cc
                        irtkConvolution_3D.cc
                        irtkConvolution.cc
                        irtkConvolutionWithGaussianDerivative2.cc
                        irtkConvolutionWithGaussianDerivative.cc
                        irtkConvolutionWithPadding_1D.cc
                        irtkConvolutionWithPadding_2D.cc
                        irtkConvolutionWithPadding_3D.cc
                        irtkDensity.cc
                        irtkDifferenceOfCompositionLieBracketImageFilter3D.cc
                        irtkDilation.cc
                        irtkDisplacementToVelocityFieldBCH.cc
                        irtkDownsampling.cc
                        irtkErosion.cc
                        irtkEuclideanDistanceTransform.cc
                        irtkExtrapolateImageFunction.cc
                        irtkFileANALYZEToImage.cc
                        irtkFileGIPLToImage.cc
                        irtkFileNIFTIToImage.cc
                        irtkFileOpenCVToImage.cc
                        irtkFilePGMToImage.cc
                        irtkFileToImage.cc
                        irtkGaussianBlurring2D.cc
                        irtkGaussianBlurring4D.cc
                        irtkGaussianBlurring.cc
                        irtkGaussianBlurringWithPadding2D.cc
                        irtkGaussianBlurringWithPadding.cc
                        irtkGaussianNoise.cc
                        irtkGaussianNoiseWithPadding.cc
                        irtkGaussianPyramidFilter.cc
                        irtkGenericImage.cc
                        irtkGradientImage.cc
                        irtkGradientImageFilter.cc
                        irtkGradientImageX.cc
                        irtkGradientImageY.cc
                        irtkGradientImageZ.cc
                        irtkHessianImageFilter.cc
                        irtkHistogram_1D.cc
                        irtkHistogram_2D.cc
                        irtkImageAttributes.cc
                        irtkImageFunction.cc
                        irtkImageGradientFunction.cc
                        irtkImageHistogram_1D.cc
                        irtkImageToFileANALYZE.cc
                        irtkImageToFile.cc
                        irtkImageToFileGIPL.cc
                        irtkImageToFileNIFTI.cc
                        irtkImageToFilePGM.cc
                        irtkImageToFilePNG.cc
                        irtkImageToImage2.cc
                        irtkImageToImage.cc
                        irtkImageToOpenCv.cc
                        irtkInterpolateImageFunction.cc
                        irtkLargestConnectedComponent.cc
                        irtkLargestConnectedComponentIterative.cc
                        irtkMedianFilter.cc
                        irtkModeFilter.cc
                        irtkNoise.cc
                        irtkNonLocalMedianFilter.cc
                        irtkNormalizeNyul.cc
                        irtkRegionFilter.cc
                        irtkResampling.cc
                        irtkResamplingWithPadding.cc
                        irtkRicianNoise.cc
                        irtkRicianNoiseWithPadding.cc
                        irtkScalarFunctionToImage.cc
                        irtkScalingAndSquaring.cc
                        irtkShapeBasedInterpolateImageFunction.cc
                        irtkUniformNoise.cc
                        irtkUniformNoiseWithPadding.cc
                        irtkVelocityToDisplacementField.cc
                        irtkVelocityToDisplacementFieldEuler.cc
                        irtkVelocityToDisplacementFieldSS.cc
                        irtkVoxel.cc)

if(WITH_VTK)
  list(APPEND IRTK_MODULE_SOURCES irtkFileVTKToImage.cc
                                  irtkImageToFileVTK.cc)
endif()

add_library(${IRTK_MODULE_NAME} ${IRTK_MODULE_SOURCES})

target_link_libraries(${IRTK_MODULE_NAME} Common Geometry)

set_target_properties(${IRTK_MODULE_NAME} PROPERTIES
                      VERSION ${IRTK_VERSION}
                      SOVERSION ${IRTK_SOVERSION}
                      OUTPUT_NAME ${IRTK_MODULE_OUTPUT_NAME})

install(TARGETS ${IRTK_MODULE_NAME}
        DESTINATION "${IRTK_INSTALL_LIBDIR}"
        INCLUDES DESTINATION "${IRTK_MODULE_INSTALL_INCLUDEDIR}"
        COMPONENT lib)