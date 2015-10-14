/* The Image Registration Toolkit (IRTK)
 *
 * Copyright 2008-2015 Imperial College London
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. */

#ifndef _IRTKGAUSSIANBLURRINGWITHPADDING_H
#define _IRTKGAUSSIANBLURRINGWITHPADDING_H

#include <irtkGaussianBlurring.h>


/**
 * Class for Gaussian blurring of padded images
 *
 * This class defines and implements the Gaussian blurring of padded images.
 * The blurring is implemented by three successive 1D convolutions with a 1D
 * Gaussian kernel. If more than 50% of the voxels used for the convolution
 * have intensities smaller or equal to the padding value, the blurred voxel
 * will be filled with the padding value.
 */

template <class VoxelType>
class irtkGaussianBlurringWithPadding : public irtkGaussianBlurring<VoxelType>
{
  irtkInPlaceImageFilterMacro(irtkGaussianBlurringWithPadding);

  /// Padding value
  irtkAttributeMacro(VoxelType, PaddingValue);

public:

  /// Constructor
  irtkGaussianBlurringWithPadding(double, VoxelType);

  /// Constructor
  irtkGaussianBlurringWithPadding(double, double, VoxelType);

  /// Constructor
  irtkGaussianBlurringWithPadding(double, double, double, VoxelType);

  /// Constructor
  irtkGaussianBlurringWithPadding(double, double, double, double, VoxelType);

  /// Run Gaussian blurring
  virtual void Run();

};


#endif
