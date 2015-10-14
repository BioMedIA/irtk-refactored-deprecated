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

#ifndef _IRTKGAUSSIANBLURRINGWITHPADDING2D_H
#define _IRTKGAUSSIANBLURRINGWITHPADDING2D_H

#include <irtkGaussianBlurringWithPadding.h>


/**
 * Class for Gaussian blurring of padded images
 *
 * This class defines and implements the Gaussian blurring of padded images.
 * It takes 2D and 3D images but blurres only in the x and y direction.
 * The blurring is implemented by two successive 1D convolutions with a 1D
 * Gaussian kernel. If more than 50% of the voxels used for the convolution
 * have intensities smaller or equal to the padding value, the blurred voxel
 * will be filled with the padding value.
 */

template <class VoxelType>
class irtkGaussianBlurringWithPadding2D : public irtkGaussianBlurringWithPadding<VoxelType>
{
  irtkInPlaceImageFilterMacro(irtkGaussianBlurringWithPadding2D);

public:

  /// Constructor
  irtkGaussianBlurringWithPadding2D(double = 1.0, VoxelType = -1);

  /// Constructor
  irtkGaussianBlurringWithPadding2D(double, double, VoxelType);

};


#endif
