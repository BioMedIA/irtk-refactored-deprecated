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

#ifndef _IRTKGAUSSIANBLURRING2D_H
#define _IRTKGAUSSIANBLURRING2D_H

#include <irtkGaussianBlurring.h>


/**
 * Class for Gaussian blurring of images
 *
 * This class defines and implements the Gaussian blurring of images. It takes
 * 2D and 3D images but blurs only in the x and y direction. The
 * blurring is implemented by two successive 1D convolutions with a 1D
 * Gaussian kernel.
 */

template <class VoxelType>
class irtkGaussianBlurring2D : public irtkGaussianBlurring<VoxelType>
{
  irtkInPlaceImageFilterMacro(irtkGaussianBlurring2D);

public:

  /// Constructor
  irtkGaussianBlurring2D(double = 1.0);

  /// Constructor
  irtkGaussianBlurring2D(double, double);

  /// Destructor
  ~irtkGaussianBlurring2D();

};


#endif
