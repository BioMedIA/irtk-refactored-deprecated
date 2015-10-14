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

#ifndef _IRTKGAUSSIANBLURRING4D_H
#define _IRTKGAUSSIANBLURRING4D_H

#include <irtkGaussianBlurring.h>


/**
 * Class for Gaussian blurring of image sequences
 *
 * This class defines and implements the Gaussian blurring of image sequences.
 * The blurring is implemented by four successive 1D convolutions with a 1D
 * Gaussian kernel.
 */

template <class VoxelType>
class irtkGaussianBlurring4D : public irtkGaussianBlurring<VoxelType>
{
  irtkInPlaceImageFilterMacro(irtkGaussianBlurring4D);

public:

  /// Constructor
  irtkGaussianBlurring4D(double = 1.0);

  /// Constructor
  irtkGaussianBlurring4D(double, double, double, double);

  /// Destructor
  ~irtkGaussianBlurring4D();

};


#endif
