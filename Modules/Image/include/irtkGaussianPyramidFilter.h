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

#ifndef _IRTKGAUSSIANPYRAMIDFILTER_H
#define _IRTKGAUSSIANPYRAMIDFILTER_H

#include <irtkImageToImage.h>
#include <irtkImageFunction.h>


/**
 * Filter for down-/upsampling of a Gaussian pyramid image from one level to another
 */

template <class VoxelType>
class irtkGaussianPyramidFilter : public irtkImageToImage<VoxelType>
{
  irtkImageFilterMacro(irtkGaussianPyramidFilter);

  /// Level of the input image
  irtkPublicAttributeMacro(int, InputLevel);

  /// Level of the output image
  irtkPublicAttributeMacro(int, OutputLevel);

protected:

  /// Initialize the filter
  virtual void Initialize();

  /// Downsample image by one level
  virtual void Downsample();

  /// Upsample image by one level
  virtual void Upsample();

public:

  /// Constructor
  irtkGaussianPyramidFilter(int = 1);

  /// Constructor
  irtkGaussianPyramidFilter(int, int);

  /// Run the downsampling filter
  virtual void Run();

};


#endif
