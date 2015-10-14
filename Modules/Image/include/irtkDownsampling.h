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

#ifndef _IRTKDOWNSAMPLING_H

#define _IRTKDOWNSAMPLING_H

#include <irtkImageToImage.h>
#include <irtkImageFunction.h>


/**
 * Filter for downsampling of images by a power of two
 *
 * By default, a Gaussian kernel is used to average the voxel intensities within
 * the neighborhood corresponding to each downsampled voxel. This can be disabled
 * by specifying a zero sigma value for the Gaussian kernel.
 */

template <class VoxelType>
class irtkDownsampling : public irtkImageToImage<VoxelType>
{
  irtkImageFilterMacro(irtkDownsampling);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Downsampling factor in x dimension
  irtkPublicAttributeMacro(int, DownsampleFactorX);

  /// Downsampling factor in y dimension
  irtkPublicAttributeMacro(int, DownsampleFactorY);

  /// Downsampling factor in z dimension
  irtkPublicAttributeMacro(int, DownsampleFactorZ);

  /// Custom kernel to be applied at each downsampled voxel or NULL
  irtkScalarFunction *_Kernel;

  /// Extend of discretized kernel
  int _KernelSize;

  /// Whether to normalize the discretized kernel weights
  bool _NormalizeKernel;

protected:

  /// Initialize the filter
  virtual void Initialize();

  /// Discretize continious kernel
  void DiscreteKernel(irtkGenericImage<irtkRealPixel> &, double);

public:

  /// Constructor
  irtkDownsampling(int = 2);

  /// Constructor
  irtkDownsampling(int, int, int = 1);

  /// Set downsampling kernel and desired radius of discretized kernel
  void Kernel(irtkScalarFunction *, int, bool = true);

  /// Run the downsampling filter
  virtual void Run();

};


#endif
