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

#ifndef _IRTKHESSIANIMAGEFILTER_H

#define _IRTKHESSIANIMAGEFILTER_H

/**
 * Class for calculating the second order gradient of an image.
 *
 * The class provides an interface to calculating the second order gradient in the
 * x-x-, x-y-, x-z-, y-y-, y-z- and z-z- directions.
 */

#include <irtkImageToImage.h>

template <class VoxelType>
class irtkHessianImageFilter : public irtkImageToImage<VoxelType>
{
  irtkImageFilterMacro(irtkHessianImageFilter);

  /// Whether to return gradient in mm or voxel units
  irtkPublicAttributeMacro(bool, UseVoxelSize);

  /// Whether to return gradient in world orientation
  irtkPublicAttributeMacro(bool, UseOrientation);

protected:

  /// Type of gradient
  int _type;

  /// Padding
  int _Padding;

  /** Initialize the filter. This function must be called by any derived
   *  filter class to perform some initialize tasks. */
  virtual void Initialize();

public:

  // Type of gradient vector to compute
  const static int HESSIAN_XX = 0;
  const static int HESSIAN_XY = 1;
  const static int HESSIAN_XZ = 2;
  const static int HESSIAN_YY = 3;
  const static int HESSIAN_YZ = 4;
  const static int HESSIAN_ZZ = 5;
  const static int HESSIAN_VECTOR = 6;
  const static int HESSIAN_MATRIX = 7;

  /// Constructor
  irtkHessianImageFilter(int type = irtkHessianImageFilter::HESSIAN_MATRIX);

  /// Run the convolution filter
  virtual void Run();

  /// Set Padding
  virtual SetMacro(Padding,VoxelType);
};

#endif
