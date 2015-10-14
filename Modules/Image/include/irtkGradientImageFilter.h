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

#ifndef _IRTKGRADIENTIMAGEFILTER_H

#define _IRTKGRADIENTIMAGEFILTER_H

/**
 * Class for calculating the gradient of an image.
 *
 * The class provides an interface to calculating the gradient in the
 * x- , y- and z- directions.
 */

#include <irtkImageToImage.h>

template <class VoxelType>
class irtkGradientImageFilter : public irtkImageToImage<VoxelType>
{
  irtkImageFilterMacro(irtkGradientImageFilter);

  /// Whether to return gradient in mm or voxel units
  irtkPublicAttributeMacro(bool, UseVoxelSize);

  /// Whether to return gradient in world orientation
  irtkPublicAttributeMacro(bool, UseOrientation);

protected:

  /// Padding value
  VoxelType _Padding;

  /// Type of gradient
  int _type;

  /** Initialize the filter. This function must be called by any derived
   *  filter class to perform some initialize tasks. */
  virtual void Initialize();

public:

  // Type of gradient vector to compute
  const static int GRADIENT_X          = 0;
  const static int GRADIENT_Y          = 1;
  const static int GRADIENT_Z          = 2;
  const static int GRADIENT_MAGNITUDE  = 3;
  const static int GRADIENT_VECTOR     = 4;
  const static int NORMALISED_GRADIENT_VECTOR = 5;

  /// Constructor
  irtkGradientImageFilter(int type = irtkGradientImageFilter::GRADIENT_MAGNITUDE);

  /// Run the convolution filter
  virtual void Run();

  irtkSetMacro(Padding, VoxelType);
  irtkGetMacro(Padding, VoxelType);
};


#endif
