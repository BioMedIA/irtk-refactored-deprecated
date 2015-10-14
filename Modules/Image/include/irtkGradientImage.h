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

#ifndef _IRTKGRADIENTIMAGE_H

#define _IRTKGRADIENTIMAGE_H

/**
 * Class for caluclating the gradient of an image.
 * The class provides an iterface to subclasses which calculate the gradient in
 * x- , y- and z- directions.
 */

#include <irtkImageToImage.h>

template <class VoxelType>
class irtkGradientImage : public irtkImageToImage<VoxelType>
{
  irtkImageFilterMacro(irtkGradientImage);

protected:

  /// padding value
  VoxelType _Padding;

  /// Runs the filter on a single voxel.
  virtual double Run(int, int, int, int);

public:
  /// Constructor
  irtkGradientImage();
  /// Set Padding
  virtual SetMacro(Padding,VoxelType);

  /// Run the convolution filter
  virtual void Run();

};

#include <irtkGradientImageX.h>
#include <irtkGradientImageY.h>
#include <irtkGradientImageZ.h>

#endif
