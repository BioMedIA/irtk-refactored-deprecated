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

#ifndef _IRTKRESAMPLING_H

#define _IRTKRESAMPLING_H

#include <irtkImageToImage.h>

#include <irtkImageFunction.h>


/**
 * Class for resampling of images
 *
 * This class defines and implements the resampling of images with arbitrary
 * voxel dimensions.  The new image intensity of the voxels is calculated by
 * interpolation of the old image intensities. Possible interpolation schemes
 * are nearest neighbor, linear, cubic spline and B-spline interpolation.
 */

template <class VoxelType>
class irtkResampling : public irtkImageToImage<VoxelType>
{
  irtkImageFilterMacro(irtkResampling);

protected:

  /// Size of output after resampling
  int _X;
  int _Y;
  int _Z;

  /// Voxel size of output after resampling
  double _XSize;
  double _YSize;
  double _ZSize;

  /// Image function used to interpolate/extrapolate input
  irtkImageFunction *_Interpolator;

  /// Initialize the filter
  virtual void Initialize();

  /// Initialize filter output
  virtual void InitializeOutput();

public:

  /// Constructor
  irtkResampling(double, double, double);

  /// Constructor
  irtkResampling(int, int, int);

  /// Constructor
  irtkResampling(int, int, int, double, double, double);

  /// Set image size after resampling
  SetMacro(X, int);

  /// Get image size after resampling
  GetMacro(X, int);

  /// Set image size after resampling
  SetMacro(Y, int);

  /// Get image size after resampling
  GetMacro(Y, int);

  /// Set image size after resampling
  SetMacro(Z, int);

  /// Get image size after resampling
  GetMacro(Z, int);

  /// Set voxel size after resampling
  SetMacro(XSize, double);

  /// Get voxel size after resampling
  GetMacro(XSize, double);

  /// Set voxel size after resampling
  SetMacro(YSize, double);

  /// Get voxel size after resampling
  GetMacro(YSize, double);

  /// Set voxel size after resampling
  SetMacro(ZSize, double);

  /// Get voxel size after resampling
  GetMacro(ZSize, double);

  /// Set interpolator
  SetMacro(Interpolator, irtkImageFunction *);

  /// Get interpolator
  GetMacro(Interpolator, irtkImageFunction *);

  /// Run the resampling filter
  virtual void Run();

};

#include <irtkResamplingWithPadding.h>

#endif
