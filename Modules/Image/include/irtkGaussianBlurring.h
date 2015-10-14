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

#ifndef _IRTKGAUSSIANBLURRING_H
#define _IRTKGAUSSIANBLURRING_H

#include <irtkImageToImage.h>


/**
 * Class for Gaussian blurring of images
 *
 * This class defines and implements the Gaussian blurring of images. The
 * blurring is implemented by three successive 1D convolutions with a 1D
 * Gaussian kernel.
 *
 * By default, if one isotropic Gaussian standard deviation is specified,
 * the first three dimensions of the image are blurred only. Otherwise,
 * the 1D convolution with a 1D Gaussian kernel is performed only for
 * dimensions of more than one voxel size and for which a non-zero
 * standard deviation for the Gaussian kernel has been set.
 */

template <class VoxelType>
class irtkGaussianBlurring : public irtkImageToImage<VoxelType>
{
  irtkInPlaceImageFilterMacro(irtkGaussianBlurring);

  /// Standard deviation of Gaussian kernel in x
  irtkAttributeMacro(double, SigmaX);

  /// Standard deviation of Gaussian kernel in y
  irtkAttributeMacro(double, SigmaY);

  /// Standard deviation of Gaussian kernel in z
  irtkAttributeMacro(double, SigmaZ);

  /// Standard deviation of Gaussian kernel in t
  irtkAttributeMacro(double, SigmaT);

protected:

  /// Gaussian convolution kernel
  irtkGenericImage<irtkRealPixel> *_Kernel;

  /// Initialize filter
  virtual void Initialize();

  /// Initialize 1D Gaussian kernel with sigma given in voxel units
  virtual void InitializeKernel(double);

  /// Finalize filter
  virtual void Finalize();

public:

  /// Constructor
  irtkGaussianBlurring(double = 1.0);

  /// Constructor
  irtkGaussianBlurring(double, double, double = .0, double = .0);

  /// Destructor
  ~irtkGaussianBlurring();

  /// Run Gaussian blurring
  virtual void Run();

  /// Run Gaussian blurring along x only
  virtual void RunX();

  /// Run Gaussian blurring along y only
  virtual void RunY();

  /// Run Gaussian blurring along z only
  virtual void RunZ();

  /// Run Gaussian blurring along t only
  virtual void RunT();

  /// Set sigma
  virtual void SetSigma(double);

  /// Set sigma
  virtual void SetSigma(double, double, double = .0, double = .0);

  /// Kernel size used for a given sigma (divided by voxel size)
  static int KernelSize(double);

};

#include <irtkGaussianBlurringWithPadding.h>


#endif
