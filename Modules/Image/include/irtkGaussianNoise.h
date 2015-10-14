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

#ifndef _IRTKGAUSSIANNOISE_H

#define _IRTKGAUSSIANNOISE_H

/**
 * Class for adding Gaussian noise to images
 *
 * This class implements a filter for adding Gaussian noise to images.
 *
 */

template <class VoxelType>
class irtkGaussianNoise : public irtkNoise<VoxelType>
{
  irtkInPlaceImageFilterMacro(irtkGaussianNoise);

protected:

  /// The mean value of the noise distribution.
  double _Mean;

  /** The standard deviation of the noise distribution. This is the same as
      the amplitude of the noise.
  */
  double _Sigma;

  /// Minimum voxel value
  VoxelType _MinVal;

  /// Maximum voxel value
  VoxelType _MaxVal;

public:

  // Default constructor
  irtkGaussianNoise();

  /** Constructor.  Sets mean value and standard deviation of the noise
   *  distribution, and range of values permitted for each voxel in the image.
   */
  irtkGaussianNoise(double mean, double sigma, VoxelType min, VoxelType max);

  /// Destructor (empty).
  ~irtkGaussianNoise() {};

  // Run gaussian noise filter
  virtual double Run(int, int, int, int);

  /// Set mean of Gaussian noise
  SetMacro(Mean, double);

  /// Get mean of Gaussian noise
  GetMacro(Mean, double);

  /// Set standard deviation of Gaussian noise
  SetMacro(Sigma, double);

  /// Get standard deviation of Gaussian noise
  GetMacro(Sigma, double);

  /// Set minimum value for Gaussian noise
  SetMacro(MinVal, VoxelType);

  /// Get minimum value for Gaussian noise
  GetMacro(MinVal, VoxelType);

  /// Set maximum value for Gaussian noise
  SetMacro(MaxVal, VoxelType);

  /// Get maximum value for Gaussian noise
  GetMacro(MaxVal, VoxelType);
};

#endif
