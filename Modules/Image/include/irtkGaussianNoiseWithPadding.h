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

#ifndef _IRTKGAUSSIANNOISEWITHPADDING_H

#define _IRTKGAUSSIANNOISEWITHPADDING_H

#include <irtkGaussianNoise.h>

/**
 * Class for adding Gaussian noise to images which are padded
 *
 * This class implements a filter for adding Gaussian noise to images which
 * are padded.
 *
 */

template <class VoxelType>
class irtkGaussianNoiseWithPadding : public irtkGaussianNoise<VoxelType>
{
  irtkInPlaceImageFilterMacro(irtkGaussianNoiseWithPadding);

protected:

  /// The padding value.
  VoxelType _PaddingValue;

public:

  // Default constructor
  irtkGaussianNoiseWithPadding();

  /** Constructor.  Sets mean value and standard deviation of the noise
   *  distribution, the padding value, and the image range.
   */
  irtkGaussianNoiseWithPadding(double Mean, double Sigma, VoxelType PaddingValue, VoxelType MinVal, VoxelType MaxVal);

  /// Destructor (empty).
  ~irtkGaussianNoiseWithPadding() {};

  /// Set padding value
  SetMacro(PaddingValue, VoxelType);

  /// Get padding value
  GetMacro(PaddingValue, VoxelType);

  /// Adds noise to input image.
  virtual double Run(int, int, int, int);

};

#endif
