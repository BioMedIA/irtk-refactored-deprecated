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

#ifndef _IRTKUNIFORMNOISEWITHPADDING_H

#define _IRTKUNIFORMNOISEWITHPADDING_H

#include <irtkUniformNoise.h>

/**
 * Class for adding uniform noise to images which are padded
 *
 * This class implements a filter for adding uniform noise to images which
 * are padded.
 *
 */

template <class VoxelType> class irtkUniformNoiseWithPadding : public irtkUniformNoise<VoxelType>
{
  irtkInPlaceImageFilterMacro(irtkUniformNoiseWithPadding);

protected:

  /// The padding value.
  VoxelType _PaddingValue;

public:

  // Default constructor
  irtkUniformNoiseWithPadding();

  /** Constructor. Sets standard deviation of the noise distribution and
   *  the padding value.
   */
  irtkUniformNoiseWithPadding(double Amplitude, VoxelType PaddingValue);

  /// Destructor (empty).
  ~irtkUniformNoiseWithPadding() {};

  /// Set padding value
  SetMacro(PaddingValue, VoxelType);

  /// Get padding value
  GetMacro(PaddingValue, VoxelType);

  /// Run uniform noise filter
  virtual double Run(int, int, int, int);

};

#endif
