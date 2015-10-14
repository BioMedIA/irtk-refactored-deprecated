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

#ifndef _IRTKNOISE_H

#define _IRTKNOISE_H

// irtk includes
#include <irtkImageToImage.h>

/**
 * Class for adding noise to images
 *
 * This class defines an abstract filter for adding noise to images.
 *
 */
template <class VoxelType>
class irtkNoise : public irtkImageToImage<VoxelType>
{
  irtkInPlaceImageFilterMacro(irtkNoise);

protected:

  /// Initialization for the noise filter
  long int _Init;

  /// Amplitude for the noise filter
  double _Amplitude;

public:

  /// Constructor
  irtkNoise(double amplitude = 1);

  /// Destructor (empty).
  ~irtkNoise() {};

  /// Set amplitude
  SetMacro(Amplitude, double);

  /// Get amplitude
  GetMacro(Amplitude, double);

};

#include <irtkUniformNoise.h>
#include <irtkUniformNoiseWithPadding.h>
#include <irtkGaussianNoise.h>
#include <irtkGaussianNoiseWithPadding.h>
#include <irtkRicianNoise.h>
#include <irtkRicianNoiseWithPadding.h>

#endif
