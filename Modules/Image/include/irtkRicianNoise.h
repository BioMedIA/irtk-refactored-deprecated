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

#ifndef _IRTKRICIANNOISE_H

#define _IRTKRICIANNOISE_H

/**
 * Class for adding Rician noise to images
 *
 * This class implements a filter for adding Rician noise to images. The
 * Rician distribution is approiximately Gaussian for high intensity
 * signals, and Rayleigh distributed for low intensities. The Rayleigh
 * intensity distribution can be expressed as:
 * $P_M(M)=\frac{m}{\sigma^2}exp\left(\frac{-M^2}{2\sigma^2}\right)
 * where $M$ is the actual intensity and $\sigma$ is the standard
 * deviation of Gaussian noise. For more information, see Holden et al.,
 * IEEE-TMI 19(2) 2000.
 *
 */

template <class VoxelType>
class irtkRicianNoise : public irtkNoise<VoxelType>
{
  irtkInPlaceImageFilterMacro(irtkRicianNoise);

public:

  /// Constructor
  irtkRicianNoise(double amplitude = 1);

  /// Destructor (empty).
  ~irtkRicianNoise() {};

  /// Run rician noise filter
  virtual double Run(int, int, int, int);

};

#endif

