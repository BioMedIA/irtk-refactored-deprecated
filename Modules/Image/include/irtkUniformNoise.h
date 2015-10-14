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

#ifndef _IRTKUNIFORMNOISE_H

#define _IRTKUNIFORMNOISE_H

/**
 * Class for adding uniform noise to images
 *
 * This class implements a filter for adding uniformly distributed noise to
 * images.
 *
 */

template <class VoxelType>
class irtkUniformNoise : public irtkNoise<VoxelType>
{
  irtkInPlaceImageFilterMacro(irtkUniformNoise);

public:

  /// Constructor
  irtkUniformNoise(double amplitude = 1);

  /// Destructor (empty).
  ~irtkUniformNoise() {};

  /// Run uniform noise filter
  virtual double Run(int, int, int, int);

};

#endif
