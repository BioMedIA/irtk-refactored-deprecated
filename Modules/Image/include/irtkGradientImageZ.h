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

#ifndef _IRTKGRADIENTIMAGE_Z_H

#define _IRTKGRADIENTIMAGE_Z_H

#include <irtkGradientImage.h>
/**
 * Class for caluclating the gradient of an image in the z-direction.
 */

template <class VoxelType>
class irtkGradientImageZ : public irtkGradientImage<VoxelType>
{
  irtkImageFilterMacro(irtkGradientImageZ);

protected:

  // Calculate the gradient on a single voxel.
  virtual double Run(int, int, int, int);

public:

  /// Run the convolution filter
  virtual void Run();
};


#endif
