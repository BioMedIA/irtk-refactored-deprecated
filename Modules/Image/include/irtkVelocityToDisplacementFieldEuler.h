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

#ifndef _IRTKVELOCITYTODISPLACEMENTFIELDEULER_H

#define _IRTKVELOCITYTODISPLACEMENTFIELDEULER_H

#include <irtkVelocityToDisplacementField.h>
#include <irtkImageFunction.h>


/**
 * Computes a displacement field from a stationary velocity field.
 *
 * This class implements an image filter which computes the exponential map
 * of a stationary velocity field using the forward Euler integration method.
 * The result is a diffeomorphic displacement field.
 */

template <class VoxelType>
class irtkVelocityToDisplacementFieldEuler : public irtkVelocityToDisplacementField<VoxelType>
{
  irtkImageFilterMacro(irtkVelocityToDisplacementFieldEuler);

protected:

  /// Image function for interpolation of velocities
  irtkInterpolateImageFunction *_VelocityInterpolator;

  /// Initialize filter
  virtual void Initialize();

public:

  /// Constructor
  irtkVelocityToDisplacementFieldEuler();

  /// Destructor
  virtual ~irtkVelocityToDisplacementFieldEuler();

  /// Compute output = log(input)
  virtual void Run();

};

#endif
