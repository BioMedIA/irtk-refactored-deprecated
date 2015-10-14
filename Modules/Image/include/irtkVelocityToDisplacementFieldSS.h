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

#ifndef _IRTKVELOCITYTODISPLACEMENTFIELDSS_H

#define _IRTKVELOCITYTODISPLACEMENTFIELDSS_H


#include <irtkVelocityToDisplacementField.h>
#include <irtkImageFunction.h>


/**
 * Computes a displacement field from a stationary velocity field.
 *
 * This class implements an image filter which computes the exponential map
 * of a stationary velocity field using the scaling and squaring method.
 * The result is a diffeomorphic displacement field.
 */

template <class VoxelType>
class irtkVelocityToDisplacementFieldSS : public irtkVelocityToDisplacementField<VoxelType>
{
  irtkImageFilterMacro(irtkVelocityToDisplacementFieldSS);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Number of squaring steps
  irtkPublicAttributeMacro(int, NumberOfSquaringSteps);

  /// Whether to upsample the input vector field before
  irtkPublicAttributeMacro(bool, Upsample);

  /// Whether to use a Gaussian pyramid downsampler when downsampling the
  /// previously upsample vector field to obey Shannon's sampling theorem.
  /// Otherwise a simple interpolation without smoothing kernel is used.
  irtkPublicAttributeMacro(bool, SmoothBeforeDownsampling);

  /// Maximum velocity after scaling
  ///
  /// Set to zero in order to scale velocities by exactly pow(2.0, _NumberOfSquaringSteps).
  /// Otherwise, the number of squaring steps is increased in order to ensure that
  /// the maximum velocity in each dimension is less or equal the specified value.
  irtkPublicAttributeMacro(VoxelType, MaxScaledVelocity);

  /// External memory that can be used for intermedate displacement field
  irtkPublicAggregateMacro(irtkGenericImage<VoxelType>, ExternalCache);

  /// Intermediate displacement field
  irtkAggregateMacro(irtkGenericImage<VoxelType>, Displacement);

  /// Interpolator of intermediate displacement field
  irtkAggregateMacro(irtkInterpolateImageFunction, Interpolator);

protected:

  /// Returns whether the filter requires buffering
  virtual bool RequiresBuffering();

  /// Initialize filter
  virtual void Initialize();

  /// Finalize filter
  virtual void Finalize();

public:

  /// Constructor
  irtkVelocityToDisplacementFieldSS();

  /// Destructor
  virtual ~irtkVelocityToDisplacementFieldSS();

  /// Compute output = log(input)
  virtual void Run();

};


#endif
