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

#ifndef _IRTKVELOCITYTODISPLACEMENTFIELD_H

#define _IRTKVELOCITYTODISPLACEMENTFIELD_H

#include <irtkImage.h>
#include <irtkImageToImage.h>

/**
 * Computes a displacement field from a stationary velocity field.
 *
 * Base class for image filters which compute the exponential map of a
 * stationary velocity field. The filter output is a diffeomorphic
 * displacement field.
 */

template <class VoxelType>
class irtkVelocityToDisplacementField : public irtkImageToImage<VoxelType>
{
  irtkImageFilterMacro(irtkVelocityToDisplacementField);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Vector field interpolation mode
  irtkPublicAttributeMacro(irtkInterpolationMode, Interpolation);
  
  /// Vector field extrapolation mode
  irtkPublicAttributeMacro(irtkExtrapolationMode, Extrapolation);
  
  /// Whether to compute interpolation coefficients from the given input
  /// or if the input images contain the coefficients already
  irtkPublicAttributeMacro(bool, ComputeInterpolationCoefficients);

  /// Whether filter should compute the inverse displacement field
  ///
  /// Note that this is equivalent to simply changing the sign of \c _T.
  /// It is in particular used by the irtkDisplacementToVelocityField filter.
  irtkPublicAttributeMacro(bool, ComputeInverse);

  /// Number of integration steps
  irtkPublicAttributeMacro(int, NumberOfSteps);

  /// Upper integration limit (negative <=> inverse)
  irtkPublicAttributeMacro(double, T);

protected:

  /// Optional input displacement field which is composed with the
  /// exponential of the velocity field, i.e., exp(v) ° phi
  irtkGenericImage<VoxelType> *_InputDisplacementField;

  /// Initialize filter
  virtual void Initialize();

  /// Finalize filter
  virtual void Finalize();

  /// Constructor
  irtkVelocityToDisplacementField();

public:  

  /// Destructor
  virtual ~irtkVelocityToDisplacementField();

  // Do not overwrite base class methods
  using irtkImageToImage<VoxelType>::SetInput;
  using irtkImageToImage<VoxelType>::GetInput;

  /// Set n-th input (0: velocity field, 1: input displacement field, optional)
  virtual void SetInput(int, irtkGenericImage<VoxelType> *);

  /// Get n-th input (0: velocity field, 1: input displacement field, optional)
  virtual irtkGenericImage<VoxelType> *GetInput(int);

};

///////////////////////////////////////////////////////////////////////////////
// Actual implementations
///////////////////////////////////////////////////////////////////////////////

#include <irtkVelocityToDisplacementFieldEuler.h>
#include <irtkVelocityToDisplacementFieldSS.h>

///////////////////////////////////////////////////////////////////////////////
// exp function for vector fields
///////////////////////////////////////////////////////////////////////////////

template <class VoxelType>
void exp(irtkGenericImage<VoxelType> *v)
{
  irtkVelocityToDisplacementFieldSS<VoxelType> vtod;
  vtod.SetInput (v);
  vtod.SetOutput(v);
  vtod.Run();
}


#endif
