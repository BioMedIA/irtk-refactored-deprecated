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

#ifndef _IRTKDISPLACEMENTTOVELOCITYFIELDBCH_H

#define _IRTKDISPLACEMENTTOVELOCITYFIELDBCH_H


#include <irtkVelocityToDisplacementField.h>


#ifdef USE_CUDA
template <class VoxelType> class irtkCUGenericImage;
#endif


/**
 * Computes a stationary velocity field from a given displacement field.
 *
 * This class implements an image filter which computes the group logarithm
 * of a (diffeomorphic) displacement field using the Baker-Campbell-Hausdorff (BCH)
 * formula. The result is a stationary velocity field. Integrating this velocity
 * field from 0 to _T, e.g., using N forward Euler integration steps or the
 * scaling and squaring (SS) method, yields the original displacement field
 * (with an approximation error).
 *
 * M. Bossa and S. Olmos, "A new algorithm for the computation of the group
 * logarithm of diffeomorphisms", MFCA'08
 */

template <class VoxelType>
class irtkDisplacementToVelocityFieldBCH : public irtkDisplacementToVelocityField<VoxelType>
{
  irtkObjectMacro(irtkDisplacementToVelocityFieldBCH);

public:

#ifdef USE_CUDA
  typedef irtkCUGenericImage<VoxelType> ImageType;
#else
  typedef irtkGenericImage<VoxelType>   ImageType;
#endif

private:

  /// Result of composition: exp(-v) ° phi
  ImageType *_dv;

  /// Result of Lie bracket: [v, dv]
  ImageType *_l1;

  /// Result of Lie bracket: [v, [v, dv]]
  ImageType *_l2;

  /// Result of Lie bracket: [dv, [v, dv]]
  ImageType *_l3;

  /// Result of Lie bracket: [dv, [v, [v, dv]]]
  ImageType *_l4;

protected:

  /// Image filter for computation of exponential map of inverse velocities
  irtkVelocityToDisplacementField<VoxelType> *_ExponentialFilter;

  /// Whether exponential image filter was set by user
  bool _CustomExponentialFilter;

  /// Number of update steps
  int _NumberOfIterations;

  /// Number of terms to use
  int _NumberOfTerms;

  /// Whether to use the Jacobian for the Lie bracket computation
  bool _UseJacobian;

  /// Whether to smooth the velocity fields to stabilize the computation
  bool _SmoothVelocities;

  /// Initialize filter
  virtual void Initialize();

  /// Finalize filter
  virtual void Finalize();

public:

  /// Constructor
  irtkDisplacementToVelocityFieldBCH();

  /// Destructor
  virtual ~irtkDisplacementToVelocityFieldBCH();

  /// Compute output = log(input)
  virtual void Run();

  /// Set image filter for computation of exponential map of inverse velocities
  void SetExponentialFilter(irtkVelocityToDisplacementField<VoxelType> *);

  /// Get image filter for computation of exponential map of inverse velocities
  irtkVelocityToDisplacementField<VoxelType> *GetExponentialFilter();

  /// Get image filter for computation of exponential map of inverse velocities
  const irtkVelocityToDisplacementField<VoxelType> *GetExponentialFilter() const;

  /// Set upper integration limit of exponential filter
  void SetT(double);

  /// Get upper integration limit of exponential filter
  double GetT();

  /// Set number of integration steps
  void SetNumberOfSteps(int);

  /// Get number of integration steps
  int GetNumberOfSteps();

  /// Set number of update steps
  SetMacro(NumberOfIterations, int);

  /// Get number of update steps
  GetMacro(NumberOfIterations, int);

  /// Set number of BCH terms
  SetMacro(NumberOfTerms, int);

  /// Get number of BCH terms
  GetMacro(NumberOfTerms, int);

  /// Set whether to use the Jacobian for the Lie bracket computation
  SetMacro(UseJacobian, bool);
  
  /// Get whether to use the Jacobian for the Lie bracket computation
  GetMacro(UseJacobian, bool);

  /// Set whether to smooth the velocity field to stabilize the computation
  SetMacro(SmoothVelocities, bool);

  /// Get whether to smooth the velocity field to stabilize the computation
  GetMacro(SmoothVelocities, bool);

};


#endif
