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

#ifndef _IRTKSCALINGANDSQUARING_H
#define _IRTKSCALINGANDSQUARING_H

#include <irtkImage.h>
#include <irtkInterpolateImageFunction.h>


/**
 * Computes the exponential map of a SVF and its derivatives
 *
 * This class implements an image filter which computes the exponential map
 * of a stationary velocity field using the scaling and squaring method.
 * The result is a diffeomorphic displacement field. Additionally, this filter
 * computes also the derivatives of the exponential map either with respect to
 * the spatial coordinate or the stationary velocity field itself.
 */

template <class TReal>
class irtkScalingAndSquaring : public irtkObject
{
  irtkObjectMacro(irtkScalingAndSquaring);

  // ---------------------------------------------------------------------------
  // Types

public:

  /// Image type of input, output, and intermediate images
  typedef irtkGenericImage<TReal> ImageType;

  /// Type of displacement field interpolator
  typedef irtkGenericInterpolateImageFunction<ImageType> DisplacementField;

  /// Type of Jacobian field interpolator
  typedef irtkGenericInterpolateImageFunction<ImageType> JacobianField;

  /// Type of input velocity field interpolator
  typedef irtkGenericFastCubicBSplineInterpolateImageFunction<ImageType>  VelocityField;

private:

  // ---------------------------------------------------------------------------
  // Input

  /// Input velocity field
  irtkPublicAggregateMacro(const ImageType, InputVelocity);

  /// Input displacement field
  irtkPublicAggregateMacro(const ImageType, InputDisplacement);

  /// Input deformation field
  irtkPublicAggregateMacro(const ImageType, InputDeformation);

  // ---------------------------------------------------------------------------
  // Interim

  /// Intermediate displacement field
  irtkComponentMacro(ImageType, InterimDisplacement);

  /// Intermediate Jacobian w.r.t. x
  irtkComponentMacro(ImageType, InterimJacobian);

  /// Intermediate determinant of Jacobian w.r.t. x
  irtkComponentMacro(ImageType, InterimDetJacobian);

  /// Intermediate log of determinant of Jacobian w.r.t. x
  irtkComponentMacro(ImageType, InterimLogJacobian);

  /// Intermediate Jacobian w.r.t. v
  irtkComponentMacro(ImageType, InterimJacobianDOFs);

  /// Interpolator of intermediate displacement field
  irtkComponentMacro(DisplacementField, Displacement);

  /// Interpolator of intermediate Jacobian w.r.t. x
  irtkComponentMacro(JacobianField, Jacobian);

  /// Interpolator of intermediate determinant of Jacobian w.r.t. x
  irtkComponentMacro(JacobianField, DetJacobian);

  /// Interpolator of intermediate log of determinant of Jacobian w.r.t. x
  irtkComponentMacro(JacobianField, LogJacobian);

  /// Interpolator of intermediate Jacobian w.r.t. v
  irtkComponentMacro(JacobianField, JacobianDOFs);

  // ---------------------------------------------------------------------------
  // Output

  /// Output displacement field
  irtkPublicAggregateMacro(ImageType, OutputDisplacement);

  /// Output deformation field
  irtkPublicAggregateMacro(ImageType, OutputDeformation);

  /// Jacobian of output deformation field w.r.t. x
  irtkPublicAggregateMacro(ImageType, OutputJacobian);

  /// Determinant of Jacobian of output deformation field w.r.t. x
  irtkPublicAggregateMacro(ImageType, OutputDetJacobian);

  /// Log of determinant of Jacobian of output deformation field w.r.t. x
  irtkPublicAggregateMacro(ImageType, OutputLogJacobian);

  /// Jacobian of output deformation (or displacement) field w.r.t. v
  irtkPublicAggregateMacro(ImageType, OutputJacobianDOFs);

  // ---------------------------------------------------------------------------
  // Settings

  /// Attributes of intermediate images, defaults to output attributes
  irtkPublicAttributeMacro(irtkImageAttributes, InterimAttributes);

  /// Attributes of output images, defaults to input attributes
  irtkPublicAttributeMacro(irtkImageAttributes, OutputAttributes);

  /// Interpolation used for each squaring step
  irtkPublicAttributeMacro(irtkInterpolationMode, Interpolation);

  /// Whether to compute interpolation coefficients from the given input.
  /// Set to \c false if the input is the coefficient image.
  irtkPublicAttributeMacro(bool, ComputeInterpolationCoefficients);

  /// Whether filter should invert the velocities
  /// This is equivalent to changing the sign of the upper integration limit.
  irtkPublicAttributeMacro(bool, ComputeInverse);

  /// Upper integration limit (negative <=> inverse)
  irtkPublicAttributeMacro(double, IntegrationLimit);

  /// Number of integration steps
  irtkPublicAttributeMacro(int, NumberOfSteps);

  /// Number of squaring steps
  irtkPublicAttributeMacro(int, NumberOfSquaringSteps);

  /// Maximum velocity after scaling
  ///
  /// Set to zero in order to scale velocities by exactly
  /// pow(2.0, _NumberOfSquaringSteps). Otherwise, the number of squaring steps
  /// is increased in order to ensure that the maximum velocity in each
  /// dimension is less or equal the specified value.
  irtkPublicAttributeMacro(double, MaxScaledVelocity);

  /// Whether to upsample the input velocity field
  irtkPublicAttributeMacro(bool, Upsample);

  /// Whether to use a Gaussian pyramid downsampler when downsampling the
  /// previously upsampled output fields to obey Shannon's sampling theorem.
  /// Otherwise a simple interpolation without smoothing kernel is used.
  irtkPublicAttributeMacro(bool, SmoothBeforeDownsampling);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Free allocated memory
  void Clear();

public:

  /// Constructor
  irtkScalingAndSquaring();

  /// Destructor
  virtual ~irtkScalingAndSquaring();

  // ---------------------------------------------------------------------------
  // Evaluation

protected:

  /// Initialize filter
  virtual void Initialize();

  /// Finalize filter
  virtual void Finalize();

  /// Resample intermediate filter output
  template <class TInterpolator> void Resample(ImageType *, TInterpolator *);

public:

  /// Compute output deformation/displacement field and its derivatives
  virtual void Run();

};


#endif
