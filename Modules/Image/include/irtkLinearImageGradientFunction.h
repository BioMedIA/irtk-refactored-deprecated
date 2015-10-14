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

#ifndef _IRTKLINEARIMAGEGRADIENTFUNCTION_H
#define _IRTKLINEARIMAGEGRADIENTFUNCTION_H

#include <irtkImageGradientFunction.h>
#include <irtkLinearInterpolateImageFunction.h>


/**
 * Linear interpolation of generic image gradient
 */

template <class TImage>
class irtkGenericLinearImageGradientFunction
: public irtkGenericImageGradientFunction<TImage>
{
  irtkGenericGradientInterpolatorMacro(
    irtkGenericLinearImageGradientFunction,
    Interpolation_Linear
  );

protected:

  /// Image intensity interpolator
  irtkGenericLinearInterpolateImageFunction<TImage> _ContinuousImage;

public:

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  irtkGenericLinearImageGradientFunction();

  /// Destructor
  virtual ~irtkGenericLinearImageGradientFunction();

  /// Initialize interpolation function
  virtual void Initialize(bool = false);

  // ---------------------------------------------------------------------------
  // Domain checks

  /// Returns interval of discrete image indices whose values are needed for
  /// interpolation of the image value at a given continuous coordinate
  virtual void BoundingInterval(double, int &, int &) const;

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Get gradient of given 2D image without handling boundary conditions
  GradientType GetInside2D(double, double, double = 0, double = 0) const;

  /// Get gradient of given 2D image without handling boundary conditions
  ///
  /// This function is used to only interpolate foreground image values.
  /// If fully outside the foreground region, the _DefaultValue is returned.
  GradientType GetWithPaddingInside2D(double, double, double = 0, double = 0) const;

  /// Get gradient of given 3D image without handling boundary conditions
  GradientType GetInside3D(double, double, double = 0, double = 0) const;

  /// Get gradient of given 3D image without handling boundary conditions
  ///
  /// This function is used to only interpolate foreground image values.
  /// If fully outside the foreground region, the _DefaultValue is returned.
  GradientType GetWithPaddingInside3D(double, double, double = 0, double = 0) const;

  /// Evaluate generic image gradient without handling boundary conditions
  ///
  /// This version is faster than EvaluateOutside, but is only defined inside
  /// the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual GradientType GetInside(double, double, double = 0, double = 0) const;

  /// Evaluate generic image gradient at an arbitrary location (in pixels)
  /// \returns Always returns _DefaultValue.
  virtual GradientType GetOutside(double, double, double = 0, double = 0) const;

  /// Evaluate generic image gradient without handling boundary conditions
  virtual GradientType GetWithPaddingInside(double, double, double = 0, double = 0) const;

  /// Evaluate generic image gradient at an arbitrary location (in pixels)
  /// \returns Always returns _DefaultValue.
  virtual GradientType GetWithPaddingOutside(double, double, double = 0, double = 0) const;

};

/**
 * Linear interpolation of any scalar image gradient
 */

class irtkLinearImageGradientFunction
: public irtkGenericLinearImageGradientFunction<irtkBaseImage>
{
  irtkObjectMacro(irtkLinearImageGradientFunction);

public:

  /// Constructor
  irtkLinearImageGradientFunction() {}

};


#endif
