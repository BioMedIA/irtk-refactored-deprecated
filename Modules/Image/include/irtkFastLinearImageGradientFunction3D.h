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

#ifndef _IRTKFASTLINEARIMAGEGRADIENTFUNCTION3D_H
#define _IRTKFASTLINEARIMAGEGRADIENTFUNCTION3D_H

#include <irtkFastLinearImageGradientFunction.h>


/**
 * Fast linear interpolation of generic 3D image gradient
 */

template <class TImage>
class irtkGenericFastLinearImageGradientFunction3D
: public irtkGenericFastLinearImageGradientFunction<TImage>
{
  irtkObjectMacro(irtkGenericFastLinearImageGradientFunction3D);
  irtkGenericGradientInterpolatorTypes(irtkGenericFastLinearImageGradientFunction);

public:

  /// Default constructor
  irtkGenericFastLinearImageGradientFunction3D();

  /// Get gradient of given image at arbitrary location (in pixels)
  ///
  /// This function is used to interpolate the image value at arbitrary
  /// locations when no extrapolator was set.
  GradientType Get(double, double, double = 0, double = 0) const;

  /// Get gradient of given image at arbitrary location (in pixels)
  ///
  /// This function is used to only interpolate foreground image values.
  /// If fully outside the foreground region, the _DefaultValue is returned.
  GradientType GetWithPadding(double, double, double = 0, double = 0) const;

  /// Get gradient of given image at arbitrary location (in pixels)
  ///
  /// If the location is inside the finite domain of the image, an actual image
  /// instance can be passed as first argument directly such as an instance of
  /// irtkGenericImage. Otherwise, an image function which extends the finite
  /// image domain to an infinite lattice is needed, i.e., an instance of a
  /// subclass of irtkExtrapolateImageFunction.
  template <class TOtherImage>
  GradientType Get(const TOtherImage *, double, double, double = 0, double = 0) const;

  /// Get gradient of given image at arbitrary location (in pixels)
  ///
  /// This function is used to only interpolate foreground image values.
  /// If fully outside the foreground region, the _DefaultValue is returned.
  ///
  /// If the location is inside the finite domain of the image, an actual image
  /// instance can be passed as first argument directly such as an instance of
  /// irtkGenericImage. Otherwise, an image function which extends the finite
  /// image domain to an infinite lattice is needed, i.e., an instance of a
  /// subclass of irtkExtrapolateImageFunction.
  template <class TOtherImage>
  GradientType GetWithPadding(const TOtherImage *, double, double, double = 0, double = 0) const;

  /// Evaluate image gradient without handling boundary conditions
  ///
  /// This version is faster than EvaluateOutside, but is only defined inside
  /// the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual GradientType GetInside(double, double, double = 0, double = 0) const;

  /// Evaluate image gradient at an arbitrary location (in pixels)
  virtual GradientType GetOutside(double, double, double = 0, double = 0) const;

  /// Evaluate image gradient without handling boundary conditions
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  ///
  /// This version is faster than GetWithPaddingOutside, but is only defined
  /// inside the domain for which all image values required for interpolation
  /// are defined and thus require no extrapolation of the finite image.
  virtual GradientType GetWithPaddingInside(double, double, double = 0, double = 0) const;

  /// Evaluate image gradient at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  virtual GradientType GetWithPaddingOutside(double, double, double = 0, double = 0) const;

};

/**
 * Fast linear interpolation of any 3D image gradient
 */

class irtkFastLinearImageGradientFunction3D
: public irtkGenericFastLinearImageGradientFunction3D<irtkBaseImage>
{
  irtkObjectMacro(irtkFastLinearImageGradientFunction3D);

public:

  /// Constructor
  irtkFastLinearImageGradientFunction3D() {}
  
};


#endif
