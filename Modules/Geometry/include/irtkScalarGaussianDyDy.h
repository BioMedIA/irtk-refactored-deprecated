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

#ifndef _IRTKSCALARGAUSSIANDYDY_H
#define _IRTKSCALARGAUSSIANDYDY_H

#include <irtkScalarGaussian.h>


/**
 * Second derivative w.r.t y then y of (multi-dimensional) Gaussian function
 */

class irtkScalarGaussianDyDy : public irtkScalarGaussian
{
  irtkObjectMacro(irtkScalarGaussianDyDy);

public:

  /// Construct Gaussian with isotropic sigma of 1 and center at origin
  irtkScalarGaussianDyDy();

  /// Construct Gaussian with isotropic sigma and center at origin
  irtkScalarGaussianDyDy(double);

  /// Construct Gaussian with isotropic sigma and specific center
  irtkScalarGaussianDyDy(double, double, double, double = .0, double = .0);

  /// Construct 3D Gaussian with anisotropic sigma and specific center
  irtkScalarGaussianDyDy(double, double, double,
                         double, double, double);

  /// Construct 4D Gaussian with anisotropic sigma and specific center
  irtkScalarGaussianDyDy(double, double, double, double,
                         double, double, double, double);

  /// Destructor
  ~irtkScalarGaussianDyDy();

  /// Evaluate derivative of 1D Gaussian function w.r.t y then y
  double Evaluate(double);

  /// Evaluate derivative of 2D Gaussian function w.r.t y then y
  double Evaluate(double, double);

  /// Evaluate derivative of 3D Gaussian function w.r.t y then y
  double Evaluate(double, double, double);

  /// Evaluate derivative of 4D Gaussian function w.r.t y then y
  double Evaluate(double, double, double, double);

};


#endif
