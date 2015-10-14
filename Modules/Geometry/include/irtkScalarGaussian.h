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

#ifndef _IRTKSCALARGAUSSIAN_H
#define _IRTKSCALARGAUSSIAN_H


#include <irtkVector4D.h>


/**
 * Scalar Gaussian function class.
 */

class irtkScalarGaussian : public irtkScalarFunction
{
  irtkObjectMacro(irtkScalarGaussian);

protected:

  /// Variance of the Gaussian
  irtkVector4D<double> _Variance;

  /// Center of Gaussian
  irtkVector4D<double> _Center;

  /// Normalization of ND Gaussian function
  double _Norm[4];

  /// Pre-compute normalization factor(s)
  void ComputeNorm();

public:

  /// Construct Gaussian with isotropic sigma of 1 and center at origin
  irtkScalarGaussian();

  /// Construct Gaussian with isotropic sigma and center at origin
  irtkScalarGaussian(double);

  /// Construct Gaussian with isotropic sigma and specific center
  irtkScalarGaussian(double, double, double, double = .0, double = .0);

  /// Construct 3D Gaussian with anisotropic sigma and specific center
  irtkScalarGaussian(double, double, double,
                     double, double, double);

  /// Construct 4D Gaussian with anisotropic sigma and specific center
  irtkScalarGaussian(double, double, double, double,
                     double, double, double, double);

  /// Destructor
  virtual ~irtkScalarGaussian();

  /// Evaluate 1D Gaussian
  virtual double Evaluate(double);

  /// Evaluate 2D Gaussian
  virtual double Evaluate(double, double);

  /// Evaluate 3D Gaussian
  virtual double Evaluate(double, double, double);

  /// Evaluate 4D Gaussian
  virtual double Evaluate(double, double, double, double);
};

#endif
