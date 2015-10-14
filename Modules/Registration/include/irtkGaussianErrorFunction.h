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

#ifndef IRTKGAUSSIANERRORFUNCTION_H
#define IRTKGAUSSIANERRORFUNCTION_H

#include <irtkRadialErrorFunction.h>


/**
 * Gaussian weighted registration error function
 */
class irtkGaussianErrorFunction : public irtkRadialErrorFunction
{
  irtkObjectMacro(irtkGaussianErrorFunction);

  /// Variance of Gaussian function
  irtkPublicAttributeMacro(double, Variance);
 
public:

  /// Constructor
  ///
  /// \param[in] sigma Standard deviation of Gaussian function.
  irtkGaussianErrorFunction(double sigma = 1.0)
  :
    _Variance(sigma * sigma)
  {}

  /// Copy constructor
  irtkGaussianErrorFunction(const irtkGaussianErrorFunction &other)
  :
    _Variance(other._Variance)
  {}

  /// Destructor
  ~irtkGaussianErrorFunction() {}

  /// Copy construct a new instance
  virtual irtkRadialErrorFunction *NewInstance() const
  {
    return new irtkGaussianErrorFunction(*this);
  }

  /// Type enumeration value
  virtual TypeId Type() const
  {
    return Gaussian;
  }

  /// Set parameter value from string
  virtual bool Set(const char *name, const char *value)
  {
    if (strcmp(name, "sigma")              == 0 ||
        strcmp(name, "standard deviation") == 0) {
      double sigma = .0;
      FromString(value, sigma);
      if (sigma == .0) return false;
      _Variance = sigma * sigma;
      return true;
    } else if (strcmp(name, "variance") == 0) {
      return FromString(value, _Variance) && _Variance > .0;
    }
    return false;
  }

  // Import other overloads
  using irtkRadialErrorFunction::Parameter;

  /// Get parameter key/value as string map
  virtual irtkParameterList Parameter() const
  {
    irtkParameterList params;
    Insert(params, "sigma", ToString(sqrt(_Variance)));
    return params;
  }

  /// Evaluate radial registration error
  virtual double Value(double d) const
  {
    return exp(-.5 * d / _Variance);
  }

  /// Evaluate derivative of radial registration error
  virtual double Derivative(double d) const
  {
    return exp(-.5 * d / _Variance) / _Variance;
  }

};


#endif
