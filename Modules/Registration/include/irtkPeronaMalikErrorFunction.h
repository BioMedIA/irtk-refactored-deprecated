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

#ifndef IRTKPERONAMALIKERRORFUNCTION_H
#define IRTKPERONAMALIKERRORFUNCTION_H

#include <irtkRadialErrorFunction.h>


/**
 * Perona and Malik fiducial registration error function
 */
class irtkPeronaMalikErrorFunction : public irtkRadialErrorFunction
{
  irtkObjectMacro(irtkPeronaMalikErrorFunction);

  /// Squared fiducial registration error threshold
  irtkPublicAttributeMacro(double, SquaredThreshold);

public:

  /// Constructor
  irtkPeronaMalikErrorFunction(double threshold = 1.0)
  :
    _SquaredThreshold(threshold * threshold)
  {}

  /// Copy constructor
  irtkPeronaMalikErrorFunction(const irtkPeronaMalikErrorFunction &other)
  :
    _SquaredThreshold(other._SquaredThreshold)
  {}

  /// Destructor
  ~irtkPeronaMalikErrorFunction() {}

  /// Copy construct a new instance
  virtual irtkRadialErrorFunction *NewInstance() const
  {
    return new irtkPeronaMalikErrorFunction(*this);
  }

  /// Type enumeration value
  virtual TypeId Type() const
  {
    return PeronaMalik;
  }

  /// Set parameter value from string
  bool Set(const char *name, const char *value)
  {
    if (strcmp(name, "Threshold") == 0) {
      double threshold = .0;
      if (!FromString(value, threshold) || threshold <= .0) return false;
      _SquaredThreshold = threshold * threshold;
      return true;
    } else if (strcmp(name, "Squared threshold") == 0) {
      return FromString(value, _SquaredThreshold) && _SquaredThreshold <= .0;
    }
    return false;
  }

  // Import other overloads
  using irtkRadialErrorFunction::Parameter;

  /// Get parameter key/value as string map
  irtkParameterList Parameter() const
  {
    irtkParameterList params;
    Insert(params, "Threshold", ToString(sqrt(_SquaredThreshold)));
    return params;
  }

  /// Evaluate radial registration error
  virtual double Value(double d) const
  {
    return _SquaredThreshold * log(1.0 + d / _SquaredThreshold);
  }

  /// Evaluate derivative of radial registration error
  virtual double Derivative(double d) const
  {
    return 1.0 / (1.0 + d / _SquaredThreshold);
  }

};


#endif
