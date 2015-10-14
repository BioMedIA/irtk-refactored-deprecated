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

#ifndef IRTKSQUAREDERRORFUNCTION_H
#define IRTKSQUAREDERRORFUNCTION_H

#include <irtkRadialErrorFunction.h>


/**
 * Squared Euclidean distance registration error function
 */
class irtkSquaredErrorFunction : public irtkRadialErrorFunction
{
  irtkObjectMacro(irtkSquaredErrorFunction);

public:

  /// Constructor
  irtkSquaredErrorFunction() {}

  /// Copy constructor
  irtkSquaredErrorFunction(const irtkSquaredErrorFunction &) {}

  /// Destructor
  ~irtkSquaredErrorFunction() {}

  /// Copy construct a new instance
  virtual irtkRadialErrorFunction *NewInstance() const
  {
    return new irtkSquaredErrorFunction(*this);
  }

  /// Type enumeration value
  virtual TypeId Type() const
  {
    return Squared;
  }

  /// Evaluate radial registration error
  virtual double Value(double d) const
  {
    return d;
  }

  /// Evaluate derivative of radial registration error
  virtual double Derivative(double d) const
  {
    return 1.0;
  }

};


#endif
