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

#ifndef _IRTKEULERMETHODWITHMOMENTUM_H
#define _IRTKEULERMETHODWITHMOMENTUM_H

#include <irtkEulerMethod.h>


/**
 * Minimizes deformable surface model using Euler's method with momentum
 *
 * This method has for example been used in the following paper on which this
 * implementation is also based on:
 *
 *   Park et al. (2001), A non-self-intersecting adaptive deformable surface for
 *   complex boundary extraction from volumetric images,
 *   Computer & Graphics, 25, 421–440
 */
class irtkEulerMethodWithMomentum : public irtkEulerMethod
{
  irtkOptimizerMacro(irtkEulerMethodWithMomentum, EulerMethodWithMomentum);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Common mass of all nodes
  irtkPublicAttributeMacro(double, BodyMass);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  irtkEulerMethodWithMomentum(irtkObjectiveFunction * = NULL);

  /// Copy constructor
  irtkEulerMethodWithMomentum(const irtkEulerMethodWithMomentum &);

  /// Assignment operator
  irtkEulerMethodWithMomentum &operator =(const irtkEulerMethodWithMomentum &);

  /// Destructor
  virtual ~irtkEulerMethodWithMomentum();

  // ---------------------------------------------------------------------------
  // Parameters
  using irtkLocalOptimizer::Parameter;

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  /// Get parameters as key/value as string map
  virtual irtkParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Execution

  /// Compute node velocity given its last velocity and the sum of node forces
  virtual void ComputeVelocity(double *, const double *) const;

};


#endif
