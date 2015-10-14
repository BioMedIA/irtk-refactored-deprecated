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

#ifndef _IRTKTOPOLOGYPRESERVATIONCONSTRAINT_H
#define _IRTKTOPOLOGYPRESERVATIONCONSTRAINT_H

#include <irtkTransformationConstraint.h>


/**
 * Topology preservation constraint for deformable image registration
 *
 */
class irtkTopologyPreservationConstraint : public irtkTransformationConstraint
{
  irtkObjectMacro(irtkTopologyPreservationConstraint);

public:

  /// Constructor
  irtkTopologyPreservationConstraint(const char * = "");

protected:

  /// Compute penalty for current transformation estimate
  double Evaluate();

  /// Compute gradient of penalty term w.r.t transformation parameters
  void EvaluateGradient(double *, double, double);

};


#endif
