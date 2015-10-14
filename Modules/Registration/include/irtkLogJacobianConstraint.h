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

#ifndef _IRTKLOGJACOBIANCONSTRAINT_H
#define _IRTKLOGJACOBIANCONSTRAINT_H

#include <irtkJacobianConstraint.h>


/**
 * Constrains log of Jacobian determinant of FFD parameterization
 *
 * This constraint prevents folding of the transformation parameterization,
 * i.e., either of the control point displacements or velocities. It preserves
 * volume in case of a classical FFD model and is in this case equivalent to
 * the irtkVolumePreservationConstraint.
 *
 * Torsten Rohlﬁng and Calvin R. Maurer, Jr., Intensity-Based Non-rigid
 * Registration Using Adaptive Multilevel Free-Form Deformation with an
 * Incompressibility Constraint, MICCAI 2001.
 *
 * Modat et al., Log-Euclidean free-form deformation, SPIE 2011.
 *
 * @sa irtkVolumePreservationConstraint
 */
class irtkLogJacobianConstraint : public irtkJacobianConstraint
{
  irtkObjectMacro(irtkLogJacobianConstraint);

public:

  /// Constructor
  irtkLogJacobianConstraint(const char * = "");

  /// Destructor
  virtual ~irtkLogJacobianConstraint();

  /// Compute determinant and adjugate of Jacobian of transformation
  virtual double Jacobian(const irtkFreeFormTransformation *ffd,
                          double x, double y, double z, double t,
                          irtkMatrix &adj) const
  {
    double det;
    ffd->FFDJacobianWorld(adj, x, y, z, t, t);
    adj.Adjugate(det);
    if (det < 1e-7) det = 1e-7;
    return det;
  }

protected:

  /// Compute penalty for current transformation estimate
  virtual double Evaluate();

  /// Compute gradient of penalty term w.r.t transformation parameters
  virtual void EvaluateGradient(double *, double, double);

};


#endif
