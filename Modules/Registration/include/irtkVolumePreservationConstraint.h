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

#ifndef _IRTKVOLUMEPRESERVATIONCONSTRAINT_H
#define _IRTKVOLUMEPRESERVATIONCONSTRAINT_H

#include <irtkLogJacobianConstraint.h>


/**
 * Volume preservation constraint for deformable image registration
 *
 * The volume preservation constraint is based on a penalty imposed by the
 * squared logarithm of the Jacobian determinant. Unlike the base class
 * implementation which computes the Jacobian w.r.t the spline function,
 * this constraint computes the Jacobian always w.r.t the transformed
 * coordinates. These two derivative computations are identical for a
 * FFD parameterized by control point displacements, but differ in case of
 * a stationary velocity field parameterization.
 *
 * Torsten Rohlﬁng and Calvin R. Maurer, Jr., Intensity-Based Non-rigid
 * Registration Using Adaptive Multilevel Free-Form Deformation with an
 * Incompressibility Constraint, MICCAI 2001.
 */
class irtkVolumePreservationConstraint : public irtkLogJacobianConstraint
{
  irtkObjectMacro(irtkVolumePreservationConstraint);

public:

  /// Constructor
  irtkVolumePreservationConstraint(const char * name = "")
  :
    irtkLogJacobianConstraint(name)
  {}

  /// Destructor
  virtual ~irtkVolumePreservationConstraint() {}

  /// Compute determinant and adjugate of Jacobian of transformation
  virtual double Jacobian(const irtkFreeFormTransformation *ffd,
                          double x, double y, double z, double t,
                          irtkMatrix &adj) const
  {
    double det;
    ffd->LocalJacobian(adj, x, y, z, t, t);
    adj.Adjugate(det);
    if (det < 1e-7) det = 1e-7;
    return det;
  }

};


#endif
