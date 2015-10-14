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

#ifndef _IRTKPOINTSETCONSTRAINT_H
#define _IRTKPOINTSETCONSTRAINT_H

#include <irtkPointSetForce.h>


/**
 * Base class for a penalty term imposed on a transformed point set
 *
 * Subclasses represent in particular internal forces on deformable simplicial
 * complexes such as elasticity, strain, curvature, and non-self-intersection.
 *
 * The penalty is minimized by the registration using an instance of
 * irtkRegistrationEnergy with set data similarity and regularization terms.
 * Higher penalties lead to stronger enforcement of the constraint.
 */
class irtkPointSetConstraint : public irtkPointSetForce
{
  irtkAbstractMacro(irtkPointSetConstraint);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Constructor
  irtkPointSetConstraint(const char * = "", double = 1.0);

  /// Copy constructor
  irtkPointSetConstraint(const irtkPointSetConstraint &);

  /// Assignment operator
  irtkPointSetConstraint &operator =(const irtkPointSetConstraint &);

public:

  /// Instantiate new constraint representing specified internal forces
  static irtkPointSetConstraint *New(irtkPointSetConstraintMeasure,
                                     const char * = "", double = 1.0);

  /// Destructor
  virtual ~irtkPointSetConstraint();

};


#endif
