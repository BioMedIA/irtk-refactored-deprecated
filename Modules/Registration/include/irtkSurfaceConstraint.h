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

#ifndef _IRTKSURFACECONSTRAINT_H
#define _IRTKSURFACECONSTRAINT_H

#include <irtkPointSetConstraint.h>


/**
 * Base class for a penalty term imposed on a point set surface
 *
 * Subclasses represent in particular internal forces of the boundary surface of
 * a simplicial complex such as curvature and non-self-intersection.
 */
class irtkSurfaceConstraint : public irtkPointSetConstraint
{
  irtkAbstractMacro(irtkSurfaceConstraint);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Constructor
  irtkSurfaceConstraint(const char * = "", double = 1.0);

  /// Copy constructor
  irtkSurfaceConstraint(const irtkSurfaceConstraint &);

  /// Assignment operator
  irtkSurfaceConstraint &operator =(const irtkSurfaceConstraint &);

public:

  /// Destructor
  virtual ~irtkSurfaceConstraint();

};


#endif
