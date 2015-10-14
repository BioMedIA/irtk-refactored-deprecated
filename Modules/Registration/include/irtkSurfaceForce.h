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

#ifndef _IRTKSURFACEFORCE_H
#define _IRTKSURFACEFORCE_H

#include <irtkExternalForce.h>


/**
 * Base class for external surface force terms
 *
 * Subclasses implement in particular external forces for deformable surface
 * models such as inflation/balloon forces and intensity edge forces. In case
 * of a tetrahedral mesh, these forces only apply to the corresponding
 * triangular surface of the simplicial complex boundary.
 */
class irtkSurfaceForce : public irtkExternalForce
{
  irtkAbstractMacro(irtkSurfaceForce);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
protected:

  /// Constructor
  irtkSurfaceForce(const char * = "", double = 1.0);

  /// Copy constructor
  irtkSurfaceForce(const irtkSurfaceForce &);

  /// Assignment operator
  irtkSurfaceForce &operator =(const irtkSurfaceForce &);

  /// Copy attributes of this class from another instance
  void Copy(const irtkSurfaceForce &);

public:

  /// Destructor
  virtual ~irtkSurfaceForce();

};


#endif
