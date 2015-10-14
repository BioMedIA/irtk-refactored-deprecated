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

#ifndef _IRTKEXTERNALFORCE_H
#define _IRTKEXTERNALFORCE_H

#include <irtkPointSetForce.h>
#include <irtkRegisteredImage.h>


/**
 * Base class for external point set force terms
 *
 * Subclasses implement in particular external forces for deformable surface
 * models such as inflation/balloon forces and intensity edge forces.
 */
class irtkExternalForce : public irtkPointSetForce
{
  irtkAbstractMacro(irtkExternalForce);

  // ---------------------------------------------------------------------------
  // Attributes
protected:

  /// (Transformed) reference image
  irtkPublicAggregateMacro(irtkRegisteredImage, Image);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
protected:

  /// Constructor
  irtkExternalForce(const char * = "", double = 1.0);

  /// Copy constructor
  irtkExternalForce(const irtkExternalForce &);

  /// Assignment operator
  irtkExternalForce &operator =(const irtkExternalForce &);

  /// Copy attributes of this class from another instance
  void Copy(const irtkExternalForce &);

public:

  /// Instantiate specified external force
  static irtkExternalForce *New(irtkPointSetForceMeasure);

  /// Destructor
  virtual ~irtkExternalForce();

  // ---------------------------------------------------------------------------
  // Initialization
public:

  /// Initialize external force once input and parameters have been set
  virtual void Initialize();

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Update moving input points and internal state of external force term
  virtual void Update(bool = true);

protected:

  /// Evaluate external force term
  virtual double Evaluate();

};


#endif
