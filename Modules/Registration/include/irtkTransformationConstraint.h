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

#ifndef _IRTKTRANSFORMATIONCONSTRAINT_H
#define _IRTKTRANSFORMATIONCONSTRAINT_H

#include <irtkEnergyTerm.h>


/**
 * Base class for a penalty term imposed on a transformation
 *
 * The penalty is minimized by the registration using an instance of
 * irtkRegistrationEnergy with set data similarity and regularization terms.
 * Higher penalties lead to stronger enforcement of the constraint.
 */
class irtkTransformationConstraint : public irtkEnergyTerm
{
  irtkAbstractMacro(irtkTransformationConstraint);

  /// Image domain on which penalty is applied
  irtkPublicAttributeMacro(irtkImageAttributes, Domain);

  /// Whether to apply constraint also at passive DoFs (control points)
  irtkPublicAttributeMacro(bool, ConstrainPassiveDoFs);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
protected:

  /// Constructor
  irtkTransformationConstraint(const char * = "", double = 1.0);

  /// Copy constructor
  irtkTransformationConstraint(const irtkTransformationConstraint &);

  /// Assignment operator
  irtkTransformationConstraint &operator =(const irtkTransformationConstraint &);

public:

  /// Instantiate new transformation constraint of given kind
  static irtkTransformationConstraint *New(irtkConstraintMeasure);

  /// Destructor
  virtual ~irtkTransformationConstraint();

  // ---------------------------------------------------------------------------
  // Configuration

  // Import other overloads
  using irtkEnergyTerm::Parameter;

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  /// Get parameter name/value pairs
  virtual irtkParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Subclass helper
protected:

  /// Get transformation as specific type or NULL if dynamic cast fails
  template <class TransformationType>
  const TransformationType *TransformationAs() const;

  /// Get transformation as free-form deformation or NULL if it is none
  const irtkFreeFormTransformation *FFD() const;

  /// Get transformation as multi-level transformation or NULL if it is none
  const irtkMultiLevelTransformation *MFFD() const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Subclass helper
// =============================================================================

// -----------------------------------------------------------------------------
template <class TransformationType>
inline const TransformationType *irtkTransformationConstraint::TransformationAs() const
{
  return dynamic_cast<const TransformationType *>(Transformation());
}

// -----------------------------------------------------------------------------
inline const irtkFreeFormTransformation *irtkTransformationConstraint::FFD() const
{
  return TransformationAs<irtkFreeFormTransformation>();
}

// -----------------------------------------------------------------------------
inline const irtkMultiLevelTransformation *irtkTransformationConstraint::MFFD() const
{
  return TransformationAs<irtkMultiLevelTransformation>();
}


#endif
