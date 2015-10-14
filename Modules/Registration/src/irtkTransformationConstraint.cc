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

#include <irtkTransformationConstraint.h>

#include <irtkSmoothnessConstraint.h>
#include <irtkMinJacobianConstraint.h>
#include <irtkLogJacobianConstraint.h>
#include <irtkVolumePreservationConstraint.h>
#include <irtkTopologyPreservationConstraint.h>
#include <irtkSparsityConstraint.h>


// =============================================================================
// Factory method
// =============================================================================

// -----------------------------------------------------------------------------
irtkTransformationConstraint *irtkTransformationConstraint::New(irtkConstraintMeasure c)
{
  switch (c) {
    case CM_Smoothness:           return new irtkSmoothnessConstraint();
    case CM_MinDetJac:            return new irtkMinJacobianConstraint();
    case CM_SqLogDetJac:          return new irtkLogJacobianConstraint();
    case CM_VolumePreservation:   return new irtkVolumePreservationConstraint();
    case CM_TopologyPreservation: return new irtkTopologyPreservationConstraint();
    case CM_Sparsity:             return new irtkSparsityConstraint();
    case CM_BendingEnergy:        return new irtkSmoothnessConstraint();
    default: {
      cerr << "irtkTransformationConstraint::New: Unknown transformation constraint: " << c << endl;
      exit(1);
    }
  }
}

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkTransformationConstraint::irtkTransformationConstraint(const char *name, double weight)
:
  irtkEnergyTerm(name, weight),
  _ConstrainPassiveDoFs(false)
{
}

// -----------------------------------------------------------------------------
irtkTransformationConstraint::irtkTransformationConstraint(const irtkTransformationConstraint &other)
:
  irtkEnergyTerm(other),
  _ConstrainPassiveDoFs(other._ConstrainPassiveDoFs)
{
}

// -----------------------------------------------------------------------------
irtkTransformationConstraint &irtkTransformationConstraint::operator =(const irtkTransformationConstraint &other)
{
  irtkEnergyTerm::operator =(other);
  _ConstrainPassiveDoFs = other._ConstrainPassiveDoFs;
  return *this;
}

// -----------------------------------------------------------------------------
irtkTransformationConstraint::~irtkTransformationConstraint()
{
}

// =============================================================================
// Configuration
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkTransformationConstraint::Set(const char *param, const char *value)
{
  const string name = ParameterNameWithoutPrefix(param);
  if (strcmp(param, "Constrain passive DoFs")           == 0 ||
      strcmp(param, "Constrain passive CPs")            == 0 ||
      strcmp(param, "Constrain passive control points") == 0 ||
      strcmp(param, "Constrain passive parameters")     == 0 ||
      name == "Of passive DoFs" || name == "Of passive parameters" ||
      name == "At passive control points") {
    return FromString(value, _ConstrainPassiveDoFs);
  }
  return irtkEnergyTerm::Set(param, value);
}

// -----------------------------------------------------------------------------
irtkParameterList irtkTransformationConstraint::Parameter() const
{
  irtkParameterList params = irtkEnergyTerm::Parameter();
  if (!_Name.empty()) {
    Insert(params, _Name + " at passive control points", ToString(_ConstrainPassiveDoFs));
  }
  return params;
}
