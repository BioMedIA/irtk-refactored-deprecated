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

#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkStructuredGrid.h>

#include <irtkPointSetConstraint.h>

#include <irtkEdgeLengthConstraint.h>
#include <irtkCurvatureConstraint.h>
#include <irtkNonSelfIntersectionConstraint.h>
#include <irtkInflationForce.h>


// =============================================================================
// Factory method
// =============================================================================

// -----------------------------------------------------------------------------
irtkPointSetConstraint *irtkPointSetConstraint::New(irtkPointSetConstraintMeasure c, const char *name, double weight)
{
  switch (c) {
    case PCM_Stretching:          return new irtkEdgeLengthConstraint(name, weight);
    case PCM_Curvature:           return new irtkCurvatureConstraint(name, weight);
    case PCM_NonSelfIntersection: return new irtkNonSelfIntersectionConstraint(name, weight);
    case PCM_Inflation:           return new irtkInflationForce(name, weight);
    default: {
      cerr << "irtkPolyDataConstraint::New: Unknown polydata constraint: " << c << endl;
      exit(1);
    }
  }
}

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkPointSetConstraint::irtkPointSetConstraint(const char *name, double weight)
:
  irtkPointSetForce(name, weight)
{
}

// -----------------------------------------------------------------------------
irtkPointSetConstraint::irtkPointSetConstraint(const irtkPointSetConstraint &other)
:
  irtkPointSetForce(other)
{
}

// -----------------------------------------------------------------------------
irtkPointSetConstraint &irtkPointSetConstraint::operator =(const irtkPointSetConstraint &other)
{
  irtkPointSetForce::operator =(other);
  return *this;
}

// -----------------------------------------------------------------------------
irtkPointSetConstraint::~irtkPointSetConstraint()
{
}
