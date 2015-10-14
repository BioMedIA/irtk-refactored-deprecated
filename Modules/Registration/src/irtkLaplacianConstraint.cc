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

#include <irtkLaplacianConstraint.h>

#include <irtkEdgeTable.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkIdList.h>
#include <vtkMath.h>

using namespace irtk::polydata;


// =============================================================================
// Auxiliaries
// =============================================================================

namespace irtkLaplacianConstraintUtils {


// -----------------------------------------------------------------------------
/// Evaluate internal force, i.e., apply umbrella operator
struct EvaluateGradient
{
  typedef irtkLaplacianConstraint::GradientType Force;

  vtkPolyData   *_DataSet;
  irtkEdgeTable *_EdgeTable;
  Force         *_Gradient;

  void operator ()(const blocked_range<int> &re) const
  {
    int npts;
    const int *pts;
    double p[3], c[3];

    for (int ptId = re.begin(); ptId != re.end(); ++ptId) {
      _DataSet->GetPoint(ptId, c);
      _EdgeTable->GetAdjacentPoints(ptId, npts, pts);
      if (npts > 0) {
        for (int i = 0; i < npts; ++i) {
          _DataSet->GetPoint(pts[i], p);
          _Gradient[ptId] += Force(c[0] - p[0], c[1] - p[1], c[2] - p[2]);
        }
        _Gradient[ptId] /= npts;
      }
    }
  }
};


} // namespace irtkLaplacianConstraintUtils

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkLaplacianConstraint::irtkLaplacianConstraint(const char *name, double weight)
:
  irtkSurfaceConstraint(name, weight)
{
  _ParameterPrefix.push_back("Surface curvature ");
  _ParameterPrefix.push_back("Surface bending ");
  _ParameterPrefix.push_back("Mesh curvature ");
  _ParameterPrefix.push_back("Mesh bending ");
  _ParameterPrefix.push_back("Surface mesh curvature ");
  _ParameterPrefix.push_back("Surface mesh bending ");
}

// -----------------------------------------------------------------------------
irtkLaplacianConstraint::irtkLaplacianConstraint(const irtkLaplacianConstraint &other)
:
  irtkSurfaceConstraint(other),
  _EdgeTable(other._EdgeTable)
{
}

// -----------------------------------------------------------------------------
irtkLaplacianConstraint &irtkLaplacianConstraint::operator =(const irtkLaplacianConstraint &other)
{
  irtkSurfaceConstraint::operator =(other);
  _EdgeTable = other._EdgeTable;
  return *this;
}

// -----------------------------------------------------------------------------
irtkLaplacianConstraint::~irtkLaplacianConstraint()
{
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void irtkLaplacianConstraint::Initialize()
{
  // Initialize base class
  irtkSurfaceConstraint::Initialize();

  // Initialize this class
  irtkLaplacianConstraint::Init();
}

// -----------------------------------------------------------------------------
void irtkLaplacianConstraint::Reinitialize()
{
  // Reinitialize base class
  irtkSurfaceConstraint::Reinitialize();

  // Reinitialize this class
  irtkLaplacianConstraint::Init();
}

// -----------------------------------------------------------------------------
void irtkLaplacianConstraint::Init()
{
  _EdgeTable.Initialize(_PointSet->InputSurface());
}

// -----------------------------------------------------------------------------
double irtkLaplacianConstraint::Evaluate()
{
  // Use irtkCurvatureConstraint if you want to minimize a well-defined energy
  return numeric_limits<double>::infinity();
}

// -----------------------------------------------------------------------------
void irtkLaplacianConstraint::EvaluateGradient(double *gradient, double step, double weight)
{
  if (_NumberOfPoints == 0) return;

  memset(_Gradient, 0, _NumberOfPoints * sizeof(GradientType));

  irtkLaplacianConstraintUtils::EvaluateGradient eval;
  eval._DataSet   = _PointSet->Surface();
  eval._EdgeTable = &_EdgeTable;
  eval._Gradient  = _Gradient;
  parallel_for(blocked_range<int>(0, _NumberOfPoints), eval);

  irtkSurfaceConstraint::EvaluateGradient(gradient, step, weight / _NumberOfPoints);
}
