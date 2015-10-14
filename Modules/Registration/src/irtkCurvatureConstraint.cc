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

#include <irtkCurvatureConstraint.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkIdList.h>
#include <vtkMath.h>

using namespace irtk::polydata;


// =============================================================================
// Auxiliaries
// =============================================================================

namespace irtkCurvatureConstraintUtils {


// -----------------------------------------------------------------------------
/// Compute centroids of adjacent nodes
struct ComputeCentroids
{
  vtkPoints           *_Points;
  const irtkEdgeTable *_EdgeTable;
  vtkPoints           *_Centroids;

  void operator ()(const blocked_range<int> &re) const
  {
    double     c[3], p[3];
    const int *adjPtIds;
    int        numAdjPts;

    for (int ptId = re.begin(); ptId != re.end(); ++ptId) {
      _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPtIds);
      if (numAdjPts > 0) {
        c[0] = c[1] = c[2] = .0;
        for (int i = 0; i < numAdjPts; ++i) {
          _Points->GetPoint(adjPtIds[i], p);
          c[0] += p[0], c[1] += p[1], c[2] += p[2];
        }
        c[0] /= numAdjPts, c[1] /= numAdjPts, c[2] /= numAdjPts;
      } else {
        _Points->GetPoint(ptId, c);
      }
      _Centroids->SetPoint(ptId, c);
    }
  }
};

// -----------------------------------------------------------------------------
/// Evaluate bending penalty
struct Evaluate
{
  vtkPoints           *_Points;
  const irtkEdgeTable *_EdgeTable;
  vtkPoints           *_Centroids;
  vtkDataArray        *_Normals;
  double               _ConcavityWeight;
  double               _ConvexityWeight;
  double               _Sum;

  Evaluate() : _Sum(.0) {}

  Evaluate(const Evaluate &other, split)
  :
    _Points(other._Points),
    _EdgeTable(other._EdgeTable),
    _Centroids(other._Centroids),
    _Normals(other._Normals),
    _ConcavityWeight(other._ConcavityWeight),
    _ConvexityWeight(other._ConvexityWeight),
    _Sum(.0)
  {}

  void join(const Evaluate &other)
  {
    _Sum += other._Sum;
  }

  void operator ()(const blocked_range<int> &re)
  {
    double w, p[3], n[3], c[3];

    // Either smooth concave and convex regions equally strong
    if (fequal(_ConvexityWeight, _ConcavityWeight)) {

      w = _ConvexityWeight;
      for (int ptId = re.begin(); ptId != re.end(); ++ptId) {
        _Points->GetPoint(ptId, p);
        _Centroids->GetPoint(ptId, c);
        _Sum += w * vtkMath::Distance2BetweenPoints(c, p);
      }

    // or smooth either concave or convex regions more strongly
    } else {

      for (int ptId = re.begin(); ptId != re.end(); ++ptId) {
        _Points->GetPoint(ptId, p);
        _Normals->GetTuple(ptId, n);
        _Centroids->GetPoint(ptId, c);
        vtkMath::Subtract(c, p, p);     // p = c - p
        w = (vtkMath::Dot(n, p) < .0 ? _ConvexityWeight : _ConcavityWeight);
        _Sum += w * vtkMath::Dot(p, p); // |c - p|^2
      }

    }
  }
};

// -----------------------------------------------------------------------------
/// Evaluate gradient of bending penalty (i.e., negative internal bending force)
struct EvaluateGradient
{
  typedef irtkCurvatureConstraint::GradientType Force;

  vtkPoints           *_Points;
  vtkPoints           *_Centroids;
  vtkDataArray        *_Normals;
  const irtkEdgeTable *_EdgeTable;
  Force               *_Gradient;
  double               _ConvexityWeight;
  double               _ConcavityWeight;

  void operator ()(const blocked_range<int> &re) const
  {
    const int *adjPtIds;
    int        numAdjPts;
    double     p[3], c[3], n[3], w, W;

    // Either smooth concave and convex regions equally strong
    if (fequal(_ConvexityWeight, _ConcavityWeight)) {

      for (int ptId = re.begin(); ptId != re.end(); ++ptId) {
        // Derivative of sum terms of adjacent points
        _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPtIds);
        if (numAdjPts > 0) {
          for (int i = 0; i < numAdjPts; ++i) {
            _Points->GetPoint(adjPtIds[i], p);
            _Centroids->GetPoint(adjPtIds[i], c);
            w = 1.0 / _EdgeTable->NumberOfAdjacentPoints(adjPtIds[i]);
            _Gradient[ptId] += w * Force(c[0] - p[0], c[1] - p[1], c[2] - p[2]);
          }
        }
        // Derivative of sum term of this point
        _Points->GetPoint(ptId, p);
        _Centroids->GetPoint(ptId, c);
        _Gradient[ptId] -= Force(c[0] - p[0], c[1] - p[1], c[2] - p[2]);
      }

    // or smooth either concave or convex regions more strongly
    } else {

      for (int ptId = re.begin(); ptId != re.end(); ++ptId) {
        // Derivative of sum terms of adjacent points
        _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPtIds);
        if (numAdjPts > 0) {
          W = .0;
          for (int i = 0; i < numAdjPts; ++i) {
            _Points->GetPoint(adjPtIds[i], p);
            _Normals->GetTuple(adjPtIds[i], n);
            _Centroids->GetPoint(adjPtIds[i], c);
            vtkMath::Subtract(c, p, p);
            W += (w = (vtkMath::Dot(n, p) < .0 ? _ConvexityWeight : _ConcavityWeight));
            _Gradient[ptId] += w * Force(p[0], p[1], p[2]);
          }
          _Gradient[ptId] /= W;
        }
        // Derivative of sum term of this point
        _Points->GetPoint(ptId, p);
        _Normals->GetTuple(ptId, n);
        _Centroids->GetPoint(ptId, c);
        vtkMath::Subtract(c, p, p);
        w = (vtkMath::Dot(n, p) < .0) ? _ConvexityWeight : _ConcavityWeight;
        _Gradient[ptId] -= w * Force(p[0], p[1], p[2]);
      }

    }
  }
};


} // namespace irtkCurvatureConstraintUtils

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
void irtkCurvatureConstraint::Copy(const irtkCurvatureConstraint &other)
{
  if (other._Centroids) {
    if (!_Centroids) _Centroids = vtkSmartPointer<vtkPoints>::New();
    _Centroids->DeepCopy(other._Centroids);
  } else {
    _Centroids = NULL;
  }
}

// -----------------------------------------------------------------------------
irtkCurvatureConstraint::irtkCurvatureConstraint(const char *name, double weight)
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
irtkCurvatureConstraint::irtkCurvatureConstraint(const irtkCurvatureConstraint &other)
:
  irtkSurfaceConstraint(other)
{
  Copy(other);
}

// -----------------------------------------------------------------------------
irtkCurvatureConstraint &irtkCurvatureConstraint::operator =(const irtkCurvatureConstraint &other)
{
  if (this != &other) {
    irtkSurfaceConstraint::operator =(other);
    Copy(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
irtkCurvatureConstraint::~irtkCurvatureConstraint()
{
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void irtkCurvatureConstraint::Initialize()
{
  // Initialize base class
  irtkSurfaceConstraint::Initialize();

  // Initialize this class
  irtkCurvatureConstraint::Init();
}

// -----------------------------------------------------------------------------
void irtkCurvatureConstraint::Reinitialize()
{
  // Reinitialize base class
  irtkSurfaceConstraint::Reinitialize();

  // Reinitialize this class
  irtkCurvatureConstraint::Init();
}

// -----------------------------------------------------------------------------
void irtkCurvatureConstraint::Init()
{
  if (_Centroids == NULL) _Centroids = vtkSmartPointer<vtkPoints>::New();
  _Centroids->SetNumberOfPoints(_NumberOfPoints);
}

// -----------------------------------------------------------------------------
void irtkCurvatureConstraint::Update(bool gradient)
{
  // Update base class
  irtkSurfaceConstraint::Update(gradient);

  // Update centroids
  IRTK_START_TIMING();
  irtkCurvatureConstraintUtils::ComputeCentroids eval;
  eval._Points    = _PointSet->SurfacePoints();
  eval._EdgeTable = _PointSet->SurfaceEdges();
  eval._Centroids = _Centroids;
  parallel_for(blocked_range<int>(0, _NumberOfPoints), eval);
  IRTK_DEBUG_TIMING(3, "update of curvature centroids");
}

// -----------------------------------------------------------------------------
double irtkCurvatureConstraint::Evaluate()
{
  if (_NumberOfPoints == 0) return .0;
  IRTK_START_TIMING();
  irtkCurvatureConstraintUtils::Evaluate eval;
  eval._Points          = _PointSet->SurfacePoints();
  eval._EdgeTable       = _PointSet->SurfaceEdges();
  eval._Centroids       = _Centroids;
  eval._ConvexityWeight = 1.0;
  eval._ConcavityWeight = 1.0;
  parallel_reduce(blocked_range<int>(0, _NumberOfPoints), eval);
  IRTK_DEBUG_TIMING(3, "evaluation of curvature penalty");
  return eval._Sum / _NumberOfPoints;
}

// -----------------------------------------------------------------------------
void irtkCurvatureConstraint::EvaluateGradient(double *gradient, double step, double weight)
{
  if (_NumberOfPoints == 0) return;

  IRTK_START_TIMING();
  memset(_Gradient, 0, _NumberOfPoints * sizeof(GradientType));

  irtkCurvatureConstraintUtils::EvaluateGradient eval;
  eval._Points          = _PointSet->SurfacePoints();
  eval._EdgeTable       = _PointSet->SurfaceEdges();
  eval._Centroids       = _Centroids;
  eval._ConvexityWeight = 1.0;
  eval._ConcavityWeight = 1.0;
  eval._Gradient        = _Gradient;
  parallel_for(blocked_range<int>(0, _NumberOfPoints), eval);

  irtkSurfaceConstraint::EvaluateGradient(gradient, step, weight / _NumberOfPoints);
  IRTK_DEBUG_TIMING(3, "evaluation of curvature force");
}
