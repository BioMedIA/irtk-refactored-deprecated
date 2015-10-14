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

#include <irtkInflationForce.h>

#include <vtkPolyData.h>
#include <vtkMath.h>

using namespace irtk::polydata;


// =============================================================================
// Auxiliaries
// =============================================================================

namespace irtkInflationForceUtils {


// -----------------------------------------------------------------------------
/// Evaluate penalty
struct Evaluate
{
  vtkPoints           *_Points;
  vtkDataArray        *_Normals;
  const irtkEdgeTable *_EdgeTable;
  double               _ConcavityWeight;
  double               _ConvexityWeight;
  vtkDataArray        *_Norm;
  double               _Sum;

  Evaluate() : _Sum(.0) {}

  Evaluate(const Evaluate &other, split)
  :
    _Points(other._Points),
    _Normals(other._Normals),
    _EdgeTable(other._EdgeTable),
    _ConcavityWeight(other._ConcavityWeight),
    _ConvexityWeight(other._ConvexityWeight),
    _Norm(other._Norm),
    _Sum(.0)
  {}

  void join(const Evaluate &other)
  {
    _Sum += other._Sum;
  }

  void operator ()(const blocked_range<int> &re)
  {
    double     p1[3], p2[3], s;
    int        numAdjPts;
    const int *adjPtIds;

    if (fequal(_ConvexityWeight, _ConcavityWeight)) {

      for (int ptId = re.begin(); ptId != re.end(); ++ptId) {
        _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPtIds);
        if (numAdjPts > 0) {
          s = .0;
          _Points->GetPoint(ptId, p1);
          for (int i = 0; i < numAdjPts; ++i) {
            _Points->GetPoint(adjPtIds[i], p2);
            s += vtkMath::Distance2BetweenPoints(p1, p2);
          }
          _Sum += s / numAdjPts;
        }
      }

    } else {

      double n1[3], w, W;

      for (int ptId = re.begin(); ptId != re.end(); ++ptId) {
        s = W = .0;
        _Points ->GetPoint(ptId, p1);
        _Normals->GetTuple(ptId, n1);
        _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPtIds);
        for (int i = 0; i < numAdjPts; ++i) {
          _Points ->GetPoint(adjPtIds[i], p2);
          vtkMath::Subtract(p2, p1, p2);
          w = vtkMath::Dot(p2, n1);
          w = (w < .0 ? _ConvexityWeight * -w : _ConcavityWeight * w);
          s += w * vtkMath::Dot(p2, p2);
          W += w;
        }
        if (W != .0) _Sum += s / W;
        _Norm->SetComponent(ptId, 0, W);
      }

    }
  }
};

// -----------------------------------------------------------------------------
/// Evaluate gradient of penalty term (i.e., negative inflation force)
struct EvaluateGradient
{
  typedef irtkInflationForce::GradientType Force;

  vtkPoints           *_Points;
  vtkDataArray        *_Normals;
  const irtkEdgeTable *_EdgeTable;
  Force               *_Gradient;
  double               _ConvexityWeight;
  double               _ConcavityWeight;
  vtkDataArray        *_Norm;

  void operator ()(const blocked_range<int> &re) const
  {
    const int *adjPtIds;
    int        numAdjPts;
    double     p1[3], p2[3];
    Force      f1, f2;

    if (fequal(_ConvexityWeight, _ConcavityWeight)) {

      for (int ptId = re.begin(); ptId != re.end(); ++ptId) {
        _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPtIds);
        if (numAdjPts > 0) {
          _Points->GetPoint(ptId, p1);
          f1 = f2 = .0;
          for (int i = 0; i < numAdjPts; ++i) {
            _Points->GetPoint(adjPtIds[i], p2);
            vtkMath::Subtract(p2, p1, p2);
            f1 -= Force(p2[0], p2[1], p2[2]);
            f2 -= Force(p2[0], p2[1], p2[2]) / _EdgeTable->NumberOfAdjacentPoints(adjPtIds[i]);
          }
          f1 /= numAdjPts;
          _Gradient[ptId]  = f1;
          _Gradient[ptId] += f2;
          _Gradient[ptId] *= 2.0;
        }
      }

    } else {

      double w1, w2, W1, W2, n1[3], n2[3];
      Force  f1, f2;

      for (int ptId = re.begin(); ptId != re.end(); ++ptId) {
        f1 = f2 = .0;
        W1 = _Norm->GetComponent(ptId, 0);
        _Points ->GetPoint(ptId, p1);
        _Normals->GetTuple(ptId, n1);
        _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPtIds);
        for (int i = 0; i < numAdjPts; ++i) {
          _Points ->GetPoint(adjPtIds[i], p2);
          _Normals->GetTuple(adjPtIds[i], n2);
          vtkMath::Subtract(p2, p1, p2);
          w1 =  vtkMath::Dot(p2, n1); // = <p2 - p1, n1>
          w2 = -vtkMath::Dot(p2, n2); // = <p1 - p2, n2>
          w1 = (w1 < .0 ? _ConvexityWeight * -w1 : _ConcavityWeight * w1);
          w2 = (w2 < .0 ? _ConvexityWeight * -w2 : _ConcavityWeight * w2);
          W2 = _Norm->GetComponent(adjPtIds[i], 0);
          if (W1 > .0) f1 -= (w1 / W1) * Force(p2[0], p2[1], p2[2]);
          if (W2 > .0) f2 -= (w2 / W2) * Force(p2[0], p2[1], p2[2]);
        }
        _Gradient[ptId]  = f1;
        _Gradient[ptId] += f2;
        _Gradient[ptId] *= 2.0;
      }

//      double min, max, e[3];
//
//      for (int ptId = re.begin(); ptId != re.end(); ++ptId) {
//        min = max = .0;
//        _Points ->GetPoint(ptId, p1);
//        _Normals->GetTuple(ptId, n1);
//        _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPtIds);
//        for (int i = 0; i < numAdjPts; ++i) {
//          _Points ->GetPoint(adjPtIds[i], p2);
////          _Normals->GetTuple(adjPtIds[i], n2);
//          vtkMath::Subtract(p2, p1, e);
//          vtkMath::Normalize(e);
//          double dp = vtkMath::Dot(e, n1);
//          if (min > dp) min = dp;
//          if (max < dp) max = dp;
//        }
//        if (min < -.6 && max > .3) {
//          _Gradient[ptId] = - max * Force(n1[0], n1[1], n1[2]);
//        }
////        _Gradient[ptId]  = f1;
////        _Gradient[ptId] += f2;
////        _Gradient[ptId] *= 2.0;
//      }

    }
  }
};


} // namespace irtkInflationForceUtils

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkInflationForce::irtkInflationForce(const char *name, double weight)
:
  irtkSurfaceConstraint(name, weight),
  _ConvexityWeight(1.0),
  _ConcavityWeight(1.0),
  _EnergyValue(numeric_limits<double>::quiet_NaN())
{
  _ParameterPrefix.push_back("Surface inflation ");
  _ParameterPrefix.push_back("Mesh inflation ");
  _ParameterPrefix.push_back("Surface mesh inflation ");
}

// -----------------------------------------------------------------------------
void irtkInflationForce::Copy(const irtkInflationForce &other)
{
  _ConvexityWeight = other._ConvexityWeight;
  _ConcavityWeight = other._ConcavityWeight;
  if (other._Norm) {
    if (!_Norm) _Norm = other._Norm->NewInstance();
    _Norm->DeepCopy(other._Norm);
  } else {
    _Norm = NULL;
  }
  _EnergyValue = other._EnergyValue;
}

// -----------------------------------------------------------------------------
irtkInflationForce::irtkInflationForce(const irtkInflationForce &other)
:
  irtkSurfaceConstraint(other)
{
  Copy(other);
}

// -----------------------------------------------------------------------------
irtkInflationForce &irtkInflationForce::operator =(const irtkInflationForce &other)
{
  if (this != &other) {
    irtkSurfaceConstraint::operator =(other);
    Copy(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
irtkInflationForce::~irtkInflationForce()
{
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void irtkInflationForce::Initialize()
{
  // Initialize base class
  irtkSurfaceConstraint::Initialize();

  // Initialize this class
  irtkInflationForce::Init();
}

// -----------------------------------------------------------------------------
void irtkInflationForce::Reinitialize()
{
  // Reinitialize base class
  irtkSurfaceConstraint::Reinitialize();

  // Reinitialize this class
  irtkInflationForce::Init();
}

// -----------------------------------------------------------------------------
void irtkInflationForce::Init()
{
  if (fequal(_ConvexityWeight, _ConcavityWeight)) {
    _Norm = NULL;
  } else {
    if (_Norm == NULL) {
      _Norm = vtkSmartPointer<vtkFloatArray>::New();
      _Norm->SetNumberOfComponents(1);
    }
    _Norm->SetNumberOfTuples(_NumberOfPoints);
  }
  _EnergyValue = numeric_limits<double>::quiet_NaN();
}

// -----------------------------------------------------------------------------
void irtkInflationForce::Update(bool gradient)
{
  // Update base class
  irtkSurfaceConstraint::Update(gradient);

  // Reset cached energy value
  _EnergyValue = numeric_limits<double>::quiet_NaN();
  if (_NumberOfPoints == 0) return;

  // Compute sum of convexity/concavity weights
  //
  // This computation is a side effect of the evaluation of the energy of the
  // inflation force term, which is thus cached here to avoid its recomputation.
  if (!fequal(_ConvexityWeight, _ConcavityWeight)) this->Evaluate();
}

// -----------------------------------------------------------------------------
double irtkInflationForce::Evaluate()
{
  if (_NumberOfPoints == 0) return .0;
  if (IsNaN(_EnergyValue)) {
    IRTK_START_TIMING();
    irtkInflationForceUtils::Evaluate eval;
    eval._Points          = _PointSet->SurfacePoints();
    eval._Normals         = _PointSet->SurfaceNormals();
    eval._EdgeTable       = _PointSet->SurfaceEdges();
    eval._ConvexityWeight = _ConvexityWeight;
    eval._ConcavityWeight = _ConcavityWeight;
    eval._Norm            = _Norm;
    parallel_reduce(blocked_range<int>(0, _NumberOfPoints), eval);
    IRTK_DEBUG_TIMING(3, "evaluation of inflation energy");
    _EnergyValue = eval._Sum / _NumberOfPoints;
  }
  return _EnergyValue;
}

// -----------------------------------------------------------------------------
void irtkInflationForce::EvaluateGradient(double *gradient, double step, double weight)
{
  if (_NumberOfPoints == 0) return;

  IRTK_START_TIMING();
  memset(_Gradient, 0, _NumberOfPoints * sizeof(GradientType));

  irtkInflationForceUtils::EvaluateGradient eval;
  eval._Points          = _PointSet->SurfacePoints();
  eval._Normals         = _PointSet->SurfaceNormals();
  eval._EdgeTable       = _PointSet->SurfaceEdges();
  eval._ConvexityWeight = _ConvexityWeight;
  eval._ConcavityWeight = _ConcavityWeight;
  eval._Gradient        = _Gradient;
  eval._Norm            = _Norm;
  parallel_for(blocked_range<int>(0, _NumberOfPoints), eval);

  irtkSurfaceConstraint::EvaluateGradient(gradient, step, weight / _NumberOfPoints);
  IRTK_DEBUG_TIMING(3, "evaluation of inflation force");
}
