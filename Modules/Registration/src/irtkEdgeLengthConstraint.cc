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

#include <irtkEdgeLengthConstraint.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkExtractEdges.h>
#include <vtkIdList.h>
#include <vtkMath.h>

#include <irtkPolyDataUtils.h>
using namespace irtk::polydata;


// =============================================================================
// Auxiliaries
// =============================================================================

namespace irtkEdgeLengthConstraintUtils {


// -----------------------------------------------------------------------------
/// Evaluate stretching penalty
struct Evaluate
{
  vtkPoints     *_Points;
  vtkPoints     *_InitialPoints;
  irtkEdgeTable *_EdgeTable;
  double         _RestLength;
  double         _Sum;

  Evaluate() : _Sum(.0) {}

  Evaluate(const Evaluate &other)
  :
    _Points(other._Points),
    _InitialPoints(other._InitialPoints),
    _EdgeTable(other._EdgeTable),
    _RestLength(other._RestLength),
    _Sum(.0)
  {}

  Evaluate(const Evaluate &other, split)
  :
    _Points(other._Points),
    _InitialPoints(other._InitialPoints),
    _EdgeTable(other._EdgeTable),
    _RestLength(other._RestLength),
    _Sum(.0)
  {}

  void join(const Evaluate &other)
  {
    _Sum += other._Sum;
  }

  void operator ()(const blocked_range<int> &re)
  {
    int    ptId1, ptId2, edgeId;
    double p1[3], p2[3], d, d0 = _RestLength;

    irtkEdgeIterator it(*_EdgeTable);
    for (it.InitTraversal(re); (edgeId = it.GetNextEdge(ptId1, ptId2) != -1);) {
      if (_InitialPoints) {
        _InitialPoints->GetPoint(ptId1, p1);
        _InitialPoints->GetPoint(ptId2, p2);
        d0 = sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
      }
      _Points->GetPoint(ptId1, p1);
      _Points->GetPoint(ptId2, p2);
      d = sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
      _Sum += pow(d - d0, 2);
    }
  }
};

// -----------------------------------------------------------------------------
/// Evaluate gradient of penalty (i.e., negative internal stretching force)
struct EvaluateGradient
{
  typedef irtkEdgeLengthConstraint::GradientType Force;

  vtkPoints     *_Points;
  vtkPoints     *_InitialPoints;
  irtkEdgeTable *_EdgeTable;
  double         _RestLength;
  Force         *_Gradient;

  void operator ()(const blocked_range<int> &re) const
  {
    double     p0[3], p1[3], p2[3], e[3], w, d, d0 = _RestLength;
    const int *adjPts;
    int        numAdjPts;

    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId) {
      if (_InitialPoints) _InitialPoints->GetPoint(ptId, p0);
      _Points->GetPoint(ptId, p1);
      _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPts);
      for (int i = 0; i < numAdjPts; ++i) {
        if (_InitialPoints) {
          _InitialPoints->GetPoint(adjPts[i], p2);
          d0 = sqrt(vtkMath::Distance2BetweenPoints(p0, p2));
        }
        _Points->GetPoint(adjPts[i], p2);
        vtkMath::Subtract(p1, p2, e);
        d = vtkMath::Norm(e);
        w = 2.0 * (d - d0) / d;
        _Gradient[ptId] += w * Force(e[0], e[1], e[2]);
      }
    }
  }
};


} // namespace irtkEdgeLengthConstraintUtils

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkEdgeLengthConstraint::irtkEdgeLengthConstraint(const char *name, double weight)
:
  irtkPointSetConstraint(name, weight),
  _RestLength(-2.0)
{
  _ParameterPrefix.push_back("Stretching ");
  _ParameterPrefix.push_back("Distortion ");
  _ParameterPrefix.push_back("Edge stretching ");
  _ParameterPrefix.push_back("Surface stretching ");
  _ParameterPrefix.push_back("Surface distortion ");
  _ParameterPrefix.push_back("Mesh stretching ");
  _ParameterPrefix.push_back("Mesh distortion ");
  _ParameterPrefix.push_back("Surface mesh stretching ");
  _ParameterPrefix.push_back("Surface mesh distortion ");
  _ParameterPrefix.push_back("Metric distortion ");
}

// -----------------------------------------------------------------------------
irtkEdgeLengthConstraint::irtkEdgeLengthConstraint(const irtkEdgeLengthConstraint &other)
:
  irtkPointSetConstraint(other),
  _RestLength   (other._RestLength),
  _InitialPoints(other._InitialPoints),
  _EdgeTable    (other._EdgeTable)
{
}

// -----------------------------------------------------------------------------
irtkEdgeLengthConstraint &irtkEdgeLengthConstraint::operator =(const irtkEdgeLengthConstraint &other)
{
  irtkPointSetConstraint::operator =(other);
  _RestLength    = other._RestLength;
  _InitialPoints = other._InitialPoints;
  _EdgeTable     = other._EdgeTable;
  return *this;
}

// -----------------------------------------------------------------------------
irtkEdgeLengthConstraint::~irtkEdgeLengthConstraint()
{
}

// =============================================================================
// Configuration
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkEdgeLengthConstraint::Set(const char *param, const char *value)
{
  const string name = ParameterNameWithoutPrefix(param);
  if (name == "Rest length"         ||
      name == "Edge length"         ||
      name == "Rest edge length"    ||
      name == "Average edge length" ||
      name == "Average length") {
    string lvalue = value;
    transform(lvalue.begin(), lvalue.end(), lvalue.begin(), ::tolower);
    if (lvalue.size() > 12 && lvalue.compare(lvalue.size() - 12, 12, " edge length") == 0) {
      lvalue = lvalue.substr(0, lvalue.size() - 12);
    }
    if (lvalue == "average" || lvalue == "avg" || lvalue == "mean") {
      _RestLength = -1;
      return true;
    } else if (lvalue == "initial" || lvalue == "undeformed" || lvalue == "original") {
      _RestLength = -2;
      return true;
    } else {
      return FromString(value, _RestLength);
    }
  }
  return irtkPointSetConstraint::Set(param, value);
}

// -----------------------------------------------------------------------------
irtkParameterList irtkEdgeLengthConstraint::Parameter() const
{
  irtkParameterList params = irtkPointSetConstraint::Parameter();
  if (!_Name.empty()) {
    if (round(_RestLength) > -2) {
      Insert(params, _Name + " edge length", "Average edge length");
    } else if (_RestLength < 0) {
      Insert(params, _Name + " edge length", "Initial edge length");
    } else {
      Insert(params, _Name + " edge length", _RestLength);
    }
  }
  return params;
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void irtkEdgeLengthConstraint::Initialize()
{
  // Clear previous edge table
  _EdgeTable.Clear();
  _InitialPoints = NULL;

  // Initialize base class
  irtkPointSetConstraint::Initialize();
  if (_NumberOfPoints == 0) return;

  // Compute average edge length of undeformed mesh (after global transformation)
  if (_RestLength < .0) {
    _EdgeTable.Initialize(_PointSet->InputPointSet());
    if (_EdgeTable.NumberOfEdges() > 0) {
      if (round(_RestLength) > -2) {
        _RestLength = AverageEdgeLength(GetInitialPoints(), _EdgeTable);
      } else {
        _InitialPoints = GetInitialPoints();
      }
    } else {
      _RestLength = .0;
    }
    _EdgeTable.Clear();
  }

  // Initialize this class
  irtkEdgeLengthConstraint::Init();
}

// -----------------------------------------------------------------------------
void irtkEdgeLengthConstraint::Reinitialize()
{
  // Reinitialize base class
  irtkPointSetConstraint::Reinitialize();

  // Reinitialize this class
  irtkEdgeLengthConstraint::Init();
}

// -----------------------------------------------------------------------------
void irtkEdgeLengthConstraint::Init()
{
  // Update edge table
  _EdgeTable.Initialize(_PointSet->PointSet());

  // Update initial point locations
  if (_InitialPoints) _InitialPoints = GetInitialPoints();
}

// -----------------------------------------------------------------------------
double irtkEdgeLengthConstraint::Evaluate()
{
  if (_EdgeTable.NumberOfEdges() == 0) return .0;
  IRTK_START_TIMING();
  irtkEdgeLengthConstraintUtils::Evaluate eval;
  eval._Points        =  _PointSet->Points();
  eval._InitialPoints = _InitialPoints;
  eval._EdgeTable     = &_EdgeTable;
  eval._RestLength    =  _RestLength;
  parallel_reduce(blocked_range<int>(0, _EdgeTable.NumberOfEdges()), eval);
  IRTK_DEBUG_TIMING(3, "evaluation of distortion");
  return eval._Sum / _EdgeTable.NumberOfEdges();
}

// -----------------------------------------------------------------------------
void irtkEdgeLengthConstraint::EvaluateGradient(double *gradient, double step, double weight)
{
  if (_EdgeTable.NumberOfEdges() == 0) return;

  IRTK_START_TIMING();
  memset(_Gradient, 0, _NumberOfPoints * sizeof(GradientType));

  irtkEdgeLengthConstraintUtils::EvaluateGradient eval;
  eval._Points        =  _PointSet->Points();
  eval._InitialPoints = _InitialPoints;
  eval._EdgeTable     = &_EdgeTable;
  eval._RestLength    =  _RestLength;
  eval._Gradient      =  _Gradient;
  parallel_for(blocked_range<int>(0, _NumberOfPoints), eval);

  irtkPointSetConstraint::EvaluateGradient(gradient, step, weight / _EdgeTable.NumberOfEdges());
  IRTK_DEBUG_TIMING(3, "evaluation of distortion force");
}
