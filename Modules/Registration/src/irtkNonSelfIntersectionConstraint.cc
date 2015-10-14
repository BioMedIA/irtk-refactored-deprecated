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

#include <irtkNonSelfIntersectionConstraint.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkGenericCell.h>
#include <vtkIdList.h>
#include <vtkMath.h>
#include <vtkLine.h>
#include <vtkAbstractPointLocator.h>
#include <vtkOctreePointLocator.h>
#include <vtkPolyDataNormals.h>

#include <irtkEdgeTable.h>

using namespace irtk::polydata;


#define APPROXIMATE_NSI 1

// =============================================================================
// Auxiliaries
// =============================================================================

namespace irtkNonSelfIntersectionConstraintUtils {

typedef irtkNonSelfIntersectionConstraint::GradientType  Force;


// -----------------------------------------------------------------------------
/// Calculate sum of edge lengths
struct SumEdgeLengths
{
  vtkPoints     *_Points;
  irtkEdgeTable *_EdgeTable;
  double         _Sum;

  SumEdgeLengths() : _Sum(.0) {}

  SumEdgeLengths(const SumEdgeLengths &other, split)
  :
    _Points(other._Points),
    _EdgeTable(other._EdgeTable),
    _Sum(.0)
  {}

  void join(const SumEdgeLengths &other)
  {
    _Sum += other._Sum;
  }

  void operator ()(const blocked_range<int> &re)
  {
    vtkIdType ptId1, ptId2;
    double    p1[3], p2[3];

    irtkEdgeIterator it(*_EdgeTable);
    for (it.InitTraversal(re); it.GetNextEdge(ptId1, ptId2) != -1;) {
      _Points->GetPoint(ptId1, p1);
      _Points->GetPoint(ptId2, p2);
      _Sum += sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
    }
  }
};

// -----------------------------------------------------------------------------
/// Get center point and radius of bounding sphere
inline double GetBoundingSphereRadius(vtkGenericCell *cell, double c[3])
{
  double p[3];
  // Get center of bounding sphere
  c[0] = .0, c[1] = .0, c[2] = .0;
  double *weights = new double[cell->GetNumberOfPoints()];
  int subId = cell->GetParametricCenter(p);
  cell->EvaluateLocation(subId, p, c, weights);
  delete[] weights;
  // Get radius of bounding sphere
  double r = .0;
  vtkPoints *points = cell->GetPoints();
  for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i) {
    points->GetPoint(i, p);
    r = max(r, vtkMath::Distance2BetweenPoints(c, p));
  }
  return sqrt(r);
}

#ifdef APPROXIMATE_NSI

// -----------------------------------------------------------------------------
/// Evaluate non-self-intersection penalty
struct Evaluate
{
  vtkPolyData             *_DataSet;
  vtkAbstractPointLocator *_Locator;
  double                   _Distance;
  double                   _Sum;
  int                      _Num;

  Evaluate() : _Sum(.0), _Num(0) {}

  Evaluate(const Evaluate &other, split)
  :
    _DataSet (other._DataSet),
    _Locator (other._Locator),
    _Distance(other._Distance),
    _Sum(.0), _Num(0)
  {}

  void join(const Evaluate &other)
  {
    _Sum += other._Sum;
    _Num += other._Num;
  }

  void operator ()(const blocked_range<vtkIdType> &re)
  {
    int       subId;
    vtkIdType ptId;
    double    p1[3], p2[3], c[3], r, pcoords[3], dist2;
    double   *weights = new double[_DataSet->GetMaxCellSize()];

    vtkSmartPointer<vtkGenericCell> cell      = vtkSmartPointer<vtkGenericCell>::New();
    vtkSmartPointer<vtkIdList>      cellPtIds = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList>      ptIds     = vtkSmartPointer<vtkIdList>::New();

    const double maxdist2 = _Distance * _Distance;

    for (vtkIdType cellId = re.begin(); cellId != re.end(); ++cellId) {
      _DataSet->GetCell(cellId, cell);
      _DataSet->GetCellPoints(cellId, cellPtIds);
      r = GetBoundingSphereRadius(cell, c);
      _Locator->FindPointsWithinRadius(3 * r + _Distance, c, ptIds);
      for (vtkIdType i = 0; i < ptIds->GetNumberOfIds(); ++i) {
        ptId = ptIds->GetId(i);
        if (cellPtIds->IsId(ptId) == -1) {
          _DataSet->GetPoint(ptId, p1);
          if (cell->EvaluatePosition(p1, p2, subId, pcoords, dist2, weights) == 1 && dist2 < maxdist2) {
            _Sum += pow(sqrt(dist2) - _Distance, 2);
            ++_Num;
          }
        }
      }
    }

    delete[] weights;
  }
};

// -----------------------------------------------------------------------------
struct EvaluateGradient
{
  vtkPolyData             *_DataSet;
  vtkAbstractPointLocator *_Locator;
  double                   _Distance;
  Force                   *_Gradient;
  int                     *_Count;
  int                      _Num;

  EvaluateGradient() : _Gradient(NULL), _Count(NULL), _Num(0) {}

  EvaluateGradient(const EvaluateGradient &other, split)
  :
    _DataSet (other._DataSet),
    _Locator (other._Locator),
    _Distance(other._Distance)
  {
    _Gradient = CAllocate<Force>(_DataSet->GetNumberOfPoints());
    _Count    = CAllocate<int>  (_DataSet->GetNumberOfPoints());
    _Num      = 0;
  }

  void join(EvaluateGradient &other)
  {
    for (vtkIdType i = 0; i < _DataSet->GetNumberOfPoints(); ++i) {
      _Gradient[i] += other._Gradient[i];
      _Count   [i] += other._Count   [i];
    }
    Deallocate(other._Gradient);
    Deallocate(other._Count);
    _Num += other._Num;
  }

  void operator ()(const blocked_range<vtkIdType> &re)
  {
    int       subId;
    vtkIdType ptId;
    double    p1[3], p2[3], c[3], r, pcoords[3], dist2, d, w;
    double   *weights = new double[_DataSet->GetMaxCellSize()];

    vtkSmartPointer<vtkGenericCell> cell      = vtkSmartPointer<vtkGenericCell>::New();
    vtkSmartPointer<vtkIdList>      cellPtIds = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList>      ptIds     = vtkSmartPointer<vtkIdList>::New();

    const double maxdist2 = _Distance * _Distance;

    for (vtkIdType cellId = re.begin(); cellId != re.end(); ++cellId) {
      _DataSet->GetCell(cellId, cell);
      _DataSet->GetCellPoints(cellId, cellPtIds);
      r = GetBoundingSphereRadius(cell, c);
      _Locator->FindPointsWithinRadius(3 * r + _Distance, c, ptIds);
      for (vtkIdType i = 0; i < ptIds->GetNumberOfIds(); ++i) {
        ptId = ptIds->GetId(i);
        if (cellPtIds->IsId(ptId) == -1) {
          _DataSet->GetPoint(ptId, p1);
          if (cell->EvaluatePosition(p1, p2, subId, pcoords, dist2, weights) == 1 && dist2 < maxdist2) {
            d = sqrt(dist2);
            w = 2.0 * (d - _Distance) / (d + 1e-6);
            _Gradient[ptId] += Force(w * (p1[0] - p2[0]),
                                     w * (p1[1] - p2[1]),
                                     w * (p1[2] - p2[2]));
            ++_Count[ptId];
            ++_Num;
          }
        }
      }
    }

    delete[] weights;
  }
};

#else

/// Structure storing pre-computed information of candidate pairs of cells
typedef irtkNonSelfIntersectionConstraint::CellPair CellPair;

// -----------------------------------------------------------------------------
/// Determine whether two triangles are at least min distance apart
inline bool NotCloserThan(double dist, vtkGenericCell *cell1, vtkGenericCell *cell2)
{
  double c1[3], c2[3], r1, r2;
  r1 = GetBoundingSphereRadius(cell1, c1);
  r2 = GetBoundingSphereRadius(cell2, c2);
  return sqrt(vtkMath::Distance2BetweenPoints(c1, c2)) >= r1 + r2 + dist;
}

// -----------------------------------------------------------------------------
inline void CopyPoint(double *dst, double *src)
{
  memcpy(dst, src, 3 * sizeof(double));
}

// -----------------------------------------------------------------------------
/// Find closest points of two triangles and their distance
///
/// TODO: Use code available freely of the following paper:
/// "A fast triangle to triangle intersection test for collision detection"
/// Oren Tropp, Ayellet Tal, Ilan Shimshoni
/// Computer Animation and Virtual Worlds 17(5) 2006, pp 527-535
inline double DistanceBetweenTriangles(vtkGenericCell *cell1, vtkGenericCell *cell2, double closestPoint1[3], double closestPoint2[3])
{
  irtkAssert(cell1->GetNumberOfPoints() == 3, "cell1 is a triangles");
  irtkAssert(cell2->GetNumberOfPoints() == 3, "cell2 is a triangles");

  double dist2, min_dist2 = numeric_limits<double>::infinity();
  double p1[3], p2[3], w[3], pcoords[3], t1, t2;
  int    subId;

  // Get corners of each triangle
  vtkPoints *points1 = cell1->GetPoints();
  vtkPoints *points2 = cell2->GetPoints();

  // Point in cell2 closest to any of the corners of cell1
  if (cell2->EvaluatePosition(points1->GetPoint(0), p2, subId, pcoords, dist2, w) != -1 && dist2 < min_dist2) {
    CopyPoint(closestPoint1, points1->GetPoint(0));
    CopyPoint(closestPoint2, p2);
    min_dist2 = dist2;
  }
  if (cell2->EvaluatePosition(points1->GetPoint(1), p2, subId, pcoords, dist2, w) != -1 && dist2 < min_dist2) {
    CopyPoint(closestPoint1, points1->GetPoint(1));
    CopyPoint(closestPoint2, p2);
    min_dist2 = dist2;
  }
  if (cell2->EvaluatePosition(points1->GetPoint(2), p2, subId, pcoords, dist2, w) != -1 && dist2 < min_dist2) {
    CopyPoint(closestPoint1, points1->GetPoint(2));
    CopyPoint(closestPoint2, p2);
    min_dist2 = dist2;
  }

  // Point in cell1 closest to any of the corners of cell2
  if (cell1->EvaluatePosition(points2->GetPoint(0), p1, subId, pcoords, dist2, w) != -1 && dist2 < min_dist2) {
    CopyPoint(closestPoint1, p1);
    CopyPoint(closestPoint2, points2->GetPoint(0));
    min_dist2 = dist2;
  }
  if (cell1->EvaluatePosition(points2->GetPoint(1), p1, subId, pcoords, dist2, w) != -1 && dist2 < min_dist2) {
    CopyPoint(closestPoint1, p1);
    CopyPoint(closestPoint2, points2->GetPoint(1));
    min_dist2 = dist2;
  }
  if (cell1->EvaluatePosition(points2->GetPoint(2), p1, subId, pcoords, dist2, w) != -1 && dist2 < min_dist2) {
    CopyPoint(closestPoint1, p1);
    CopyPoint(closestPoint2, points2->GetPoint(2));
    min_dist2 = dist2;
  }

  // Distance between each pair of edges
  dist2 = vtkLine::DistanceBetweenLineSegments(points1->GetPoint(0), points1->GetPoint(1),
                                               points2->GetPoint(0), points2->GetPoint(1),
                                               p1, p2, t1, t2);
  if (dist2 < min_dist2) {
    CopyPoint(closestPoint1, p1);
    CopyPoint(closestPoint2, p2);
    min_dist2 = dist2;
  }
  dist2 = vtkLine::DistanceBetweenLineSegments(points1->GetPoint(0), points1->GetPoint(1),
                                               points2->GetPoint(1), points2->GetPoint(2),
                                               p1, p2, t1, t2);
  if (dist2 < min_dist2) {
    CopyPoint(closestPoint1, p1);
    CopyPoint(closestPoint2, p2);
    min_dist2 = dist2;
  }
  dist2 = vtkLine::DistanceBetweenLineSegments(points1->GetPoint(0), points1->GetPoint(1),
                                               points2->GetPoint(2), points2->GetPoint(0),
                                               p1, p2, t1, t2);
  if (dist2 < min_dist2) {
    CopyPoint(closestPoint1, p1);
    CopyPoint(closestPoint2, p2);
    min_dist2 = dist2;
  }
  dist2 = vtkLine::DistanceBetweenLineSegments(points1->GetPoint(1), points1->GetPoint(2),
                                               points2->GetPoint(0), points2->GetPoint(1),
                                               p1, p2, t1, t2);
  if (dist2 < min_dist2) {
    CopyPoint(closestPoint1, p1);
    CopyPoint(closestPoint2, p2);
    min_dist2 = dist2;
  }
  dist2 = vtkLine::DistanceBetweenLineSegments(points1->GetPoint(1), points1->GetPoint(2),
                                               points2->GetPoint(1), points2->GetPoint(2),
                                               p1, p2, t1, t2);
  if (dist2 < min_dist2) {
    CopyPoint(closestPoint1, p1);
    CopyPoint(closestPoint2, p2);
    min_dist2 = dist2;
  }
  dist2 = vtkLine::DistanceBetweenLineSegments(points1->GetPoint(1), points1->GetPoint(2),
                                               points2->GetPoint(2), points2->GetPoint(0),
                                               p1, p2, t1, t2);
  if (dist2 < min_dist2) {
    CopyPoint(closestPoint1, p1);
    CopyPoint(closestPoint2, p2);
    min_dist2 = dist2;
  }
  dist2 = vtkLine::DistanceBetweenLineSegments(points1->GetPoint(2), points1->GetPoint(0),
                                               points2->GetPoint(0), points2->GetPoint(1),
                                               p1, p2, t1, t2);
  if (dist2 < min_dist2) {
    CopyPoint(closestPoint1, p1);
    CopyPoint(closestPoint2, p2);
    min_dist2 = dist2;
  }
  dist2 = vtkLine::DistanceBetweenLineSegments(points1->GetPoint(2), points1->GetPoint(0),
                                               points2->GetPoint(1), points2->GetPoint(2),
                                               p1, p2, t1, t2);
  if (dist2 < min_dist2) {
    CopyPoint(closestPoint1, p1);
    CopyPoint(closestPoint2, p2);
    min_dist2 = dist2;
  }
  dist2 = vtkLine::DistanceBetweenLineSegments(points1->GetPoint(2), points1->GetPoint(0),
                                               points2->GetPoint(2), points2->GetPoint(0),
                                               p1, p2, t1, t2);
  if (dist2 < min_dist2) {
    CopyPoint(closestPoint1, p1);
    CopyPoint(closestPoint2, p2);
    min_dist2 = dist2;
  }

  return sqrt(min_dist2);
}

// -----------------------------------------------------------------------------
/// Find pairs of cells (triangles) which are too close to each other
struct FindCandidates
{
  vtkOctreePointLocator *_Locator;
  vector<CellPair>       _Candidates;
  double                 _MinDistance;

  FindCandidates() {}

  FindCandidates(const FindCandidates &other, split)
  :
    _Locator    (other._Locator),
    _MinDistance(other._MinDistance)
  {}

  void join(const FindCandidates &other)
  {
//    _Candidates.reserve(_Candidates.size() + other._Candidates.size());
//    for (vector<CellPair>::const_iterator i = other._Candidates.begin(); i != other._Candidates.end(); ++i) {
//      bool exists = false;
//      for (vector<CellPair>::const_iterator j = _Candidates.begin(); j != _Candidates.end(); ++j) {
//        if (i->_CellId1 == j->_CellId1 && i->_CellId2 == j->_CellId2) {
//          exists = true;
//          break;
//        }
//      }
//      if (!exists) _Candidates.push_back(*i);
//    }
    _Candidates.insert(_Candidates.end(), other._Candidates.begin(), other._Candidates.end());
  }

  void operator ()(const blocked_range<vtkIdType> &re)
  {
    vtkSmartPointer<vtkGenericCell> cell1       = vtkSmartPointer<vtkGenericCell>::New();
    vtkSmartPointer<vtkGenericCell> cell2       = vtkSmartPointer<vtkGenericCell>::New();
    vtkSmartPointer<vtkIdList>      nearCellIds = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList>      nearPtIds   = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList>      cellIds     = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList>      ptIds       = vtkSmartPointer<vtkIdList>::New();

    vtkDataSet *dataset = _Locator->GetDataSet();

    CellPair pair;
    double c1[3], r1;
    double c2[3], r2;
    double dist2;

    for (pair._CellId1 = re.begin(); pair._CellId1 != re.end(); ++pair._CellId1) {
      // Get first cell
      dataset->GetCell(pair._CellId1, cell1);
      // Get bounding sphere
      r1 = GetBoundingSphereRadius(cell1, c1);
      // Consider only cells with points within a certain radius
      // Note: This assumes that all cells have roughly the same bounding sphere
      //       radius. To account for difference, we use a factor > 2.
      _Locator->FindPointsWithinRadius(3.0 * r1 + _MinDistance, c1, nearPtIds);
      cellIds->Reset();
      for (vtkIdType i = 0; i < nearPtIds->GetNumberOfIds(); ++i) {
        dataset->GetPointCells(nearPtIds->GetId(i), nearCellIds);
        for (vtkIdType j = 0; j < nearCellIds->GetNumberOfIds(); ++j) {
          pair._CellId2 = nearCellIds->GetId(j);
          if (pair._CellId1 < pair._CellId2) { // unique pairs
            dataset->GetCellPoints(pair._CellId2, ptIds);
            ptIds->IntersectWith(cell1->GetPointIds());
            if (ptIds->GetNumberOfIds() == 0) { // non-neighboring pairs
              cellIds->InsertUniqueId(pair._CellId2);
            }
          }
        }
      }
      // Determine distance of cell pairs
      for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {
        pair._CellId2 = cellIds->GetId(i);
        dataset->GetCell(pair._CellId2, cell2);
        r2 = GetBoundingSphereRadius(cell2, c2);
        dist2 = r1 + r2 + _MinDistance;
        dist2 *= dist2;
        if (vtkMath::Distance2BetweenPoints(c1, c2) >= dist2) continue;
        pair._Distance = DistanceBetweenTriangles(cell1, cell2, pair._Point1, pair._Point2);
        if (pair._Distance < _MinDistance) {
          _Candidates.push_back(pair);
        }
      }
    }
  }

  static vector<CellPair> Run(vtkPolyData *polydata, double dist)
  {
    if (polydata->GetNumberOfCells() == 0) return vector<CellPair>();
    // Build point locator
    vtkSmartPointer<vtkOctreePointLocator> octree = vtkSmartPointer<vtkOctreePointLocator>::New();
    octree->SetDataSet(polydata);
    octree->BuildLocator();
    // Find candidate pairs
    FindCandidates body;
    body._Locator     = octree;
    body._MinDistance = dist;
    parallel_reduce(blocked_range<vtkIdType>(0, polydata->GetNumberOfCells()), body);
    return body._Candidates;
  }
};

// -----------------------------------------------------------------------------
/// Update distances of candidate pairs
struct UpdateDistances
{
  vtkPolyData      *_DataSet;
  vector<CellPair> *_Candidates;

  void operator ()(const blocked_range<size_t> &re) const
  {
    vtkSmartPointer<vtkGenericCell> cell1 = vtkSmartPointer<vtkGenericCell>::New();
    vtkSmartPointer<vtkGenericCell> cell2 = vtkSmartPointer<vtkGenericCell>::New();
    for (size_t i = re.begin(); i != re.end(); ++i) {
      CellPair &pair = (*_Candidates)[i];
      _DataSet->GetCell(pair._CellId1, cell1);
      _DataSet->GetCell(pair._CellId2, cell2);
      pair._Distance = DistanceBetweenTriangles(cell1, cell2, pair._Point1, pair._Point2);
    }
  }
};

// -----------------------------------------------------------------------------
/// Evaluate non-self-intersection penalty
struct Evaluate
{
  const vector<CellPair> *_Candidates;
  double                  _MinDistance;
  double                  _Sum;

  Evaluate() : _Sum(.0) {}

  Evaluate(const Evaluate &other, split)
  :
    _Candidates (other._Candidates),
    _MinDistance(other._MinDistance),
    _Sum(.0)
  {}

  void join(const Evaluate &other)
  {
    _Sum += other._Sum;
  }

  void operator ()(const blocked_range<size_t> &re)
  {
    for (size_t i = re.begin(); i != re.end(); ++i) {
      const CellPair &pair = (*_Candidates)[i];
      _Sum += _MinDistance - pair._Distance;
    }
  }
};

// -----------------------------------------------------------------------------
// Evaluate non-self-intersection force
struct EvaluateGradient
{
  typedef irtkNonSelfIntersectionConstraint::Force Force;

  vtkSmartPointer<vtkPolyData> _DataSet;
  const vector<CellPair>      *_Candidates;
  double                       _Distance;
  Force                       *_Gradient;
  int                         *_Count;

  EvaluateGradient() : _Gradient(NULL), _Count(NULL) {}

  EvaluateGradient(const EvaluateGradient &other, split)
  :
    _DataSet   (other._DataSet),
    _Candidates(other._Candidates),
    _Distance  (other._Distance)
  {
    _Gradient = CAllocate<Force>(_DataSet->GetNumberOfPoints());
    _Count    = CAllocate<int>  (_DataSet->GetNumberOfPoints());
  }

  void join(EvaluateGradient &other)
  {
    for (vtkIdType i = 0; i < _DataSet->GetNumberOfPoints(); ++i) {
      _Gradient[i] += other._Gradient[i];
      _Count   [i] += other._Count   [i];
    }
    Deallocate(other._Gradient);
    Deallocate(other._Count);
  }

  void operator ()(const blocked_range<size_t> &re)
  {
    double p[3], d, dist_weight, w;
    vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
    vtkIdType ptId;
    for (size_t i = re.begin(); i != re.end(); ++i) {
      const CellPair &pair = (*_Candidates)[i];
      dist_weight = abs(_Distance - pair._Distance) / pair._Distance;
      _DataSet->GetCellPoints(pair._CellId1, ptIds);
      for (vtkIdType i = 0; i < ptIds->GetNumberOfIds(); ++i) {
        ptId = ptIds->GetId(i);
        _DataSet->GetPoint(ptId, p);
        d = vtkMath::Distance2BetweenPoints(p, pair._Point1);
        w = (d > .0 ? dist_weight / d : dist_weight);
        _Gradient[ptId] = Force(w * (pair._Point1[0] - pair._Point2[0]),
                                w * (pair._Point1[1] - pair._Point2[1]),
                                w * (pair._Point1[2] - pair._Point2[2]));
        ++_Count[ptId];
      }
      _DataSet->GetCellPoints(pair._CellId2, ptIds);
      for (vtkIdType i = 0; i < ptIds->GetNumberOfIds(); ++i) {
        ptId = ptIds->GetId(i);
        _DataSet->GetPoint(ptId, p);
        d = vtkMath::Distance2BetweenPoints(p, pair._Point2);
        w = (d > .0 ? dist_weight / d : dist_weight);
        _Gradient[ptId] = Force(w * (pair._Point2[0] - pair._Point1[0]),
                                w * (pair._Point2[1] - pair._Point1[1]),
                                w * (pair._Point2[2] - pair._Point1[2]));
        ++_Count[ptId];
      }
    }
  }
};

#endif // APPROXIMATE_NSI


} // namespace irtkNonSelfIntersectionConstraintUtils
using namespace irtkNonSelfIntersectionConstraintUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkNonSelfIntersectionConstraint
::irtkNonSelfIntersectionConstraint(const char *name, double weight)
:
  irtkSurfaceConstraint(name, weight),
  _MinDistance(-1.0)
{
  _ParameterPrefix.push_back("Non-self-intersection ");
  _ParameterPrefix.push_back("Non-self-intersection force ");
}

// -----------------------------------------------------------------------------
irtkNonSelfIntersectionConstraint
::irtkNonSelfIntersectionConstraint(const irtkNonSelfIntersectionConstraint &other)
:
  irtkSurfaceConstraint(other),
  _MinDistance(other._MinDistance),
  _Candidates (other._Candidates)
{
}

// -----------------------------------------------------------------------------
irtkNonSelfIntersectionConstraint &irtkNonSelfIntersectionConstraint
::operator =(const irtkNonSelfIntersectionConstraint &other)
{
  irtkSurfaceConstraint::operator =(other);
  _MinDistance = other._MinDistance;
  _Candidates  = other._Candidates;
  return *this;
}

// -----------------------------------------------------------------------------
irtkNonSelfIntersectionConstraint::~irtkNonSelfIntersectionConstraint()
{
}

// =============================================================================
// Configuration
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkNonSelfIntersectionConstraint::Set(const char *param, const char *value)
{
  const string name = ParameterNameWithoutPrefix(param);
  if (name == "Distance"           ||
      name == "Distance threshold" ||
      name == "Minimum distance"   ||
      name == "Minimum surface distance") {
    return FromString(value, _MinDistance);
  }
  return irtkSurfaceConstraint::Set(param, value);
}

// -----------------------------------------------------------------------------
irtkParameterList irtkNonSelfIntersectionConstraint::Parameter() const
{
  irtkParameterList params = irtkSurfaceConstraint::Parameter();
  if (!_Name.empty()) {
    Insert(params, _Name + " distance threshold", ToString(_MinDistance));
  }
  return params;
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void irtkNonSelfIntersectionConstraint::Init()
{
  AllocateCount(_NumberOfPoints);
}

// -----------------------------------------------------------------------------
void irtkNonSelfIntersectionConstraint::Initialize()
{
  irtkSurfaceConstraint::Initialize();
  irtkNonSelfIntersectionConstraint::Init();

  // By default, set minimum distance relative to average edge length of undeformed mesh
  if (_MinDistance < .0) {
    irtkEdgeTable edgeTable(_PointSet->InputSurface());
    if (edgeTable.NumberOfEdges() > 0) {
      irtkNonSelfIntersectionConstraintUtils::SumEdgeLengths eval;
      vtkSmartPointer<vtkPoints> points = GetInitialPoints();
      eval._Points    = points;
      eval._EdgeTable = &edgeTable;
      parallel_reduce(blocked_range<int>(0, edgeTable.NumberOfEdges()), eval);
      _MinDistance = 10.0 * (eval._Sum / edgeTable.NumberOfEdges());
    }
  }
}

// -----------------------------------------------------------------------------
void irtkNonSelfIntersectionConstraint::Reinitialize()
{
  irtkSurfaceConstraint::Reinitialize();
  irtkNonSelfIntersectionConstraint::Init();
}

// -----------------------------------------------------------------------------
void irtkNonSelfIntersectionConstraint::Update(bool gradient)
{
#if APPROXIMATE_NSI
#else
//  if (gradient) {
    _Candidates = FindCandidates::Run(_PointSet->Surface(), _MinDistance);
//  } else {
//    UpdateDistances update;
//    update._DataSet    = _PointSet->Surface();
//    update._Candidates = &_Candidates;
//    parallel_for(blocked_range<size_t>(0, _Candidates.size()), update);
//  }
#endif
}

// -----------------------------------------------------------------------------
double irtkNonSelfIntersectionConstraint::Evaluate()
{
#if APPROXIMATE_NSI
  if (_NumberOfPoints == 0) return .0;

  // TODO: Build Octree in Update function and use for both Evaluate and EvaluateGradient
  vtkSmartPointer<vtkOctreePointLocator> octree = vtkSmartPointer<vtkOctreePointLocator>::New();
  octree->SetDataSet(_PointSet->Surface());
  octree->BuildLocator();

  irtkNonSelfIntersectionConstraintUtils::Evaluate eval;
  eval._DataSet  = _PointSet->Surface();
  eval._Locator  = octree;
  eval._Distance = _MinDistance;
  parallel_reduce(blocked_range<vtkIdType>(0, _PointSet->NumberOfSurfaceCells()), eval);

  return (eval._Num > 0 ? eval._Sum / eval._Num : .0);
#else
  if (_Candidates.empty()) return .0;

  irtkNonSelfIntersectionConstraintUtils::Evaluate eval;
  eval._Candidates  = &_Candidates;
  eval._MinDistance = _MinDistance;
  parallel_reduce(blocked_range<size_t>(0, _Candidates.size()), eval);

  return eval._Sum / _Candidates.size();
#endif
}

// -----------------------------------------------------------------------------
void irtkNonSelfIntersectionConstraint::EvaluateGradient(double *gradient, double step, double weight)
{
#if APPROXIMATE_NSI
  if (_NumberOfPoints == 0) return;

  vtkSmartPointer<vtkOctreePointLocator> octree = vtkSmartPointer<vtkOctreePointLocator>::New();
  octree->SetDataSet(_PointSet->Surface());
  octree->BuildLocator();

  irtkNonSelfIntersectionConstraintUtils::EvaluateGradient eval;
  eval._DataSet  = _PointSet->Surface();
  eval._Distance = _MinDistance;
  eval._Locator  = octree;
  eval._Gradient = _Gradient;
  eval._Count    = _Count;
  parallel_reduce(blocked_range<vtkIdType>(0, _PointSet->NumberOfSurfaceCells()), eval);

  for (int i = 0; i < _NumberOfPoints; ++i) {
    if (_Count[i] > 0) _Gradient[i] /= _Count[i];
  }

  irtkSurfaceConstraint::EvaluateGradient(gradient, step, weight / eval._Num);
#else
  if (_Candidates.empty()) return;

  memset(_Gradient, 0, _NumberOfPoints * sizeof(GradientType));
  memset(_Count,    0, _NumberOfPoints * sizeof(int));

  irtkNonSelfIntersectionConstraintUtils::EvaluateGradient eval;
  eval._DataSet    = _PointSet->Surface();
  eval._Candidates = &_Candidates;
  eval._Distance   = _MinDistance;
  eval._Gradient   = _Gradient;
  eval._Count      = _Count;
  parallel_reduce(blocked_range<size_t>(0, _Candidates.size()), eval);

  for (int i = 0; i < _NumberOfPoints; ++i) {
    if (_Count[i] > 0) _Gradient[i] /= _Count[i];
  }

  irtkSurfaceConstraint::EvaluateGradient(gradient, step, weight / _Candidates.size());
#endif
}
