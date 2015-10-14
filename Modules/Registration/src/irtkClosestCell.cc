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

#include <irtkClosestCell.h>

#include <vtkSmartPointer.h>
#include <vtkAbstractCellLocator.h>
#include <vtkCellLocator.h>
#include <vtkCellTreeLocator.h>
#include <vtkModifiedBSPTree.h>
#include <vtkOBBTree.h>

#include <irtkParallel.h>
#include <irtkPoint.h>
#include <irtkTransformation.h>


// =============================================================================
// Auxiliary functors
// =============================================================================

namespace irtkClosestCellUtils {


// Maximum point coordinate difference to consider two points equal
// TODO: This should be relative to the "resolution" of the input data sets.
const double TOL = 1e-3;

// -----------------------------------------------------------------------------
typedef enum irtkClosestCell::LocatorType LocatorType;

// -----------------------------------------------------------------------------
class UpdateCorrespondences
{
private:

  vtkPointSet              *_Target;
  const vector<int>        *_Sample;
  vtkAbstractCellLocator   *_Source;
  irtkPointSet             *_Points;
  vector<double>           *_Distance;
  const irtkTransformation *_Transformation;
  bool                      _Changed;

public:

  UpdateCorrespondences() : _Changed(false) {}

  UpdateCorrespondences(const UpdateCorrespondences &lhs, split)
  :
    _Target(lhs._Target),
    _Sample(lhs._Sample),
    _Source(lhs._Source),
    _Points(lhs._Points),
    _Distance(lhs._Distance),
    _Transformation(lhs._Transformation),
    _Changed(false)
  {}

  void join(const UpdateCorrespondences &rhs)
  {
    if (rhs._Changed) _Changed = true;
  }

  void operator()(const blocked_range<int> &re)
  {
    vtkIdType cellId;
    int       subId;
    double    p[3];

    vector<double> &distance = *_Distance;

    for (int k = re.begin(); k != re.end(); ++k) {
      _Target->GetPoint(irtkPointCorrespondence::GetPointIndex(_Target, _Sample, k), p);
      _Source->FindClosestPoint(p, p, cellId, subId, distance[k]);
      distance[k] = sqrt(distance[k]);
      if (_Transformation) {
        if (!_Transformation->Inverse(p[0], p[1], p[2])) {
          // TODO: Use cell index and parametric coordinates of closest point
          //       to find (approximate) corresponding point in input data set.
        }
      }
      const irtkPoint &old = _Points->GetPoint(k);
      if (!fequal(old._x, p[0], TOL) ||
          !fequal(old._y, p[1], TOL) ||
          !fequal(old._z, p[2], TOL)) {
        _Points->SetPoint(k, p);
        _Changed = true;
      }
    }
  }

  static bool Run(const irtkRegisteredPointSet *target,
                  const vector<int>            *sample,
                  const irtkRegisteredPointSet *source,
                  LocatorType                   type,
                  int                           cells_per_node,
                  irtkPointSet                 &points,
                  vector<double>               &distance,
                  const irtkTransformation     *transformation = NULL)
  {
    const int n = irtkPointCorrespondence::GetNumberOfPoints(target, sample);
    points  .Reserve(n);
    points  .Resize (n);
    distance.resize (n);
    if (n == 0) return false;
    // Initialize cell locator
    vtkSmartPointer<vtkAbstractCellLocator> locator;
    switch (type) {
      case LocatorType::Default:  locator = vtkSmartPointer<vtkCellLocator>    ::New(); break;
      case LocatorType::CellTree: locator = vtkSmartPointer<vtkCellTreeLocator>::New(); break;
      case LocatorType::BSPTree:  locator = vtkSmartPointer<vtkModifiedBSPTree>::New(); break;
      case LocatorType::OBBTree:  locator = vtkSmartPointer<vtkOBBTree>        ::New(); break;
      default:
        cerr << "irtkClosestCell::Initialize: Invalid locator type = " << type << endl;
        exit(1);
    }
    locator->SetNumberOfCellsPerNode(cells_per_node);
    locator->SetDataSet(source->PointSet());
    locator->BuildLocator();
    // Find closest cell points
    UpdateCorrespondences body;
    body._Target         = target->PointSet();
    body._Sample         = sample;
    body._Source         = locator;
    body._Points         = &points;
    body._Distance       = &distance;
    body._Transformation = transformation;
    blocked_range<int> range(0, n);
    // TODO: vtkAbstractCellLocator's are not thread-safe!
    //       Why not use separate cell locator instances for each thread?
    //       Use concurrent_bounded_queue to store pointers to these.
    //parallel_reduction(range, body);
    body(range);
    return body._Changed;
  }
};


} // namespace irtkClosestCellUtils
using namespace irtkClosestCellUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkClosestCell::irtkClosestCell()
:
  _LocatorType(Default),
  _NumberOfCellsPerNode(10),
  _Sigma(numeric_limits<double>::quiet_NaN()),
  _MaxDistance(numeric_limits<double>::infinity())
{
}

// -----------------------------------------------------------------------------
irtkClosestCell
::irtkClosestCell(const irtkClosestCell &other)
:
  irtkPointCorrespondence(other),
  _LocatorType         (other._LocatorType),
  _NumberOfCellsPerNode(other._NumberOfCellsPerNode),
  _Sigma               (other._Sigma),
  _MaxDistance         (other._MaxDistance)
{
}

// -----------------------------------------------------------------------------
irtkPointCorrespondence *irtkClosestCell::NewInstance() const
{
  return new irtkClosestCell(*this);
}

// -----------------------------------------------------------------------------
irtkClosestCell::~irtkClosestCell()
{
}

// -----------------------------------------------------------------------------
irtkClosestCell::TypeId irtkClosestCell::Type() const
{
  return ClosestCell;
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkClosestCell::Set(const char *name, const char *value)
{
  if (strcmp(name, "Locator") == 0 || strcmp(name, "Locator type") == 0) {
    if (strcmp(value, "Default") == 0) {
      _LocatorType = Default;
      return true;
    }
    if (strcmp(value, "CellTree") == 0) {
      _LocatorType = CellTree;
      return true;
    }
    if (strcmp(value, "BSPTree")  == 0 ||
        strcmp(value, "BSP Tree") == 0 ||
        strcmp(value, "BSP")      == 0) {
      _LocatorType = BSPTree;
      return true;
    }
    if (strcmp(value, "OBBTree")  == 0 ||
        strcmp(value, "OBB Tree") == 0 ||
        strcmp(value, "OBB")      == 0) {
      _LocatorType = OBBTree;
      return true;
    }
    return false;
  }
  if (strcmp(name, "No. of cells per node") == 0 ||
      strcmp(name, "No. of cells")          == 0) {
    return FromString(value, _NumberOfCellsPerNode) && _NumberOfCellsPerNode > 0;
  }
  if (strcmp(name, "Sigma") == 0) {
    return FromString(value, _Sigma);
  }
  if (strcmp(name, "Maximum distance") == 0) {
    return FromString(value, _MaxDistance);
  }
  return false;
}

// -----------------------------------------------------------------------------
irtkParameterList irtkClosestCell::Parameter() const
{
  irtkParameterList params;
  switch (_LocatorType) {
    case Default:  Insert(params, "Locator type", "Default" );  break;
    case CellTree: Insert(params, "Locator type", "CellTree"); break;
    case BSPTree:  Insert(params, "Locator type", "BSPTree" );  break;
    case OBBTree:  Insert(params, "Locator type", "OBBTree" );  break;
  }
  Insert(params, "No. of cells per node", _NumberOfCellsPerNode);
  if (_Sigma >= .0) {
    Insert(params, "Sigma", _Sigma);
  } else if (_MaxDistance > .0 && !IsInf(_MaxDistance)) {
    Insert(params, "Maximum distance", _MaxDistance);
  }
  return params;
}

// =============================================================================
// Correspondences
// =============================================================================

// -----------------------------------------------------------------------------
void irtkClosestCell::Initialize()
{
  // Initialize base class
  irtkPointCorrespondence::Initialize();

  // TODO: Implement cell locator for N-D points
  for (size_t i = 0; i < _TargetFeatures.size(); ++i) {
    if (_TargetFeatures[i]._Index != -1 && _TargetFeatures[i]._Slope != .0) {
      cerr << "irtkClosestCell::Initialize: Cannot use this correspondence map with extra features" << endl;
      exit(1);
    }
  }
  for (size_t i = 0; i < _SourceFeatures.size(); ++i) {
    if (_SourceFeatures[i]._Index != -1 && _SourceFeatures[i]._Slope != .0) {
      cerr << "irtkClosestCell::Initialize: Cannot use this correspondence map with extra features" << endl;
      exit(1);
    }
  }

  // Ensure that point sets contain cells
  if (_FromTargetToSource && _Source->NumberOfCells() == 0) {
    cerr << "irtkClosestCell::Initialize: Source dataset has no cells!" << endl;
    exit(1);
  }
  if (_FromSourceToTarget && _Target->NumberOfCells() == 0) {
    cerr << "irtkClosestCell::Initialize: Target dataset has no cells!" << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
void irtkClosestCell::Update()
{
  this->Upgrade();
}

// -----------------------------------------------------------------------------
bool irtkClosestCell::Upgrade()
{
  bool changed = false;
  // Update target to source correspondences
  if (_FromTargetToSource) {
    if (UpdateCorrespondences::Run(_Target, _TargetSample, _Source,
                                   _LocatorType, _NumberOfCellsPerNode,
                                   _SourcePoints, _SourceDistance,
                                   _Source->Transformation())) {
      changed = true;
    }
  } else {
    _SourcePoints  .Clear();
    _SourceDistance.clear();
  }
  // Update source to target correspondences
  if (_FromSourceToTarget) {
    if (UpdateCorrespondences::Run(_Source, _SourceSample, _Target,
                                   _LocatorType, _NumberOfCellsPerNode,
                                   _TargetPoints, _TargetDistance,
                                   _Target->Transformation())) {
      changed = true;
    }
  } else {
    _TargetPoints  .Clear();
    _TargetDistance.clear();
  }
  // Update maximum distance threshold
  if (changed && _Sigma >= .0) {
    vector<double>::const_iterator d;
    // Robust evaluation of variance of corresponding point distances
    int    m = 0;
    double mean = .0, var = .0, delta;
    for (d = _TargetDistance.begin(); d != _TargetDistance.end(); ++d) {
      ++m;
      delta  = (*d) - mean;
      mean  += delta / m;
      var   += delta * ((*d) - mean);
    }
    for (d = _SourceDistance.begin(); d != _SourceDistance.end(); ++d) {
      ++m;
      delta  = (*d) - mean;
      mean  += delta / m;
      var   += delta * ((*d) - mean);
    }
    if (m > 1) var /= m - 1;
    else       var  = .0;
    // Set maximum distance to mean plus _Sigma times standard deviation
    _MaxDistance = mean + _Sigma * sqrt(var);
  }
  return changed;
}

// -----------------------------------------------------------------------------
bool irtkClosestCell::GetInputTargetPoint(int i, irtkPoint &p) const
{
  _TargetPoints.GetPoint(i, p);
  return _TargetDistance[i] <= _MaxDistance;
}

// -----------------------------------------------------------------------------
bool irtkClosestCell::GetInputSourcePoint(int i, irtkPoint &p) const
{
  _SourcePoints.GetPoint(i, p);
  return _SourceDistance[i] <= _MaxDistance;
}

// -----------------------------------------------------------------------------
bool irtkClosestCell::GetTargetPoint(int i, irtkPoint &p) const
{
  _TargetPoints.GetPoint(i, p);
  if (_Target->Transformation()) _Target->Transformation()->Transform(p);
  return _TargetDistance[i] <= _MaxDistance;
}

// -----------------------------------------------------------------------------
bool irtkClosestCell::GetSourcePoint(int i, irtkPoint &p) const
{
  _SourcePoints.GetPoint(i, p);
  if (_Source->Transformation()) _Source->Transformation()->Transform(p);
  return _SourceDistance[i] <= _MaxDistance;
}
