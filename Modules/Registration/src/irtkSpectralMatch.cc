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

#include <irtkSpectralMatch.h>

#include <irtkEdgeTable.h>

using namespace irtk::polydata;


// =============================================================================
// Auxiliary functors
// =============================================================================

namespace irtkSpectralMatchUtils {

// -----------------------------------------------------------------------------
struct FindCorrespondingPoints
{
  typedef irtkSpectralMatch::FeatureList FeatureList;

  vtkPointSet       *_Target;
  const vector<int> *_TargetSample;
  const FeatureList *_TargetFeatures;
  vtkPointSet       *_Source;
  const vector<int> *_SourceSample;
  const FeatureList *_SourceFeatures;
  int                _NumberOfNeighbors;
  int                _NumberOfEigenmodes;
  double             _Sigma2;
  irtkPointSet      *_Points;

  const irtkEdgeTable *_SourceEdgeTable;
  irtkPointLocator    *_SourceLocator;

  void operator ()(const blocked_range<vtkIdType> &ids) const
  {
    double w, W, x[3], y[3];
    double *u = Allocate<double>(_NumberOfEigenmodes);

    if (_NumberOfNeighbors > 0) {
      vector<int>    index;
      vector<double> dist2;
      double         denom;
      for (vtkIdType i = ids.begin(); i != ids.end(); ++i) {
        // Find k nearest neighboring source vertices in the spectral domain
        irtkPointCorrespondence::GetPoint(u, _Target, _TargetSample, i, _TargetFeatures);
        index = _SourceLocator->FindClosestNPoints(_NumberOfNeighbors, u, &dist2);
        if (_SourceSample && _SourceSample->size() > 0) {
          for (size_t j = 0; j < index.size(); ++j) {
            index[j] = (*_SourceSample)[ index[j] ];
          }
        }
        denom = 2.0 * _Sigma2 * dist2.back();
        _Source->GetPoint(index[0], x);
        W = exp(- dist2[0] / denom);
        x[0] *= W, x[1] *= W, x[2] *= W;
        for (size_t j = 1; j < index.size(); ++j) {
          _Source->GetPoint(index[j], y);
          w = exp(- dist2[j] / denom);
          x[0] += w * y[0], x[1] += w * y[1], x[2] += w * y[2];
          W += w;
        }
        // Normalize by sum of weights
        x[0] /= W, x[1] /= W, x[2] /= W;
        _Points->SetPoint(i, x);
      }
    } else {
      double *v = Allocate<double>(_NumberOfEigenmodes);
      const int *n, *end;
      for (vtkIdType i = ids.begin(); i != ids.end(); ++i) {
        // Find source vertex with most similar spectral coordinates
        irtkPointCorrespondence::GetPoint(u, _Target, _TargetSample, i, _TargetFeatures);
        vtkIdType j = _SourceLocator->FindClosestPoint(u, &w);
        W = w = exp(-.5 * w / _Sigma2);
        _Source->GetPoint(j, x);
        x[0] *= w, x[1] *= w, x[2] *= w;
        // Add contribution of vertices adjacent to closest vertex
        for (_SourceEdgeTable->GetAdjacentPoints(j, n, end); n != end; ++n) {
          _Source->GetPoint(*n, y);
          irtkPointCorrespondence::GetPoint(v, _Source, _SourceSample, *n, _SourceFeatures);
          w = exp(-.5 * irtkPointCorrespondence::Distance2BetweenPoints(u, v, _NumberOfEigenmodes) / _Sigma2);
          x[0] += w * y[0], x[1] += w * y[1], x[2] += w * y[2];
          W += w;
        }
        // Normalize by sum of weights
        x[0] /= W, x[1] /= W, x[2] /= W;
        _Points->SetPoint(i, x);
      }
      Deallocate(v);
    }
    Deallocate(u);
  }
};


} // namespace irtkSpectralMatchUtils
using namespace irtkSpectralMatchUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkSpectralMatch::irtkSpectralMatch()
:
  _NumberOfNeighbors(10),
  _Sigma            (1.),
  _TargetLocator    (NULL),
  _SourceLocator    (NULL)
{
}

// -----------------------------------------------------------------------------
irtkSpectralMatch::irtkSpectralMatch(const irtkSpectralMatch &other)
:
  irtkPointCorrespondence(other),
  _NumberOfNeighbors(other._NumberOfNeighbors),
  _Sigma            (other._Sigma),
  _TargetPoints     (other._TargetPoints),
  _SourcePoints     (other._SourcePoints),
  _TargetLocator    (NULL),
  _SourceLocator    (NULL)
{
}

// -----------------------------------------------------------------------------
irtkPointCorrespondence *irtkSpectralMatch::NewInstance() const
{
  return new irtkSpectralMatch(*this);
}

// -----------------------------------------------------------------------------
irtkSpectralMatch::~irtkSpectralMatch()
{
  Delete(_TargetLocator);
  Delete(_SourceLocator);
}

// -----------------------------------------------------------------------------
irtkSpectralMatch::TypeId irtkSpectralMatch::Type() const
{
  return SpectralMatch;
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkSpectralMatch::Set(const char *name, const char *value)
{
  if (strcmp(name, "No. of spectral neighbors")    == 0 ||
      strcmp(name, "Number of spectral neighbors") == 0 ||
      strcmp(name, "No. of neighbors")             == 0 ||
      strcmp(name, "Number of neighbors")          == 0) {
    return FromString(value, _NumberOfNeighbors);
  }
  if (strcmp(name, "Sigma") == 0) {
    return FromString(value, _Sigma);
  }
  return irtkPointCorrespondence::Set(name, value);
}

// -----------------------------------------------------------------------------
irtkParameterList irtkSpectralMatch::Parameter() const
{
  irtkParameterList params = irtkPointCorrespondence::Parameter();
  Insert(params, "No. of spectral neighbors", ToString(_NumberOfNeighbors));
  Insert(params, "Sigma",                     ToString(_Sigma));
  return params;
}

// =============================================================================
// Correspondences
// =============================================================================

// -----------------------------------------------------------------------------
void irtkSpectralMatch::Initialize()
{
  // Check that all inputs are set
  if (_Target == NULL) {
    cerr << "irtkPointCorrespondence::Initialize: Target data set not set" << endl;
    exit(1);
  }
  if (_Source == NULL) {
    cerr << "irtkPointCorrespondence::Initialize: Source data set not set" << endl;
    exit(1);
  }

  // Delete possible previously instantiated locator objects
  Delete(_TargetLocator);
  Delete(_SourceLocator);

  // Select spectral point coordinates as features
  // (before base class is initialized which may compute them for us)
  size_t i;
  for (i = 0; i < _TargetFeatures.size(); ++i) {
    if (_TargetFeatures[i]._Name == "joint_eigenmodes") break;
  }
  if (i == _TargetFeatures.size()) {
    for (i = 0; i < _TargetFeatures.size(); ++i) {
      if (_TargetFeatures[i]._Name == "eigenmodes") break;
    }
    if (i == _TargetFeatures.size()) {
      if (GetPointDataIndexByCaseInsensitiveName(_Target->PointSet()->GetPointData(), "eigenmodes") >= 0) {
        _TargetFeatures.push_back(FeatureInfo("eigenmodes", -2));
      } else {
        // If not present, it is up to the superclass to compute these features
        _TargetFeatures.push_back(FeatureInfo("joint_eigenmodes", -2));
      }
    }
  }
  if (i != 0) swap(_TargetFeatures[0], _TargetFeatures[i]);
  _TargetFeatures.resize(1);

  for (i = 0; i < _SourceFeatures.size(); ++i) {
    if (_SourceFeatures[i]._Name == "joint_eigenmodes") break;
  }
  if (i == _SourceFeatures.size()) {
    for (i = 0; i < _SourceFeatures.size(); ++i) {
      if (_SourceFeatures[i]._Name == "eigenmodes") break;
    }
    if (i == _SourceFeatures.size()) {
      if (GetPointDataIndexByCaseInsensitiveName(_Source->PointSet()->GetPointData(), "eigenmodes") >= 0) {
        _SourceFeatures.push_back(FeatureInfo("eigenmodes", -2));
      } else {
        // If not present, may be computed by superclass
        _SourceFeatures.push_back(FeatureInfo("joint_eigenmodes", -2));
      }
    }
  }
  if (i != 0) swap(_SourceFeatures[0], _SourceFeatures[i]);
  _SourceFeatures.resize(1);

  // Initialize base class
  irtkPointCorrespondence::Initialize();

  // Ensure that datasets contain cells
  if (_NumberOfNeighbors <= 0) {
    if (_FromTargetToSource && _Source->NumberOfCells() == 0) {
      cerr << "irtkSpectralMatch::Initialize: Source has no cells!" << endl;
      exit(1);
    }
    if (_FromSourceToTarget && _Target->NumberOfCells() == 0) {
      cerr << "irtkSpectralMatch::Initialize: Target has no cells!" << endl;
      exit(1);
    }
  }

  // Build spectral point locators
  if (_FromSourceToTarget) {
    _TargetLocator = irtkPointLocator::New(_Target->PointSet(), _TargetSample, &_TargetFeatures);
  }
  if (_FromTargetToSource) {
    _SourceLocator = irtkPointLocator::New(_Source->PointSet(), _SourceSample, &_SourceFeatures);
  }
}

// -----------------------------------------------------------------------------
void irtkSpectralMatch::Update()
{
  if (_FromTargetToSource) {
    IRTK_START_TIMING();
    irtkEdgeTable edgeTable;
    if (_NumberOfNeighbors <= 0) edgeTable.Initialize(_Source->PointSet());
    _SourcePoints.Resize(_Target->NumberOfPoints());
    FindCorrespondingPoints search;
    search._Target             = _Target->PointSet();
    search._TargetSample       = _TargetSample;
    search._TargetFeatures     = &_TargetFeatures;
    search._Source             = _Source->InputPointSet();
    search._SourceSample       = _SourceSample;
    search._SourceFeatures     = &_SourceFeatures;
    search._SourceLocator      = _SourceLocator;
    search._SourceEdgeTable    = &edgeTable;
    search._NumberOfEigenmodes = _NumberOfFeatures;
    search._NumberOfNeighbors  = _NumberOfNeighbors;
    search._Sigma2             = _Sigma * _Sigma;
    search._Points             = &_SourcePoints;
    parallel_for(blocked_range<vtkIdType>(0, _SourcePoints.Size()), search);
    IRTK_DEBUG_TIMING(7, "finding corresponding source points");
  }
  if (_FromSourceToTarget) {
    IRTK_START_TIMING();
    irtkEdgeTable edgeTable;
    if (_NumberOfNeighbors <= 0) edgeTable.Initialize(_Target->PointSet());
    _TargetPoints.Resize(_Source->NumberOfPoints());
    FindCorrespondingPoints search;
    search._Target             = _Source->PointSet();
    search._TargetSample       = _SourceSample;
    search._TargetFeatures     = &_SourceFeatures;
    search._Source             = _Target->InputPointSet();
    search._SourceSample       = _TargetSample;
    search._SourceFeatures     = &_TargetFeatures;
    search._SourceLocator      = _TargetLocator;
    search._SourceEdgeTable    = &edgeTable;
    search._NumberOfEigenmodes = _NumberOfFeatures;
    search._NumberOfNeighbors  = _NumberOfNeighbors;
    search._Sigma2             = _Sigma * _Sigma;
    search._Points             = &_TargetPoints;
    parallel_for(blocked_range<vtkIdType>(0, _TargetPoints.Size()), search);
    IRTK_DEBUG_TIMING(7, "finding corresponding target points");
  }
}

// -----------------------------------------------------------------------------
bool irtkSpectralMatch::Upgrade()
{
  return false;
}

/// -----------------------------------------------------------------------------
bool irtkSpectralMatch::GetInputTargetPoint(int i, irtkPoint &p) const
{
  _TargetPoints.GetPoint(i, p);
  return true;
}

// -----------------------------------------------------------------------------
bool irtkSpectralMatch::GetInputSourcePoint(int i, irtkPoint &p) const
{
  _SourcePoints.GetPoint(i, p);
  return true;
}

// -----------------------------------------------------------------------------
bool irtkSpectralMatch::GetTargetPoint(int i, irtkPoint &p) const
{
  _TargetPoints.GetPoint(i, p);
  if (_Target->Transformation()) _Target->Transformation()->Transform(p);
  return true;
}

// -----------------------------------------------------------------------------
bool irtkSpectralMatch::GetSourcePoint(int i, irtkPoint &p) const
{
  _SourcePoints.GetPoint(i, p);
  if (_Source->Transformation()) _Source->Transformation()->Transform(p);
  return true;
}
