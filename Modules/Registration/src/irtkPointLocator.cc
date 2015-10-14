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

#include <irtkPointLocator.h>

#include <vtkIdList.h>
#include <vtkPolyData.h>

// =============================================================================
// FLANN helpers
// =============================================================================
#ifdef HAS_FLANN

// -----------------------------------------------------------------------------
flann::Matrix<irtkPointLocator::FlannType>
irtkPointLocator::FlannMatrix(vtkPointSet       *dataset,
                              const vector<int> *sample,
                              const FeatureList *features)
{
  const int n = GetNumberOfPoints(dataset, sample);
  double *point = Allocate<double>(_PointDimension);
  flann::Matrix<FlannType> matrix(Allocate<FlannType>(n * _PointDimension), n, _PointDimension);
  for (int i = 0; i < n; ++i) {
    GetPoint(point, dataset, sample, i, features);
    FlannType *v = matrix[i];
    for (int j = 0; j < _PointDimension; ++j, ++v) {
      (*v) = static_cast<FlannType>(point[j]);
    }
  }
  Deallocate(point);
  return matrix;
}

// -----------------------------------------------------------------------------
flann::Matrix<irtkPointLocator::FlannType> irtkPointLocator::FlannMatrix(double *point)
{
  float *m = Allocate<FlannType>(_PointDimension);
  flann::Matrix<FlannType> matrix(m, 1, _PointDimension);
  float *v = matrix[0];
  for (int j = 0; j < _PointDimension; ++j, ++v) {
    (*v) = static_cast<FlannType>(point[j]);
  }
  return matrix;
}

#endif
// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkPointLocator::irtkPointLocator()
:
  _DataSet(NULL), _Sample(NULL), _NumberOfPoints(0), _PointDimension(0)
#ifdef HAS_FLANN
  , _FlannTree(flann::KDTreeSingleIndexParams(10))
#endif
{
}

// -----------------------------------------------------------------------------
irtkPointLocator::~irtkPointLocator()
{
}

// -----------------------------------------------------------------------------
void irtkPointLocator::Initialize()
{
  // Check inputs
  if (!_DataSet) {
    cerr << "irtkPointLocator: Missing dataset!" << endl;
    exit(1);
  }
  _NumberOfPoints = GetNumberOfPoints(_DataSet, _Sample);
  if (_NumberOfPoints == 0) {
    cerr << "irtkPointLocator: No points in search tree!" << endl;
    exit(1);
  }
  _PointDimension = GetPointDimension(_DataSet, &_Features);
  if (_PointDimension == 0) {
    cerr << "irtkPointLocator: Point feature vector size is zero!" << endl;
    exit(1);
  }
  // Build VTK locator for 3-D feature vectors
  if (_PointDimension <= 3) {
    // Dataset of 2/3-D sample feature points
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(_NumberOfPoints);
    double point[3] = {0};
    for (int i = 0; i < _NumberOfPoints; ++i) {
      GetPoint(point, _DataSet, _Sample, i, &_Features);
      points->SetPoint(i, point);
    }
    vtkSmartPointer<vtkPolyData> dataset = vtkSmartPointer<vtkPolyData>::New();
    dataset->SetPoints(points);
    // Note: vtkOctreeLocator preferred over vtkKdTree because the
    //       latter is not thread safe even after BuildLocator was called!
    _VtkLocator = vtkSmartPointer<vtkOctreePointLocator>::New();
    _VtkLocator->SetDataSet(dataset);
    _VtkLocator->BuildLocator();
  } else {
    _VtkLocator = NULL;
    // Build FLANN tree for N-D feature vectors
#ifdef HAS_FLANN
    flann::Matrix<FlannType> points = FlannMatrix(_DataSet, _Sample, &_Features);
    _FlannTree.buildIndex(points);
    delete[] points.ptr();
    // Build ITK tree for N-D feature vectors
#elif HAS_ITK
    // TODO: Build ITK tree for N-D feature vectors
#endif
  }
}

// -----------------------------------------------------------------------------
irtkPointLocator *irtkPointLocator::New(vtkPointSet       *dataset,
                                        const vector<int> *sample,
                                        const FeatureList *features)
{
  irtkPointLocator *locator = new irtkPointLocator();
  locator->DataSet(dataset);
  locator->Sample(sample);
  if (features) locator->Features(*features);
  locator->Initialize();
  return locator;
}

// =============================================================================
// Closest point
// =============================================================================

// -----------------------------------------------------------------------------
int irtkPointLocator::FindClosestPoint(double *point, double *dist2)
{
  // ---------------------------------------------------------------------------
  // Using VTK
  if (_VtkLocator) {
    vtkIdType j = _VtkLocator->FindClosestPoint(point);
    if (dist2) {
      double p[3] = {0};
      _DataSet->GetPoint(j, p);
      *dist2 = Distance2BetweenPoints(p, point);
    }
    return static_cast<int>(j);
  }

  // ---------------------------------------------------------------------------
  // Using FLANN
#ifdef HAS_FLANN
  vector<vector<int> > indices;
  vector<vector<FlannType> > dists;
  flann::Matrix<FlannType> queries = FlannMatrix(point);
  _FlannTree.knnSearch(queries, indices, dists, 1,
                       flann::SearchParams(flann::FLANN_CHECKS_UNLIMITED));
  delete[] queries.ptr();
  if (dist2) (*dist2) = dists[0][0];
  return indices[0][0];

  // ---------------------------------------------------------------------------
  // Brute force
#else
  double *p = new double[_PointDimension];
  int    idx, index = -1;
  double d2, mind2 = numeric_limits<double>::infinity();
  for (int i = 0; i < _NumberOfPoints; ++i) {
    idx = GetPointIndex(_DataSet, _Sample, i);
    GetPoint(p, _DataSet, idx, &_Features);
    d2 = Distance2BetweenPoints(p, point, _PointDimension);
    if (d2 < mind2) {
      index = idx;
      mind2 = d2;
    }
  }
  delete[] p;
  if (dist2) (*dist2) = mind2;
  return index;
#endif
}

// -----------------------------------------------------------------------------
namespace irtkPointLocatorUtils {
struct FindClosestPoint
{
  vtkPointSet                         *_DataSet;
  const vector<int>                   *_Sample;
  const irtkPointLocator::FeatureList *_Features;
  irtkPointLocator                    *_Locator;
  vector<int>                         *_Index;
  vector<double>                      *_Dist2;

  void operator ()(const blocked_range<int> &idx) const
  {
    if (_Dist2) {
      for (int i = idx.begin(); i != idx.end(); ++i) {
        (*_Index)[i] = _Locator->FindClosestPoint(_DataSet, _Sample, i, _Features, &((*_Dist2)[i]));
      }
    } else {
      for (int i = idx.begin(); i != idx.end(); ++i) {
        (*_Index)[i] = _Locator->FindClosestPoint(_DataSet, _Sample, i, _Features);
      }
    }
  }
};
} // namespace irtkPointLocatorUtils

// -----------------------------------------------------------------------------
vector<int> irtkPointLocator
::FindClosestPoint(vtkPointSet *dataset, const vector<int> *sample,
                   const FeatureList *features, vector<double> *dist2)
{
#ifdef HAS_FLANN
  if (!_VtkLocator) {
    vector<vector<int   > >  indices;
    vector<vector<double> > *dists;
    dists   = (dist2 ? new vector<vector<double> >() : NULL);
    indices = FindClosestNPoints(1, dataset, sample, features, dists);
    vector<int> index(indices.size());
    for (size_t i = 0; i < indices.size(); ++i) {
      index[i] = indices[i][0];
    }
    if (dist2) {
      dist2->resize(dists->size());
      for (size_t i = 0; i < dists->size(); ++i) {
        (*dist2)[i] = (*dists)[i][0];
      }
    }
    delete dists;
    return index;
  }
#endif

  vector<int> index(GetNumberOfPoints(dataset, sample));
  if (dist2) dist2->resize(index.size());
  irtkPointLocatorUtils::FindClosestPoint query;
  query._DataSet  = dataset;
  query._Sample   = sample;
  query._Features = features;
  query._Locator  = this;
  query._Index    = &index;
  query._Dist2    = dist2;
  parallel_for(blocked_range<int>(0, index.size()), query);
  return index;
}

// =============================================================================
// Nearest neighbors
// =============================================================================

// -----------------------------------------------------------------------------
vector<int> irtkPointLocator::FindClosestNPoints(int k, double *point, vector<double> *dist2)
{
  // ---------------------------------------------------------------------------
  // Using VTK
  if (_VtkLocator) {
    double p[3];
    vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
    _VtkLocator->FindClosestNPoints(k, point, ids);
    vector<int> indices(ids->GetNumberOfIds());
    if (dist2) dist2->resize(ids->GetNumberOfIds());
    for (vtkIdType i = 0; i < ids->GetNumberOfIds(); ++i) {
      indices[i] = ids->GetId(i);
      if (dist2) {
        _VtkLocator->GetDataSet()->GetPoint(indices[i], p);
        (*dist2)[i] = Distance2BetweenPoints(point, p);
      }
    }
    return indices;
  }

  // ---------------------------------------------------------------------------
  // Using FLANN
#ifdef HAS_FLANN
  vector<vector<int> > indices;
  vector<vector<FlannType> > dists;
  flann::Matrix<FlannType> queries = FlannMatrix(point);
  _FlannTree.knnSearch(queries, indices, dists, k,
                       flann::SearchParams(flann::FLANN_CHECKS_UNLIMITED));
  delete[] queries.ptr();
  if (dist2) {
    dist2->resize(dists[0].size());
    for (size_t i = 0; i < dists[0].size(); ++i) {
      (*dist2)[i] = dists[0][i];
    }
  }
  return indices[0];
#else

  // ---------------------------------------------------------------------------
  // Brute force
  struct Comp
  {
    bool operator()(const pair<double, int> &a, const pair<double, int> &b)
    {
      return a.first > b.first;
    }
  } comp;

  double *p = Allocate<double>(_PointDimension);
  vector<pair<double, int> > dists(_NumberOfPoints);
  for (int i = 0; i < _NumberOfPoints; ++i) {
    int idx = GetPointIndex(_DataSet, _Sample, i);
    GetPoint(p, _DataSet, idx, &_Features);
    dists[i] = make_pair(Distance2BetweenPoints(point, p, _PointDimension), idx);
  }
  make_heap(dists.begin(), dists.end(), comp);
  vector<int> indices(k);
  for (int i = 0; i < k; ++i) {
    pair<double, int> &min = dists.front();
    if (dist2) (*dist2)[i] = min.first;
    indices[i] = min.second;
    pop_heap(dists.begin(), dists.end(), comp);
    dists.pop_back();
  }
  Deallocate(p);
  return indices;
#endif
}

// -----------------------------------------------------------------------------
namespace irtkPointLocatorUtils {
struct FindClosestNPoints
{
  vtkPointSet                         *_DataSet;
  const vector<int>                   *_Sample;
  const irtkPointLocator::FeatureList *_Features;
  irtkPointLocator                    *_Locator;
  vector<vector<int> >                *_Indices;
  vector<vector<double> >             *_Dist2;
  int                                  _K;

  void operator ()(const blocked_range<int> &idx) const
  {
    if (_Dist2) {
      for (int i = idx.begin(); i != idx.end(); ++i) {
        (*_Indices)[i] = _Locator->FindClosestNPoints(_K, _DataSet, _Sample, i, _Features, &((*_Dist2)[i]));
      }
    } else {
      for (int i = idx.begin(); i != idx.end(); ++i) {
        (*_Indices)[i] = _Locator->FindClosestNPoints(_K, _DataSet, _Sample, i, _Features);
      }
    }
  }
};
} // namespace irtkPointLocatorUtils

// -----------------------------------------------------------------------------
vector<vector<int> > irtkPointLocator
::FindClosestNPoints(int k, vtkPointSet *dataset, const vector<int> *sample,
                     const FeatureList *features, vector<vector<double> > *dist2)
{
  irtkAssert(GetPointDimension(dataset, features) == _PointDimension,
             "Query points must have same dimension as feature points");

  if (k > _NumberOfPoints) {
    cerr << "irtkPointLocator::FindClosestNPoints: Cannot find more points than there are in total" << endl;
    exit(1);
  }

#ifdef HAS_FLANN
  if (!_VtkLocator) {
    vector<vector<int> > indices;
    vector<vector<FlannType> > dists;
    flann::Matrix<FlannType> queries = FlannMatrix(dataset, sample, features);
    _FlannTree.knnSearch(queries, indices, dists, k,
                         flann::SearchParams(flann::FLANN_CHECKS_UNLIMITED));
    delete[] queries.ptr();
    if (dist2) {
      dist2->resize(dists.size());
      for (size_t i = 0; i < dists.size(); ++i) {
        vector<double> &row = (*dist2)[i];
        row.resize(dists[i].size());
        for (size_t j = 0; j < dists[i].size(); ++j) {
          row[j] = static_cast<double>(dists[i][j]);
        }
      }
    }
    return indices;
  }
#endif

  vector<vector<int> > indices(GetNumberOfPoints(dataset, sample));
  if (dist2) dist2->resize(indices.size());
  irtkPointLocatorUtils::FindClosestNPoints query;
  query._DataSet  = dataset;
  query._Sample   = sample;
  query._Features = features;
  query._Locator  = this;
  query._Indices  = &indices;
  query._Dist2    = dist2;
  query._K        = k;
  parallel_for(blocked_range<int>(0, indices.size()), query);
  return indices;
}

// =============================================================================
// Radius search
// =============================================================================

// -----------------------------------------------------------------------------
vector<int> irtkPointLocator::FindPointsWithinRadius(double radius, double *point, vector<double> *dist2)
{
  // ---------------------------------------------------------------------------
  // Using VTK
  if (_VtkLocator) {
    double p[3] = {0};
    vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
    _VtkLocator->FindPointsWithinRadius(radius, point, ids);
    vector<int> indices(ids->GetNumberOfIds());
    if (dist2) dist2->resize(ids->GetNumberOfIds());
    for (vtkIdType i = 0; i < ids->GetNumberOfIds(); ++i) {
      indices[i] = ids->GetId(i);
      if (dist2) {
        _VtkLocator->GetDataSet()->GetPoint(indices[i], p);
        (*dist2)[i] = Distance2BetweenPoints(point, p);
      }
    }
    return indices;
  }

  // ---------------------------------------------------------------------------
  // Using FLANN
#ifdef HAS_FLANN

  vector<vector<int> > indices;
  vector<vector<FlannType> > dists;
  flann::Matrix<FlannType> queries = FlannMatrix(point);
  _FlannTree.radiusSearch(queries, indices, dists, radius * radius,
                          flann::SearchParams(flann::FLANN_CHECKS_UNLIMITED));
  delete[] queries.ptr();
  if (dist2) {
    dist2->resize(dists[0].size());
    for (size_t i = 0; i < dist2->size(); ++i) {
      (*dist2)[i] = dists[0][i];
    }
  }
  return indices[0];

  // ---------------------------------------------------------------------------
  // Brute force
#else
  struct Comp
  {
    bool operator()(const pair<double, int> &a, const pair<double, int> &b)
    {
      return a.first > b.first;
    }
  } comp;

  const double maxdist2 = radius * radius;
  double *p = Allocate<double>(_PointDimension);
  vector<pair<double, int> > dists(_NumberOfPoints);
  int idx, k = 0;
  for (int i = 0; i < _NumberOfPoints; ++i) {
    idx = GetPointIndex(_DataSet, _Sample, i);
    GetPoint(p, _DataSet, idx, &_Features);
    dists[i] = make_pair(Distance2BetweenPoints(point, p, _PointDimension), idx);
    if (dists[i].first <= maxdist2) ++k;
  }
  make_heap(dists.begin(), dists.end(), comp);
  vector<int> indices(k);
  if (dist2) dist2->resize(k);
  for (int i = 0; i < k; ++i) {
    pair<double, int> &min = dists.front();
    if (dist2) (*dist2)[i] = min.first;
    indices[i] = min.second;
    pop_heap(dists.begin(), dists.end(), comp);
    dists.pop_back();
  }
  Deallocate(p);
  return indices;
#endif
}

// -----------------------------------------------------------------------------
namespace irtkPointLocatorUtils {
struct FindPointsWithinRadius
{
  vtkPointSet                         *_DataSet;
  const vector<int>                   *_Sample;
  const irtkPointLocator::FeatureList *_Features;
  irtkPointLocator                    *_Locator;
  vector<vector<int> >                *_Indices;
  vector<vector<double> >             *_Dist2;
  double                               _Radius;

  void operator ()(const blocked_range<int> &idx) const
  {
    if (_Dist2) {
      for (int i = idx.begin(); i != idx.end(); ++i) {
        (*_Indices)[i] = _Locator->FindPointsWithinRadius(_Radius, _DataSet, _Sample, i, _Features, &((*_Dist2)[i]));
      }
    } else {
      for (int i = idx.begin(); i != idx.end(); ++i) {
        (*_Indices)[i] = _Locator->FindPointsWithinRadius(_Radius, _DataSet, _Sample, i, _Features);
      }
    }
  }
};
} // namespace irtkPointLocatorUtils

// -----------------------------------------------------------------------------
vector<vector<int> > irtkPointLocator
::FindPointsWithinRadius(double radius, vtkPointSet             *dataset,
                                        const vector<int>       *sample,
                                        const FeatureList       *features,
                                        vector<vector<double> > *dist2)
{
  irtkAssert(GetPointDimension(dataset, features) == _PointDimension,
             "Query points must have same dimension as feature points");

#ifdef HAS_FLANN
  if (!_VtkLocator) {
    vector<vector<int> > indices;
    vector<vector<FlannType> > dists;
    flann::Matrix<FlannType> queries = FlannMatrix(dataset, sample, features);
    _FlannTree.radiusSearch(queries, indices, dists, radius,
                            flann::SearchParams(flann::FLANN_CHECKS_UNLIMITED));
    delete[] queries.ptr();
    if (dist2) {
      dist2->resize(dists.size());
      for (size_t i = 0; i < dists.size(); ++i) {
        vector<double> &row = (*dist2)[i];
        row.resize(dists[0].size());
        for (size_t j = 0; j < row.size(); ++j) {
          row[j] = static_cast<double>(dists[i][j]);
        }
      }
    }
    return indices;
  }
#endif

  vector<vector<int> > indices(GetNumberOfPoints(dataset, sample));
  if (dist2) dist2->resize(indices.size());
  irtkPointLocatorUtils::FindPointsWithinRadius query;
  query._DataSet  = dataset;
  query._Sample   = sample;
  query._Features = features;
  query._Locator  = this;
  query._Indices  = &indices;
  query._Dist2    = dist2;
  query._Radius   = radius;
  parallel_for(blocked_range<int>(0, indices.size()), query);
  return indices;
}
