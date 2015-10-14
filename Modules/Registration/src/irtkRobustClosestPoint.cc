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

#include <irtkRobustClosestPoint.h>


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkRobustClosestPoint::irtkRobustClosestPoint()
:
  _Sigma(3.0)
{
}

// -----------------------------------------------------------------------------
irtkRobustClosestPoint::irtkRobustClosestPoint(const irtkRegisteredPointSet *target,
                                               const irtkRegisteredPointSet *source)
:
  _Sigma(3.0)
{
  Target(target);
  Source(source);
  Initialize();
}

// -----------------------------------------------------------------------------
irtkRobustClosestPoint::irtkRobustClosestPoint(const irtkRobustClosestPoint &other)
:
  irtkFuzzyCorrespondence(other),
  _Sigma(other._Sigma)
{
}

// -----------------------------------------------------------------------------
irtkPointCorrespondence *irtkRobustClosestPoint::NewInstance() const
{
  return new irtkRobustClosestPoint(*this);
}

// -----------------------------------------------------------------------------
irtkRobustClosestPoint::~irtkRobustClosestPoint()
{
}

// -----------------------------------------------------------------------------
irtkRobustClosestPoint::TypeId irtkRobustClosestPoint::Type() const
{
  return RobustClosestPoint;
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkRobustClosestPoint::Set(const char *name, const char *value)
{
  // Standard deviation of outlier rejection
  if (strcmp(name, "Sigma") == 0) {
    return FromString(value, _Sigma);
  }
  return irtkFuzzyCorrespondence::Set(name, value);
}

// -----------------------------------------------------------------------------
irtkParameterList irtkRobustClosestPoint::Parameter() const
{
  irtkParameterList params = irtkFuzzyCorrespondence::Parameter();
  Insert(params, "Sigma", _Sigma);
  return params;
}

// =============================================================================
// Correspondences
// =============================================================================

// -----------------------------------------------------------------------------
void irtkRobustClosestPoint::CalculateWeights()
{
  IRTK_START_TIMING();

  // Initialize weight matrix
  _Weight.Initialize(_M, _N);

  // Find closest points
  vector<int   > corr12, corr21;
  vector<double> dist12, dist21;

  irtkPointLocator *source_locator;
  source_locator = irtkPointLocator::New(_Source->PointSet(), _SourceSample, &_SourceFeatures);
  corr12 = source_locator->FindClosestPoint(_Target->PointSet(), _TargetSample, &_TargetFeatures, &dist12);
  delete source_locator;

  irtkPointLocator *target_locator;
  target_locator = irtkPointLocator::New(_Target->PointSet(), _TargetSample, &_TargetFeatures);
  corr21 = target_locator->FindClosestPoint(_Source->PointSet(), _SourceSample, &_SourceFeatures, &dist21);
  delete target_locator;

  // Allocate lists for non-zero weight entries
  const int nentries = (_Weight.Layout() == WeightMatrix::CRS ? _M : _N);
  WeightMatrix::Entries *entries = new WeightMatrix::Entries[nentries];
  for (int i = 0; i < nentries; ++i) entries[i].reserve(2);

  // Symmetric closest point matching with outlier rejection
  if (_Sigma > .0) {

    double max_dist, dist_mean = .0, dist_mean2 = .0;
    for (vector<double>::iterator dist = dist12.begin(); dist != dist12.end(); ++dist) {
      dist_mean  += (*dist);
      dist_mean2 += (*dist) * (*dist);
    }
    for (vector<double>::iterator dist = dist21.begin(); dist != dist21.end(); ++dist) {
      dist_mean  += (*dist);
      dist_mean2 += (*dist) * (*dist);
    }
    dist_mean  /= (dist12.size() + dist21.size());
    dist_mean2 /= (dist12.size() + dist21.size());

    max_dist = dist_mean + _Sigma * sqrt(dist_mean2 - dist_mean * dist_mean);

    if (_Weight.Layout() == WeightMatrix::CRS) {
      for (size_t r = 0; r < corr12.size(); ++r) {
        if (dist12[r] <= max_dist) entries[r].push_back(make_pair(corr12[r], 1.0));
      }
      for (size_t c = 0; c < corr21.size(); ++c) {
        if (dist21[c] <= max_dist) entries[corr21[c]].push_back(make_pair(c, 1.0));
      }
    } else {
      for (size_t r = 0; r < corr12.size(); ++r) {
        if (dist12[r] <= max_dist) entries[corr12[r]].push_back(make_pair(r, 1.0));
      }
      for (size_t c = 0; c < corr21.size(); ++c) {
        if (dist21[c] <= max_dist) entries[c].push_back(make_pair(corr21[c], 1.0));
      }
    }

  // Symmetric closest point matching without outlier rejection
  } else {

    if (_Weight.Layout() == WeightMatrix::CRS) {
      for (size_t r = 0; r < corr12.size(); ++r) {
        entries[r].push_back(make_pair(corr12[r], 1.0));
      }
      for (size_t c = 0; c < corr21.size(); ++c) {
        entries[corr21[c]].push_back(make_pair(c, 1.0));
      }
    } else {
      for (size_t r = 0; r < corr12.size(); ++r) {
        entries[corr12[r]].push_back(make_pair(r, 1.0));
      }
      for (size_t c = 0; c < corr21.size(); ++c) {
        entries[c].push_back(make_pair(corr21[c], 1.0));
      }
    }

  }

  // Initialize weight matrix
  _Weight.Initialize(_M, _N, entries);

  // Normalize weights to sum up to 1
  _Weight *= 0.5;

  IRTK_DEBUG_TIMING(6, "calculating correspondence weights");
}
