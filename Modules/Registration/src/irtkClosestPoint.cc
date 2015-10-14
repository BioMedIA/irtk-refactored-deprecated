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

#include <irtkClosestPoint.h>


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkClosestPoint::irtkClosestPoint()
:
  _Sigma(-1.0),
  _MaxDistance(numeric_limits<double>::infinity()),
  _MaxSquaredDistance(numeric_limits<double>::infinity())
{
}

// -----------------------------------------------------------------------------
irtkClosestPoint::irtkClosestPoint(const irtkClosestPoint &other)
:
  irtkPointCorrespondence(other),
  _Sigma(other._Sigma),
  _MaxDistance(other._MaxDistance),
  _MaxSquaredDistance(other._MaxSquaredDistance),
  _TargetIndex(other._TargetIndex),
  _SourceIndex(other._SourceIndex)
{
}

// -----------------------------------------------------------------------------
irtkPointCorrespondence *irtkClosestPoint::NewInstance() const
{
  return new irtkClosestPoint(*this);
}

// -----------------------------------------------------------------------------
irtkClosestPoint::~irtkClosestPoint()
{
}

// -----------------------------------------------------------------------------
irtkClosestPoint::TypeId irtkClosestPoint::Type() const
{
  return ClosestPoint;
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkClosestPoint::Set(const char *name, const char *value)
{
  if (strcmp(name, "Sigma") == 0) {
    return FromString(value, _Sigma);
  }
  if (strcmp(name, "Maximum distance") == 0) {
    return FromString(value, _MaxDistance);
  }
  return irtkPointCorrespondence::Set(name, value);
}

// -----------------------------------------------------------------------------
irtkParameterList irtkClosestPoint::Parameter() const
{
  irtkParameterList params = irtkPointCorrespondence::Parameter();
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
void irtkClosestPoint::Initialize()
{
  // Initialize base class
  irtkPointCorrespondence::Initialize();

  // Set maximum squared distance threshold
  if (_Sigma < .0) {
    if (_MaxDistance > .0 && !IsInf(_MaxDistance)) {
      _MaxSquaredDistance = _MaxDistance * _MaxDistance;
    } else {
      _MaxSquaredDistance = numeric_limits<double>::infinity();
    }
  }
}

// -----------------------------------------------------------------------------
void irtkClosestPoint::Update()
{
  // Update point features
  irtkPointCorrespondence::Update();

  // Find closest points
  if (_FromTargetToSource) {
    irtkPointLocator *locator = irtkPointLocator::New(_Source->PointSet(), _SourceSample, &_SourceFeatures);
    _SourceIndex = locator->FindClosestPoint(_Target->PointSet(), _TargetSample, &_TargetFeatures, &_SourceDistance);
    delete locator;
  } else {
    _SourceIndex   .clear();
    _SourceDistance.clear();
  }

  if (_FromSourceToTarget) {
    irtkPointLocator *locator = irtkPointLocator::New(_Target->PointSet(), _TargetSample, &_TargetFeatures);
    _TargetIndex = locator->FindClosestPoint(_Source->PointSet(), _SourceSample, &_SourceFeatures, &_TargetDistance);
    delete locator;
  } else {
    _TargetIndex   .clear();
    _TargetDistance.clear();
  }

  // Update maximum distance threshold
  if (_Sigma >= .0) {
    vector<double>::const_iterator d;
    // Robust evaluation of variance of corresponding squared point distances
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
    _MaxSquaredDistance = mean + _Sigma * sqrt(var);
  }
}

// -----------------------------------------------------------------------------
bool irtkClosestPoint::Upgrade()
{
  vector<int> target_index(_TargetIndex), source_index(_SourceIndex);
  this->Update();
  return _TargetIndex != target_index || _SourceIndex != source_index;
}

// -----------------------------------------------------------------------------
bool irtkClosestPoint::GetInputTargetPoint(int i, irtkPoint &p) const
{
  if (_TargetIndex[i] < 0) return false;
  _Target->GetInputPoint(_TargetIndex[i], p);
  return _TargetDistance[i] <= _MaxSquaredDistance;
}

// -----------------------------------------------------------------------------
bool irtkClosestPoint::GetInputSourcePoint(int i, irtkPoint &p) const
{
  if (_SourceIndex[i] < 0) return false;
  _Source->GetInputPoint(_SourceIndex[i], p);
  return _SourceDistance[i] <= _MaxSquaredDistance;
}

// -----------------------------------------------------------------------------
bool irtkClosestPoint::GetTargetPoint(int i, irtkPoint &p) const
{
  if (_TargetIndex[i] < 0) return false;
  _Target->GetPoint(_TargetIndex[i], p);
  return _TargetDistance[i] <= _MaxSquaredDistance;
}

// -----------------------------------------------------------------------------
bool irtkClosestPoint::GetSourcePoint(int i, irtkPoint &p) const
{
  if (_SourceIndex[i] < 0) return false;
  _Source->GetPoint(_SourceIndex[i], p);
  return _SourceDistance[i] <= _MaxSquaredDistance;
}

// -----------------------------------------------------------------------------
int irtkClosestPoint::GetTargetIndex(int i) const
{
  return _TargetIndex[i];
}

// -----------------------------------------------------------------------------
int irtkClosestPoint::GetSourceIndex(int i) const
{
  return _SourceIndex[i];
}
