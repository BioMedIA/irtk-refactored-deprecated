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

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

#include <irtkFuzzyCorrespondence.h>

#include <vtkSmartPointer.h>
#include <irtkParallel.h>
#include <irtkPoint.h>
#include <irtkTransformation.h>


// =============================================================================
// Utilities
// =============================================================================

namespace irtkFuzzyCorrespondenceUtils {


// -----------------------------------------------------------------------------
class UpdateClusters
{
  typedef WeightMatrix::Entries   WeightEntries;
  typedef WeightEntries::iterator WeightIterator;

  const irtkRegisteredPointSet *_DataSet;
  const vector<int>            *_Sample;
  const WeightMatrix           *_Weight;
  bool                          _Transpose;
  irtkPointSet                 *_Cluster;
  vector<bool>                 *_Outlier;
  bool                          _Native;

public:

  void operator()(const blocked_range<int> &re) const
  {
    double         wout, wsum;
    irtkPoint      p, pavg;
    WeightEntries  weight;
    WeightIterator wend;

    const int npts = irtkPointCorrespondence::GetNumberOfPoints(_DataSet, _Sample);

    for (int i = re.begin(); i != re.end(); ++i) {
      if (_Transpose) _Weight->GetCol(i, weight);
      else            _Weight->GetRow(i, weight);
      wend = weight.end();
      wout = .0;
      while (wend != weight.begin()) {
        --wend;
        if (wend->first < npts) {
          ++wend;
          break;
        }
        wout += wend->second;
      }
      wsum = pavg._x = pavg._y = pavg._z = .0;
      for (WeightIterator w = weight.begin(); w != wend; ++w) {
        if (_Sample) w->first = (*_Sample)[w->first];
        if (_Native) _DataSet->GetInputPoint(w->first, p);
        else         _DataSet->GetPoint     (w->first, p);
        pavg += p * w->second;
        wsum +=     w->second;
      }
      if (wsum == .0) {
        if (_Outlier) (*_Outlier)[i] = true;
      } else {
        pavg /= wsum;
        if (_Outlier) (*_Outlier)[i] = (wout / (wout + wsum) > .75);
      }
      _Cluster->SetPoint(i, pavg);
    }
  }

  static void Run(const irtkRegisteredPointSet *dataset,
                  const vector<int>            *sample,
                  const WeightMatrix           &weight,
                  bool                          transpose,
                  irtkPointSet                 &cluster,
                  vector<bool>                 *outlier = NULL,
                  bool                          native  = false)
  {
    UpdateClusters body;
    body._DataSet   = dataset;
    body._Sample    = sample;
    body._Weight    = &weight;
    body._Transpose = transpose;
    body._Cluster   = &cluster;
    body._Outlier   = outlier;
    if (outlier) outlier->resize(cluster.Size());
    blocked_range<int> range(0, cluster.Size());
    parallel_for(range, body);
  }
};


} // namespace irtkFuzzyCorrespondenceUtils
using namespace irtkFuzzyCorrespondenceUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkFuzzyCorrespondence::irtkFuzzyCorrespondence()
:
  _MinWeight    (0.001),
  _GaussianNoise(-1.0) // -1: Set to default in Initialize
                       //  0: Do not add any Gaussian noise
{
}

// -----------------------------------------------------------------------------
irtkFuzzyCorrespondence
::irtkFuzzyCorrespondence(const irtkFuzzyCorrespondence &other)
:
  irtkPointCorrespondence(other),
  _MinWeight          (other._MinWeight),
  _Weight             (other._Weight),
  _InputTargetClusters(other._InputTargetClusters),
  _InputSourceClusters(other._InputSourceClusters),
  _TargetClusters     (other._TargetClusters),
  _SourceClusters     (other._SourceClusters),
  _TargetOutlier      (other._TargetOutlier),
  _SourceOutlier      (other._SourceOutlier),
  _GaussianNoise      (other._GaussianNoise)
{
}

// -----------------------------------------------------------------------------
irtkFuzzyCorrespondence::~irtkFuzzyCorrespondence()
{
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkFuzzyCorrespondence::Set(const char *name, const char *value)
{
  // Standard deviation of Gaussian noise
  if (strcmp(name, "Correspondence noise") == 0) {
    return FromString(value, _GaussianNoise);
  }
  return irtkPointCorrespondence::Set(name, value);
}

// -----------------------------------------------------------------------------
irtkParameterList irtkFuzzyCorrespondence::Parameter() const
{
  irtkParameterList params = irtkPointCorrespondence::Parameter();
  Insert(params, "Correspondence noise", ToString(_GaussianNoise));
  return params;
}

// =============================================================================
// Correspondences
// =============================================================================

// -----------------------------------------------------------------------------
void irtkFuzzyCorrespondence::Initialize()
{
  // Initialize base class
  irtkPointCorrespondence::Initialize();

  // Set default standard deviation of Gaussian noise
  if (_GaussianNoise < .0) _GaussianNoise = 2.0 * _MinWeight / (_M + _N);

  // Initialize this class
  irtkFuzzyCorrespondence::Init();
}

// -----------------------------------------------------------------------------
void irtkFuzzyCorrespondence::Reinitialize()
{
  // Reinitialize base class
  irtkPointCorrespondence::Reinitialize();

  // Reinitialize this class
  irtkFuzzyCorrespondence::Init();
}

// -----------------------------------------------------------------------------
void irtkFuzzyCorrespondence::Init()
{
  // Allocate memory for estimated point clusters
  if (_FromTargetToSource) {
    if (_Source->Transformation()) {
      _InputSourceClusters.Reserve(_M);
      _InputSourceClusters.Resize (_M);
    } else {
      _InputSourceClusters.Clear();
    }
    _SourceClusters.Reserve(_M);
    _SourceClusters.Resize (_M);
    _SourceOutlier .reserve(_M);
    _SourceOutlier .resize (_M, false);
  }
  if (_FromSourceToTarget) {
    if (_Target->Transformation()) {
      _InputTargetClusters.Reserve(_N);
      _InputTargetClusters.Resize (_N);
    } else {
      _InputTargetClusters.Clear();
    }
    _TargetClusters.Reserve(_N);
    _TargetClusters.Resize (_N);
    _TargetOutlier .reserve(_N);
    _TargetOutlier .resize (_N, false);
  }
}

// -----------------------------------------------------------------------------
void irtkFuzzyCorrespondence::Update()
{
  IRTK_START_TIMING();

  // Update point features
  irtkPointCorrespondence::Update();

  // Recalculate correspondence matrix
  this->CalculateWeights();

  // Iterative column and row normalization
  this->NormalizeWeights();

  // Add some normally distributed noise to break symmetry
  this->AddGaussianNoise();

  // Update cluster centers
  this->CalculateClusters();

  IRTK_DEBUG_TIMING(5, "update of correspondences");
}

// -----------------------------------------------------------------------------
bool irtkFuzzyCorrespondence::Upgrade()
{
  // Check current ratio of outliers
  vector<bool>::const_iterator it;
  double target_outlier_ratio = .0;
  double source_outlier_ratio = .0;
  int    noutlier;

  if (!_TargetOutlier.empty()) {
    noutlier = 0;
    for (it = _TargetOutlier.begin(); it != _TargetOutlier.end(); ++it) {
      if (*it) ++noutlier;
    }
    target_outlier_ratio = double(noutlier) / _TargetOutlier.size();
  }

  if (!_SourceOutlier.empty()) {
    noutlier = 0;
    for (it = _SourceOutlier.begin(); it != _SourceOutlier.end(); ++it) {
      if (*it) ++noutlier;
    }
    source_outlier_ratio = double(noutlier) / _SourceOutlier.size();
  }

  double outlier_ratio = target_outlier_ratio + source_outlier_ratio;
  if (!_TargetOutlier.empty() && !_SourceOutlier.empty()) outlier_ratio *= 0.5;

  // Interrupt annealing process if number of outliers is too high
  // This happens in particular when the search range (i.e., temperature)
  // dropped too low already such that almost no points are within this range.
  if (outlier_ratio > 0.2) {
    IRTK_LOG_EVENT("Number of outliers exceeds 20% of the total number of samples\n");
    return false;
  }

  return true;
}

// -----------------------------------------------------------------------------
void irtkFuzzyCorrespondence::NormalizeWeights()
{
  if (_Weight.Rows() == _M && _Weight.Cols() == _N) return;
  IRTK_START_TIMING();
  irtkFuzzyCorrespondenceUtils::NormalizeWeights(_Weight, _M, _N);
  IRTK_DEBUG_TIMING(6, "normalization of correspondence weights (" << _M << "x" << _N << ")");
}

// -----------------------------------------------------------------------------
void irtkFuzzyCorrespondence::AddGaussianNoise()
{
  if (_GaussianNoise <= .0) return;
  IRTK_START_TIMING();

  typedef boost::normal_distribution<double> NormalDistribution;

  static boost::mt19937 rng;
  NormalDistribution nd(0.0, _GaussianNoise);
  boost::variate_generator<boost::mt19937&, NormalDistribution> noise(rng, nd);

  WeightMatrix::EntryType *w = _Weight.RawPointer();
  for (int i = 0; i < _Weight.NNZ(); ++i, ++w) (*w) += noise();

  IRTK_DEBUG_TIMING(6, "addition of Gaussian noise (NNZ=" << _Weight.NNZ() << ")");
}

// -----------------------------------------------------------------------------
void irtkFuzzyCorrespondence::CalculateClusters()
{
  if (!_FromTargetToSource && !_FromSourceToTarget) return;
  IRTK_START_TIMING();

  // Precompute column/row indices
  _Weight.Index();

  // Calculate source clusters
  if (_FromTargetToSource) {
    UpdateClusters::Run(_Source, _SourceSample, _Weight, false,
                        _SourceClusters, &_SourceOutlier);
    if (_InputSourceClusters.Size() > 0) {
      UpdateClusters::Run(_Source, _SourceSample, _Weight, false,
                          _InputSourceClusters, NULL, true);
    }
  }

  // Calculate target clusters
  if (_FromSourceToTarget) {
    UpdateClusters::Run(_Target, _TargetSample, _Weight, true,
                        _TargetClusters, &_TargetOutlier);
    if (_InputTargetClusters.Size() > 0) {
      UpdateClusters::Run(_Target, _TargetSample, _Weight, true,
                          _InputTargetClusters, NULL, true);
    }
  }

  IRTK_DEBUG_TIMING(6, "calculation of clusters");
}

// -----------------------------------------------------------------------------
bool irtkFuzzyCorrespondence::GetInputTargetPoint(int i, irtkPoint &p) const
{
  if (_InputTargetClusters.Size() > 0) _InputTargetClusters.GetPoint(i, p);
  else                                 _TargetClusters     .GetPoint(i, p);
  return !_TargetOutlier[i];
}

// -----------------------------------------------------------------------------
bool irtkFuzzyCorrespondence::GetInputSourcePoint(int i, irtkPoint &p) const
{
  if (_InputSourceClusters.Size() > 0) _InputSourceClusters.GetPoint(i, p);
  else                                 _SourceClusters     .GetPoint(i, p);
  return !_SourceOutlier[i];
}

// -----------------------------------------------------------------------------
bool irtkFuzzyCorrespondence::GetTargetPoint(int i, irtkPoint &p) const
{
  _TargetClusters.GetPoint(i, p);
  return !_TargetOutlier[i];
}

// -----------------------------------------------------------------------------
bool irtkFuzzyCorrespondence::GetSourcePoint(int i, irtkPoint &p) const
{
  _SourceClusters.GetPoint(i, p);
  return !_SourceOutlier[i];
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
void irtkFuzzyCorrespondence::WriteDataSets(const char *prefix, const char *suffix, bool all) const
{
#ifdef HAVE_MATLAB
  if (all && debug >= 3) {
    const int sz = 1024;
    char      fname[sz];
    snprintf(fname, sz, "%scorrespondences%s.mat", prefix, suffix);
    _Weight.WriteMAT(fname, "W");
  }
#endif
}
