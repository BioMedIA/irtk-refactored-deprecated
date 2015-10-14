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

#include <irtkImage.h>

#include <irtkHistogram.h>

// =============================================================================
// Auxiliary functors
// =============================================================================

namespace irtkHistogram_2DUtils {


// -----------------------------------------------------------------------------
template <class T>
class MarginalizeX
{
  const T *_JointProbability;
  T       *_Probability;
  int      _Nx, _Ny;

public:

  void operator ()(const blocked_range<int> &re) const
  {
    T sum;
    const T *bin;
    for (int i = re.begin(); i != re.end(); ++i) {
      sum = T(0);
      bin = _JointProbability + i;
      for (int j = 0; j < _Ny; ++j, bin += _Nx) sum += (*bin);
      _Probability[i] = sum;
    }
  }

  static void Run(const irtkHistogram_2D<T> *hxy, irtkHistogram_1D<T> &hx)
  {
    hx.PutNumberOfBins(hxy->NumberOfBinsX());
    hx.Min(hxy->MinX());
    hx.Max(hxy->MaxX());
    hx.NumberOfSamples(hxy->NumberOfSamples());
    MarginalizeX body;
    body._JointProbability = hxy->RawPointer();
    body._Probability      = hx  .RawPointer();
    body._Nx               = hxy->NumberOfBinsX();
    body._Ny               = hxy->NumberOfBinsY();
    blocked_range<int> i(0, hx.NumberOfBins());
    parallel_for(i, body);
  }
};

// -----------------------------------------------------------------------------
template <class T>
class MarginalizeY
{
  const T *_JointProbability;
  T       *_Probability;
  int      _Nx, _Ny;

public:

  void operator ()(const blocked_range<int> &re) const
  {
    T sum;
    const T *bin = _JointProbability + re.begin() * _Nx;
    for (int j = re.begin(); j != re.end(); ++j) {
      sum = T(0);
      for (int i = 0; i < _Nx; ++i, ++bin) sum += (*bin);
      _Probability[j] = sum;
    }
  }

  static void Run(const irtkHistogram_2D<T> *hxy, irtkHistogram_1D<T> &hy)
  {
    hy.PutNumberOfBins(hxy->NumberOfBinsY());
    hy.Min(hxy->MinY());
    hy.Max(hxy->MaxY());
    hy.NumberOfSamples(hxy->NumberOfSamples());
    MarginalizeY body;
    body._JointProbability = hxy->RawPointer();
    body._Probability      = hy  .RawPointer();
    body._Nx               = hxy->NumberOfBinsX();
    body._Ny               = hxy->NumberOfBinsY();
    blocked_range<int> i(0, hy.NumberOfBins());
    parallel_for(i, body);
  }
};

// -----------------------------------------------------------------------------
template <class T>
class LogTransform
{
  T     *_Bins;
  double _Num;

public:

  void operator ()(const blocked_range<int> &re) const
  {
    double p;
    T *ptr = _Bins + re.begin();
    for (int i = re.begin(); i != re.end(); ++i, ++ptr) {
      p = static_cast<double>(*ptr);
      if (p > .0) p = log(p / _Num);
      else        p = .0;
      (*ptr) = static_cast<T>(p);
    }
  }

  static void Run(irtkHistogram_2D<T> *hxy)
  {
    LogTransform body;
    body._Bins = hxy->RawPointer();
    body._Num  = hxy->NumberOfSamples();
    blocked_range<int> i(0, hxy->NumberOfBins());
    parallel_for(i, body);
  }
};


} // namespace irtkHistogram_2DUtils
using namespace irtkHistogram_2DUtils;

// =============================================================================
// irtkHistogram_2D
// =============================================================================

// -----------------------------------------------------------------------------
template <class HistogramType>
irtkHistogram_2D<HistogramType>::irtkHistogram_2D(const irtkHistogram_2D &h)
:
  irtkObject(h)
{
  _min_x   = h._min_x;
  _min_y   = h._min_y;
  _max_x   = h._max_x;
  _max_y   = h._max_y;
  _width_x = h._width_x;
  _width_y = h._width_y;
  _nbins_x = h._nbins_x;
  _nbins_y = h._nbins_y;
  _nsamp   = h._nsamp;
  if ((_nbins_x < 1) || (_nbins_y < 1)) {
    cerr << "irtkHistogram_2D<HistogramType>::irtkHistogram_2D: Should have at least one bin" << endl;
    exit(1);
  }
  Allocate(_bins, _nbins_x, _nbins_y);
  memcpy(RawPointer(), h.RawPointer(), NumberOfBins() * sizeof(HistogramType));
}

// -----------------------------------------------------------------------------
template <class HistogramType>
irtkHistogram_2D<HistogramType>::irtkHistogram_2D(int nbins_x, int nbins_y)
{
  if ((nbins_x < 1) || (nbins_y < 1)) {
    cerr << "irtkHistogram_2D<HistogramType>::irtkHistogram_2D: Should have at least one bin" << endl;
    exit(1);
  }
  _min_x   = 0;
  _min_y   = 0;
  _max_x   = nbins_x;
  _max_y   = nbins_y;
  _width_x = 1;
  _width_y = 1;
  _nbins_x = nbins_x;
  _nbins_y = nbins_y;
  _nsamp   = 0;
  CAllocate(_bins, _nbins_x, _nbins_y);
}

// -----------------------------------------------------------------------------
template <class HistogramType>
irtkHistogram_2D<HistogramType>::irtkHistogram_2D(double min_x, double max_x, double width_x,
                                                  double min_y, double max_y, double width_y)
{
  _min_x   = min_x;
  _min_y   = min_y;
  _max_x   = max_x;
  _max_y   = max_y;
  _nbins_x = round((max_x - min_x) / width_x);
  _nbins_y = round((max_y - min_y) / width_y);
  _width_x = (_max_x - _min_x) / (double)_nbins_x;
  _width_y = (_max_y - _min_y) / (double)_nbins_y;
  _nsamp = 0;
  if ((_nbins_x < 1) || (_nbins_y < 1)) {
    cerr << "irtkHistogram_2D<HistogramType>::irtkHistogram_2D: Should have at least one bin" << endl;
    exit(1);
  }
  CAllocate(_bins, _nbins_x, _nbins_y);
}

// -----------------------------------------------------------------------------
template <class HistogramType>
irtkHistogram_2D<HistogramType>::~irtkHistogram_2D()
{
  Deallocate(_bins);
  _nbins_x = 0;
  _nbins_y = 0;
  _min_x   = 0;
  _min_y   = 0;
  _max_x   = 0;
  _max_y   = 0;
  _width_x = 0;
  _width_y = 0;
  _nsamp   = 0;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_2D<HistogramType>::Reset()
{
  memset(RawPointer(), 0, NumberOfBins() * sizeof(HistogramType));
  _nsamp = 0;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_2D<HistogramType>::Reset(const irtkHistogram_2D &h)
{
  if ((_nbins_x != h._nbins_x) || (_nbins_y != h._nbins_y)) {
    Deallocate(_bins);
    Allocate(_bins, h._nbins_x, h._nbins_y);
  }
  _min_x   = h._min_x;
  _min_y   = h._min_y;
  _max_x   = h._max_x;
  _max_y   = h._max_y;
  _width_x = h._width_x;
  _width_y = h._width_y;
  _nbins_x = h._nbins_x;
  _nbins_y = h._nbins_y;
  _nsamp   = h._nsamp;
  memcpy(RawPointer(), h.RawPointer(), NumberOfBins() * sizeof(HistogramType));
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_2D<HistogramType>::PutMin(double min_x, double min_y)
{
  _min_x = min_x;
  _min_y = min_y;
  _width_x = (_max_x - _min_x) / (double)_nbins_x;
  _width_y = (_max_y - _min_y) / (double)_nbins_y;
  this->Reset();
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_2D<HistogramType>::GetMin(double *min_x, double *min_y) const
{
  *min_x = _min_x;
  *min_y = _min_y;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_2D<HistogramType>::GetMin(double &min_x, double &min_y) const
{
  min_x = _min_x;
  min_y = _min_y;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_2D<HistogramType>::PutMax(double max_x, double max_y)
{
  _max_x = max_x;
  _max_y = max_y;
  _width_x = (_max_x - _min_x) / (double)_nbins_x;
  _width_y = (_max_y - _min_y) / (double)_nbins_y;
  this->Reset();
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_2D<HistogramType>::GetMax(double *max_x, double *max_y) const
{
  *max_x = _max_x;
  *max_y = _max_y;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_2D<HistogramType>::GetMax(double &max_x, double &max_y) const
{
  max_x = _max_x;
  max_y = _max_y;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_2D<HistogramType>::PutWidth(double width_x, double width_y)
{
  if ((_nbins_x > 0) && (_nbins_y > 0)) {
    Deallocate(_bins);
  }
  const int nbins_x = round((_max_x - _min_x) / width_x);
  const int nbins_y = round((_max_y - _min_y) / width_y);
  if (_nbins_x != nbins_x || _nbins_y != nbins_y) {
    Deallocate(_bins);
    Allocate(_bins, _nbins_x, _nbins_y);
  }
  _nbins_x = nbins_x;
  _nbins_y = nbins_y;
  _width_x = (_max_x - _min_x) / (double)_nbins_x;
  _width_y = (_max_y - _min_y) / (double)_nbins_y;
  if ((_nbins_x < 1) || (_nbins_y < 1)) {
    cerr << "irtkHistogram_2D<HistogramType>::PutWidth: Should have at least one bin" << endl;
    exit(1);
  }
  this->Reset();
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_2D<HistogramType>::GetWidth(double *width_x, double *width_y) const
{
  *width_x = _width_x;
  *width_y = _width_y;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_2D<HistogramType>::GetWidth(double &width_x, double &width_y) const
{
  width_x = _width_x;
  width_y = _width_y;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_2D<HistogramType>::PutNumberOfBins(int nbins_x, int nbins_y)
{
  if (_nbins_x != nbins_x || _nbins_y != nbins_y) {
    Deallocate(_bins);
    Allocate(_bins, _nbins_x, _nbins_y);
  }
  _nbins_x = nbins_x;
  _nbins_y = nbins_y;
  _width_x = (_max_x - _min_x) / (double)_nbins_x;
  _width_y = (_max_y - _min_y) / (double)_nbins_y;
  if ((_nbins_x < 1) || (_nbins_y < 1)) {
    cerr << "irtkHistogram_2D<HistogramType>::PutWidth: Should have at least one bin" << endl;
    exit(1);
  }
  this->Reset();
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_2D<HistogramType>::GetNumberOfBins(int *nbins_x, int *nbins_y) const
{
  *nbins_x = _nbins_x;
  *nbins_y = _nbins_y;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_2D<HistogramType>::GetNumberOfBins(int &nbins_x, int &nbins_y) const
{
  nbins_x = _nbins_x;
  nbins_y = _nbins_y;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_2D<HistogramType>::AddSample(double x, double y, HistogramType n)
{
  if (x < _min_x || x > _max_x || y < _min_y || y > _max_y) return;
  int i = round(_nbins_x * (x - _min_x - 0.5*_width_x) / (_max_x - _min_x));
  int j = round(_nbins_y * (y - _min_y - 0.5*_width_y) / (_max_y - _min_y));
  if (i < 0) i = 0;
  if (j < 0) j = 0;
  if (i >= _nbins_x) i = _nbins_x - 1;
  if (j >= _nbins_y) j = _nbins_y - 1;
  _bins[j][i] += n;
  _nsamp      += n;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_2D<HistogramType>::DelSample(double x, double y, HistogramType n)
{
  if (x < _min_x || x > _max_x || y < _min_y || y > _max_y) return;
  int i = round(_nbins_x * (x - _min_x - 0.5*_width_x) / (_max_x - _min_x));
  int j = round(_nbins_y * (y - _min_y - 0.5*_width_y) / (_max_y - _min_y));
  if (i < 0) i = 0;
  if (j < 0) j = 0;
  if (i >= _nbins_x) i = _nbins_x - 1;
  if (j >= _nbins_y) j = _nbins_y - 1;
  _bins[j][i] -= n;
  _nsamp      -= n;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_2D<HistogramType>::HistogramX(irtkHistogram_1D<HistogramType> &hx) const
{
  MarginalizeX<HistogramType>::Run(this, hx);
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_2D<HistogramType>::HistogramY(irtkHistogram_1D<HistogramType> &hy) const
{
  MarginalizeY<HistogramType>::Run(this, hy);
}

// -----------------------------------------------------------------------------
template <>
void irtkHistogram_2D<double>::Log()
{
  LogTransform<double>::Run(this);
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double irtkHistogram_2D<HistogramType>::MeanX() const
{
  int i, j;
  double val, tmp;

  if (_nsamp == 0) {
    if (debug) {
      cerr << "irtkHistogram_2D<HistogramType>::MeanX: No samples in Histogram" << endl;
    }
    return 0;
  }
  val = 0;
  for (i = 0; i < _nbins_x; i++) {
    tmp = this->BinToValX(i);
    for (j = 0; j < _nbins_y; j++) {
      val += _bins[j][i] * tmp;
    }
  }
  return val / (double)_nsamp;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double irtkHistogram_2D<HistogramType>::MeanY() const
{
  int i, j;
  double val, tmp;

  if (_nsamp == 0) {
    if (debug) {
      cerr << "Histogram_2D::MeanY: No samples in Histogram" << endl;
    }
    return 0;
  }
  val = 0;
  for (j = 0; j < _nbins_y; j++) {
    tmp = this->BinToValY(j);
    for (i = 0; i < _nbins_x; i++) {
      val += _bins[j][i] * tmp;
    }
  }
  return val / (double)_nsamp;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double irtkHistogram_2D<HistogramType>::VarianceX() const
{
  int i;
  double val;

  if (_nsamp == 0) {
    if (debug) {
      cerr << "irtkHistogram_2D<HistogramType>::VarianceX: No samples in Histogram" << endl;
    }
    return 0;
  }
  val  = 0;
  for (i = 0; i < _nbins_x; i++) {
    val += this->MarginalProbabilityX(i) * pow(this->BinToValX(i), 2.0);
  }
  return val - pow(this->MeanX(), 2.0);
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double irtkHistogram_2D<HistogramType>::VarianceY() const
{
  int i;
  double val;

  if (_nsamp == 0) {
    if (debug) {
      cerr << "irtkHistogram_2D<HistogramType>::VarianceY: No samples in Histogram" << endl;
    }
    return 0;
  }
  val  = 0;
  for (i = 0; i < _nbins_y; i++) {
    val += this->MarginalProbabilityY(i) * pow(this->BinToValY(i), 2.0);
  }
  return val - pow(this->MeanY(), 2.0);
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double irtkHistogram_2D<HistogramType>::StandardDeviationX() const
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "irtkHistogram_2D<HistogramType>::StandardDeviationX: No samples in Histogram" << endl;
    }
    return 0;
  }
  return sqrt(this->VarianceX());
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double irtkHistogram_2D<HistogramType>::StandardDeviationY() const
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "irtkHistogram_2D<HistogramType>::StandardDeviationY: No samples in Histogram" << endl;
    }
    return 0;
  }
  return sqrt(this->VarianceY());
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double irtkHistogram_2D<HistogramType>::Covariance() const
{
  int i, j;
  double val, mean_x, mean_y;

  if (_nsamp == 0) {
    if (debug) {
      cerr << "irtkHistogram_2D<HistogramType>::Covariance: No samples in Histogram" << endl;
    }
    return 0;
  }
  val  = 0;
  mean_x = this->MeanX();
  mean_y = this->MeanY();
  for (j = 0; j < _nbins_y; j++) {
    for (i = 0; i < _nbins_x; i++) {
      val += _bins[j][i] *
             (this->BinToValX(i) - mean_x) *
             (this->BinToValY(j) - mean_y);
    }
  }
  return val / (double)_nsamp;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double irtkHistogram_2D<HistogramType>::EntropyX() const
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "irtkHistogram_2D<HistogramType>::EntropyX: No samples in Histogram" << endl;
    }
    return 0;
  }
  return HistogramX().Entropy();
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double irtkHistogram_2D<HistogramType>::EntropyY() const
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "irtkHistogram_2D<HistogramType>::EntropyY: No samples in Histogram" << endl;
    }
    return 0;
  }
  return HistogramY().Entropy();
}

template <class HistogramType>
double irtkHistogram_2D<HistogramType>::JointEntropy() const
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "irtkHistogram_2D<HistogramType>::JointEntropy: No samples in Histogram" << endl;
    }
    return 0;
  }
  // Attention: Parallel summation yielded slightly different results each
  //            time this function was executed. This might be caused by a
  //            different summation of values, which causes different numerical
  //            cancelations. Qhen used for NMI gradient computation, the
  //            registration result could differ from run to run!
  double p, sum = .0;
  const HistogramType *bin = _bins[0];
  const int nbins = _nbins_x * _nbins_y;
  for (int i = 0; i != nbins; ++i, ++bin) {
    p = static_cast<double>(*bin);
    if (p > .0) sum += p * log(p);
  }
  // H = - sum (p/n) log(p/n) = log(n) - sum(p log p) / n
  return log(_nsamp) - sum / _nsamp;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double irtkHistogram_2D<HistogramType>::ConditionalMeanXY(int i) const
{
  int j;
  double p, m;

  if ((i < 0) || (i > _nbins_y-1)) {
    cerr << "irtkHistogram_2D<HistogramType>::ConditionalMeanXY: No such bin " << i << endl;
    return 0;
  }
  m = 0;
  p = 0;
  for (j = 0; j < _nbins_x; j++) {
    m += this->JointProbability(j, i) * this->BinToValX(j);
    p += this->JointProbability(j, i);
  }
  if (p > 0) {
    return m / p;
  } else {
    return 0;
  }
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double irtkHistogram_2D<HistogramType>::ConditionalMeanYX(int i) const
{
  int j;
  double p, m;

  if ((i < 0) || (i > _nbins_x-1)) {
    cerr << "irtkHistogram_2D<HistogramType>::ConditionalMeanYX: No such bin " << i << endl;
    return 0;
  }
  m = 0;
  p = 0;
  for (j = 0; j < _nbins_y; j++) {
    m += this->JointProbability(i, j) * this->BinToValY(j);
    p += this->JointProbability(i, j);
  }
  if (p > 0) {
    return m / p;
  } else {
    return 0;
  }
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double irtkHistogram_2D<HistogramType>::CorrelationRatioXY() const
{
  int i;
  double c, m;

  if (_nsamp == 0) {
    if (debug) {
      cerr << "irtkHistogram_2D<HistogramType>::CorrelationRatioXY: No samples in Histogram" << endl;
    }
    return 0;
  }

  c = 0;
  m = this->MeanX();
  for (i = 0; i < _nbins_y; i++) {
    c += this->MarginalProbabilityY(i) * pow(this->ConditionalMeanXY(i)-m, 2.0);
  }
  return (c / this->VarianceX());
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double irtkHistogram_2D<HistogramType>::CorrelationRatioYX() const
{
  int i;
  double c, m;

  if (_nsamp == 0) {
    if (debug) {
      cerr << "irtkHistogram_2D<HistogramType>::CorrelationRatioYX: No samples in Histogram" << endl;
    }
    return 0;
  }

  c = 0;
  m = this->MeanY();
  for (i = 0; i < _nbins_x; i++) {
    c += this->MarginalProbabilityX(i) * pow(this->ConditionalMeanYX(i)-m, 2.0);
  }
  return (c / this->VarianceY());
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double irtkHistogram_2D<HistogramType>::MutualInformation() const
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "irtkHistogram_2D<HistogramType>::MutualInformation: No samples in Histogram" << endl;
    }
    return 0;
  }
  return this->EntropyX() + this->EntropyY() - this->JointEntropy();
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double irtkHistogram_2D<HistogramType>::NormalizedMutualInformation() const
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "irtkHistogram_2D<HistogramType>::NormalizedMutualInformation: No samples in Histogram" << endl;
    }
    return 0;
  }
  return (this->EntropyX() + this->EntropyY()) / this->JointEntropy();
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double irtkHistogram_2D<HistogramType>::CrossCorrelation() const
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "irtkHistogram_2D<HistogramType>::CrossCorrelation: No samples in Histogram" << endl;
    }
    return 0;
  }
  return fabs(this->Covariance() / (sqrt(this->VarianceX()) *
                                    sqrt(this->VarianceY())));
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double irtkHistogram_2D<HistogramType>::SumsOfSquaredDifferences() const
{
  int i, j;
  double val_x, val_y, ssd;

  if (_nsamp == 0) {
    if (debug) {
      cerr << "irtkHistogram_2D<HistogramType>::SumsOfSquaredDifferences: ";
      cerr << "No samples in Histogram" << endl;
    }
    return 0;
  }

  ssd   = 0;
  val_x = this->BinToValX(0);
  val_y = this->BinToValY(0);
  for (j = 0; j < _nbins_y; j++) {
    val_x = this->BinToValX(0);
    for (i = 0; i < _nbins_x; i++) {
      ssd   += _bins[j][i] * (val_x - val_y) * (val_x - val_y);
      val_x += _width_x;
    }
    val_y += _width_y;
  }
  return ssd;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double irtkHistogram_2D<HistogramType>::LabelConsistency() const
{
  int i;
  HistogramType n;

  if (_nsamp == 0) {
    if (debug) {
      cerr << "irtkHistogram_2D<HistogramType>::LabelConsistency: No samples in Histogram" << endl;
    }
    return 0;
  }
  if (_nbins_x != _nbins_y) {
    cerr << "irtkHistogram_2D<HistogramType>::LabelConsistency: Histogram must have equal number of bins in X and Y" << endl;
    return 0;
  }

  n = 0;
  for (i = 0; i < _nbins_x; i++) {
    n += _bins[i][i];
  }
  return n / (double)_nsamp;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double irtkHistogram_2D<HistogramType>::Kappa() const
{
  int i, j;

  if (_nsamp == 0) {
    if (debug) {
      cerr << "irtkHistogram_2D<HistogramType>::Kappa: No samples in Histogram" << endl;
    }
    return 0;
  }
  if (_nbins_x != _nbins_y) {
    cerr << "irtkHistogram_2D<HistogramType>::Kappa: Histogram must have equal number of bins in X and Y" << endl;
    return 0;
  }

  HistogramType *col_sum = new HistogramType[_nbins_x] ;
  HistogramType *row_sum = new HistogramType[_nbins_x] ;
  for (j = 0; j < _nbins_x; j++) {
    col_sum[j] = 0;
    row_sum[j] = 0;
    for (i = 0; i < _nbins_x; i++) {
      col_sum[j] += _bins[i][j];
      row_sum[j] += _bins[j][i];
    }
  }

  double po = 0.0;
  double pe = 0.0;
  for (j = 0; j < _nbins_x; j++) {
    po += _bins[j][j];
    pe += col_sum[j] * (double)row_sum[j];
  }
  po /= _nsamp;
  pe /= _nsamp * (double)_nsamp;

  delete []row_sum;
  delete []col_sum;
  return ((po-pe)/(1-pe));
}

// -----------------------------------------------------------------------------
template <>
void irtkHistogram_2D<double>::Smooth()
{
  int i, j, k;
  double **tmp, value;

  if (_nsamp == 0) {
    if (debug) {
      cerr << "irtkHistogram_2D<HistogramType>::Smooth: No samples in Histogram" << endl;
    }
    return;
  }

  // Smoothing kernel
  double kernel[3] = { 1.0/6.0, 2.0/3.0, 1.0/6.0 };

  // Allocate temporary memory
  Allocate(tmp, _nbins_x, _nbins_y);

  // Smooth along the x-axis
  for (j = 0; j < _nbins_y; j++) {
    for (i = 0; i < _nbins_x; i++) {
      value = 0;
      for (k = 0; k < 3; k++) {
        if ((i-1+k >= 0) && (i-1+k < _nbins_x)) {
          value += kernel[k] * _bins[j][i-1+k];
        }
      }
      tmp[j][i] = value;
    }
  }

  // Smooth along the y-axis
  _nsamp = 0;
  for (i = 0; i < _nbins_x; i++) {
    for (j = 0; j < _nbins_y; j++) {
      value = 0;
      for (k = 0; k < 3; k++) {
        if ((j-1+k >= 0) && (j-1+k < _nbins_y)) {
          value += kernel[k] * tmp[j-1+k][i];
        }
      }
      _bins[j][i] = value;
      _nsamp += value;
    }
  }

  // Free tmp memory
  Deallocate(tmp);
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_2D<HistogramType>::Read(const char *filename)
{
  int i, j;
  char buffer[255];

  ifstream from(filename);
  if (!from) {
    cerr << "irtkHistogram_2D<HistogramType>::Read: Can't open file " << filename << endl;
    exit(1);
  }
  if ((_nbins_x > 0) && (_nbins_y > 0)) {
    Deallocate(_bins);
    _nbins_x = 0;
    _nbins_y = 0;
  }

  from >> buffer;
  if (strcmp(buffer, "irtkHistogram_2D") != 0) {
    cerr << "irtkHistogram_2D<HistogramType>::Read: Invalid format" << endl;
    exit(1);
  }

  // Read no. of bins
  from >> _nbins_x;
  from >> _nbins_y;

  // Read no. of samples
  from >> _nsamp;

  // Read min and max of bins
  from >> _min_x;
  from >> _max_x;
  from >> _min_y;
  from >> _max_y;

  // Read width of bins
  from >> _width_x;
  from >> _width_y;

  Allocate(_bins, _nbins_x, _nbins_y);
  for (j = 0; j < _nbins_y; j++) {
    for (i = 0; i < _nbins_x; i++) {
      from >> _bins[j][i];
    }
  }
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_2D<HistogramType>::Write(const char *filename) const
{
  int i, j;

  ofstream to(filename);
  if (!to) {
    cerr << "irtkHistogram_2D<HistogramType>::Write: Can't open file " << filename << endl;
    exit(1);
  }
  to << "irtkHistogram_2D\n";
  to << _nbins_x << " "
     << _nbins_y << " "
     << _nsamp << " "
     << _min_x << " "
     << _max_x << " "
     << _min_y << " "
     << _max_y << " "
     << _width_x << " "
     << _width_y << endl;
  for (j = 0; j < _nbins_y; j++) {
    for (i = 0; i < _nbins_x; i++) {
      to << _bins[j][i] << " ";
    }
    to << endl;
  }
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_2D<HistogramType>::WriteAsImage(const char *filename) const
{
  int i, j;
  irtkGenericImage<HistogramType> image(_nbins_x, _nbins_y, 1);

  for (j = 0; j < _nbins_y; j++) {
    for (i = 0; i < _nbins_x; i++) {
      image(i, j, 0) = _bins[j][i];
    }
  }
  image.Write(filename);
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_2D<HistogramType>::Print() const
{
  int i, j;

  cout << _nbins_x << " "
  << _nbins_y << " "
  << _nsamp << " "
  << _min_x << " "
  << _max_x << " "
  << _min_y << " "
  << _max_y << " "
  << _width_x << " "
  << _width_y << endl;
  for (j = 0; j < _nbins_y; j++) {
    for (i = 0; i < _nbins_x; i++) {
      cout << _bins[j][i] << " ";
    }
    cout << endl;
  }
}

// =============================================================================
// Explicit instantiations
// =============================================================================

template class irtkHistogram_2D<int>;
template class irtkHistogram_2D<double>;
