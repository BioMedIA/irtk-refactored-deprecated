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
// irtkHistogram_1D
// =============================================================================

// -----------------------------------------------------------------------------
template <class HistogramType>
irtkHistogram_1D<HistogramType>::irtkHistogram_1D(const irtkHistogram_1D &h)
:
  irtkObject(h), _bins(NULL)
{
  _min   = h._min;
  _max   = h._max;
  _width = h._width;
  _nbins = h._nbins;
  _nsamp = h._nsamp;
  if (_nbins > 0) {
    Allocate(_bins, _nbins);
    memcpy(RawPointer(), h.RawPointer(), _nbins * sizeof(HistogramType));
  }
}

// -----------------------------------------------------------------------------
template <class HistogramType>
irtkHistogram_1D<HistogramType>::irtkHistogram_1D(int nbins)
{
  _min   = 0;
  _max   = nbins;
  _width = 1;
  _nbins = nbins;
  _nsamp = 0;
  CAllocate(_bins, _nbins);
}

// -----------------------------------------------------------------------------
template <class HistogramType>
irtkHistogram_1D<HistogramType>::irtkHistogram_1D(double min, double max, double width)
{
  _min   = min;
  _max   = max;
  _nbins = int((_max - _min) / width);
  _width = (_max - _min) / (double)_nbins;
  _nsamp = 0;
  if (_nbins < 1) {
    cerr << "irtkHistogram_1D<HistogramType>::irtkHistogram_1D: Should have at least one bin << endl";
    exit(1);
  }
  CAllocate(_bins, _nbins);
}

// -----------------------------------------------------------------------------
template <class HistogramType>
irtkHistogram_1D<HistogramType>::irtkHistogram_1D(const char *filename)
{
  ifstream from(filename);
  if (!from) {
    cerr << "irtkHistogram_1D<HistogramType>::Read: Can't open file " << filename << endl;
    exit(1);
  }
  char buffer[255];
  from >> buffer;
  if (strcmp(buffer, "irtkHistogram_1D") != 0) {
    cerr << "irtkHistogram_1D<HistogramType>::Read: Invalid format" << endl;
    exit(1);
  }
  from >> _nbins >> _nsamp >> _min >> _max >> _width;
  Allocate(_bins, _nbins);
  for (int i = 0; i < _nbins; ++i) from >> _bins[i];
}

// -----------------------------------------------------------------------------
template <class HistogramType>
irtkHistogram_1D<HistogramType>::~irtkHistogram_1D()
{
  Deallocate(_bins);
  _nbins = 0;
  _nsamp = 0;
  _min   = 0;
  _max   = 0;
  _width = 0;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_1D<HistogramType>::Reset()
{
  memset(_bins, 0, _nbins * sizeof(HistogramType));
  _nsamp = 0;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_1D<HistogramType>::PutMin(double min)
{
  Min(min);
  this->Reset();
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_1D<HistogramType>::PutMax(double max)
{
  Max(max);
  this->Reset();
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_1D<HistogramType>::PutWidth(double width)
{
  const int nbins = round((_max - _min) / width);
  if (nbins < 1) {
    cerr << "irtkHistogram_1D<HistogramType>::PutWidth: Should have at least one bin";
    exit(1);
  }
  if (_nbins != nbins) {
    Deallocate(_bins);
    Allocate(_bins, nbins);
    _nbins = nbins;
  }
  _width = (_max - _min) / _nbins;
  this->Reset();
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_1D<HistogramType>::PutNumberOfBins(int nbins)
{
  if (nbins < 1) {
    cerr << "irtkHistogram_1D<HistogramType>::PutNumberOfBins: Should have at least one bin";
    exit(1);
  }
  NumberOfBins(nbins);
  this->Reset();
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double irtkHistogram_1D<HistogramType>::CDFToVal(double p) const
{
  if ((p < 0) || (p > 1)) {
    cerr << "irtkHistogram_1D<HistogramType>::CDFToVal: Must be between 0 and 1" << endl;
    exit(1);
  }
  int i;
  HistogramType sum = 0;
  for (i = 0; i < _nbins; ++i) {
    sum += _bins[i];
    if (sum / _nsamp >= p) {
      break;
    }
  }
  return BinToVal(i);
}

// -----------------------------------------------------------------------------
template <>
void irtkHistogram_1D<double>::Log()
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "irtkHistogram_1D<HistogramType>::Log: No samples in irtkHistogram" << endl;
    }
    return;
  }
  for (int i = 0; i < _nbins; ++i) {
    if (_bins[i] > 0) {
      _bins[i] = log(static_cast<double>(_bins[i]/(double)_nsamp));
    } else {
      _bins[i] = 0;
    }
  }
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double irtkHistogram_1D<HistogramType>::Mean() const
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "irtkHistogram_1D<HistogramType>::Mean: No samples in irtkHistogram" << endl;
    }
    return 0;
  }
  double val = 0;
  for (int i = 0; i < _nbins; ++i) {
    val += _bins[i] * this->BinToVal(i);
  }
  return val / _nsamp;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double irtkHistogram_1D<HistogramType>::Variance() const
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "irtkHistogram_1D<HistogramType>::Variance: No samples in irtkHistogram" << endl;
    }
    return 0;
  }
  double val = 0, mean = this->Mean();
  for (int i = 0; i < _nbins; ++i) {
    val += _bins[i] * (this->BinToVal(i) - mean) * (this->BinToVal(i) - mean);
  }
  return val / _nsamp;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double irtkHistogram_1D<HistogramType>::StandardDeviation() const
{
  return sqrt(this->Variance());
}

// -----------------------------------------------------------------------------
template <class HistogramType>
double irtkHistogram_1D<HistogramType>::Entropy() const
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "irtkHistogram_1D<HistogramType>::Entropy: No samples in irtkHistogram" << endl;
    }
    return 0;
  }
  // Attention: Parallel summation yielded slightly different results each
  //            time this function was executed. This might be caused by a
  //            different summation of values, which causes different numerical
  //            cancelations. Qhen used for NMI gradient computation, the
  //            registration result could differ from run to run!
  double p, sum = .0;
  const HistogramType *bin = _bins;
  for (int i = 0; i != _nbins; ++i, ++bin) {
    p = static_cast<double>(*bin);
    if (p > .0) sum += p * log(p);
  }
  // H = - sum (p/n) log(p/n) = log(n) - sum(p log p) / n
  return log(_nsamp) - sum / _nsamp;
}

// -----------------------------------------------------------------------------
template <>
void irtkHistogram_1D<double>::Smooth()
{
  if (_nsamp == 0) {
    if (debug) {
      cerr << "irtkHistogram_1D<HistogramType>::Smooth: No samples in Histogram" << endl;
    }
    return;
  }
  if (_nbins == 0) return;

  // Smoothing kernel
  const double kernel[3] = { 1.0/6.0, 2.0/3.0, 1.0/6.0 };

  // Allocate new histogram
  double *bins = Allocate<double>(_nbins);

  // Smooth histogram
  bins[0]                  = kernel[1] * _bins[0];
  if (_nbins > 1) bins[0] += kernel[0] * _bins[1];
  for (int i = 1; i < _nbins - 1; ++i) {
    bins[i] = kernel[2] * _bins[i - 1]
            + kernel[1] * _bins[i    ]
            + kernel[0] * _bins[i + 1];
  }
  if (_nbins > 1) {
    bins[_nbins-1] = kernel[1] * _bins[_nbins-1]
                   + kernel[2] * _bins[_nbins-2];
  }

  // Replace histogram by smoothed version
  Deallocate(_bins);
  _bins = bins;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_1D<HistogramType>::Read(const char *filename)
{
  ifstream from(filename);
  if (!from) {
    cerr << "irtkHistogram_1D<HistogramType>::Read: Can't open file " << filename << "\n";
    exit(1);
  }

  char buffer[255];
  from >> buffer;
  if (strcmp(buffer, "irtkHistogram_1D") != 0) {
    cerr << "irtkHistogram_1D<HistogramType>::Read: Invalid format" << endl;
    exit(1);
  }

  int nbins;
  from >> nbins >> _nsamp >> _min >> _max >> _width;
  if (_nbins != nbins) {
    Deallocate(_bins);
    Allocate(_bins, nbins);
    _nbins = nbins;
  }
  for (int i = 0; i < _nbins; ++i) from >> _bins[i];
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_1D<HistogramType>::Write(const char *filename) const
{
  ofstream to(filename);
  if (!to) {
    cerr << "irtkHistogram_1D<HistogramType>::Write: Can't open file " << filename << "\n";
    exit(1);
  }
  to << "irtkHistogram_1D\n";
  to << _nbins << " " << _nsamp << " " << _min << " " << _max << " " << _width << endl;
  for (int i = 0; i < _nbins; ++i) to << _bins[i] << endl;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
void irtkHistogram_1D<HistogramType>::Print() const
{
  cout << _nbins << " " << _nsamp << " " << _min << " " << _max << " " << _width << endl;
  for (int i = 0; i < _nbins; ++i) cout << _bins[i] << endl;
}

// =============================================================================
// Explicit instantiations
// =============================================================================

template class irtkHistogram_1D<int>;
template class irtkHistogram_1D<double>;
template class irtkHistogram_1D<unsigned char>;
template class irtkHistogram_1D<short>;
template class irtkHistogram_1D<unsigned short>;
template class irtkHistogram_1D<float>;
