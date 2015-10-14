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

#ifndef _IRTKHISTOGRAM_1D_H

#define _IRTKHISTOGRAM_1D_H

/** Class for 2D histograms.
 *
 *  This class defines and implements 2D histograms.
 */

template <class HistogramType> class irtkHistogram_1D : public irtkObject
{
  irtkObjectMacro(irtkHistogram_1D);

protected:

  /// Number of bins
  int _nbins;

  /// Number of samples
  HistogramType _nsamp;

  /// Min. value for samples, everything below is ignored
  double _min;

  /// Max. value for samples, everything below is ignored
  double _max;

  /// Width of bins
  double _width;

  /// Dynamic memory for bins
  HistogramType *_bins;

public:

  /// Construct a histogram from another histogram
  irtkHistogram_1D(const irtkHistogram_1D &);

  /// Construct a histogram with 256 bins and samples ranging from 0 to 255
  irtkHistogram_1D(int nbins = 256);

  /// Construct a histogram for samples ranging from min to max and width
  irtkHistogram_1D(double min, double max, double width);

  /// Read constructor
  irtkHistogram_1D(const char *);

  /// Destructor
  ~irtkHistogram_1D();

  /// Get raw pointer to histogram bins
  HistogramType *RawPointer();

  /// Get raw pointer to histogram bins
  const HistogramType *RawPointer() const;

  /// Clear and reset histogram
  void Reset();

  /// Get number of bins in histogram
  int NumberOfBins() const;

  /// Put number of bins in histogram
  /// \note Unlike PutNumberOfBins, this function does not reset the histogram.
  void NumberOfBins(int);

  /// Put number of bins in histogram
  void PutNumberOfBins(int);

  /// Get minimum value in histogram
  double Min() const;

  /// Set minimum value in histogram
  /// \note Unlike PutMin, this function does not reset the histogram.
  void Min(double);

  /// Get maximum value in histogram
  double Max() const;

  /// Set maximum value in histogram
  /// \note Unlike PutMax, this function does not reset the histogram.
  void Max(double);

  /// Get minimum value in histogram
  double GetMin() const;

  /// Put minimum value in histogram
  void   PutMin(double);

  /// Get maximum value in histogram
  double GetMax() const;

  /// Put maximum value in histogram
  void   PutMax(double);

  /// Get width of bins in histogram
  double GetWidth() const;

  /// Put width of bins in histogram
  void   PutWidth(double);

  /// Get number of samples in histogram
  HistogramType NumberOfSamples() const;

  /// Set number of samples in histogram
  void NumberOfSamples(HistogramType);

  /// Get number of samples in bin(i)
  HistogramType &operator()(int);

  /// Get number of samples in bin(i)
  const HistogramType &operator()(int) const;

  /// Add counts to bin
  void Add(int, HistogramType = 1);

  /// Delete counts from bin
  void Delete(int, HistogramType = 1);

  /// Add sample to bin
  void AddSample(double, HistogramType = 1);

  /// Delete sample from bin
  void DelSample(double, HistogramType = 1);

  /// Convert sample value to continuous bin index
  double ValToRange(double val) const;

  /// Convert sample value to bin index
  int ValToBin(double val) const;

  /// Convert bin index to sample value
  double BinToVal(int bin) const;

  /// Convert bin into probability density distributions
  double BinToPDF(int bin) const;

  /// Convert sample value into probability density distributions
  double ValToPDF(double val) const;

  /// Convert bin into cumulative density distributions
  double BinToCDF(int bin) const;

  /// Convert sample value into cumulative  density distributions
  double ValToCDF(double val) const;

  /// Convert cumulative density distributions to bin value
  double CDFToBin(double p) const;

  /// Convert cumulative density distributions to sample value
  double CDFToVal(double p) const;

  /// Log transform histogram
  void Log();

  /// Smooth histogram
  void Smooth();

  /// Calculate mean
  double Mean() const;

  /// Calculate variance
  double Variance() const;

  /// Calculate standard deviation
  double StandardDeviation() const;

  /// Calculate entropy
  double Entropy() const;

  /// Read histogram
  void Read(const char *);

  /// Wrirte histogram
  void Write(const char *) const;

  /// Print histogram
  void Print() const;
};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
template <class HistogramType>
inline HistogramType *irtkHistogram_1D<HistogramType>::RawPointer()
{
  return _bins;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline const HistogramType *irtkHistogram_1D<HistogramType>::RawPointer() const
{
  return _bins;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline int irtkHistogram_1D<HistogramType>::NumberOfBins() const
{
  return _nbins;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline void irtkHistogram_1D<HistogramType>::NumberOfBins(int nbins)
{
  if (_nbins != nbins) {
    Deallocate(_bins);
    Allocate(_bins, nbins);
    _nbins = nbins;
    _width = (_max - _min) / _nbins;
    _nsamp = .0;
  }
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double irtkHistogram_1D<HistogramType>::Min() const
{
  return _min;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline void irtkHistogram_1D<HistogramType>::Min(double min)
{
  _min   = min;
  _width = (_max - _min) / _nbins;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double irtkHistogram_1D<HistogramType>::GetMin() const
{
  return _min;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double irtkHistogram_1D<HistogramType>::Max() const
{
  return _max;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline void irtkHistogram_1D<HistogramType>::Max(double max)
{
  _max   = max;
  _width = (_max - _min) / _nbins;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double irtkHistogram_1D<HistogramType>::GetMax() const
{
  return _max;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double irtkHistogram_1D<HistogramType>::GetWidth() const
{
  return _width;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline HistogramType irtkHistogram_1D<HistogramType>::NumberOfSamples() const
{
  return _nsamp;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline void irtkHistogram_1D<HistogramType>::NumberOfSamples(HistogramType n)
{
  _nsamp = n;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline HistogramType &irtkHistogram_1D<HistogramType>::operator()(int i)
{
#ifndef NO_BOUNDS
  if ((i < 0) || (i >= _nbins)) {
    cerr << "irtkHistogram_1D<HistogramType>::operator(): No such bin" << endl;
    exit(1);
  }
#endif
  return _bins[i];
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline const HistogramType &irtkHistogram_1D<HistogramType>::operator()(int i) const
{
#ifndef NO_BOUNDS
  if ((i < 0) || (i >= _nbins)) {
    cerr << "irtkHistogram_1D<HistogramType>::operator(): No such bin" << endl;
    exit(1);
  }
#endif
  return _bins[i];
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline void irtkHistogram_1D<HistogramType>::Add(int i, HistogramType n)
{
#ifndef NO_BOUNDS
  if ((i < 0) || (i >= _nbins)) {
    cerr << "irtkHistogram_1D<HistogramType>::Add: No such bin" << endl;
    exit(1);
  }
#endif
  _bins[i] += n;
  _nsamp   += n;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline void irtkHistogram_1D<HistogramType>::Delete(int i, HistogramType n)
{
#ifndef NO_BOUNDS
  if ((i < 0) || (i >= _nbins)) {
    cerr << "irtkHistogram_1D<HistogramType>::Delete: No such bin" << endl;
    exit(1);
  }
#endif
  _bins[i] -= n;
  _nsamp   -= n;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline void irtkHistogram_1D<HistogramType>::AddSample(double x, HistogramType n)
{
  if (x < _min || x > _max) return;
  int index = round(_nbins * (x - _min - 0.5*_width) / (_max - _min));
  if (index <  0     ) index = 0;
  if (index >= _nbins) index = _nbins - 1;
  _bins[index] += n;
  _nsamp       += n;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline void irtkHistogram_1D<HistogramType>::DelSample(double x, HistogramType n)
{
  if (x < _min || x > _max) return;
  int index = round(_nbins * (x - _min - 0.5*_width) / (_max - _min));
  if (index <  0     ) index = 0;
  if (index >= _nbins) index = _nbins - 1;
  _bins[index] -= n;
  _nsamp       -= n;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double irtkHistogram_1D<HistogramType>::ValToRange(double val) const
{
#ifndef NO_BOUNDS
  if ((val < _min) || (val > _max)) {
    cerr << "irtkHistogram_1D<HistogramType>::ValToRange: Must be between " << _min << " and " << _max << endl;
    exit(1);
  }
#endif
  return _nbins * (val - _min - 0.5 * _width) / (_max - _min);
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline int irtkHistogram_1D<HistogramType>::ValToBin(double val) const
{
#ifndef NO_BOUNDS
  if ((val < _min) || (val > _max)) {
    cerr << "irtkHistogram_1D<HistogramType>::ValToBin: Must be between " << _min << " and " << _max << endl;
    exit(1);
  }
#endif
  const int index = round(ValToRange(val));
  return (index < 0 ? 0 : (index >= _nbins ? _nbins - 1 : index));
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double irtkHistogram_1D<HistogramType>::BinToVal(int i) const
{
#ifndef NO_BOUNDS
  if ((i < 0) || (i > _nbins)) {
    cerr << "irtkHistogram_1D<HistogramType>::BinToVal: Must be between 0 and " << _nbins << endl;
    exit(1);
  }
#endif
  return (i*(_max - _min)/_nbins + _min) + _width/2.0;
}

// -----------------------------------------------------------------------------
template <class HistogramType>
inline double irtkHistogram_1D<HistogramType>::BinToPDF(int i) const
{
#ifndef NO_BOUNDS
  if ((i < 0) || (i >= _nbins)) {
    cerr << "irtkHistogram_1D<HistogramType>::BinToPDF: No such bin" << endl;
    exit(1);
  }
#endif
  return ((_nsamp == .0) ? .0 : _bins[i] / _nsamp);
}


#include <irtkImageHistogram_1D.h>

#endif
