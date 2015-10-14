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

#include <irtkProbabilisticImageSimilarity.h>

#include <irtkParallel.h>
#include <irtkVoxelFunction.h>


// =============================================================================
// Auxiliary functor
// =============================================================================

namespace irtkProbabilisticImageSimilarityUtils {


// -----------------------------------------------------------------------------
/// Add rescaled samples to joint histogram
class FillHistogram_v1
{
  const irtkProbabilisticImageSimilarity               *_Similarity;
  irtkProbabilisticImageSimilarity::JointHistogramType *_Histogram;
  irtkProbabilisticImageSimilarity::JointHistogramType *_Output;

public:

  FillHistogram_v1(const irtkProbabilisticImageSimilarity               *sim,
                   irtkProbabilisticImageSimilarity::JointHistogramType *hist)
  :
    _Similarity(sim), _Histogram(hist), _Output(hist)
  {}

  FillHistogram_v1(const FillHistogram_v1 &lhs, split)
  :
    _Similarity(lhs._Similarity), _Histogram(NULL), _Output(lhs._Output)
  {
    _Histogram = new irtkProbabilisticImageSimilarity::JointHistogramType(_Output->NumberOfBinsX(),
                                                                          _Output->NumberOfBinsY());
  }

  ~FillHistogram_v1()
  {
    if (_Histogram != _Output) Delete(_Histogram);
  }

  void join(const FillHistogram_v1 &rhs)
  {
    const int nbins = _Histogram->NumberOfBins();
    irtkProbabilisticImageSimilarity::JointHistogramType::BinType *l = _Histogram->RawPointer();
    irtkProbabilisticImageSimilarity::JointHistogramType::BinType *r = rhs._Histogram->RawPointer();
    for (int i = 0; i < nbins; ++i, ++l, ++r) (*l) += (*r);
    _Histogram->NumberOfSamples(_Histogram->NumberOfSamples() + rhs._Histogram->NumberOfSamples());
  }

  void operator ()(const blocked_range<int> &re)
  {
    const irtkRegisteredImage::VoxelType *tgt = _Similarity->Target()->Data(re.begin());
    const irtkRegisteredImage::VoxelType *src = _Similarity->Source()->Data(re.begin());
    for (int idx = re.begin(); idx != re.end(); ++idx, ++tgt, ++src) {
      if (_Similarity->IsForeground(idx)) {
        _Histogram->Add(*tgt, round(*src));
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Add samples to joint histogram (no rescaling required)
class FillHistogram_v2
{
  const irtkProbabilisticImageSimilarity               *_Similarity;
  irtkProbabilisticImageSimilarity::JointHistogramType *_Histogram;
  irtkProbabilisticImageSimilarity::JointHistogramType *_Output;

public:

  FillHistogram_v2(const irtkProbabilisticImageSimilarity               *sim,
                   irtkProbabilisticImageSimilarity::JointHistogramType *hist)
  :
    _Similarity(sim), _Histogram(hist), _Output(hist)
  {}

  FillHistogram_v2(const FillHistogram_v2 &lhs, split)
  :
    _Similarity(lhs._Similarity), _Histogram(NULL), _Output(lhs._Output)
  {
    double xmin, ymin, xmax, ymax, xwidth, ywidth;
    _Output->GetMin  (&xmin,   &ymin);
    _Output->GetMax  (&xmax,   &ymax);
    _Output->GetWidth(&xwidth, &ywidth);
    _Histogram = new irtkProbabilisticImageSimilarity::JointHistogramType(xmin, xmax, xwidth,
                                                                          ymin, ymax, ywidth);
    if (_Histogram->NumberOfBinsX() != _Output->NumberOfBinsX() ||
        _Histogram->NumberOfBinsY() != _Output->NumberOfBinsY()) {
      _Histogram->PutNumberOfBins(_Output->NumberOfBinsX(), _Output->NumberOfBinsY());
    }
  }

  ~FillHistogram_v2()
  {
    if (_Histogram != _Output) Delete(_Histogram);
  }

  void join(const FillHistogram_v2 &rhs)
  {
    const int nbins = _Histogram->NumberOfBins();
    irtkProbabilisticImageSimilarity::JointHistogramType::BinType *l = _Histogram->RawPointer();
    irtkProbabilisticImageSimilarity::JointHistogramType::BinType *r = rhs._Histogram->RawPointer();
    for (int i = 0; i < nbins; ++i, ++l, ++r) (*l) += (*r);
    _Histogram->NumberOfSamples(_Histogram->NumberOfSamples() + rhs._Histogram->NumberOfSamples());
  }

  void operator ()(const blocked_range<int> &re)
  {
    const irtkRegisteredImage::VoxelType *tgt = _Similarity->Target()->Data(re.begin());
    const irtkRegisteredImage::VoxelType *src = _Similarity->Source()->Data(re.begin());
    for (int idx = re.begin(); idx != re.end(); ++idx, ++tgt, ++src) {
      if (_Similarity->IsForeground(idx)) {
        _Histogram->Add(_Histogram->ValToBinX(*tgt), _Histogram->ValToBinY(*src));
      }
    }
  }
};


} // namespace irtkProbabilisticImageSimilarityUtils
using namespace irtkProbabilisticImageSimilarityUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkProbabilisticImageSimilarity::irtkProbabilisticImageSimilarity(const char *name, double weight)
:
  irtkImageSimilarity(name, weight),
  _Samples           (NULL),
  _Histogram         (NULL),
  _NumberOfTargetBins(0),
  _NumberOfSourceBins(0)
{
}

// -----------------------------------------------------------------------------
irtkProbabilisticImageSimilarity::irtkProbabilisticImageSimilarity(const irtkProbabilisticImageSimilarity &other)
:
  irtkImageSimilarity(other),
  _Samples           (other._Samples   ? new JointHistogramType(*other._Samples)   : NULL),
  _Histogram         (other._Histogram ? new JointHistogramType(*other._Histogram) : NULL),
  _NumberOfTargetBins(other._NumberOfTargetBins),
  _NumberOfSourceBins(other._NumberOfSourceBins)
{
}

// -----------------------------------------------------------------------------
irtkProbabilisticImageSimilarity &irtkProbabilisticImageSimilarity::operator =(const irtkProbabilisticImageSimilarity &other)
{
  irtkImageSimilarity::operator =(other);
  _Samples            = other._Samples   ? new JointHistogramType(*other._Samples)   : NULL;
  _Histogram          = other._Histogram ? new JointHistogramType(*other._Histogram) : NULL;
  _NumberOfTargetBins = other._NumberOfTargetBins;
  _NumberOfSourceBins = other._NumberOfSourceBins;
  return *this;
}

// -----------------------------------------------------------------------------
irtkProbabilisticImageSimilarity::~irtkProbabilisticImageSimilarity()
{
  Delete(_Samples);
  Delete(_Histogram);
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkProbabilisticImageSimilarity::Set(const char *param, const char *value)
{
  if (strcmp(param, "No. of bins") == 0) {
    if (!FromString(value, _NumberOfTargetBins) && _NumberOfTargetBins < 1) return false;
    _NumberOfSourceBins = _NumberOfTargetBins;
    return true;
  }
  if (strcmp(param, "No. of target bins") == 0) {
    return FromString(value, _NumberOfTargetBins) && _NumberOfTargetBins > 0;
  }
  if (strcmp(param, "No. of source bins") == 0) {
    return FromString(value, _NumberOfSourceBins) && _NumberOfSourceBins > 0;
  }
  return irtkImageSimilarity::Set(param, value);
}

// -----------------------------------------------------------------------------
irtkParameterList irtkProbabilisticImageSimilarity::Parameter() const
{
  irtkParameterList params = irtkImageSimilarity::Parameter();
  if (_NumberOfTargetBins == _NumberOfSourceBins) {
    Insert(params, "No. of bins", ToString(_NumberOfTargetBins));
  } else {
    Insert(params, "No. of target bins", ToString(_NumberOfTargetBins));
    Insert(params, "No. of source bins", ToString(_NumberOfSourceBins));
  }
  return params;
}

// =============================================================================
// Initialization/Update
// =============================================================================

// -----------------------------------------------------------------------------
int DefaultNumberOfBins(const irtkImage *image, double min_intensity, double max_intensity)
{
  int nbins = min(round((max_intensity - min_intensity) / 5.0),
                  round(image->NumberOfVoxels() / 1000.0));
  if      (nbins < 16) nbins = 16;
  else if (nbins > 64) nbins = 64;
  return nbins;
}

// -----------------------------------------------------------------------------
void irtkProbabilisticImageSimilarity::Initialize()
{
  double tmin = numeric_limits<double>::quiet_NaN(), tmax;
  double smin = numeric_limits<double>::quiet_NaN(), smax;

  // Initialize base class
  irtkImageSimilarity::Initialize();

  // Set default number of bins
  if (_NumberOfTargetBins <= 0) {
    Target()->InputImage()->GetMinMaxAsDouble(&tmin, &tmax);
    _NumberOfTargetBins = DefaultNumberOfBins(Target()->InputImage(), tmin, tmax);
  }
  if (_NumberOfSourceBins <= 0) {
    Source()->InputImage()->GetMinMaxAsDouble(&smin, &smax);
    _NumberOfSourceBins = DefaultNumberOfBins(Source()->InputImage(), smin, smax);
  }

  // Initialize container for raw joint histogram samples
  if (version >= irtkVersion(2, 2)) {
    Delete(_Samples);
    if (IsNaN(tmin)) Target()->InputImage()->GetMinMaxAsDouble(&tmin, &tmax);
    if (IsNaN(smin)) Source()->InputImage()->GetMinMaxAsDouble(&smin, &smax);
    if (fequal(tmin, tmax)) {
      cerr << this->NameOfClass() << "::Initialize(): Input target image has homogeneous intensity values only" << endl;
      exit(1);
    }
    if (fequal(smin, smax)) {
      cerr << this->NameOfClass() << "::Initialize(): Input source image has homogeneous intensity values only" << endl;
      exit(1);
    }
    double twidth = (tmax - tmin) / _NumberOfTargetBins;
    double swidth = (smax - smin) / _NumberOfSourceBins;
    _Samples = new JointHistogramType(tmin, tmax, twidth,
                                      smin, smax, swidth);
  } else {
    if (_Samples) _Samples->PutNumberOfBins        (_NumberOfTargetBins, _NumberOfSourceBins);
    else          _Samples = new JointHistogramType(_NumberOfTargetBins, _NumberOfSourceBins);
  }

  // Initialize joint histogram
  if (!_Histogram) _Histogram = new JointHistogramType(*_Samples);
}

// -----------------------------------------------------------------------------
void irtkProbabilisticImageSimilarity::Update(bool gradient)
{
  // Update base class and moving image(s)
  irtkImageSimilarity::Update(gradient);

  IRTK_START_TIMING();

  // Reset histogram
  _Samples->Reset();

  // Add histogram samples
  blocked_range<int> voxels(0, _NumberOfVoxels, _NumberOfVoxels / 8);
  if (version >= irtkVersion(2, 2)) {
    FillHistogram_v2 add(this, _Samples);
    parallel_reduce(voxels, add);
  } else {
    FillHistogram_v1 add(this, _Samples);
    parallel_reduce(voxels, add);
  }

  // Smooth histogram
  //
  // Note that the _Samples cannot be smoothed directly because of the
  // Include/Exclude functions needed for the (optional) finite difference
  // approximation of the gradient.
  _Histogram->Reset(*_Samples);
  _Histogram->Smooth();

  IRTK_DEBUG_TIMING(2, "update of joint histogram");
}

// -----------------------------------------------------------------------------
void irtkProbabilisticImageSimilarity::Exclude(const blocked_range3d<int> &region)
{
  if (version >= irtkVersion(2, 2)) {
    for (int k = region.pages().begin(); k < region.pages().end(); ++k)
    for (int j = region.rows ().begin(); j < region.rows ().end(); ++j)
    for (int i = region.cols ().begin(); i < region.cols ().end(); ++i) {
      if (IsForeground(i, j, k)) {
        _Samples->Delete(_Samples->ValToBinX(_Target->Get(i, j, k)),
                         _Samples->ValToBinY(_Source->Get(i, j, k)));
      }
    }
  } else {
    for (int k = region.pages().begin(); k < region.pages().end(); ++k)
    for (int j = region.rows ().begin(); j < region.rows ().end(); ++j)
    for (int i = region.cols ().begin(); i < region.cols ().end(); ++i) {
      if (IsForeground(i, j, k)) {
        _Samples->DelSample(_Target->Get(i, j, k), round(_Source->Get(i, j, k)));
      }
    }
  }
}

// -----------------------------------------------------------------------------
void irtkProbabilisticImageSimilarity::Include(const blocked_range3d<int> &region)
{
  bool changed = false;
  if (version >= irtkVersion(2, 2)) {
    for (int k = region.pages().begin(); k < region.pages().end(); ++k)
    for (int j = region.rows ().begin(); j < region.rows ().end(); ++j)
    for (int i = region.cols ().begin(); i < region.cols ().end(); ++i) {
      if (IsForeground(i, j, k)) {
        _Samples->Add(_Samples->ValToBinX(_Target->Get(i, j, k)),
                      _Samples->ValToBinY(_Source->Get(i, j, k)));
        changed = true;
      }
    }
  } else {
    for (int k = region.pages().begin(); k < region.pages().end(); ++k)
    for (int j = region.rows ().begin(); j < region.rows ().end(); ++j)
    for (int i = region.cols ().begin(); i < region.cols ().end(); ++i) {
      if (IsForeground(i, j, k)) {
        _Samples->Add(_Target->Get(i, j, k), round(_Source->Get(i, j, k)));
        changed = true;
      }
    }
  }
  if (changed) {
    _Histogram->Reset(*_Samples);
    _Histogram->Smooth();
  }
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
void irtkProbabilisticImageSimilarity::Print(irtkIndent indent) const
{
  irtkImageSimilarity::Print(indent);

  double xmin, xmax, ymin, ymax, xwidth, ywidth;
  _Samples->GetMin  (&xmin,   &ymin);
  _Samples->GetMax  (&xmax,   &ymax);
  _Samples->GetWidth(&xwidth, &ywidth);

  cout << indent << "Intensity range: [" << xmin << ", " << xmax << "] x [" << endl
                                         << ymin << ", " << ymax << "]" << endl;
  cout << indent << "No. of bins:     " << _Samples->NumberOfBinsX() << " x "
                                        << _Samples->NumberOfBinsY() << endl;
  cout << indent << "Bin size:        " << xwidth << " x " << ywidth << endl;
  cout << indent << "No. of samples:  " << _Samples->NumberOfSamples() << endl;
}

// -----------------------------------------------------------------------------
void irtkProbabilisticImageSimilarity::WriteDataSets(const char *p, const char *suffix, bool all) const
{
  irtkImageSimilarity::WriteDataSets(p, suffix, all);

  const int   sz = 1024;
  char        fname[sz];
  string      _prefix = Prefix(p);
  const char *prefix  = _prefix.c_str();

#ifdef HAS_NIFTI
  snprintf(fname, sz, "%sjoint_histogram%s.nii.gz", prefix, suffix);
#else
  snprintf(fname, sz, "%sjoint_histogram%s.gipl", prefix, suffix);
#endif
  if (_Histogram) _Histogram->WriteAsImage(fname);
  else            _Samples  ->WriteAsImage(fname);
}
