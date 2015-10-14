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

#include <irtkNormalizedMutualImageInformation.h>


// =============================================================================
// Auxiliary functors
// =============================================================================

namespace irtkNormalizedMutualImageInformationUtils {


// -----------------------------------------------------------------------------
inline double ValToRange(const irtkHistogram_1D<double> &hist, double val)
{
  if (version >= irtkVersion(2, 2)) {
    if (val < hist.Min()) return .0;
    if (val > hist.Max()) return static_cast<double>(hist.NumberOfBins() - 1);
    return hist.ValToRange(val);
  }
  return val;
}

// -----------------------------------------------------------------------------
class CalculateGradient : public irtkVoxelFunction
{
  const irtkNormalizedMutualImageInformation *_This;
  const irtkHistogram_2D<double>             &_LogJointHistogram;
  const irtkHistogram_1D<double>             &_LogMarginalXHistogram;
  const irtkHistogram_1D<double>             &_LogMarginalYHistogram;
  double                                      _JointEntropy;
  double                                      _NormalizedMutualInformation;

public:

  CalculateGradient(const irtkNormalizedMutualImageInformation *_this,
                    const irtkHistogram_2D<double> &logJointHistogram,
                    const irtkHistogram_1D<double> &logMarginalXHistogram,
                    const irtkHistogram_1D<double> &logMarginalYHistogram,
                    double je, double nmi)
  :
    _This(_this),
    _LogJointHistogram(logJointHistogram),
    _LogMarginalXHistogram(logMarginalXHistogram),
    _LogMarginalYHistogram(logMarginalYHistogram),
    _JointEntropy(je),
    _NormalizedMutualInformation(nmi)
  {}

  template <class TImage, class TIntensity, class TGradient>
  void operator ()(const TImage &, int idx, const TIntensity *tgt, const TIntensity *src, TGradient *deriv)
  {
    if (_This->IsForeground(idx)) {
      const int    target_nbins = _LogMarginalXHistogram.NumberOfBins();
      const int    source_nbins = _LogMarginalYHistogram.NumberOfBins();
      const double target_value = ValToRange(_LogMarginalXHistogram, *tgt);
      const double source_value = ValToRange(_LogMarginalYHistogram, *src);

      int t1 = static_cast<int>(     target_value ) - 1;
      int t2 = static_cast<int>(ceil(target_value)) + 1;
      int s1 = static_cast<int>(     source_value ) - 1;
      int s2 = static_cast<int>(ceil(source_value)) + 1;

      if (t1 <  0           ) t1 = 0;
      if (t2 >= target_nbins) t2 = target_nbins - 1;
      if (s1 <  0           ) s1 = 0;
      if (s2 >= source_nbins) s2 = source_nbins - 1;

      double jointEntropyGrad  = .0;
      double targetEntropyGrad = .0;
      double sourceEntropyGrad = .0;
      double w;

      for (int t = t1; t <= t2; ++t)
      for (int s = s1; s <= s2; ++s) {
        w = irtkBSpline<double>::B  (static_cast<double>(t) - target_value) *
            irtkBSpline<double>::B_I(static_cast<double>(s) - source_value);
        jointEntropyGrad  += w * _LogJointHistogram(t, s);
        targetEntropyGrad += w * _LogMarginalXHistogram(t);
        sourceEntropyGrad += w * _LogMarginalYHistogram(s);
      }

      (*deriv) = (targetEntropyGrad + sourceEntropyGrad - _NormalizedMutualInformation * jointEntropyGrad) / _JointEntropy;
    }
  }
};


} // namespace irtkNormalizedMutualImageInformationUtils
using namespace irtkNormalizedMutualImageInformationUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkNormalizedMutualImageInformation
::irtkNormalizedMutualImageInformation(const char *name)
:
  irtkProbabilisticImageSimilarity(name)
{
  if (version < irtkVersion(3, 1)) this->Weight(-1.0);
}

// -----------------------------------------------------------------------------
irtkNormalizedMutualImageInformation
::irtkNormalizedMutualImageInformation(const irtkNormalizedMutualImageInformation &other)
:
  irtkProbabilisticImageSimilarity(other)
{
}

// -----------------------------------------------------------------------------
irtkNormalizedMutualImageInformation::~irtkNormalizedMutualImageInformation()
{
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
double irtkNormalizedMutualImageInformation::Evaluate()
{
  const double nmi = _Histogram->NormalizedMutualInformation();
  if (version < irtkVersion(3, 1)) return       nmi;
  else                             return 2.0 - nmi;
}

// -----------------------------------------------------------------------------
double irtkNormalizedMutualImageInformation::RawValue(double value) const
{
  value = irtkProbabilisticImageSimilarity::RawValue(value);
  if (version >= irtkVersion(3, 1)) value = 2.0 - value;
  return value;
}

// -----------------------------------------------------------------------------
// This code is based on an idea from Marc Modat for computing the NMI
// derivative as suggested in his niftyreg package
bool irtkNormalizedMutualImageInformation
::NonParametricGradient(const irtkRegisteredImage *image, GradientImageType *gradient)
{
  IRTK_START_TIMING();

  const double sign = ((version < irtkVersion(3, 1)) ? 1.0 : -1.0);

  // Swap target and source if similarity derived w.r.t transformed "target"
  irtkRegisteredImage *fixed = Target();
  int tbin = Histogram()->NumberOfBinsX();
  int sbin = Histogram()->NumberOfBinsY();
  double tmin, smin, tmax, smax;
  Histogram()->GetMin(&tmin, &smin);
  Histogram()->GetMax(&tmax, &smax);

  if (image == Target()) {
    fixed = Source();
    swap(tbin, sbin);
    swap(tmin, smin);
    swap(tmax, smax);
  }

  // Compute joint entropy and normalized mutual information
  const double je  = _Histogram->JointEntropy();
  const double nmi = (_Histogram->EntropyX() + _Histogram->EntropyY()) / je;

  // Log transform histograms
  irtkHistogram_2D<double>  logJointHistogram(tbin, sbin);
  irtkHistogram_1D<double>  logMarginalXHistogram(tbin);
  irtkHistogram_1D<double>  logMarginalYHistogram(sbin);

  logJointHistogram.PutMin(tmin, smin);
  logJointHistogram.PutMax(tmax, smax);
  logMarginalXHistogram.PutMin(tmin);
  logMarginalXHistogram.PutMax(tmax);
  logMarginalYHistogram.PutMin(smin);
  logMarginalYHistogram.PutMax(smax);

  for (int s = 0; s < sbin; s++)
  for (int t = 0; t < tbin; t++) {
    double num = ((image == Target()) ? (*_Histogram)(s, t) : (*_Histogram)(t, s));
    logJointHistogram.Add(t, s, num);
    logMarginalXHistogram.Add(t, num);
    logMarginalYHistogram.Add(s, num);
  }

  logJointHistogram    .Log();
  logMarginalXHistogram.Log();
  logMarginalYHistogram.Log();

  // Evaluate similarity gradient w.r.t given transformed image
  CalculateGradient eval(this, logJointHistogram,
                               logMarginalXHistogram,
                               logMarginalYHistogram,
     /* denominator = */ sign * je * _Histogram->NumberOfSamples(), nmi);
  memset(gradient->Data(), 0, _NumberOfVoxels * sizeof(GradientType));
  ParallelForEachVoxel(fixed, image, gradient, eval);

  // Apply chain rule to obtain gradient w.r.t y = T(x)
  MultiplyByImageGradient(image, gradient);

  IRTK_DEBUG_TIMING(2, "similarity gradient computation (NMI)");
  return true;
}
