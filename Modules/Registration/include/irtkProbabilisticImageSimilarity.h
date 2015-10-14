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

#ifndef _IRTKPROBABILISTICIMAGESIMILARITY_H
#define _IRTKPROBABILISTICIMAGESIMILARITY_H

#include <irtkImageSimilarity.h>

#include <irtkHistogram_1D.h>
#include <irtkHistogram_2D.h>


/**
 * Base class for probabilistic image similarity measures
 *
 * Subclasses of this intensity-based image similarity measure compute similarity
 * from the joint and marginal probabilities of the intensities in the images.
 * An estimate of the probabilities is obtained using a joint histogram and
 * cubic B-spline Parzen Windows for a continuous representation.
 */

class irtkProbabilisticImageSimilarity : public irtkImageSimilarity
{
  irtkAbstractMacro(irtkProbabilisticImageSimilarity);

  // ---------------------------------------------------------------------------
  // Types
public:

  /// Type of joint histogram
  typedef irtkHistogram_2D<double> JointHistogramType;

  // ---------------------------------------------------------------------------
  // Attributes

  /// Joint histogram of image intensities (unsmoothed samples)
  irtkComponentMacro(JointHistogramType, Samples);

  /// Joint histogram of image intensities (cubic B-spline Parzen windows)
  irtkComponentMacro(JointHistogramType, Histogram);

  /// Number of histogram bins for target image intensities
  irtkPublicAttributeMacro(int, NumberOfTargetBins);

  /// Number of histogram bins for source image intensities
  irtkPublicAttributeMacro(int, NumberOfSourceBins);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
protected:

  /// Constructor
  irtkProbabilisticImageSimilarity(const char * = "", double = 1.0);

  /// Copy constructor
  irtkProbabilisticImageSimilarity(const irtkProbabilisticImageSimilarity &);

  /// Assignment operator
  irtkProbabilisticImageSimilarity &operator =(const irtkProbabilisticImageSimilarity &);

  /// Destructor
  virtual ~irtkProbabilisticImageSimilarity();

  // ---------------------------------------------------------------------------
  // Parameters
public:

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  // Do not hide other overload
  using irtkImageSimilarity::Parameter;

  /// Get parameter key/value as string map
  virtual irtkParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Initialization

  /// Initialize similarity measure once input and parameters have been set
  virtual void Initialize();

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Update moving image and internal state of similarity measure
  virtual void Update(bool = true);

  /// Exclude region from similarity evaluation
  ///
  /// Called by ApproximateGradient \b before the registered image region of
  /// the transformed image is updated.
  virtual void Exclude(const blocked_range3d<int> &);

  /// Include region in similarity evaluation
  ///
  /// Called by ApproximateGradient \b after the registered image region of
  /// the transformed image is updated.
  virtual void Include(const blocked_range3d<int> &);

  // ---------------------------------------------------------------------------
  // Debugging

  /// Print debug information
  virtual void Print(irtkIndent = 0) const;

  /// Write input of data fidelity term
  virtual void WriteDataSets(const char *, const char *, bool = true) const;

};


#endif
