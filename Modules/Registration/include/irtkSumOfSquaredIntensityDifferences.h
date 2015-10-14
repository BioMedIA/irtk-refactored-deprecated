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

#ifndef _IRTKSUMOFSQUAREDINTENSITYDIFFERENCES_H
#define _IRTKSUMOFSQUAREDINTENSITYDIFFERENCES_H

#include <irtkImageSimilarity.h>


/**
 * Sum of squared differences image similarity measure
 */
class irtkSumOfSquaredIntensityDifferences : public irtkImageSimilarity
{
  irtkObjectMacro(irtkSumOfSquaredIntensityDifferences);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Maximum squared intensity difference used for normalization
  irtkAttributeMacro(double, MaxSqDiff);

  /// Sum of squared intensity difference value
  double _Value;

  /// Number of foreground voxels for which similarity is evaluated
  int _N;

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  irtkSumOfSquaredIntensityDifferences(const char * = "");

  /// Copy constructor
  irtkSumOfSquaredIntensityDifferences(const irtkSumOfSquaredIntensityDifferences &);

  /// Assignment operator
  irtkSumOfSquaredIntensityDifferences &operator =(const irtkSumOfSquaredIntensityDifferences &);

  /// Destructor
  ~irtkSumOfSquaredIntensityDifferences();

  // ---------------------------------------------------------------------------
  // Initialization

  /// Initialize similarity measure
  virtual void Initialize();

  /// Update moving input image(s) and internal state of similarity measure
  virtual void Update(bool = true);

  // ---------------------------------------------------------------------------
  // Evaluation
protected:

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

  /// Evaluate similarity of images
  virtual double Evaluate();

  /// Evaluate non-parametric similarity gradient w.r.t the given image
  virtual bool NonParametricGradient(const irtkRegisteredImage *, GradientImageType *);

};


#endif
