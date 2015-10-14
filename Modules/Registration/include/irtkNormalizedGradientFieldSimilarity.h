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

#ifndef _IRTKNORMALIZEDGRADIENTFIELDSIMILARITY_H
#define _IRTKNORMALIZEDGRADIENTFIELDSIMILARITY_H

#include <irtkGradientFieldSimilarity.h>

/**
 * Base class for normalized gradient image similarity measures
 *
 * Subclasses of this intensity-based image similarity measure compute similarity
 * from the normalized gradient field of the images.
 */

class irtkNormalizedGradientFieldSimilarity : public irtkGradientFieldSimilarity
{
  irtkAbstractMacro(irtkNormalizedGradientFieldSimilarity);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Noise parameter for target image
  irtkPublicAttributeMacro(double, TargetNoise);

  /// Noise parameter for source image
  irtkPublicAttributeMacro(double, SourceNoise);

  /// Normalization constant for target image
  irtkPublicAttributeMacro(double, TargetNormalization);

  /// Normalization constant for source image
  irtkPublicAttributeMacro(double, SourceNormalization);

  /// Compute normalization constant for given image and noise value
  double NormalizationFactor(irtkRegisteredImage *, double = 1.0) const;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Constructor
  irtkNormalizedGradientFieldSimilarity(const char * = "", double = 1.0);

  /// Copy constructor
  irtkNormalizedGradientFieldSimilarity(const irtkNormalizedGradientFieldSimilarity &);

  /// Destructor
  virtual ~irtkNormalizedGradientFieldSimilarity();

  // ---------------------------------------------------------------------------
  // Parameters

public:

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  // Import other overloads
  using irtkGradientFieldSimilarity::Parameter;

  /// Get parameter name/value map
  virtual irtkParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Update moving input image(s) and internal state of similarity measure
  virtual void Update(bool = true);

};


#endif
