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

#ifndef _IRTKCOSINEOFNORMALIZEDGRADIENTFIELD_H
#define _IRTKCOSINEOFNORMALIZEDGRADIENTFIELD_H

#include <irtkNormalizedGradientFieldSimilarity.h>


/**
 * Cosine of normalized gradient field similarity measure
 */
class irtkCosineOfNormalizedGradientField : public irtkNormalizedGradientFieldSimilarity
{
  irtkObjectMacro(irtkCosineOfNormalizedGradientField);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Power (exponent) of cosine function
  irtkPublicAttributeMacro(int, Power);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  irtkCosineOfNormalizedGradientField(const char * = "");

  /// Copy constructor
  irtkCosineOfNormalizedGradientField(const irtkCosineOfNormalizedGradientField &);

  /// Assignment operator
  irtkCosineOfNormalizedGradientField &operator =(const irtkCosineOfNormalizedGradientField &);

  /// Destructor
  ~irtkCosineOfNormalizedGradientField();

  // ---------------------------------------------------------------------------
  // Parameters

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  // Import other overloads
  using irtkImageSimilarity::Parameter;

  /// Get parameter key/value as string map
  virtual irtkParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Evaluation
protected:

  /// Evaluate similarity of images
  virtual double Evaluate();

  /// Evaluate non-parametric similarity gradient w.r.t the given image
  virtual bool NonParametricGradient(const irtkRegisteredImage *, GradientImageType *);

  // ---------------------------------------------------------------------------
  // Debugging
public:

  /// Return unweighted and unnormalized raw energy term value
  /// \remarks Use for progress reporting only.
  virtual double RawValue(double) const;

  /// Write input of data fidelity term
  virtual void WriteDataSets(const char *, const char *, bool = true) const;

};


#endif
