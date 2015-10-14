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

#ifndef _IRTKINTENSITYCROSSCORRELATION_H
#define _IRTKINTENSITYCROSSCORRELATION_H

#include <irtkImageSimilarity.h>


/**
 * Cross correlation image similarity measure
 */
class irtkIntensityCrossCorrelation : public irtkImageSimilarity
{
  irtkObjectMacro(irtkIntensityCrossCorrelation);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  irtkIntensityCrossCorrelation(const char * = "");

  /// Copy constructor
  irtkIntensityCrossCorrelation(const irtkIntensityCrossCorrelation &);

  /// Assignment operator
  irtkIntensityCrossCorrelation &operator =(const irtkIntensityCrossCorrelation &);

  /// Destructor
  ~irtkIntensityCrossCorrelation();

  // ---------------------------------------------------------------------------
  // Evaluation
protected:

  /// Evaluate similarity of images
  virtual double Evaluate();

  /// Evaluate non-parametric similarity gradient w.r.t the given image
  virtual bool NonParametricGradient(const irtkRegisteredImage *, GradientImageType *);

};


#endif
