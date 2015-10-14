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

#include <irtkMutualImageInformation.h>


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkMutualImageInformation
::irtkMutualImageInformation(const char *name)
:
  irtkProbabilisticImageSimilarity(name)
{
  if (version < irtkVersion(3, 1)) this->Weight(-1.0);
}

// -----------------------------------------------------------------------------
irtkMutualImageInformation
::irtkMutualImageInformation(const irtkMutualImageInformation &other)
:
  irtkProbabilisticImageSimilarity(other)
{
}

// -----------------------------------------------------------------------------
irtkMutualImageInformation::~irtkMutualImageInformation()
{
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
double irtkMutualImageInformation::Evaluate()
{
  if (version < irtkVersion(3, 1)) {
    return _Histogram->MutualInformation();
  } else {
    return 2.0 - _Histogram->MutualInformation();
  }
}

// -----------------------------------------------------------------------------
double irtkMutualImageInformation::RawValue(double value) const
{
  value = irtkProbabilisticImageSimilarity::RawValue(value);
  if (version >= irtkVersion(3, 1)) value = 2.0 - value;
  return value;
}

// -----------------------------------------------------------------------------
bool irtkMutualImageInformation
::NonParametricGradient(const irtkRegisteredImage *image, GradientImageType *gradient)
{
  // TODO Implement gradient computation of MI, using NMI as reference
  cerr << "irtkMutualImageInformation::NonParametricGradient: Not implemented" << endl;

  // Apply chain rule to obtain gradient w.r.t y = T(x)
  MultiplyByImageGradient(image, gradient);
  return true;
}
