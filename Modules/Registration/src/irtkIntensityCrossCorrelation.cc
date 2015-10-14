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

#include <irtkIntensityCrossCorrelation.h>


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkIntensityCrossCorrelation
::irtkIntensityCrossCorrelation(const char *name)
:
  irtkImageSimilarity(name)
{
}

// -----------------------------------------------------------------------------
irtkIntensityCrossCorrelation
::irtkIntensityCrossCorrelation(const irtkIntensityCrossCorrelation &other)
:
  irtkImageSimilarity(other)
{
}

// -----------------------------------------------------------------------------
irtkIntensityCrossCorrelation::~irtkIntensityCrossCorrelation()
{
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
double irtkIntensityCrossCorrelation::Evaluate()
{
  VoxelType *tgt = Target()->GetPointerToVoxels();
  VoxelType *src = Source()->GetPointerToVoxels();

  double x = .0, y = .0, xy = .0, x2 = .0, y2 = .0;
  int    n =  0;
  for (int idx = 0; idx < NumberOfVoxels(); ++idx) {
    if (IsForeground(idx)) {
      x  += *tgt;
      y  += *src;
      xy += (*tgt) * (*src);
      x2 += (*tgt) * (*tgt);
      y2 += (*tgt) * (*tgt);
      ++n;
    }
  }

  if (n == 0) return .0;
  return (xy - (x * y) / n) / (sqrt(x2 - x * x / n) * sqrt(y2 - y *y / n));
}

// -----------------------------------------------------------------------------
bool irtkIntensityCrossCorrelation
::NonParametricGradient(const irtkRegisteredImage *image, GradientImageType *gradient)
{
  cerr << "irtkIntensityCrossCorrelation::NonParametricGradient: Not implemented" << endl;

  // Apply chain rule to obtain gradient w.r.t y = T(x)
  MultiplyByImageGradient(image, gradient);
  return true;
}
