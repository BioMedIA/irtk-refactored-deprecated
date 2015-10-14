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

#include <irtkNormalizedGradientFieldSimilarity.h>


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkNormalizedGradientFieldSimilarity
::irtkNormalizedGradientFieldSimilarity(const char *name, double weight)
:
  irtkGradientFieldSimilarity(name, weight),
  _TargetNoise(1.0),
  _SourceNoise(1.0),
  _TargetNormalization (100.0),
  _SourceNormalization (100.0)
{
}

// -----------------------------------------------------------------------------
irtkNormalizedGradientFieldSimilarity::irtkNormalizedGradientFieldSimilarity(const irtkNormalizedGradientFieldSimilarity &other)
:
  irtkGradientFieldSimilarity(other),
  _TargetNoise(other._TargetNoise),
  _SourceNoise(other._SourceNoise),
  _TargetNormalization (other._TargetNormalization),
  _SourceNormalization (other._SourceNormalization)
{
}

// -----------------------------------------------------------------------------
irtkNormalizedGradientFieldSimilarity::~irtkNormalizedGradientFieldSimilarity()
{
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkNormalizedGradientFieldSimilarity::Set(const char *name, const char *value)
{
  if (strcmp(name, "Target noise") == 0) {
    return FromString(value, _TargetNoise) && _TargetNoise > 0;
  } else if (strcmp(name, "Source noise") == 0) {
    return FromString(value, _SourceNoise) && _SourceNoise > 0;
  }
  return irtkGradientFieldSimilarity::Set(name, value);
}

// -----------------------------------------------------------------------------
irtkParameterList irtkNormalizedGradientFieldSimilarity::Parameter() const
{
  irtkParameterList params = irtkGradientFieldSimilarity::Parameter();
  Insert(params, "Target noise", ToString(_TargetNoise));
  Insert(params, "Source noise", ToString(_SourceNoise));
  return params;
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
double irtkNormalizedGradientFieldSimilarity::NormalizationFactor(irtkRegisteredImage *image, double noise) const
{
  double norm = .0;
  int    n    = 0;

  VoxelType *dx = image->GetPointerToVoxels(0, 0, 0, 1);
  VoxelType *dy = image->GetPointerToVoxels(0, 0, 0, 2);
  VoxelType *dz = image->GetPointerToVoxels(0, 0, 0, 3);

  for (int idx = 0; idx < NumberOfVoxels(); ++idx, ++dx, ++dy, ++dz) {
    if (IsForeground(idx)) {
      norm += sqrt((*dx) * (*dx) + (*dy) * (*dy) + (*dz) * (*dz));
      ++n;
    }
  }

  return ((n > 0) ? (norm * noise / n) : .0);
}

// -----------------------------------------------------------------------------
void irtkNormalizedGradientFieldSimilarity::Update(bool hessian)
{
  const bool initial_update = _InitialUpdate;
  irtkGradientFieldSimilarity::Update(hessian); // sets _InitialUpdate = false
  if (initial_update || _Target->Transformation()) {
    _TargetNormalization = NormalizationFactor(_Target, _TargetNoise);
    if (_TargetNormalization == 0) {
      cerr << "irtkNormalizedGradientFieldSimilarity::Update: Target image has no structure!" << endl;
      exit(1);
    }
  }
  if (initial_update || _Source->Transformation()) {
    _SourceNormalization = NormalizationFactor(_Source, _SourceNoise);
    if (_SourceNormalization == 0) {
      cerr << "irtkNormalizedGradientFieldSimilarity::Update: Source image has no structure!" << endl;
      exit(1);
    }
  }
}
