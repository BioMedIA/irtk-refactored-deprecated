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

#include <irtkSmoothnessConstraint.h>


// -----------------------------------------------------------------------------
irtkSmoothnessConstraint::irtkSmoothnessConstraint(const char *name, double weight)
:
  irtkTransformationConstraint(name, weight),
  _WithRespectToWorld(false),
  _AnnealingRate(.0), _AnnealingWeight(1.0)
{
  _ParameterPrefix.push_back("Smoothness ");
  _ParameterPrefix.push_back("Bending energy ");
  _ParameterPrefix.push_back("Bending ");
}

// -----------------------------------------------------------------------------
bool irtkSmoothnessConstraint::Set(const char *param, const char *value)
{
  const string name = ParameterNameWithoutPrefix(param);
  if (name.find("W.r.t world"          ) == 0 ||
      name.find("W.r.t. world"         ) == 0 ||
      name.find("Wrt world"            ) == 0 ||
      name.find("With respect to world") == 0) {
    return FromString(value, _WithRespectToWorld);
  }
  if (name.find("Annealing rate") == 0) {
    return FromString(value, _AnnealingRate);
  }
  return irtkTransformationConstraint::Set(param, value);
}

// -----------------------------------------------------------------------------
irtkParameterList irtkSmoothnessConstraint::Parameter() const
{
  irtkParameterList params = irtkTransformationConstraint::Parameter();
  if (!_Name.empty()) {
    Insert(params, _Name + " w.r.t. world (experimental)",   ToString(_WithRespectToWorld));
    Insert(params, _Name + " annealing rate (experimental)", ToString(_AnnealingRate));
  }
  return params;
}

// -----------------------------------------------------------------------------
void irtkSmoothnessConstraint::Initialize()
{
  irtkTransformationConstraint::Initialize();
  _AnnealingWeight = 1.0;
}

// -----------------------------------------------------------------------------
bool irtkSmoothnessConstraint::Upgrade()
{
  irtkTransformationConstraint::Upgrade();
  if (_AnnealingRate != .0) _AnnealingWeight *= _AnnealingRate;
  return false; // whether to continue depends on similarity term, not penalty
}

// -----------------------------------------------------------------------------
double irtkSmoothnessConstraint::Evaluate()
{
  const irtkMultiLevelTransformation *mffd = NULL;
  const irtkFreeFormTransformation   *ffd  = NULL;

  (mffd = MFFD()) || (ffd = FFD());

  double bending = .0;
  if (mffd) {
    for (int l = 0; l < mffd->NumberOfLevels(); ++l) {
      if (!mffd->LocalTransformationIsActive(l)) continue;
      ffd = mffd->GetLocalTransformation(l);
      bending += ffd->BendingEnergy(_ConstrainPassiveDoFs, _WithRespectToWorld);
    }
  } else if (ffd) {
    bending = ffd->BendingEnergy(_ConstrainPassiveDoFs, _WithRespectToWorld);
  }
  return _AnnealingWeight * bending;
}

// -----------------------------------------------------------------------------
void irtkSmoothnessConstraint::EvaluateGradient(double *gradient, double, double weight)
{
  weight *= _AnnealingWeight;
  const irtkMultiLevelTransformation *mffd = NULL;
  const irtkFreeFormTransformation   *ffd  = NULL;

  (mffd = MFFD()) || (ffd = FFD());

  if (mffd) {
    for (int l = 0; l < mffd->NumberOfLevels(); ++l) {
      if (!mffd->LocalTransformationIsActive(l)) continue;
      ffd = mffd->GetLocalTransformation(l);
      ffd->BendingEnergyGradient(gradient, weight, _ConstrainPassiveDoFs, _WithRespectToWorld);
      gradient += ffd->NumberOfDOFs();
    }
  } else if (ffd) {
    ffd->BendingEnergyGradient(gradient, weight, _ConstrainPassiveDoFs, _WithRespectToWorld);
  }
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
void irtkSmoothnessConstraint
::WriteFFDGradient(const char *fname, const irtkFreeFormTransformation *ffd, const double *g) const
{
  typedef irtkFreeFormTransformation::CPValue CPValue;
  typedef irtkGenericImage<CPValue>           CPImage;
  CPValue *data = reinterpret_cast<CPValue *>(const_cast<double *>(g));
  CPImage gradient(ffd->Attributes(), data);
  gradient.Write(fname);
}

// -----------------------------------------------------------------------------
void irtkSmoothnessConstraint::WriteGradient(const char *p, const char *suffix) const
{
  const int   sz = 1024;
  char        fname[sz];
  string      _prefix = Prefix(p);
  const char *prefix  = _prefix.c_str();

  const irtkMultiLevelTransformation *mffd = NULL;
  const irtkFreeFormTransformation   *ffd  = NULL;

  (mffd = MFFD()) || (ffd = FFD());

  if (mffd) {
    for (int l = 0; l < mffd->NumberOfLevels(); ++l) {
      if (!mffd->LocalTransformationIsActive(l)) continue;
      ffd = mffd->GetLocalTransformation(l);
      double *gradient = CAllocate<double>(ffd->NumberOfDOFs());
      ffd->BendingEnergyGradient(gradient, 1.0, _ConstrainPassiveDoFs, _WithRespectToWorld);
      if (mffd->NumberOfActiveLevels() == 1) {
        snprintf(fname, sz, "%sgradient%s.nii.gz", prefix, suffix);
      } else {
        snprintf(fname, sz, "%sgradient_of_ffd_at_level_%d%s.nii.gz", prefix, l+1, suffix);
      }
      WriteFFDGradient(fname, ffd, gradient);
      Deallocate(gradient);
    }
  } else if (ffd) {
    snprintf(fname, sz, "%sgradient%s.nii.gz", prefix, suffix);
    double *gradient = CAllocate<double>(ffd->NumberOfDOFs());
    ffd->BendingEnergyGradient(gradient, 1.0, _ConstrainPassiveDoFs, _WithRespectToWorld);
    WriteFFDGradient(fname, ffd, gradient);
    Deallocate(gradient);
  }
}
