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

#include <irtkInexactLineSearch.h>


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkInexactLineSearch::irtkInexactLineSearch(irtkObjectiveFunction *f)
:
  irtkLineSearch(f),
  _MaxRejectedStreak      (_NumberOfIterations / 4),
  _ReusePreviousStepLength(true),
  _StrictStepLengthRange  (1),
  _AllowSignChange        (NULL),
  _CurrentDoFValues       (NULL),
  _ScaledDirection        (NULL)
{
  if (f) {
    Allocate(_CurrentDoFValues, f->NumberOfDOFs());
    Allocate(_ScaledDirection,  f->NumberOfDOFs());
  }
}

// -----------------------------------------------------------------------------
irtkInexactLineSearch::irtkInexactLineSearch(const irtkInexactLineSearch &other)
:
  irtkLineSearch(other),
  _MaxRejectedStreak      (other._MaxRejectedStreak),
  _ReusePreviousStepLength(other._ReusePreviousStepLength),
  _StrictStepLengthRange  (other._StrictStepLengthRange),
  _AllowSignChange        (other._AllowSignChange),
  _CurrentDoFValues       (NULL),
  _ScaledDirection        (NULL)
{
  if (Function()) {
    Allocate(_CurrentDoFValues, Function()->NumberOfDOFs());
    Allocate(_ScaledDirection,  Function()->NumberOfDOFs());
  }
}

// -----------------------------------------------------------------------------
irtkInexactLineSearch &irtkInexactLineSearch::operator =(const irtkInexactLineSearch &other)
{
  Deallocate(_CurrentDoFValues);
  Deallocate(_ScaledDirection);
  irtkLineSearch::operator =(other);
  _MaxRejectedStreak       = other._MaxRejectedStreak;
  _ReusePreviousStepLength = other._ReusePreviousStepLength;
  _StrictStepLengthRange   = other._StrictStepLengthRange;
  _AllowSignChange         = other._AllowSignChange;
  if (Function()) {
    Allocate(_CurrentDoFValues, Function()->NumberOfDOFs());
    Allocate(_ScaledDirection,  Function()->NumberOfDOFs());
  }
  return *this;
}

// -----------------------------------------------------------------------------
irtkInexactLineSearch::~irtkInexactLineSearch()
{
  Deallocate(_CurrentDoFValues);
  Deallocate(_ScaledDirection);
}

// -----------------------------------------------------------------------------
void irtkInexactLineSearch::Function(irtkObjectiveFunction *f)
{
  if (Function() != f) {
    Deallocate(_CurrentDoFValues);
    Deallocate(_ScaledDirection);
    irtkLineSearch::Function(f);
    if (f) {
      Allocate(_CurrentDoFValues, f->NumberOfDOFs());
      Allocate(_ScaledDirection,  f->NumberOfDOFs());
    }
  }
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkInexactLineSearch::Set(const char *name, const char *value)
{
  // Maximum number of consectutive rejections
  if (strcmp(name, "Maximum streak of rejected steps") == 0) {
    return FromString(value, _MaxRejectedStreak);
  // Whether to start new search using step length of previous search
  } else if (strcmp(name, "Reuse previous step length") == 0) {
    return FromString(value, _ReusePreviousStepLength);
  // Whether [min, max] step length range is strict
  } else if (strcmp(name, "Strict step length range")             == 0 ||
             strcmp(name, "Strict incremental step length range") == 0) {
    bool limit_increments;
    if (!FromString(value, limit_increments)) return false;
    _StrictStepLengthRange = limit_increments ? 1 : 0;
    return true;
  } else if (strcmp(name, "Strict total step length range")       == 0 ||
             strcmp(name, "Strict accumulated step length range") == 0) {
    bool limit_step;
    if (!FromString(value, limit_step)) return false;
    _StrictStepLengthRange = limit_step ? 2 : 0;
    return true;
  }
  return irtkLineSearch::Set(name, value);
}

// -----------------------------------------------------------------------------
irtkParameterList irtkInexactLineSearch::Parameter() const
{
  irtkParameterList params = irtkLineSearch::Parameter();
  Insert(params, "Maximum streak of rejected steps", _MaxRejectedStreak);
  Insert(params, "Reuse previous step length",       _ReusePreviousStepLength);
  if (_StrictStepLengthRange == 2) {
    Insert(params, "Strict total step length range", true);
  } else {
    Insert(params, "Strict incremental step length range", static_cast<bool>(_StrictStepLengthRange));
  }
  return params;
}

// =============================================================================
// Optimization
// =============================================================================

// -----------------------------------------------------------------------------
double irtkInexactLineSearch::Advance(double alpha)
{
  if (_StepLengthUnit == .0) return .0;
  // Backup current function parameter values
  Function()->Get(_CurrentDoFValues);
  // Compute gradient for given step length
  alpha /= _StepLengthUnit;
  if (_Descent) alpha *= -1.0;
  const int ndofs = Function()->NumberOfDOFs();
  for (int dof = 0; dof < ndofs; ++dof) {
    _ScaledDirection[dof] = alpha * _Direction[dof];
  }
  // Set scaled gradient to the negative of the current parameter value if sign
  // changes are not allowed for this parameter s.t. updated value is zero
  //
  // Note: This is used for the optimization of the registration cost
  //       function with L1-norm sparsity constraint on the multi-level
  //       free-form deformation parameters (Wenzhe et al.'s Sparse FFD).
  if (_AllowSignChange) {
    double next_value;
    for (int dof = 0; dof < ndofs; ++dof) {
      if (_AllowSignChange[dof]) continue;
      next_value = _CurrentDoFValues[dof] + _ScaledDirection[dof];
      if ((_CurrentDoFValues[dof] * next_value) <= .0) {
        _ScaledDirection[dof] = - _CurrentDoFValues[dof];
      }
    }
  }
  // Update all parameters at once to only trigger a single modified event
  return Function()->Step(_ScaledDirection);
}

// -----------------------------------------------------------------------------
void irtkInexactLineSearch::Retreat(double alpha)
{
  if (_StepLengthUnit == .0) return;
  Function()->Put(_CurrentDoFValues);
}
