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

#include <irtkAdaptiveLineSearch.h>


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkAdaptiveLineSearch::irtkAdaptiveLineSearch(irtkObjectiveFunction *f)
:
  irtkInexactLineSearch(f),
  _StepLengthRise(1.1),
  _StepLengthDrop(0.5)
{
}

// -----------------------------------------------------------------------------
irtkAdaptiveLineSearch::irtkAdaptiveLineSearch(const irtkAdaptiveLineSearch &other)
:
  irtkInexactLineSearch(other),
  _StepLengthRise(other._StepLengthRise),
  _StepLengthDrop(other._StepLengthDrop)
{
}

// -----------------------------------------------------------------------------
irtkAdaptiveLineSearch &irtkAdaptiveLineSearch::operator =(const irtkAdaptiveLineSearch &other)
{
  irtkInexactLineSearch::operator =(other);
  _StepLengthRise = other._StepLengthRise;
  _StepLengthDrop = other._StepLengthDrop;
  return *this;
}

// -----------------------------------------------------------------------------
irtkAdaptiveLineSearch::~irtkAdaptiveLineSearch()
{
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkAdaptiveLineSearch::Set(const char *name, const char *value)
{
  if (strcmp(name, "Step length rise") == 0) {
    return FromString(value, _StepLengthRise) && _StepLengthRise > .0;
  } else if (strcmp(name, "Step length drop") == 0) {
    return FromString(value, _StepLengthDrop) && _StepLengthDrop > .0;
  }
  return irtkInexactLineSearch::Set(name, value);
}

// -----------------------------------------------------------------------------
irtkParameterList irtkAdaptiveLineSearch::Parameter() const
{
  irtkParameterList params = irtkInexactLineSearch::Parameter();
  Insert(params, "Step length rise", ToString(_StepLengthRise));
  Insert(params, "Step length drop", ToString(_StepLengthDrop));
  return params;
}

// =============================================================================
// Optimization
// =============================================================================

// -----------------------------------------------------------------------------
double irtkAdaptiveLineSearch::Run()
{
  const double min_length   = _MinStepLength;
  const double max_length   = _MaxStepLength;
  const bool   strict_range = (min_length == max_length || _StrictStepLengthRange);
  const bool   strict_total = (_StrictStepLengthRange == 2);

  // Number of consecutively rejected steps after at least one accepted step
  int rejected_streak = -1;

  // Get current value of objective function
  irtkLineSearch::Run();

  // Define "aliases" for event data for prettier code below
  irtkLineSearchStep step;
  step._Info        = "incremental";
  step._Direction   = _Direction;
  step._Unit        = _StepLengthUnit;
  step._MinLength   = min_length;
  step._MaxLength   = max_length;
  step._TotalLength = .0;
  double &value     = step._Value;
  double &current   = step._Current;
  double &alpha     = step._Length;
  double &total     = step._TotalLength;
  double &delta     = step._Delta;

  // Start line search
  alpha   = max_length;
  total   = .0;            // accumulated length of accepted steps
  value   = _CurrentValue;
  current = _CurrentValue;
  delta   = .0;

  // Continue search using total step length of previous search
  if (_ReusePreviousStepLength && _StepLength > .0) alpha = _StepLength;

  // Limit incremental step length strictly to [min_length, max_length]
  if (strict_range) {
    if      (alpha > max_length) alpha = max_length;
    else if (alpha < min_length) alpha = min_length;
  }

  // Notify observers about start of line search
  Broadcast(LineSearchStartEvent, &step);

  // Walk along search direction until no further improvement
  irtkIteration iteration(0, _NumberOfIterations);
  while (iteration.Next()) {

    // Notify observers about start of iteration
    Broadcast(LineSearchIterationStartEvent, &iteration);

    // Take step along search direction
    delta = Advance(alpha);

    // Check minimum maximum change of DoFs convergence criterium
    if (delta <= _Delta) {

      Retreat(alpha);
      break;

    } else {

      // Re-evaluate objective function
      Function()->Update(false);
      value = Function()->Value();

      // Check minimum change of objective function value convergence criterium
      if (_Descent ? (value < current - _Epsilon)
                   : (value > current + _Epsilon)) {

        // Notify observers about progress
        Broadcast(AcceptedStepEvent, &step);

        // New current value
        current = value;

        // Accumulate accepted steps
        total += alpha;

        // Increase step length and continue search from new point
        // If the step length range is not strict, keep possibly
        // greater previous step length as it was accepted again
        if (strict_range || alpha < max_length) {
          alpha = min(alpha * _StepLengthRise, max_length);
          if (strict_total && total + alpha > max_length) {
            alpha = max_length - total;
            if (alpha < 1e-12) break;
          }
        }

        // Reset counter of consecutive rejections
        rejected_streak = 0;

      } else {

        // Notify observers about progress
        Broadcast(RejectedStepEvent, &step);

        // Reject step and start over at previous point
        Retreat(alpha);

        // Stop search if maximum number of consecutive rejections exceeded
        if (rejected_streak != -1) {
          ++rejected_streak;
          if (_MaxRejectedStreak >= 0 && rejected_streak > _MaxRejectedStreak) break;
        }

        // Break in case of computational error
        if (IsNaN(value)) {
          cerr << "irtkAdaptiveLineSearch::Run:" << __LINE__ << ": NaN objective function value!" << endl;
          exit(1);
        }

        // Decrease step length and continue if possible
        if (2 < version.Major() && version.Major() < 3) {
          alpha *= _StepLengthDrop;
          if (alpha < min_length) break;
        } else if (version < irtkVersion(3, 1)) {
          if (alpha == min_length) break;
          alpha = max(alpha * _StepLengthDrop, min_length);
        } else {
          if (alpha == min_length) break;
          // If previous step length was reused even though it was outside
          // the considered [min, max] step length range, start over with
          // step length limited by the original range as it would have been
          // if previous length had not been reused.
          if (alpha > max_length) alpha = max_length;
          else alpha = max(alpha * _StepLengthDrop, min_length);
        }
      }
    }

    // Notify observers about end of iteration
    Broadcast(LineSearchIterationEndEvent, &iteration);

  }

  // Notify observers about end of line search
  Broadcast(LineSearchEndEvent, &step);

  // Re-use final step length upon next line search
  _StepLength = total;

  // Return new value of objective function
  return current;
}
