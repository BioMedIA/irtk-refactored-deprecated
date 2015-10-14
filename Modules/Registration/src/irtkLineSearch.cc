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

#include <irtkLineSearch.h>

#include <irtkMaxStepLineSearch.h>
#include <irtkAdaptiveLineSearch.h>


// =============================================================================
// Factory
// =============================================================================

// -----------------------------------------------------------------------------
irtkLineSearch *irtkLineSearch::New(irtkLineSearchStrategy &strategy, irtkObjectiveFunction *f)
{
  switch (strategy) {
    case LS_None:     return new irtkMaxStepLineSearch (f);
    case LS_Adaptive: return new irtkAdaptiveLineSearch(f);
    default:
      cerr << "irtkLineSearch::New: Unknown line search strategy: " << strategy << endl;
      exit(1);
  }
  return NULL;
}

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkLineSearch::irtkLineSearch(irtkObjectiveFunction *f)
:
  irtkLocalOptimizer(f),
  _Descent           (true),
  _Direction         (NULL),
  _CurrentValue      (numeric_limits<double>::quiet_NaN()),
  _NumberOfIterations(20),
  _MinStepLength     ( 0.1),
  _MaxStepLength     (10.0),
  _StepLengthUnit    ( 1.0),
  _StepLength        ( 0.0)
{
}

// -----------------------------------------------------------------------------
irtkLineSearch::irtkLineSearch(const irtkLineSearch &other)
:
  irtkLocalOptimizer(other),
  _Descent           (other._Descent),
  _Direction         (other._Direction),
  _CurrentValue      (other._CurrentValue),
  _NumberOfIterations(other._NumberOfIterations),
  _MinStepLength     (other._MinStepLength),
  _MaxStepLength     (other._MaxStepLength),
  _StepLengthUnit    (other._StepLengthUnit),
  _StepLength        (other._StepLength)
{
}

// -----------------------------------------------------------------------------
irtkLineSearch &irtkLineSearch::operator =(const irtkLineSearch &other)
{
  irtkLocalOptimizer::operator =(other);
  _Descent            = other._Descent;
  _Direction          = other._Direction;
  _CurrentValue       = other._CurrentValue;
  _NumberOfIterations = other._NumberOfIterations;
  _MinStepLength      = other._MinStepLength;
  _MaxStepLength      = other._MaxStepLength;
  _StepLengthUnit     = other._StepLengthUnit;
  _StepLength         = other._StepLength;
  return *this;
}

// -----------------------------------------------------------------------------
irtkLineSearch::~irtkLineSearch()
{
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkLineSearch::Set(const char *name, const char *value)
{
  // Number of line search iterations
  if (strcmp(name, "Maximum no. of line search iterations")    == 0 ||
      strcmp(name, "Maximum number of line search iterations") == 0 ||
      strcmp(name, "Maximum no. of line iterations")           == 0 ||
      strcmp(name, "Maximum number of line iterations")        == 0 ||
      strcmp(name,    "No. of line search iterations")         == 0 ||
      strcmp(name, "Number of line search iterations")         == 0 ||
      strcmp(name,    "No. of line iterations")                == 0 ||
      strcmp(name, "Number of line iterations")                == 0) {
    return FromString(value, _NumberOfIterations) && _NumberOfIterations > 0;
  // Minimum length of step
  } else if (strcmp(name, "Minimum length of steps") == 0) {
    return FromString(value, _MinStepLength) && _MinStepLength > .0;
  // Maximum length of step
  } else if (strcmp(name, "Maximum length of steps") == 0) {
    return FromString(value, _MaxStepLength) && _MaxStepLength > .0;
  // Length of steps
  } else if (strcmp(name, "Length of steps") == 0) {
    if (FromString(value, _MaxStepLength) && _MaxStepLength > .0) {
      _MinStepLength = _MaxStepLength;
    }
  // Range of step lengths (e.g., "[0.1 1]" or "0.1 1")
  } else if (strcmp(name, "Range of steps")                  == 0 ||
             strcmp(name, "Range of step lengths")           == 0 ||
             strcmp(name, "Min-/Maximum length of steps")    == 0 ||
             strcmp(name, "Min-/maximum length of steps")    == 0 ||
             strcmp(name, "Minimum/Maximum length of steps") == 0 ||
             strcmp(name, "Minimum/maximum length of steps") == 0) {
    istringstream ss(value);
    char   c1, c2;
    string str;
    ss >> c1;
    if (c1 != '[') ss.putback(c1);
    ss >> str;
    if (!FromString(str.c_str(), _MinStepLength) || _MinStepLength < .0) return false;
    ss >> str;
    if (!FromString(str.c_str(), _MaxStepLength) || _MaxStepLength < _MinStepLength) return false;
    ss >> c2;
    if (c1 == '[' && c2 != ']') return false;
    return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
irtkParameterList irtkLineSearch::Parameter() const
{
  irtkParameterList params;
  Insert(params, "Maximum no. of line iterations", ToString(_NumberOfIterations));
  if (_MinStepLength == _MaxStepLength) {
    Insert(params, "Length of steps", ToString(_MaxStepLength));
  } else {
    Insert(params, "Minimum length of steps", ToString(_MinStepLength));
    Insert(params, "Maximum length of steps", ToString(_MaxStepLength));
  }
  return params;
}
