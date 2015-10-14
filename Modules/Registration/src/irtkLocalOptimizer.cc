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

#include <irtkLocalOptimizer.h>


// =============================================================================
// Factory
// =============================================================================

//#include <irtkDownhillDescent.h>
#include <irtkGradientDescent.h>
//#include <irtkConstrainedGradientDescent.h>
//#include <irtkSteepestGradientDescent.h>
#include <irtkConjugateGradientDescent.h>
#ifdef HAS_LBFGS
#include <irtkLimitedMemoryBFGSDescent.h>
#endif
#ifdef HAS_VTK
#include <irtkEulerMethod.h>
#include <irtkEulerMethodWithMomentum.h>
#endif

// -----------------------------------------------------------------------------
irtkLocalOptimizer *irtkLocalOptimizer::New(irtkOptimizationMethod m, irtkObjectiveFunction *f)
{
  switch (m) {
    //case DownhillDescent:            return new irtkDownhillDescent(f);
    case GradientDescent:            return new irtkGradientDescent(f);
    //case GradientDescentConstrained: return new irtkConstrainedGradientDescent(f);
    //case SteepestGradientDescent:    return new irtkSteepestGradientDescent(f);
    case ConjugateGradientDescent:   return new irtkConjugateGradientDescent(f);
#ifdef HAS_LBFGS
    case LBFGS:                      return new irtkLimitedMemoryBFGSDescent(f);
#endif
#ifdef HAS_VTK
    case EulerMethod:                return new irtkEulerMethod(f);
    case EulerMethodWithMomentum:    return new irtkEulerMethodWithMomentum(f);
#endif
    default:
      cerr << "irtkLocalOptimizer::New: Unknown optimization method: " << m << endl;
      exit(1);
  }
}

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkLocalOptimizer::irtkLocalOptimizer(irtkObjectiveFunction *f)
:
  irtkObservable(),
  _Function(f),
  _Epsilon (1e-4),
  _Delta   (1e-12)
{
}

// -----------------------------------------------------------------------------
irtkLocalOptimizer::irtkLocalOptimizer(const irtkLocalOptimizer &other)
:
  irtkObservable(other),
  _Function(other._Function),
  _Epsilon (other._Epsilon),
  _Delta   (other._Delta)
{
}

// -----------------------------------------------------------------------------
irtkLocalOptimizer &irtkLocalOptimizer::operator =(const irtkLocalOptimizer &other)
{
  irtkObservable::operator =(other);
  _Function = other._Function;
  _Epsilon  = other._Epsilon;
  _Delta    = other._Delta;
  return *this;
}

// -----------------------------------------------------------------------------
irtkLocalOptimizer::~irtkLocalOptimizer()
{
}

// -----------------------------------------------------------------------------
void irtkLocalOptimizer::Initialize()
{
  if (!_Function) {
    cerr << this->NameOfClass() << "::Initialize: Objective function not set" << endl;
    exit(1);
  }
  if (_Function->NumberOfDOFs() <= 0) {
    cerr << this->NameOfClass() << "::Initialize: Objective function has no free parameters (DoFs)" << endl;
    exit(1);
  }
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkLocalOptimizer::Set(const char *name, const char *value)
{
  if (strcmp(name, "Epsilon") == 0) {
    return FromString(value, _Epsilon);
  } else if (strcmp(name, "Delta") == 0) {
    return FromString(value, _Delta) && _Delta >= .0;
  }
  return false;
}

// -----------------------------------------------------------------------------
irtkParameterList irtkLocalOptimizer::Parameter() const
{
  irtkParameterList params;
  Insert(params, "Epsilon", ToString(_Epsilon));
  Insert(params, "Delta",   ToString(_Delta));
  return params;
}
