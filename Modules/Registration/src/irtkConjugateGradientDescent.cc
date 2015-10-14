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

#include <irtkConjugateGradientDescent.h>


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkConjugateGradientDescent::irtkConjugateGradientDescent(irtkObjectiveFunction *f)
:
  irtkGradientDescent(f), _UseConjugateGradient(true), _g(NULL), _h(NULL)
{
}

// -----------------------------------------------------------------------------
void irtkConjugateGradientDescent::Copy(const irtkConjugateGradientDescent &other)
{
  _UseConjugateGradient = other._UseConjugateGradient;

  Deallocate(_g);
  if (other._g && _Function) {
    Allocate(_g, _Function->NumberOfDOFs());
    memcpy(_g, other._g, _Function->NumberOfDOFs() * sizeof(double));
  }

  Deallocate(_h);
  if (other._h && _Function) {
    Allocate(_h, _Function->NumberOfDOFs());
    memcpy(_h, other._h, _Function->NumberOfDOFs() * sizeof(double));
  }
}

// -----------------------------------------------------------------------------
irtkConjugateGradientDescent::irtkConjugateGradientDescent(const irtkConjugateGradientDescent &other)
:
  irtkGradientDescent(other),
  _g(NULL), _h(NULL)
{
  Copy(other);
}

// -----------------------------------------------------------------------------
irtkConjugateGradientDescent &irtkConjugateGradientDescent::operator =(const irtkConjugateGradientDescent &other)
{
  if (this != &other) {
    irtkGradientDescent::operator =(other);
    Copy(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
irtkConjugateGradientDescent::~irtkConjugateGradientDescent()
{
  Deallocate(_g);
  Deallocate(_h);
}

// =============================================================================
// Optimization
// =============================================================================

// -----------------------------------------------------------------------------
void irtkConjugateGradientDescent::Initialize()
{
  irtkGradientDescent::Initialize();
  Deallocate(_g); Allocate(_g, _Function->NumberOfDOFs());
  Deallocate(_h); Allocate(_h, _Function->NumberOfDOFs());
  ResetConjugateGradient();
}

// -----------------------------------------------------------------------------
void irtkConjugateGradientDescent::Finalize()
{
  irtkGradientDescent::Finalize();
  Deallocate(_g);
  Deallocate(_h);
}

// -----------------------------------------------------------------------------
double irtkConjugateGradientDescent::Gradient(double *gradient, double step, bool *sgn_chg)
{
  // Compute gradient of objective function
  _Function->Gradient(gradient, step, sgn_chg);

  // Update gradient to be conjugate
  if (_UseConjugateGradient) {
    ConjugateGradient(gradient);
  } else {
    ResetConjugateGradient();
  }

  // Return max norm of gradient
  return _Function->GradientNorm(gradient);
}

// -----------------------------------------------------------------------------
void irtkConjugateGradientDescent::ConjugateGradient(double *gradient)
{
  const int ndofs = _Function->NumberOfDOFs();
  if (IsNaN(_g[0])) {
    for (int i = 0; i < ndofs; ++i) _g[i] = -gradient[i];
    memcpy(_h, _g, ndofs * sizeof(double));
  } else {
    double gg  = .0;
    double dgg = .0;
    for (int i = 0; i < ndofs; ++i) {
      gg  += _g[i] * _g[i];
      dgg += (gradient[i] + _g[i]) * gradient[i];
    }
    double gamma = max(dgg / gg, .0);
    for (int i = 0; i < ndofs; ++i) {
      _g[i] = -gradient[i];
      _h[i] = _g[i] + gamma * _h[i];
      gradient[i] = -_h[i];
    }
  }
}
