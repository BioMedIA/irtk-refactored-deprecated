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

#ifndef _IRTKCONJUGATEGRADIENTDESCENT_H
#define _IRTKCONJUGATEGRADIENTDESCENT_H

#include <irtkGradientDescent.h>


/**
 * Minimizes objective function using conjugate gradient descent
 */
class irtkConjugateGradientDescent : public irtkGradientDescent
{
  irtkOptimizerMacro(irtkConjugateGradientDescent, ConjugateGradientDescent);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Enable/disable conjugation of gradient
  ///
  /// Can be used to switch to steepest gradient descent once the conjugate
  /// gradient descent is converged. Call Run after setting ConjugateGradientOff.
  irtkPublicAttributeMacro(bool, UseConjugateGradient);

protected:

  double *_g;
  double *_h;

  /// Copy attributes of this class from another instance
  void Copy(const irtkConjugateGradientDescent &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  irtkConjugateGradientDescent(irtkObjectiveFunction * = NULL);

  /// Copy constructor
  irtkConjugateGradientDescent(const irtkConjugateGradientDescent &);

  /// Assignment operator
  irtkConjugateGradientDescent &operator =(const irtkConjugateGradientDescent &);

  /// Destructor
  virtual ~irtkConjugateGradientDescent();

  // ---------------------------------------------------------------------------
  // Optimization

  /// Reset conjugate gradient to objective function gradient
  void ResetConjugateGradient();

  /// Enable conjugation of objective function gradient
  void ConjugateGradientOn();

  /// Disable conjugation of objective function gradient
  void ConjugateGradientOff();

protected:

  /// Initialize gradient descent
  virtual void Initialize();

  /// Finalize gradient descent
  virtual void Finalize();

  /// Compute gradient of objective function and make it conjugate
  virtual double Gradient(double *, double = .0, bool * = NULL);

  /// Compute conjugate gradient
  void ConjugateGradient(double *);

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline void irtkConjugateGradientDescent::ResetConjugateGradient()
{
  // Set g[0] to NaN to indicate that vectors _g and _h are uninitialized
  if (_g) _g[0] = numeric_limits<double>::quiet_NaN();
}

// -----------------------------------------------------------------------------
inline void irtkConjugateGradientDescent::ConjugateGradientOn()
{
  this->UseConjugateGradient(true);
}

// -----------------------------------------------------------------------------
inline void irtkConjugateGradientDescent::ConjugateGradientOff()
{
  this->UseConjugateGradient(false);
}


#endif
