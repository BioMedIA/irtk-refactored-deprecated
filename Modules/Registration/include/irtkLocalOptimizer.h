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

#ifndef _IRTKLOCALOPTIMIZER_H
#define _IRTKLOCALOPTIMIZER_H

#include <irtkCommon.h>
#include <irtkObjectiveFunction.h>


/**
 * Minimizes/Maximizes objective function without guarantee of global solution
 */
class irtkLocalOptimizer : public irtkObservable
{
  irtkAbstractMacro(irtkLocalOptimizer);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Objective function
  irtkPublicAggregateMacro(irtkObjectiveFunction, Function);

  /// First convergence criterium: Required minimum change of objective function
  irtkPublicAttributeMacro(double, Epsilon);

  /// Second convergence criterium: Required maximum change of DoFs
  irtkPublicAttributeMacro(double, Delta);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
protected:

  /// Constructor
  irtkLocalOptimizer(irtkObjectiveFunction * = NULL);

  /// Copy constructor
  irtkLocalOptimizer(const irtkLocalOptimizer &);

  /// Assignment operator
  irtkLocalOptimizer &operator =(const irtkLocalOptimizer &);

public:

  /// Construct optimizer
  static irtkLocalOptimizer *New(irtkOptimizationMethod  = ConjugateGradientDescent,
                                 irtkObjectiveFunction * = NULL);

  /// Optimization method implemented by this optimizer
  virtual irtkOptimizationMethod OptimizationMethod() const = 0;

  /// Destructor
  virtual ~irtkLocalOptimizer();

  // ---------------------------------------------------------------------------
  // Parameters

  // Import other overloads
  using irtkObservable::Parameter;

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  /// Get parameters as key/value as string map
  virtual irtkParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Optimization

  /// Initialize optimization
  virtual void Initialize();

  /// Run optimization
  /// \returns Value of local minimum (maximum) of objective function
  virtual double Run() = 0;

};

////////////////////////////////////////////////////////////////////////////////
// Auxiliary macros for optimizer implementation
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
#define irtkOptimizerMacro(name, method)                                       \
  irtkObjectMacro(name);                                                       \
public:                                                                        \
  /** Optimization method implemented by this optimizer */                     \
  virtual irtkOptimizationMethod OptimizationMethod() const {return method;}   \
private:                                                                       \
  /* require developer to end macro with a semicolon */                        \
  static void _irtkOptimizerMacro_needs_trailing_semicolon()


#endif
