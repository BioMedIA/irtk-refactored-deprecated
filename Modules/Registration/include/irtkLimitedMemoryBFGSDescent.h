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

#ifndef _IRTKLIMITEDMEMORYBFGSDESCENT_H
#define _IRTKLIMITEDMEMORYBFGSDESCENT_H

#include <irtkLocalOptimizer.h>


/**
 * Minimizes objective function using L-BFGS
 */
class irtkLimitedMemoryBFGSDescent : public irtkLocalOptimizer
{
  irtkOptimizerMacro(irtkLimitedMemoryBFGSDescent, LBFGS);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Maximum number of iterations
  irtkPublicAttributeMacro(int, NumberOfIterations);

  /// Minimum length of steps
  irtkPublicAttributeMacro(double, MinStepLength);

  /// Maximum length of steps
  irtkPublicAttributeMacro(double, MaxStepLength);

public:

  /// Current line search progress
  irtkLineSearchStep _CurrentStep;

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  irtkLimitedMemoryBFGSDescent(irtkObjectiveFunction * = NULL);

  /// Copy constructor
  irtkLimitedMemoryBFGSDescent(const irtkLimitedMemoryBFGSDescent &);

  /// Assignment operator
  irtkLimitedMemoryBFGSDescent &operator =(const irtkLimitedMemoryBFGSDescent &);

  /// Destructor
  virtual ~irtkLimitedMemoryBFGSDescent();

  // ---------------------------------------------------------------------------
  // Parameters

  // Import other overloads
  using irtkLocalOptimizer::Parameter;

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  /// Get parameters as key/value as string map
  virtual irtkParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Execution

  /// Optimize objective function using gradient descent
  virtual double Run();

};


#endif
