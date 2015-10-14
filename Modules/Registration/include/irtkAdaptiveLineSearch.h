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

#ifndef _IRTKADAPTIVELINESEARCH_H
#define _IRTKADAPTIVELINESEARCH_H

#include <irtkInexactLineSearch.h>


/**
 * Searches sufficiently optimal step length along search direction
 *
 * This local optimizer implements an inexact line search with adaptive step
 * size control, increasing the step size while steps are accepted, and
 * decreasing it when a step did not yield a sufficient improvement.
 */
class irtkAdaptiveLineSearch : public irtkInexactLineSearch
{
  irtkLineSearchMacro(irtkAdaptiveLineSearch, LS_Adaptive);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Increase factor of step size if accepted
  irtkPublicAttributeMacro(double, StepLengthRise);

  /// Decrease factor of step size if rejected
  irtkPublicAttributeMacro(double, StepLengthDrop);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  irtkAdaptiveLineSearch(irtkObjectiveFunction * = NULL);

  /// Copy constructor
  irtkAdaptiveLineSearch(const irtkAdaptiveLineSearch &);

  /// Assignment operator
  irtkAdaptiveLineSearch &operator =(const irtkAdaptiveLineSearch &);

  /// Destructor
  virtual ~irtkAdaptiveLineSearch();

  // ---------------------------------------------------------------------------
  // Parameters
  using irtkLocalOptimizer::Parameter;

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  /// Get parameters as key/value as string map
  virtual irtkParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Optimization

  /// Make optimal step along search direction
  virtual double Run();

};


#endif
