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

#ifndef _IRTKLINESEARCH_H
#define _IRTKLINESEARCH_H

#include <irtkLocalOptimizer.h>


/// Enumeration of available line search strategies
enum irtkLineSearchStrategy
{
  LS_None,     ///< No line search
  // Add new enumeration values below
  LS_Adaptive, ///< Inexact line search with adaptive step length
  LS_linmin,   ///< Numerical recipes  linmin function
  LS_dlinmin,  ///< Numerical recipes dlinmin function
  // Add new enumeration values above
  LS_Last
};


/**
 * Finds optimal step length along given search direction
 */
class irtkLineSearch : public irtkLocalOptimizer
{
  irtkAbstractMacro(irtkLineSearch);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Whether to minimize or maximize the objective value
  irtkPublicAttributeMacro(bool, Descent);

  /// Ascending direction of line search
  irtkPublicAggregateMacro(double, Direction);

  /// Initial objective function value
  irtkPublicAttributeMacro(double, CurrentValue);

  /// Maximum number of search iterations
  irtkPublicAttributeMacro(int, NumberOfIterations);

  /// Minimum length of step
  irtkPublicAttributeMacro(double, MinStepLength);

  /// Maximum length of step
  irtkPublicAttributeMacro(double, MaxStepLength);

  /// Unit of step length, e.g., maximum gradient norm
  irtkPublicAttributeMacro(double, StepLengthUnit);

  /// Initial/final step length
  irtkPublicAttributeMacro(double, StepLength);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
protected:

  /// Constructor
  irtkLineSearch(irtkObjectiveFunction * = NULL);

  /// Copy constructor
  irtkLineSearch(const irtkLineSearch &);

  /// Assignment operator
  irtkLineSearch &operator =(const irtkLineSearch &);

public:

  /// Instantiate line search implementing specified strategy
  static irtkLineSearch *New(irtkLineSearchStrategy &, irtkObjectiveFunction * = NULL);

  /// Destructor
  virtual ~irtkLineSearch();

  /// Line search strategy implemented by this line search
  virtual irtkLineSearchStrategy Strategy() const = 0;

  // ---------------------------------------------------------------------------
  // Parameters

  // Import other overloads
  using irtkLocalOptimizer::Parameter;

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  /// Get parameters as key/value as string map
  virtual irtkParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Optimization

  /// Make optimal step along search direction
  /// \returns New value of objective function or previous if not step successful
  virtual double Run() = 0;

};

////////////////////////////////////////////////////////////////////////////////
// Auxiliary macros for line search implementation
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
#define irtkLineSearchMacro(name, strategy)                                    \
  irtkObjectMacro(name);                                                       \
  public:                                                                      \
    /** Optimization method implemented by this optimizer */                   \
    virtual irtkOptimizationMethod OptimizationMethod() const { return LineSearch; } \
    /** Line search strategy implemented by this line search */                \
    virtual irtkLineSearchStrategy Strategy() const { return strategy; }       \
  private:                                                                     \
  /* require developer to end macro with a semicolon */                        \
  static void _irtkLineSearchMacro_needs_trailing_semicolon()

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline string ToString(const irtkLineSearchStrategy &m)
{
  switch (m) {
    case LS_None:     return "None";
    case LS_Adaptive: return "Adaptive";
    case LS_linmin:   return "linmin";
    case LS_dlinmin:  return "dlinmin";
    default:          return "Unknown";
  }
}

// -----------------------------------------------------------------------------
inline bool FromString(const char *str, irtkLineSearchStrategy &m)
{
  if (strcmp(str, "None") == 0) {
    m = LS_None;
    return true;
  }
  m = static_cast<irtkLineSearchStrategy>(LS_Last - 1);
  while (m != LS_None) {
    if (ToString(m) == str) break;
    m = static_cast<irtkLineSearchStrategy>(m - 1);
  }
  return m != LS_None;
}

// -----------------------------------------------------------------------------
// For use by subclass Run implementation to get current function value
inline double irtkLineSearch::Run()
{
  if (IsNaN(_CurrentValue)) {
    Function()->Update(false);
    _CurrentValue = Function()->Value();
  }
  return _CurrentValue;
}


#endif
