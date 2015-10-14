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

#ifndef _IRTKENERGYTERM_H
#define _IRTKENERGYTERM_H

#include <irtkCommon.h>
#include <irtkTransformation.h>
#include <irtkImageAttributes.h>
#include <irtkRegistrationEnergy.h>


/**
 * Base class for one term of an objective function
 *
 * In particular, this is the base class for both the data similarity term
 * and transformation regularization term commonly seen in objective functions
 * used for image/surface registration.
 */
class irtkEnergyTerm : public irtkObservable
{
  irtkAbstractMacro(irtkEnergyTerm);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Name of energy term
  irtkPublicAttributeMacro(string, Name);

  /// Parameter name prefixes recognized by Set
  irtkPublicAttributeMacro(vector<string>, ParameterPrefix);

  /// Weight of energy term
  irtkPublicAttributeMacro(double, Weight);

  /// Transformation with free parameters of energy function
  irtkPublicAggregateMacro(irtkTransformation, Transformation);

  /// Whether to divide energy term by its initial value
  irtkPublicAttributeMacro(bool, DivideByInitialValue);

  /// Initial value of energy term
  double _InitialValue;

  // ---------------------------------------------------------------------------
  // Construction/Destruction
protected:

  /// Constructor
  irtkEnergyTerm(const char * = "", double = 1.0);

  /// Copy constructor
  irtkEnergyTerm(const irtkEnergyTerm &);

  /// Assignment operator
  irtkEnergyTerm &operator =(const irtkEnergyTerm &);

public:

  /// Destructor
  virtual ~irtkEnergyTerm();

  /// Default name derived from class name
  ///
  /// Example usage:
  /// \code
  /// irtkEnergyTerm *term = new MyEnergyTerm();
  /// term->Name(term->DefaultName());
  /// \endcode
  string DefaultName() const;

  // ---------------------------------------------------------------------------
  // Parameters

protected:

  /// Whether this energy term has either an explicit name or default prefix
  bool HasPrefix() const;

  /// Get default energy term name prefix (if any)
  string DefaultPrefix() const;

  /// Get name of parameter with default energy term name prefix
  string ParameterNameWithPrefix(const string &) const;

  /// Get name of parameter with default energy term name prefix
  string ParameterNameWithPrefix(const char *) const;

  /// Get name of parameter without energy term name prefix
  string ParameterNameWithoutPrefix(const char *) const;

  /// Insert parameter into name/value list with energy term name prefix
  template <class T> bool InsertWithPrefix(irtkParameterList &, string, T) const;

public:

  // Import other overloads
  using irtkObservable::Parameter;

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  /// Get parameter key/value as string map
  virtual irtkParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Initialize energy term once input and parameters have been set
  virtual void Initialize();

  /// Update internal state after change of DoFs
  ///
  /// \param[in] gradient Whether to also update internal state for evaluation
  ///                     of energy gradient. If \c false, only the internal state
  ///                     required for the energy evaluation need to be updated.
  virtual void Update(bool gradient = true);

  /// Update energy term after convergence
  virtual bool Upgrade();

  /// Reset initial value of energy term
  void ResetInitialValue();

  /// Returns initial value of energy term
  double InitialValue();

  /// Evaluate energy term
  double Value();

  /// Evaluate gradient of energy term
  ///
  /// \param[in,out] gradient Gradient to which the evaluated gradient of this
  ///                         energy term is added to with its resp. weight.
  /// \param[in]     step     Step length for finite differences.
  void Gradient(double *gradient, double step);

  /// Evaluate and normalize gradient of energy term
  ///
  /// \param[in,out] gradient Gradient to which the evaluated normalized gradient
  ///                         of this energy term is added to with its resp. weight.
  /// \param[in]     step     Step length for finite differences.
  void NormalizedGradient(double *gradient, double step);

  /// Adjust step length range
  ///
  /// \param[in]     gradient Gradient of objective function.
  /// \param[in,out] min      Minimum step length.
  /// \param[in,out] max      Maximum step length.
  virtual void GradientStep(const double *gradient, double &min, double &max) const;

protected:

  /// Evaluate unweighted energy term
  virtual double Evaluate() = 0;

  /// Evaluate and add gradient of energy term
  ///
  /// \param[in,out] gradient Gradient to which the computed gradient of the
  ///                         energy term should be added to.
  /// \param[in]     step     Step length for finite differences.
  /// \param[in]     weight   Weight to use when adding the gradient.
  virtual void EvaluateGradient(double *gradient, double step, double weight) = 0;

  // ---------------------------------------------------------------------------
  // Debugging

public:

  /// Return unweighted and unnormalized raw energy term value
  /// \remarks Use for progress reporting only.
  virtual double RawValue(double) const;

  /// Print debug information
  virtual void Print(irtkIndent = 0) const;

  /// Prefix to be used for debug output files
  string Prefix(const char * = NULL) const;

  /// Write input of data fidelity term
  virtual void WriteDataSets(const char *, const char *, bool = true) const;

  /// Write gradient of data fidelity term w.r.t each transformed input
  virtual void WriteGradient(const char *, const char *) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
template <class T>
bool irtkEnergyTerm::InsertWithPrefix(irtkParameterList &params, string name, T value) const
{
  if (name.empty() || !HasPrefix()) return false;
  Insert(params, ParameterNameWithPrefix(name), value);
  return true;
}


#endif
