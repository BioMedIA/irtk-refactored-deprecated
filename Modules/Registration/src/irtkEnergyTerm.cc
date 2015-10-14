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

#include <irtkEnergyTerm.h>


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkEnergyTerm::irtkEnergyTerm(const char *name, double weight)
:
  _Name                (name),
  _Weight              (weight),
  _Transformation      (NULL),
  _DivideByInitialValue(false),
  _InitialValue        (numeric_limits<double>::quiet_NaN())
{
}

// -----------------------------------------------------------------------------
irtkEnergyTerm::irtkEnergyTerm(const irtkEnergyTerm &other)
:
  irtkObservable(other),
  _Name                (other._Name),
  _Weight              (other._Weight),
  _Transformation      (other._Transformation),
  _DivideByInitialValue(other._DivideByInitialValue),
  _InitialValue        (other._InitialValue)
{
}

// -----------------------------------------------------------------------------
irtkEnergyTerm &irtkEnergyTerm::operator =(const irtkEnergyTerm &other)
{
  _Name                 = other._Name;
  _Weight               = other._Weight;
  _Transformation       = other._Transformation;
  _DivideByInitialValue = other._DivideByInitialValue;
  _InitialValue         = other._InitialValue;
  return *this;
}

// -----------------------------------------------------------------------------
irtkEnergyTerm::~irtkEnergyTerm()
{
}

// -----------------------------------------------------------------------------
string irtkEnergyTerm::DefaultName() const
{
  string prefix(this->NameOfClass());
  if (prefix.compare(0, 4, "irtk") == 0) prefix.erase(0, 4);
  for (string::size_type i = 1; i < prefix.length(); ++i) {
    if ('A' <= prefix[i] && prefix[i] <= 'Z') {
      prefix[i] = 'a' + (prefix[i] - 'A');
      prefix.insert(i, " ");
      ++i;
    }
  }
  return prefix;
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkEnergyTerm::HasPrefix() const
{
  return !_Name.empty() || !_ParameterPrefix.empty();
}

// -----------------------------------------------------------------------------
string irtkEnergyTerm::DefaultPrefix() const
{
  if (!_Name.empty()) return _Name + " ";
  return _ParameterPrefix.empty() ? string() : _ParameterPrefix.front();
}

// -----------------------------------------------------------------------------
string irtkEnergyTerm::ParameterNameWithPrefix(const string &str) const
{
  string name = DefaultPrefix();
  if (name.empty()) return str;
  name += ::tolower(str[0]);
  name += str.substr(1);
  return name;
}

// -----------------------------------------------------------------------------
string irtkEnergyTerm::ParameterNameWithPrefix(const char *str) const
{
  string name = DefaultPrefix();
  if (name.empty()) return str;
  name += ::tolower(str[0]);
  name += (str + 1);
  return name;
}

// -----------------------------------------------------------------------------
string irtkEnergyTerm::ParameterNameWithoutPrefix(const char *str) const
{
  string         param_name;
  vector<string> prefixes = _ParameterPrefix;
  if (!_Name.empty()) prefixes.push_back(_Name + " ");

  const size_t len = strlen(str);
  vector<string>::const_iterator prefix;
  for (prefix = prefixes.begin(); prefix != prefixes.end(); ++prefix) {
    if (len > prefix->length() &&
        strncmp(str, prefix->c_str(), prefix->length()) == 0 &&
        *(str + prefix->length()) != '\0') {
      param_name = str + prefix->length();
      if (islower(param_name[0])) param_name[0] = toupper(param_name[0]);
      break;
    }
  }

  return param_name;
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkEnergyTerm::Set(const char *param, const char *value)
{
  const string name = ParameterNameWithoutPrefix(param);

  if (strcmp(param, "Divide energy terms by initial value") == 0 ||
      name == "Relative to initial value") {
    return FromString(value, _DivideByInitialValue);
  }
  if (name == "Weight") {
    double weight;
    if (!FromString(value, weight)) return false;
    // Keep sign of weight by default, as it depends on whether the energy term
    // is being minimized or maximized; only change the relative weighting of
    // the term; a negative value can be used to flip the sign
    _Weight = copysign(1.0, _Weight) * weight;
    return true;
  }
  if (name == "Weight (signed)") {
    return FromString(value, _Weight);
  }

  return false;
}

// -----------------------------------------------------------------------------
irtkParameterList irtkEnergyTerm::Parameter() const
{
  irtkParameterList params;
  InsertWithPrefix(params, "Weight (signed)",           _Weight);
  InsertWithPrefix(params, "Relative to initial value", _DivideByInitialValue);
  return params;
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void irtkEnergyTerm::Initialize()
{
}

// -----------------------------------------------------------------------------
void irtkEnergyTerm::Update(bool)
{
}

// -----------------------------------------------------------------------------
bool irtkEnergyTerm::Upgrade()
{
  return false;
}

// -----------------------------------------------------------------------------
void irtkEnergyTerm::ResetInitialValue()
{
  _InitialValue = numeric_limits<double>::quiet_NaN();
}

// -----------------------------------------------------------------------------
double irtkEnergyTerm::InitialValue()
{
  if (IsNaN(_InitialValue)) {
    _InitialValue = ((_Weight != .0) ? this->Evaluate() : .0);
  }
  return _InitialValue;
}

// -----------------------------------------------------------------------------
double irtkEnergyTerm::Value()
{
  double value = .0;
  if (_Weight != .0) value = this->Evaluate();
  if (IsNaN(_InitialValue)) _InitialValue = value;
  if (_DivideByInitialValue && _InitialValue != .0) value /= fabs(_InitialValue);
  return _Weight * value;
}

// -----------------------------------------------------------------------------
double irtkEnergyTerm::RawValue(double value) const
{
  if (_Weight != .0) value /= _Weight;
  if (_DivideByInitialValue) {
    double abs_init_value = fabs(_InitialValue);
    if (abs_init_value > 0) value *= abs_init_value;
  }
  return value;
}

// -----------------------------------------------------------------------------
void irtkEnergyTerm::Gradient(double *gradient, double step)
{
  double weight = _Weight;
  if (_DivideByInitialValue) {
    if (IsNaN(_InitialValue)) this->InitialValue();
    if (_InitialValue != .0) weight /= fabs(_InitialValue);
  }
  if (weight != .0) this->EvaluateGradient(gradient, step, weight);
}

// -----------------------------------------------------------------------------
void irtkEnergyTerm::NormalizedGradient(double *gradient, double step)
{
  if (_Weight != .0) {
    const int ndofs = _Transformation->NumberOfDOFs();
    double *grad = CAllocate<double>(ndofs);
    this->EvaluateGradient(grad, step, 1.0);
    double norm = _Transformation->DOFGradientNorm(grad);
    if (norm > .0) norm = _Weight / norm;
    for (int dof = 0; dof < ndofs; ++dof) gradient[dof] += norm * grad[dof];
    Deallocate(grad);
  }
}

// -----------------------------------------------------------------------------
void irtkEnergyTerm::GradientStep(const double *, double &, double &) const
{
  // By default, step length range chosen by user/optimizer
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
void irtkEnergyTerm::Print(irtkIndent indent) const
{
  cout << indent << "Name:    " << _Name   << endl;
  cout << indent << "Weight:  " << _Weight << endl;
  cout << indent << "Initial: " << _InitialValue << endl;
}

// -----------------------------------------------------------------------------
string irtkEnergyTerm::Prefix(const char *prefix) const
{
  if (_Name.empty()) return string(prefix);
  string name(_Name);
  transform(name.begin(), name.end(), name.begin(), ::tolower);
  replace  (name.begin(), name.end(), ' ',  '_');
  replace  (name.begin(), name.end(), '\t', '_');
  if (prefix) name.insert(0, prefix);
  name += "_";
  return name;
}

// -----------------------------------------------------------------------------
void irtkEnergyTerm::WriteDataSets(const char *, const char *, bool) const
{
}

// -----------------------------------------------------------------------------
void irtkEnergyTerm::WriteGradient(const char *, const char *) const
{
}
