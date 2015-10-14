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

#include <irtkDataFidelity.h>


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkDataFidelity::irtkDataFidelity(const char *name, double weight)
:
  irtkEnergyTerm(name, weight)
{
  if (version >= irtkVersion(3, 2)) _DivideByInitialValue = true;
}

// -----------------------------------------------------------------------------
irtkDataFidelity::irtkDataFidelity(const irtkDataFidelity &other)
:
  irtkEnergyTerm(other)
{
}

// -----------------------------------------------------------------------------
irtkDataFidelity &irtkDataFidelity::operator =(const irtkDataFidelity &other)
{
  irtkEnergyTerm::operator =(other);
  return *this;
}

// -----------------------------------------------------------------------------
irtkDataFidelity::~irtkDataFidelity()
{
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkDataFidelity::Set(const char *param, const char *value)
{
  if (strcmp(param, "Divide data term by initial value")           == 0 ||
      strcmp(param, "Divide data terms by initial value")          == 0 ||
      strcmp(param, "Divide data fidelity term by initial value")  == 0 ||
      strcmp(param, "Divide data fidelity terms by initial value") == 0) {
    return FromString(value, _DivideByInitialValue);
  }
  return irtkEnergyTerm::Set(param, value);
}
