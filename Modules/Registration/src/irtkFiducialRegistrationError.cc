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

#include <irtkFiducialRegistrationError.h>
#include <irtkFiducialMatch.h>


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkFiducialRegistrationError::irtkFiducialRegistrationError(const char *name, double weight)
:
  irtkPointCorrespondenceDistance(name, weight, new irtkFiducialMatch)
{
  _ParameterPrefix.push_back("FRE ");
  _ParameterPrefix.push_back("Fiducial registration error ");
  _ParameterPrefix.push_back("Fiducial error ");
  _ParameterPrefix.push_back("Landmark registration error ");
  _ParameterPrefix.push_back("Landmark error ");
}

// -----------------------------------------------------------------------------
irtkFiducialRegistrationError::irtkFiducialRegistrationError(const irtkFiducialRegistrationError &other)
:
  irtkPointCorrespondenceDistance(other)
{
}

// -----------------------------------------------------------------------------
irtkFiducialRegistrationError &irtkFiducialRegistrationError::operator =(const irtkFiducialRegistrationError &other)
{
  irtkPointCorrespondenceDistance::operator =(other);
  return *this;
}

// -----------------------------------------------------------------------------
irtkFiducialRegistrationError::~irtkFiducialRegistrationError()
{
}
