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

#include <irtkRadialErrorFunction.h>

#include <irtkDistanceErrorFunction.h>
#include <irtkSquaredErrorFunction.h>
#include <irtkGaussianErrorFunction.h>
#include <irtkCharbonnierErrorFunction.h>
#include <irtkPeronaMalikErrorFunction.h>


// ----------------------------------------------------------------------------
irtkRadialErrorFunction::irtkRadialErrorFunction()
{
}

// ----------------------------------------------------------------------------
irtkRadialErrorFunction::~irtkRadialErrorFunction()
{
}

// ----------------------------------------------------------------------------
irtkRadialErrorFunction *irtkRadialErrorFunction::New(TypeId type)
{
  switch (type) {
    case Distance:    return new irtkDistanceErrorFunction();
    case Squared:     return new irtkSquaredErrorFunction();
    case Gaussian:    return new irtkGaussianErrorFunction();
    case Charbonnier: return new irtkCharbonnierErrorFunction();
    case PeronaMalik: return new irtkPeronaMalikErrorFunction();
    default:
      cerr << "irtkRadialErrorFunction::New: Unknown type = " << type << endl;
      exit(1);
  }
}

// ----------------------------------------------------------------------------
irtkRadialErrorFunction *irtkRadialErrorFunction::New(const char *type_name)
{
  TypeId type = Unknown;
  if (!FromString(type_name, type)) {
    cerr << "irtkRadialErrorFunction::New: Unknown type = " << type_name << endl;
    exit(1);
  }
  return New(type);
}

// -----------------------------------------------------------------------------
template <> bool FromString(const char *str, irtkRadialErrorFunction::TypeId &type)
{
  if (strcmp(str, "Distance") == 0) {
    type = irtkRadialErrorFunction::Distance;
  } else if (strcmp(str, "Square")           == 0 ||
             strcmp(str, "Squared")          == 0 ||
             strcmp(str, "Squared distance") == 0) {
    type = irtkRadialErrorFunction::Squared;
  } else if (strcmp(str, "Gaussian") == 0) {
    type = irtkRadialErrorFunction::Gaussian;
  } else if (strcmp(str, "Charbonnier") == 0) {
    type = irtkRadialErrorFunction::Charbonnier;
  } else if (strcmp(str, "PeronaMalik")      == 0 ||
             strcmp(str, "Perona Malik")     == 0 ||
             strcmp(str, "Perona-Malik")     == 0 ||
             strcmp(str, "Perona and Malik") == 0) {
    type = irtkRadialErrorFunction::PeronaMalik;
  } else {
    type = irtkRadialErrorFunction::Unknown;
  }
  return (type != irtkRadialErrorFunction::Unknown);
}

// -----------------------------------------------------------------------------
template <> string ToString(const irtkRadialErrorFunction::TypeId &type, int w, char c, bool left)
{
  string str;
  switch (type) {
    case irtkRadialErrorFunction::Distance:    str = "Distance";         break;
    case irtkRadialErrorFunction::Squared:     str = "Squared distance"; break;
    case irtkRadialErrorFunction::Gaussian:    str = "Gaussian";         break;
    case irtkRadialErrorFunction::Charbonnier: str = "Charbonnier";      break;
    case irtkRadialErrorFunction::PeronaMalik: str = "Perona-Malik";     break;
    default:                                   str = "Unknown";          break;
  }
  return ToString(str, w, c, left);
}
