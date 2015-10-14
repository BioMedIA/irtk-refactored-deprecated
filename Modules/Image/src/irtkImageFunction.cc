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

#include <irtkImageFunction.h>


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkImageFunction::irtkImageFunction()
:
  _DefaultValue(.0),
  _DebugFlag   (false),
  _input       (NULL)
{
}

// -----------------------------------------------------------------------------
irtkImageFunction::irtkImageFunction(const irtkImageFunction &other)
:
  irtkObject   (other),
  _DefaultValue(other._DefaultValue),
  _DebugFlag   (other._DebugFlag),
  _input       (other._input)
{
}

// -----------------------------------------------------------------------------
irtkImageFunction::~irtkImageFunction()
{
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
void irtkImageFunction::SetInput(irtkImage *image)
{
  if (image != NULL) {
    _input = image;
  } else {
    cerr << "irtkImageFunction::SetInput: Input is not an image\n";
    exit(1);
  }
}

// -----------------------------------------------------------------------------
void irtkImageFunction::Initialize()
{
  // Print debugging information
  this->Debug("irtkImageFunction::Initialize");

  // Check inputs and outputs
  if (_input == NULL) {
    cerr << this->NameOfClass() << "::Run: Filter has no input" << endl;
    exit(1);
  }
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
void irtkImageFunction::Debug(const char *message)
{
  if (_DebugFlag == true) cout << message << endl;
}
