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

#include <irtkMatlab.h>
#include <irtkCxxLib.h>


// -----------------------------------------------------------------------------
irtkMatlab irtkMatlab::_Instance;
bool irtkMatlab::_McrInitialized = false;
bool irtkMatlab::_AppInitialized = false;

// -----------------------------------------------------------------------------
irtkMatlab::irtkMatlab()
{
}

// -----------------------------------------------------------------------------
irtkMatlab::~irtkMatlab()
{
  if (_AppInitialized) mclTerminateApplication();
}

// -----------------------------------------------------------------------------
irtkMatlab &irtkMatlab::Instance()
{
  return _Instance;
}

// -----------------------------------------------------------------------------
void irtkMatlab::Initialize()
{
  if (!_McrInitialized) {
    mclmcrInitialize();
    _McrInitialized = true;
  }
}

// -----------------------------------------------------------------------------
void irtkMatlab::InitializeApplication(const char **options, int count)
{
  if (!_AppInitialized) {
    if (!mclInitializeApplication(options, count)) {
      cerr << "Failed to initialize MATLAB Compiler Runtime application" << endl;
      exit(1);
    }
    _AppInitialized = _McrInitialized = true;
  }
}
