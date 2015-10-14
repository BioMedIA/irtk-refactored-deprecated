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

#include <irtkPolyDataFilter.h>
#include <irtkProfiling.h>

namespace irtk { namespace polydata {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkPolyDataFilter::irtkPolyDataFilter()
:
  _DoublePrecision(false)
{
}

// -----------------------------------------------------------------------------
void irtkPolyDataFilter::Copy(const irtkPolyDataFilter &other)
{
  _Input = other._Input;
  if (other._Output) {
    _Output = vtkSmartPointer<vtkPolyData>::New();
    _Output->DeepCopy(other._Output);
  } else {
    _Output = NULL;
  }
  _DoublePrecision = other._DoublePrecision;
}

// -----------------------------------------------------------------------------
irtkPolyDataFilter::irtkPolyDataFilter(const irtkPolyDataFilter &other)
{
  Copy(other);
}

// -----------------------------------------------------------------------------
irtkPolyDataFilter &irtkPolyDataFilter::operator =(const irtkPolyDataFilter &other)
{
  irtkObject::operator =(other);
  Copy(other);
  return *this;
}

// -----------------------------------------------------------------------------
irtkPolyDataFilter::~irtkPolyDataFilter()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void irtkPolyDataFilter::Run()
{
  IRTK_START_TIMING();
  {
    IRTK_START_TIMING();
    this->Initialize();
    IRTK_DEBUG_TIMING(2, this->NameOfClass() << "::Initialize");
  }
  {
    IRTK_START_TIMING();
    this->Execute();
    IRTK_DEBUG_TIMING(2, this->NameOfClass() << "::Execute");
  }
  {
    IRTK_START_TIMING();
    this->Finalize();
    IRTK_DEBUG_TIMING(2, this->NameOfClass() << "::Finalize");
  }
  IRTK_DEBUG_TIMING(1, this->NameOfClass());
}

// -----------------------------------------------------------------------------
void irtkPolyDataFilter::Initialize()
{
  // Check input
  if (!_Input) {
    cerr << this->NameOfClass() << "::Initialize: Input surface mesh not set!" << endl;
    exit(1);
  }

  // By default, set output to be shallow copy of input
  _Output = _Input->NewInstance();
  _Output->ShallowCopy(_Input);
}

// -----------------------------------------------------------------------------
void irtkPolyDataFilter::Finalize()
{
}


} } // namespace irtk::polydata
