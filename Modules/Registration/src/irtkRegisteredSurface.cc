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

#include <irtkRegisteredSurface.h>


// -----------------------------------------------------------------------------
irtkRegisteredSurface::irtkRegisteredSurface(vtkPolyData              *data,
                                             const irtkTransformation *transform)
:
  irtkRegisteredPointSet(data, transform)
{
}

// -----------------------------------------------------------------------------
irtkRegisteredSurface::irtkRegisteredSurface(const irtkRegisteredSurface &other)
:
  irtkRegisteredPointSet(other)
{
}

// -----------------------------------------------------------------------------
irtkRegisteredSurface &irtkRegisteredSurface::operator =(const irtkRegisteredSurface &other)
{
  irtkRegisteredPointSet::operator =(other);
  return *this;
}

// -----------------------------------------------------------------------------
irtkRegisteredSurface::~irtkRegisteredSurface()
{
}

// -----------------------------------------------------------------------------
void irtkRegisteredSurface::Initialize()
{
  // Ensure that input dataset is of valid type
  _InputSurface = vtkPolyData::SafeDownCast(_InputPointSet);
  if (_InputSurface == NULL) {
    cerr << "irtkRegisteredSurface::Initialize: Input dataset must be a vtkPolyData" << endl;
    exit(1);
  }

  // Build cells and links
  _InputSurface->BuildLinks();

  // Initialize base class -- makes shallow copy of input surface, incl. links
  irtkRegisteredPointSet::Initialize();
}
