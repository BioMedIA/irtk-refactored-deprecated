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

#ifndef _IRTKREGISTEREDSURFACE_H
#define _IRTKREGISTEREDSURFACE_H

#include <irtkRegisteredPointSet.h>

#include <vtkPolyData.h>
#include <vtkCellArray.h>


/**
 * Registered point set/boundary surface
 */
class irtkRegisteredSurface : public irtkRegisteredPointSet
{
  irtkObjectMacro(irtkRegisteredSurface);

public:

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Constructor
  irtkRegisteredSurface(vtkPolyData * = NULL, const irtkTransformation * = NULL);

  /// Copy constructor
  irtkRegisteredSurface(const irtkRegisteredSurface &);

  /// Assignment operator
  irtkRegisteredSurface &operator =(const irtkRegisteredSurface &);

  /// Destructor
  ~irtkRegisteredSurface();

  /// Initialize dataset after input and parameters are set
  void Initialize();

  // ---------------------------------------------------------------------------
  // Polydata access

  /// Set input surface
  void InputSurface(vtkPolyData *);

  /// Get number of vertices
  int NumberOfVerts() const;

  /// Get number of lines
  int NumberOfLines() const;

  /// Get number of polygons
  int NumberOfPolys() const;

  /// Get number of triangle strips
  int NumberOfStrips() const;

  /// Get (transformed) polydata
  vtkPolyData *PolyData() const;

  /// Get the cell array defining vertices. If there are no vertices, an
  /// empty array will be returned (convenience to simplify traversal).
  vtkCellArray *Verts() const;

  /// Get the cell array defining lines. If there are no lines, an
  /// empty array will be returned (convenience to simplify traversal).
  vtkCellArray *Lines() const;

  /// Get the cell array defining polygons. If there are no polygons, an
  /// empty array will be returned (convenience to simplify traversal).
  vtkCellArray *Polys() const;

  /// Get the cell array defining triangle strips. If there are no triangle strips,
  /// an empty array will be returned (convenience to simplify traversal).
  vtkCellArray *Strips() const;

  /// Get pointer to IDs of points defining a cell
  void GetCellPoints(int, vtkIdType &, const vtkIdType *&) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Polydata access
// =============================================================================

// -----------------------------------------------------------------------------
inline void irtkRegisteredSurface::InputSurface(vtkPolyData *surface)
{
  _InputPointSet = surface;
}

// -----------------------------------------------------------------------------
inline vtkCellArray *irtkRegisteredSurface::Verts() const
{
  return _OutputSurface->GetVerts();
}

// -----------------------------------------------------------------------------
inline vtkCellArray *irtkRegisteredSurface::Lines() const
{
  return _OutputSurface->GetLines();
}

// -----------------------------------------------------------------------------
inline vtkCellArray *irtkRegisteredSurface::Polys() const
{
  return _OutputSurface->GetPolys();
}

// -----------------------------------------------------------------------------
inline vtkCellArray *irtkRegisteredSurface::Strips() const
{
  return _OutputSurface->GetStrips();
}

// -----------------------------------------------------------------------------
inline int irtkRegisteredSurface::NumberOfVerts() const
{
  return static_cast<int>(Verts()->GetNumberOfCells());
}

// -----------------------------------------------------------------------------
inline int irtkRegisteredSurface::NumberOfLines() const
{
  return static_cast<int>(Lines()->GetNumberOfCells());
}

// -----------------------------------------------------------------------------
inline int irtkRegisteredSurface::NumberOfPolys() const
{
  return static_cast<int>(Polys()->GetNumberOfCells());
}

// -----------------------------------------------------------------------------
inline int irtkRegisteredSurface::NumberOfStrips() const
{
  return static_cast<int>(Strips()->GetNumberOfCells());
}

// -----------------------------------------------------------------------------
inline void irtkRegisteredSurface::GetCellPoints(int i, vtkIdType &npts, const vtkIdType *&pts) const
{
  _OutputSurface->GetCellPoints(i, npts, const_cast<vtkIdType *&>(pts));
}


#endif
