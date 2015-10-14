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

#ifndef _IRTKPOLYHEDRON_H
#define _IRTKPOLYHEDRON_H

#include <irtkObject.h>
#include <irtkPoint.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>


namespace irtk { namespace polydata {


/**
 * Utility functions for dealing with VTK polydata representing a polyhedron.
 *
 * \note The point-in-polyhedron test based on the winding number method and
 * the polyhedron volume computation are a C++ implementation of the excellent
 * polyhedron.py Python module of Mark Dickinson found on GitHub at
 * https://github.com/mdickinson/polyhedron/blob/cd7361bcee8cbd9ef2e0aac58a7ce59bf9a52c4f/polyhedron.py
 */
class irtkPolyhedron : public irtkObject
{
  irtkObjectMacro(irtkPolyhedron);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Dataset defining the geometry and topology of this polyhedron
  irtkPublicAttributeMacro(vtkSmartPointer<vtkPolyData>, DataSet);

public:

  // ---------------------------------------------------------------------------
  // Construction/destruction

  /// Constructor
  irtkPolyhedron(vtkPolyData * = NULL);

  /// Copy constructor
  irtkPolyhedron(const irtkPolyhedron &);

  /// Assignment operator
  irtkPolyhedron &operator =(const irtkPolyhedron &);

  /// Destructor
  virtual ~irtkPolyhedron();

  // ---------------------------------------------------------------------------
  // Geometry / Topology

  /// Number of vertices
  int NumberOfPoints() const;

  /// Get vertex position
  void GetPoint(int, double &, double &, double &) const;

  /// Get vertex position
  void GetPoint(int, double [3]) const;

  /// Get vertex position
  irtkPoint GetPoint(int) const;

  // ---------------------------------------------------------------------------
  // Properties

  /// Calculate volume enclosed by polyhedron
  double Volume() const;

  /// Calculate volume enclosed by polyhedron
  static double Volume(vtkPolyData *);

  /// Compute winding number of polyhedron around given point
  int WindingNumber(double, double, double) const;

  /// Compute winding number of polyhedron around given point
  int WindingNumber(double [3]) const;

  /// Compute winding number of polyhedron around given point
  int WindingNumber(const irtkPoint &) const;

  /// Compute winding number of polyhedron around given point
  static int WindingNumber(vtkPolyData *, double, double, double);

  /// Compute winding number of polyhedron around given point
  static int WindingNumber(vtkPolyData *, double [3]);

  /// Compute winding number of polyhedron around given point
  static int WindingNumber(vtkPolyData *, const irtkPoint &);

  // ---------------------------------------------------------------------------
  // Point-in-polyhedron test

  /// Test whether point is inside the polyhedron
  bool IsInside(double, double, double) const;

  /// Test whether point is inside the polyhedron
  bool IsInside(double [3]) const;

  /// Test whether point is inside the polyhedron
  bool IsInside(const irtkPoint &) const;

  /// Test whether point is inside the polyhedron
  static bool IsInside(vtkPolyData *, double, double, double);

  /// Test whether point is inside the polyhedron
  static bool IsInside(vtkPolyData *, double [3]);

  /// Test whether point is inside the polyhedron
  static bool IsInside(vtkPolyData *, const irtkPoint &);

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Geometry / Topology
// =============================================================================

// -----------------------------------------------------------------------------
inline int irtkPolyhedron::NumberOfPoints() const
{
  return static_cast<int>(_DataSet->GetNumberOfPoints());
}

// -----------------------------------------------------------------------------
inline void irtkPolyhedron::GetPoint(int i, double p[3]) const
{
  _DataSet->GetPoint(i, p);
}

// -----------------------------------------------------------------------------
inline void irtkPolyhedron::GetPoint(int i, double &x, double &y, double &z) const
{
  double p[3];
  _DataSet->GetPoint(i, p);
  x = p[0], y = p[1], z = p[2];
}

// -----------------------------------------------------------------------------
inline irtkPoint irtkPolyhedron::GetPoint(int i) const
{
  irtkPoint p;
  GetPoint(i, p._x, p._y, p._z);
  return p;
}

// =============================================================================
// Properties
// =============================================================================

// -----------------------------------------------------------------------------
inline double irtkPolyhedron::Volume() const
{
  return Volume(_DataSet);
}

// -----------------------------------------------------------------------------
inline int irtkPolyhedron::WindingNumber(double x, double y, double z) const
{
  return WindingNumber(_DataSet, x, y, z);
}

// -----------------------------------------------------------------------------
inline int irtkPolyhedron::WindingNumber(double p[3]) const
{
  return WindingNumber(p[0], p[1], p[2]);
}

// -----------------------------------------------------------------------------
inline int irtkPolyhedron::WindingNumber(const irtkPoint &p) const
{
  return WindingNumber(p._x, p._y, p._z);
}

// -----------------------------------------------------------------------------
inline int irtkPolyhedron::WindingNumber(vtkPolyData *polydata, double p[3])
{
  return WindingNumber(polydata, p[0], p[1], p[2]);
}

// -----------------------------------------------------------------------------
inline int irtkPolyhedron::WindingNumber(vtkPolyData *polydata, const irtkPoint &p)
{
  return WindingNumber(polydata, p._x, p._y, p._z);
}

// =============================================================================
// Point-in-polyhedron test
// =============================================================================

// -----------------------------------------------------------------------------
inline bool irtkPolyhedron::IsInside(vtkPolyData *polydata, double x, double y, double z)
{
  return WindingNumber(polydata, x, y, z) != 0;
}

// -----------------------------------------------------------------------------
inline bool irtkPolyhedron::IsInside(double x, double y, double z) const
{
  return IsInside(_DataSet, x, y, z);
}

// -----------------------------------------------------------------------------
inline bool irtkPolyhedron::IsInside(double p[3]) const
{
  return IsInside(p[0], p[1], p[2]);
}

// -----------------------------------------------------------------------------
inline bool irtkPolyhedron::IsInside(const irtkPoint &p) const
{
  return IsInside(p._x, p._y, p._z);
}

// -----------------------------------------------------------------------------
inline bool irtkPolyhedron::IsInside(vtkPolyData *polydata, double p[3])
{
  return IsInside(polydata, p[0], p[1], p[2]);
}

// -----------------------------------------------------------------------------
inline bool irtkPolyhedron::IsInside(vtkPolyData *polydata, const irtkPoint &p)
{
  return IsInside(polydata, p._x, p._y, p._z);
}


} } // namespace irtk::polydata

#endif
