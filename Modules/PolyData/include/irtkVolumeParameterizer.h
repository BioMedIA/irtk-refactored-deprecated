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

#ifndef _IRTKVOLUMEPARAMETERIZER_H
#define _IRTKVOLUMEPARAMETERIZER_H

#include <irtkObject.h>
#include <irtkPoint.h>

#include <vtkSmartPointer.h>
#include <vtkPointSet.h>
#include <vtkPolyData.h>
#include <vtkAbstractCellLocator.h>
#include <vtkGenericCell.h>


namespace irtk { namespace polydata {


/**
 * Base class of filters which re-parameterize the interior of a piecewise linear complex (PLC)
 */
class irtkVolumeParameterizer : public irtkObject
{
  irtkAbstractMacro(irtkVolumeParameterizer);

  // ---------------------------------------------------------------------------
  // Types

public:

  /// Enumeration of available volume parameterization methods
  enum MapType
  {
    DefaultMap,  ///< Unknown/default method
    HarmonicMap, ///< Volumetric harmonic map
    ACAPMap      ///< As-conformal-as-possible volumetric map
  };

  // ---------------------------------------------------------------------------
  // Attributes

protected:

  /// Input point set (i.e., surface mesh or tetrahedral mesh)
  irtkPublicAttributeMacro(vtkSmartPointer<vtkPointSet>, Input);

  /// Name of input/output point data array (case insensitive)
  irtkPublicAttributeMacro(string, CoordsName);

  /// Output coordinates of tetrahedral volume mesh points
  irtkReadOnlyAttributeMacro(vtkSmartPointer<vtkDataArray>, Coords);

  /// Tetrahedralized volume mesh
  irtkReadOnlyAttributeMacro(vtkSmartPointer<vtkPointSet>, Volume);

  /// Surface of tetrahedral volume mesh
  irtkReadOnlyAttributeMacro(vtkSmartPointer<vtkPolyData>, Surface);

  /// Locates tetrahedron which is closest to a given input point
  irtkAttributeMacro(vtkSmartPointer<vtkAbstractCellLocator>, Locator);

  /// Boolean array indicating whether a point is on the boundary of the volume
  irtkReadOnlyAttributeMacro(vector<bool>, IsBoundaryPoint);

  /// Output tetrahedral mesh with modified point coordinates
  irtkReadOnlyAttributeMacro(vtkSmartPointer<vtkPointSet>, Output);

  /// Number of points
  irtkReadOnlyAttributeMacro(int, NumberOfPoints);

  /// Number of boundary points
  irtkReadOnlyAttributeMacro(int, NumberOfBoundaryPoints);

  /// Number of interior points
  irtkReadOnlyAttributeMacro(int, NumberOfInteriorPoints);

  /// Copy attributes of this class from another instance
  void Copy(const irtkVolumeParameterizer &);

public:

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  irtkVolumeParameterizer();

  /// Copy constructor
  irtkVolumeParameterizer(const irtkVolumeParameterizer &);

  /// Assignment operator
  irtkVolumeParameterizer &operator =(const irtkVolumeParameterizer &);

  /// Destructor
  virtual ~irtkVolumeParameterizer();

  /// Construct new volume parameterizer
  static irtkVolumeParameterizer *New(MapType = DefaultMap);

  /// Type of volumetric map produced by this parameterizer
  virtual MapType TypeOfMap() const = 0;

  // ---------------------------------------------------------------------------
  // Execution

public:

  /// Re-parameterize interior of input data set
  void Run();

protected:

  /// Initialize filter after input and parameters are set
  virtual void Initialize();

  /// Parameterize interior of input data set
  virtual void Parameterize() = 0;

  /// Finalize filter execution
  virtual void Finalize();

  // ---------------------------------------------------------------------------
  // Parameterization

public:

  /// Map point from input domain to output domain
  bool MapPoint(double [3]) const;

  /// Map point from input domain to output domain
  bool MapPoint(const double [3], double [3]) const;

  /// Map point from input domain to output domain
  bool MapPoint(double &, double &, double &) const;

  /// Map point from input domain to output domain
  bool MapPoint(irtkPoint &) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline bool irtkVolumeParameterizer::MapPoint(double p[3]) const
{
  vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
  vtkIdType cellId;
  double pcoords[3], weight[4];
  cellId = _Locator->FindCell(p, 1e-9, cell, pcoords, weight);
  if (cellId == -1) return false;
  vtkPoints * const points = cell->GetPoints();
  p[0]  = weight[0] * points->GetPoint(0)[0];
  p[1]  = weight[0] * points->GetPoint(0)[1];
  p[2]  = weight[0] * points->GetPoint(0)[2];
  p[0] += weight[1] * points->GetPoint(1)[0];
  p[1] += weight[1] * points->GetPoint(1)[1];
  p[2] += weight[1] * points->GetPoint(1)[2];
  p[0] += weight[2] * points->GetPoint(2)[0];
  p[1] += weight[2] * points->GetPoint(2)[1];
  p[2] += weight[2] * points->GetPoint(2)[2];
  p[0] += weight[3] * points->GetPoint(3)[0];
  p[1] += weight[3] * points->GetPoint(3)[1];
  p[2] += weight[3] * points->GetPoint(3)[2];
  return true;
}

// -----------------------------------------------------------------------------
inline bool irtkVolumeParameterizer::MapPoint(const double a[3], double b[3]) const
{
  b[0] = a[0], b[1] = a[1], b[2] = a[2];
  return MapPoint(b);
}

// -----------------------------------------------------------------------------
inline bool irtkVolumeParameterizer::MapPoint(double &x, double &y, double &z) const
{
  double p[3] = {x, y, z};
  return MapPoint(p);
}

// -----------------------------------------------------------------------------
inline bool irtkVolumeParameterizer::MapPoint(irtkPoint &p) const
{
  return MapPoint(p._x, p._y, p._z);
}


} } // namespace irtk::polydata


#endif
