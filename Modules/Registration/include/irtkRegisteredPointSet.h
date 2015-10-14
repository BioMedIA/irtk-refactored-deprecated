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

#ifndef _IRTKREGISTEREDPOINTSET_H
#define _IRTKREGISTEREDPOINTSET_H

#include <vtkSmartPointer.h>
#include <vtkPointSet.h>
#include <vtkPolyData.h>

class vtkIdTypeArray;

#include <irtkObject.h>
#include <irtkTransformation.h>
#include <irtkEdgeTable.h>


/**
 * Registered point set
 *
 * A registered point set is aligned with another data set (i.e., an image or
 * another point set) through a given transformation. If no transformation is
 * set (i.e., NULL), the point set is fixed and usually serves as target for
 * the registration. Given the pull-back character of image transformations,
 * the point set which is being transformed is generally the one defined in
 * the space of the target image (i.e., the fixed reference image). We therefore
 * refer to the transformed point set as the "target" instead, whereas the
 * "source" point set is usually not transformed. Eventually, however, both
 * target and source point sets may be transformed, e.g., in case of a
 * symmetric registration.
 *
 * Besides the geometry, i.e., point positions, a point set may also have a
 * topological structure as in case of a tetrahedralized simplicial complex or
 * a boundary surface mesh. Point set distance measures and constraints may be
 * specialized for surface meshes only, such as the surface curvature constraint,
 * for example (irtkCurvatureConstraint). These objective function terms only
 * consider the surface of a given point set which is provided by this class.
 * If the input point set itself is already a vtkPolyData object, the point set
 * and its surface are identical. Otherwise, an instance of this class extracts
 * the surface of the point set using the vtkDataSetSurfaceFilter upon
 * initialization and keeps the point positions and data of the registered point
 * set and its corresponding surface in sync.
 *
 * \note VTK does not make use of the const keyword and none of the member
 *       functions are declared const. Therefore, returning a pointer to a
 *       const VTK object is useless. Thus, this class grants anyone non-const
 *       access to the internal VTK data objects even if these are not to be
 *       modified by the caller.
 */
class irtkRegisteredPointSet : public irtkObject
{
  irtkObjectMacro(irtkRegisteredPointSet);

  // ---------------------------------------------------------------------------
  // Types
public:

  /// Feature scaling function parameters
  struct ScalingFunction
  {
    int    _InputIndex;  ///< Feature point data index in input  dataset
    int    _OutputIndex; ///< Feature point data index in output dataset
    double _Slope;       ///< Slope of linear scaling function
    double _Intercept;   ///< Intercept of linear scaling function
  };

  /// Indices and scaling function parameters of transformed point data
  typedef vector<ScalingFunction> ScalingFunctions;

  /// Adjacency matrix with edge IDs
  typedef irtk::polydata::irtkEdgeTable EdgeTable;

  // ---------------------------------------------------------------------------
  // Attributes
protected:

  /// Untransformed input point set
  irtkPublicAttributeMacro(vtkSmartPointer<vtkPointSet>, InputPointSet);

  /// Untransformed surface of input point set
  irtkReadOnlyAttributeMacro(vtkSmartPointer<vtkPolyData>, InputSurface);

  /// (Cached) Untransformed input points
  /// used by energy terms for gradient calculation
  irtkReadOnlyAttributeMacro(irtkPointSet, InputPoints);

  /// (Cached) Untransformed input surface points
  /// used by energy terms for gradient calculation
  irtkPointSet *_InputSurfacePoints;

  /// Time point of input dataset
  irtkPublicAttributeMacro(double, InputTime);

  /// Time point of (transformed) dataset
  irtkPublicAttributeMacro(double, Time);

  /// Current transformation estimate
  irtkPublicAggregateMacro(const irtkTransformation, Transformation);

  /// Indices of point data to copy and (optionally) rescale/normalize
  irtkPublicAttributeMacro(ScalingFunctions, PointDataToCopy);

  /// Whether to copy all point and cell data
  irtkPublicAttributeMacro(bool, CopyAll);

  /// Transformed output point set
  irtkReadOnlyAttributeMacro(vtkSmartPointer<vtkPointSet>, OutputPointSet);

  /// Transformed output surface
  irtkReadOnlyAttributeMacro(vtkSmartPointer<vtkPolyData>, OutputSurface);

  /// Edge table of point set (computed on demand)
  mutable EdgeTable _EdgeTable;

  /// Edge table of point set surface (computed on demand)
  mutable EdgeTable _SurfaceEdgeTable;

  /// Whether self-update is enabled
  irtkPublicAttributeMacro(bool, SelfUpdate);

  /// Whether normals of output surface need to be recomputed (on demand)
  irtkPublicAttributeMacro(bool, UpdateSurfaceNormals);

  /// Domain on which to evaluate transformation if it requires caching
  /// The obtained deformation field is interpolated linearly.
  irtkPublicAttributeMacro(irtkImageAttributes, Domain);

  /// Externally pre-computed displacements to use
  irtkPublicAggregateMacro(irtkGenericImage<double>, ExternalDisplacement);

  /// Cached displacement field evaluated at each lattice point of _Domain
  irtkComponentMacro(irtkGenericImage<double>, Displacement);

  /// Copy attributes of this class from another instance
  void Copy(const irtkRegisteredPointSet &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  irtkRegisteredPointSet(vtkPointSet * = NULL, const irtkTransformation * = NULL);

  /// Copy constructor
  irtkRegisteredPointSet(const irtkRegisteredPointSet &);

  /// Assignment operator
  irtkRegisteredPointSet &operator =(const irtkRegisteredPointSet &);

  /// Destructor
  ~irtkRegisteredPointSet();

  /// Initialize (transformed) point set after input and parameters are set
  ///
  /// \param[in] deep_copy_points Whether to force a deep copy of the input
  ///                             points. If \c false, the points of the
  ///                             output points is only made when Update
  ///                             transforms the input points.
  /// \param[in] init_edge_tables Whether to initialize the edge tables.
  ///                             If \c false, the edge tables are initialized
  ///                             on demand, i.e., when first accessed.
  void Initialize(bool deep_copy_points = false,
                  bool init_edge_tables = false);

  /// Called when the output points have been modified
  /// (cf. irtkDeformableSurfaceModel without FFD)
  void PointsChanged();

  // ---------------------------------------------------------------------------
  // Input point set

  /// Get number of points
  int NumberOfPoints() const;

  /// Get number of cells
  int NumberOfCells() const;

  /// Get untransformed input point with specified index
  void GetInputPoint(int, double &, double &, double &) const;

  /// Get untransformed input point with specified index
  void GetInputPoint(int, double *) const;

  /// Get untransformed input point with specified index
  void GetInputPoint(int, irtkPoint &) const;

  /// Get untransformed points of input data set
  void GetInputPoints(irtkPointSet &) const;

  // ---------------------------------------------------------------------------
  // Input point set surface

  /// Whether (input) point set is a polygonal surface mesh
  /// Otherwise, the surface of the point set is extracted by Initialize.
  bool IsSurface() const;

  /// Get number of surface points
  int NumberOfSurfacePoints() const;

  /// Get number of surface cells
  int NumberOfSurfaceCells() const;

  /// Get array which stores for each surface point the input point set point ID
  vtkIdTypeArray *OriginalSurfacePointIds() const;

  /// Get array which stores for each surface cell the input point set cell ID
  vtkIdTypeArray *OriginalSurfaceCellIds() const;

  /// Get untransformed input surface point with specified index
  void GetInputSurfacePoint(int, double &, double &, double &) const;

  /// Get untransformed input surface point with specified index
  void GetInputSurfacePoint(int, double *) const;

  /// Get untransformed input surface point with specified index
  void GetInputSurfacePoint(int, irtkPoint &) const;

  /// Get untransformed points of input point set surface
  void GetInputSurfacePoints(irtkPointSet &) const;

  /// Untransformed points of input point set surface
  const irtkPointSet &InputSurfacePoints() const;

  // ---------------------------------------------------------------------------
  // Point set access

  /// Implicit conversion to vtkDataSet pointer
  operator vtkDataSet *() const;

  /// Implicit conversion to vtkPointSet pointer
  operator vtkPointSet *() const;

  /// Get (transformed) point set
  vtkPointSet *PointSet() const;

  /// Get points of point set
  vtkPoints *Points() const;

  /// Get edge table of point set mesh
  ///
  /// \attention Not thread-safe unless \c _EdgeTable is already initialized,
  ///            i.e., when this function was called by main thread after Initialize.
  const EdgeTable *Edges() const;

  /// Get point with specified index
  void GetPoint(int, double &, double &, double &) const;

  /// Get point with specified index
  void GetPoint(int, double *) const;

  /// Get point with specified index
  void GetPoint(int, irtkPoint &) const;

  /// Get points of point set
  void GetPoints(irtkPointSet &) const;

  // ---------------------------------------------------------------------------
  // Point set surface access

  /// Get output surface
  vtkPolyData *Surface() const;

  /// Get points of point set surface
  vtkPoints *SurfacePoints() const;

  /// Get output surface (point) normals
  ///
  /// \attention Not thread-safe unless \c _UpdateSurfaceNormals is \c false,
  ///            i.e., when this function was called by main thread after Update.
  vtkDataArray *SurfaceNormals() const;

  /// Get edge table of point set surface mesh
  ///
  /// \attention Not thread-safe unless \c _SurfaceEdgeTable (or \c _EdgeTable
  ///            if input point set is a surface mesh) is already initialized,
  ///            i.e., when this function was called by main thread after Initialize.
  const EdgeTable *SurfaceEdges() const;

  /// Get point with specified index
  void GetSurfacePoint(int, double &, double &, double &) const;

  /// Get point with specified index
  void GetSurfacePoint(int, double *) const;

  /// Get point with specified index
  void GetSurfacePoint(int, irtkPoint &) const;

  /// Get point set surface points
  void GetSurfacePoints(irtkPointSet &) const;

  // ---------------------------------------------------------------------------
  // Update

  /// Update (transformed) dataset
  ///
  /// This function only updates the output points if the self-update attribute
  /// is enabled and only if a (changing) transformation is set. If the dataset
  /// is not transformed or only mapped by a fixed transformation, this function
  /// does nothing. Use \c force=true to initialize this point set even if it
  /// does not change over the course of the registration.
  ///
  /// \param[in] force Force update in any case.
  void Update(bool force = false);

  // ---------------------------------------------------------------------------
  // Debugging

  /// Default file name extension
  const char *DefaultExtension() const;

  /// Write transformed dataset to file
  void Write(const char *, vtkAbstractArray * = NULL, vtkAbstractArray * = NULL) const;

  /// Write transformed dataset to file
  void Write(const char *, vtkAbstractArray **, int, vtkAbstractArray ** = NULL, int = 0) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Input point set
// =============================================================================

// -----------------------------------------------------------------------------
inline int irtkRegisteredPointSet::NumberOfPoints() const
{
  return const_cast<vtkPointSet *>(_InputPointSet.GetPointer())->GetNumberOfPoints();
}

// -----------------------------------------------------------------------------
inline int irtkRegisteredPointSet::NumberOfCells() const
{
  return const_cast<vtkPointSet *>(_InputPointSet.GetPointer())->GetNumberOfCells();
}

// -----------------------------------------------------------------------------
inline void irtkRegisteredPointSet::GetInputPoint(int i, double &x, double &y, double &z) const
{
  double p[3];
  _InputPointSet->GetPoints()->GetPoint(i, p);
  x = p[0], y = p[1], z = p[2];
}

// -----------------------------------------------------------------------------
inline void irtkRegisteredPointSet::GetInputPoint(int i, double *p) const
{
  _InputPointSet->GetPoints()->GetPoint(i, p);
}

// -----------------------------------------------------------------------------
inline void irtkRegisteredPointSet::GetInputPoint(int i, irtkPoint &pt) const
{
  double p[3];
  _InputPointSet->GetPoints()->GetPoint(i, p);
  pt._x = p[0], pt._y = p[1], pt._z = p[2];
}

// =============================================================================
// Input point set surface
// =============================================================================

// -----------------------------------------------------------------------------
inline bool irtkRegisteredPointSet::IsSurface() const
{
  return vtkPolyData::SafeDownCast(_InputPointSet) != NULL;
}

// -----------------------------------------------------------------------------
inline int irtkRegisteredPointSet::NumberOfSurfacePoints() const
{
  return const_cast<vtkPolyData *>(_InputSurface.GetPointer())->GetNumberOfPoints();
}

// -----------------------------------------------------------------------------
inline int irtkRegisteredPointSet::NumberOfSurfaceCells() const
{
  return const_cast<vtkPolyData *>(_InputSurface.GetPointer())->GetNumberOfCells();
}

// -----------------------------------------------------------------------------
inline void irtkRegisteredPointSet::GetInputSurfacePoint(int i, double &x, double &y, double &z) const
{
  double p[3];
  _InputSurface->GetPoints()->GetPoint(i, p);
  x = p[0], y = p[1], z = p[2];
}

// -----------------------------------------------------------------------------
inline void irtkRegisteredPointSet::GetInputSurfacePoint(int i, double *p) const
{
  _InputSurface->GetPoints()->GetPoint(i, p);
}

// -----------------------------------------------------------------------------
inline void irtkRegisteredPointSet::GetInputSurfacePoint(int i, irtkPoint &pt) const
{
  double p[3];
  _InputSurface->GetPoints()->GetPoint(i, p);
  pt._x = p[0], pt._y = p[1], pt._z = p[2];
}

// -----------------------------------------------------------------------------
inline const irtkPointSet &irtkRegisteredPointSet::InputSurfacePoints() const
{
  return *_InputSurfacePoints;
}

// =============================================================================
// Point set
// =============================================================================

// -----------------------------------------------------------------------------
inline irtkRegisteredPointSet::operator vtkDataSet *() const
{
  return _OutputPointSet.GetPointer();
}

// -----------------------------------------------------------------------------
inline irtkRegisteredPointSet::operator vtkPointSet *() const
{
  return _OutputPointSet.GetPointer();
}

// -----------------------------------------------------------------------------
inline vtkPointSet *irtkRegisteredPointSet::PointSet() const
{
  return _OutputPointSet.GetPointer();
}

// -----------------------------------------------------------------------------
inline vtkPoints *irtkRegisteredPointSet::Points() const
{
  return _OutputPointSet->GetPoints();
}

// -----------------------------------------------------------------------------
inline void irtkRegisteredPointSet::GetPoint(int i, double &x, double &y, double &z) const
{
  double p[3];
  _OutputPointSet->GetPoints()->GetPoint(i, p);
  x = p[0], y = p[1], z = p[2];
}

// -----------------------------------------------------------------------------
inline void irtkRegisteredPointSet::GetPoint(int i, double *p) const
{
  _OutputPointSet->GetPoints()->GetPoint(i, p);
}

// -----------------------------------------------------------------------------
inline void irtkRegisteredPointSet::GetPoint(int i, irtkPoint &pt) const
{
  double p[3];
  _OutputPointSet->GetPoints()->GetPoint(i, p);
  pt._x = p[0], pt._y = p[1], pt._z = p[2];
}

// =============================================================================
// Point set surface
// =============================================================================

// -----------------------------------------------------------------------------
inline vtkPolyData *irtkRegisteredPointSet::Surface() const
{
  return _OutputSurface.GetPointer();
}

// -----------------------------------------------------------------------------
inline vtkPoints *irtkRegisteredPointSet::SurfacePoints() const
{
  return _OutputSurface->GetPoints();
}

// -----------------------------------------------------------------------------
inline void irtkRegisteredPointSet::GetSurfacePoint(int i, double &x, double &y, double &z) const
{
  double p[3];
  _OutputSurface->GetPoints()->GetPoint(i, p);
  x = p[0], y = p[1], z = p[2];
}

// -----------------------------------------------------------------------------
inline void irtkRegisteredPointSet::GetSurfacePoint(int i, double *p) const
{
  _OutputSurface->GetPoints()->GetPoint(i, p);
}

// -----------------------------------------------------------------------------
inline void irtkRegisteredPointSet::GetSurfacePoint(int i, irtkPoint &pt) const
{
  double p[3];
  _OutputSurface->GetPoints()->GetPoint(i, p);
  pt._x = p[0], pt._y = p[1], pt._z = p[2];
}


#endif
