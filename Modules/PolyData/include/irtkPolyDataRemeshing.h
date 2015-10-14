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

#ifndef _IRTKPOLYDATAREMESHING_H
#define _IRTKPOLYDATAREMESHING_H

#include <irtkPolyDataFilter.h>
#include <irtkPoint.h>
class irtkTransformation;

#include <vtkPriorityQueue.h>


namespace irtk { namespace polydata {


/**
 * Adaptive local remeshing of triangulated surface mesh
 *
 * Park et al., A non-self-intersecting adaptive deformable surface for
 * complex boundary extraction from volumetric images, 25, 421–440 (2001).
 */
class irtkPolyDataRemeshing : public irtkPolyDataFilter
{
  irtkObjectMacro(irtkPolyDataRemeshing);

  // ---------------------------------------------------------------------------
  // Types

public:

  /// Enumeration of cell order in which melting is performed
  enum Order
  {
    INDEX,         ///< Cell index
    AREA,          ///< Cell area
    SHORTEST_EDGE  ///< Length of shortest edge
  };

protected:

  /// Structure storing information about cell in priority queue
  struct CellInfo {
    vtkIdType cellId;
    double    priority;

    bool operator <(const CellInfo &other) const
    {
      return priority < other.priority;
    }
  };

  /// Structure storing information about edge in priority queue
  struct EdgeInfo {
    vtkIdType ptId1;
    vtkIdType ptId2;
    double    priority;

    bool operator <(const EdgeInfo &other) const
    {
      return priority < other.priority;
    }
  };

  /// Order array of cell information structures used by Melting pass
  typedef vector<CellInfo>          CellQueue;
  typedef CellQueue::const_iterator CellIterator;

  /// Order array of cell information structures used by Melting pass
  typedef vector<EdgeInfo>          EdgeQueue;
  typedef EdgeQueue::const_iterator EdgeIterator;

  // ---------------------------------------------------------------------------
  // Attributes

  /// Optional input transformation used to determine edge length and triangle area
  irtkPublicAggregateMacro(const irtkTransformation, Transformation);

  /// Output point labels
  irtkAttributeMacro(vtkSmartPointer<vtkDataArray>, OutputPointLabels);

  /// Whether to skip triangulation of input mesh before processing it
  /// This may only be skipped if the input mesh is triangulated already
  irtkPublicAttributeMacro(bool, SkipTriangulation);

  /// Minimum angle between edge end point normals to consider the edge as
  /// an important feature edge which is excluded from any melting operation.
  irtkPublicAttributeMacro(double, MinFeatureAngle);
  double _MinFeatureAngleCos; ///< 1 - cos(_MinFeatureAngle)

  /// If edge end point normals make up an angle greater than this maximum
  /// feature angle, the respective edge is subdivided even if the edge is
  /// shorter than the _MaxEdgeLength if both edges resulting from splitting
  /// the edge in half are at least _MinEdgeLength long.
  irtkPublicAttributeMacro(double, MaxFeatureAngle);
  double _MaxFeatureAngleCos; ///< 1 - cos(_MaxFeatureAngle)

  /// Minimum edge length
  irtkPublicAttributeMacro(double, MinEdgeLength);
  double _MinEdgeLengthSquared;

  /// Maximum edge length
  irtkPublicAttributeMacro(double, MaxEdgeLength);
  double _MaxEdgeLengthSquared;

  /// Name of point data array used to adapt the edge length range for each node
  irtkPublicAttributeMacro(string, AdaptiveEdgeLengthArrayName);

  /// Per-node squared minimum edge length
  irtkAttributeMacro(vtkSmartPointer<vtkDataArray>, MinEdgeLengthArray);

  /// Per-node squared maximum edge length
  irtkAttributeMacro(vtkSmartPointer<vtkDataArray>, MaxEdgeLengthArray);

  /// Define in which order to process the cells in the melting pass
  irtkPublicAttributeMacro(Order, MeltingOrder);

  /// Whether to melt nodes with connectivity three by merging the adjacent triangles
  irtkPublicAttributeMacro(bool, MeltNodes);

  /// Whether to melt entire triangles if all three edges are below threshold
  irtkPublicAttributeMacro(bool, MeltTriangles);

  /// Number of melted nodes with connectivity 3
  irtkReadOnlyAttributeMacro(vtkIdType, NumberOfMeltedNodes);

  /// Number of melted edges
  irtkReadOnlyAttributeMacro(vtkIdType, NumberOfMeltedEdges);

  /// Number of melted triangles
  irtkReadOnlyAttributeMacro(vtkIdType, NumberOfMeltedCells);

  /// Number of edge inversions
  irtkReadOnlyAttributeMacro(vtkIdType, NumberOfInversions);

  /// Number of bisections
  irtkReadOnlyAttributeMacro(vtkIdType, NumberOfBisections);

  /// Number of trisections
  irtkReadOnlyAttributeMacro(vtkIdType, NumberOfTrisections);

  /// Number of quadsections
  irtkReadOnlyAttributeMacro(vtkIdType, NumberOfQuadsections);

  /// Copy attributes of this class from another instance
  void Copy(const irtkPolyDataRemeshing &);

public:

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  irtkPolyDataRemeshing();

  /// Copy constructor
  irtkPolyDataRemeshing(const irtkPolyDataRemeshing &);

  /// Assignment operator
  irtkPolyDataRemeshing &operator =(const irtkPolyDataRemeshing &);

  /// Destructor
  virtual ~irtkPolyDataRemeshing();

  // ---------------------------------------------------------------------------
  // Auxiliary functions

protected:

  /// Get (transformed) surface point
  void GetPoint(vtkIdType, double [3]) const;

  /// Get (transformed) surface point normal
  void GetNormal(vtkIdType, double [3]) const;

  /// Calculate area of (transformed) triangle
  double ComputeArea(vtkIdType) const;

  /// Get (non-transformed) edge middle point
  void MiddlePoint(vtkIdType, vtkIdType, double [3]) const;

  /// Get (non-transformed) edge middle point
  irtkPoint MiddlePoint(vtkIdType, vtkIdType) const;

  /// Get node connectivity
  int NodeConnectivity(vtkIdType) const;

  /// Mark cell as deleted and remove all references to it
  void DeleteCell(vtkIdType);

  /// Mark cells as deleted and remove all references to them
  void DeleteCells(vtkIdList *);

  /// Replace cell point, also updating point to cell references
  void ReplaceCellPoint(vtkIdType, vtkIdType, vtkIdType);

  /// Get neighboring triangle sharing an edge with specified triangle
  vtkIdType GetCellEdgeNeighbor(vtkIdType, vtkIdType, vtkIdType) const;

  /// Get other vertices of cell edge neighbors
  void GetCellPointNeighbors(vtkIdType, vtkIdType, vtkIdList *) const;

  /// Check connectivity of edge neighbors
  ///
  /// If other vertex of edge neighbor triangle has connectivity three,
  /// the three triangles adjacent to it are replaced by the union of these
  /// triangle and the melting operation can be performed. Otherwise, if
  /// more than one node with connectivity three is found which is adjacent
  /// to the edge corners, melting the (triangle) edge would cause degeneration
  /// of adjacent triangles.
  ///
  /// @return ID of other vertex of edge neighbor or -1 if not unique or if
  ///         its connectivity prohibits a melting operation.
  vtkIdType GetCellEdgeNeighborPoint(vtkIdType, vtkIdType, vtkIdType, bool = false);

  /// Queue mesh cells by area, smallest first
  CellQueue QueueCellsByArea() const;

  /// Queue mesh cells by length of shortest edge, smallest first
  CellQueue QueueCellsByShortestEdge() const;

  /// Queue mesh edges by length, smallest first
  EdgeQueue QueueEdgesByLength() const;

  /// Interpolate point attributes when subdividing edge
  void InterpolatePointData(vtkIdType, vtkIdType, vtkIdType);

  /// Interpolate point attributes when melting triangle
  void InterpolatePointData(vtkIdType, vtkIdList *, double *);

  /// Get squared minimum length for specified edge
  double SquaredMinEdgeLength(vtkIdType, vtkIdType) const;

  /// Get squared maximum length for specified edge
  double SquaredMaxEdgeLength(vtkIdType, vtkIdType) const;

  // ---------------------------------------------------------------------------
  // Local remeshing operations

  /// Replace triangles adjacent to node with connectivity three by single triangle
  void MeltTriplets();

  /// Collapse single short edge of two adjacent triangles
  void MeltEdge(vtkIdType, vtkIdType, vtkIdType, vtkIdList *);

  /// Collapse entire triangle with more than one too short edges
  void MeltTriangle(vtkIdType, vtkIdList *);

  /// Melt edges or triangle if one or more edge is too short
  void Melting(vtkIdType, vtkIdList *);

  /// Melt edge if it still exists and is too short
  void Melting(vtkIdType, vtkIdList *, vtkIdType, vtkIdList *);

  /// Bisect triangle
  void Bisect(vtkIdType, vtkIdType, vtkIdType, vtkIdType, vtkPoints *, vtkCellArray *);

  /// Trisect triangle
  void Trisect(vtkIdType, vtkIdType, vtkIdType, vtkIdType, vtkPoints *, vtkCellArray *);

  /// Quadsect triangle
  void Quadsect(vtkIdType, vtkIdType, vtkIdType, vtkIdType, vtkPoints *, vtkCellArray *);

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Initialize filter after input and parameters are set
  virtual void Initialize();

  /// Initialize edge length range for each node
  void InitializeEdgeLengthRange(vtkPolyData *);

  /// Perform local remeshing passes
  virtual void Execute();

  /// Perform first pass: melt edges or triangles if one or more edge is too short
  void Melting();

  /// Perform second pass: invert of triangles sharing a long edge
  void Inversion();

  /// Perform third pass: subdivide triangles with remaining long edges
  void Subdivision();

  /// Finalize filter execution
  virtual void Finalize();

  // ---------------------------------------------------------------------------
  // Alternative VTK-like API

public:

  /// Disable/enable triangulation of input mesh before processing it
  irtkOnOffMacro(SkipTriangulation);

  /// Enable/disable melting of nodes with connectivity three
  irtkOnOffMacro(MeltNodes);

  /// Enable/disable melting of triangles when all edges are too short
  irtkOnOffMacro(MeltTriangles);

};


} } // namespace irtk::polydata

#endif
