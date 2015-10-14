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

#ifndef IRTK_REMESHER_H_
#define IRTK_REMESHER_H_

#include <irtkObject.h>

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPriorityQueue.h>

/**
 * Remeshes a given surface mesh to enforce a given edge length
 *
 * \author Robert Wright
 */
class irtkRemesher : public irtkObject
{
  irtkObjectMacro(irtkRemesher);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Maximum edge length
  irtkPublicAttributeMacro(double, MaxEdgeLength);

  /// Minimum edge length
  irtkPublicAttributeMacro(double, MinEdgeLength);

  // ---------------------------------------------------------------------------
  // Auxiliary functions

private:

  /// Get length of line segment
  double GetEdgeLength(double [], double []);

  /// Get midpoint of line segment
  void Midpoint(double [], double [],double []);

  /// Clean poly data, merging points
  vtkSmartPointer<vtkPolyData> CleanPolyData(vtkPolyData *);

  /// Preprocess the input surface mesh
  vtkSmartPointer<vtkPolyData> ProcessInput(vtkPolyData *);

  /// Construct internal edge data structure
  vtkSmartPointer<vtkPolyData> BuildEdgeStructure(vtkPolyData *);

  /// Compute lengths of all edges
  vtkSmartPointer<vtkDoubleArray> ComputeEdgeLengths(vtkPolyData *);

  /// Initialize priority queue
  vtkSmartPointer<vtkPriorityQueue> QueueEdgesByLength(vtkPolyData *, double, double);

  /// Get IDs of cells adjacent to an edge
  vtkSmartPointer<vtkIdList> GetAdjacentCellIds(int, int, vtkPolyData *);

  /// Delete collapsed edges
  void DeleteEdgesAfterCollapse(int, int, vtkPolyData *);

  /// Create new cell after edge split
  void CreateNewCells(int, int, int, int, int, vtkCellArray *,vtkPolyData *);

  /// Combine cells
  void CombineCells(vtkPolyData *, vtkCellArray *, vtkIntArray *);

  /// Process edges in the order of their priority (i.e., edge length)
  vtkSmartPointer<vtkPolyData> ProcessEdgeQueue(vtkPriorityQueue *, vtkPolyData *, vtkPolyData *, vtkPolyData *);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  irtkRemesher();

  /// Destructor
  ~irtkRemesher();

  // ---------------------------------------------------------------------------
  // Execution

  /// Remesh the given surface mesh
  vtkSmartPointer<vtkPolyData> Remesh(vtkPolyData *, int * = NULL);

  // ---------------------------------------------------------------------------
  // Deprecated

  /// \deprecated Use MinEdgeLength() setter instead.
  void SetMinEdgeLength(double l) { MinEdgeLength(l); }

  /// \deprecated Use MaxEdgeLength() setter instead.
  void SetMaxEdgeLength(double l) { MaxEdgeLength(l); }

};


#endif
