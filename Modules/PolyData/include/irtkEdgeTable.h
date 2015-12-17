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

#ifndef _IRTKEDGETABLE_H
#define _IRTKEDGETABLE_H

#include <irtkSparseMatrix.h>
#include <irtkParallel.h> // blocked_range

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>


namespace irtk { namespace polydata {


/**
 * Edge table
 *
 * This class represents the adjacency matrix of point set nodes. It provides
 * efficient access to the set of nodes adjacent to a given point. The non-zero
 * entries stored in the sparse matrix are the one-based edge IDs such that
 * sparse matrix entries are non-zero. To efficiently iterate all edges or a
 * subset of these, use the thread-safe irtkEdgeIterator.
 */
class irtkEdgeTable : public irtkGenericSparseMatrix<int>
{
  irtkObjectMacro(irtkEdgeTable);

  /// Number of undirected edges
  irtkReadOnlyAttributeMacro(int, NumberOfEdges);

public:

  /// Construct edge table for given dataset
  irtkEdgeTable(vtkDataSet * = NULL);

  /// Copy constructor
  irtkEdgeTable(const irtkEdgeTable &);

  /// Assignment operator
  irtkEdgeTable &operator =(const irtkEdgeTable &);

  /// Destructor
  virtual ~irtkEdgeTable();

  /// Initialize edge table from given dataset
  void Initialize(vtkDataSet *);

  /// Number of nodes
  int NumberOfPoints() const;

  /// Determine whether two nodes are connected by an edge
  bool IsEdge(int, int) const;

  /// Get ID of undirected edge connecting two nodes
  ///
  /// \return Zero-based edge ID or -1 if edge does not exist.
  int EdgeId(int, int) const;

  /// Get IDs of edge nodes (ptId1 < ptId2)
  bool GetEdge(int, int &ptId1, int &ptId2) const;

  /// Get number of adjacent points
  int NumberOfAdjacentPoints(int) const;

  /// Access list of adjacent nodes (thread-safe)
  void GetAdjacentPoints(int, int &, const int *&) const;

  /// Get start and end pointer into list of adjacent nodes (thread-safe)
  void GetAdjacentPoints(int, const int *&, const int *&) const;

  /// Get list of adjacent points for a point of the given surface mesh
  static void GetAdjacentPoints(vtkPolyData *, int, vtkIdList *);
};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline int irtkEdgeTable::NumberOfPoints() const
{
  return Rows();
}

// -----------------------------------------------------------------------------
inline bool irtkEdgeTable::IsEdge(int ptId1, int ptId2) const
{
  return static_cast<bool>(Get(ptId1, ptId2));
}

// -----------------------------------------------------------------------------
inline int irtkEdgeTable::EdgeId(int ptId1, int ptId2) const
{
  return Get(ptId1, ptId2) - 1;
}

// -----------------------------------------------------------------------------
inline bool irtkEdgeTable::GetEdge(int edgeId, int &ptId1, int &ptId2) const
{
  int i, j;
  const int *row = _Row;
  const int *col = _Col;
  if (_Layout == CRS) swap(row, col);

  for (ptId2 = 0; ptId2 < NumberOfPoints(); ++ptId2) {
    for (i = col[ptId2], j = col[ptId2 + 1]; i < j; ++i) {
      ptId1 = row[i];
      if (ptId1 > ptId2) break;
      if (_Data[i] == edgeId) return true;
    }
  }

  ptId1 = ptId2 = -1;
  return false;
}

// -----------------------------------------------------------------------------
inline int irtkEdgeTable::NumberOfAdjacentPoints(int ptId) const
{
  if (_Layout == CCS) {
    return _Col[ptId+1] - _Col[ptId];
  } else {
    return _Row[ptId+1] - _Row[ptId];
  }
}

// -----------------------------------------------------------------------------
inline void irtkEdgeTable::GetAdjacentPoints(int ptId, int &numAdjPts, const int *&adjPtIds) const
{
  if (_Layout == CCS) {
    numAdjPts = _Col[ptId+1] - _Col[ptId];
    adjPtIds  = _Row + _Col[ptId];
  } else {
    numAdjPts = _Row[ptId+1] - _Row[ptId];
    adjPtIds  = _Col + _Row[ptId];
  }
}

// -----------------------------------------------------------------------------
inline void irtkEdgeTable::GetAdjacentPoints(int ptId, const int *&begin, const int *&end) const
{
  if (_Layout == CCS) {
    begin = _Row + _Col[ptId];
    end   = _Row + _Col[ptId + 1];
  } else {
    begin = _Col + _Row[ptId];
    end   = _Col + _Row[ptId + 1];
  }
}

////////////////////////////////////////////////////////////////////////////////
// irtkEdgeIterator
////////////////////////////////////////////////////////////////////////////////

/**
 * Thread-safe helper class for iteration of edges
 */
class irtkEdgeIterator
{
  const irtkEdgeTable &_Table;
  int                  _EdgeId;
  int                  _EndId;
  const int           *_PointId1;
  const int           *_ListEnd;
  int                  _PointId2;

public:

  /// Constructor
  irtkEdgeIterator(const irtkEdgeTable &table)
  :
    _Table(table), _EdgeId(-1), _EndId(-1),
    _PointId1(NULL), _ListEnd(NULL), _PointId2(-1)
  {}

  /// Initialize traversal of edges
  ///
  /// @param[in] begin ID of first edge to iterate.
  /// @param[in] end   ID one behind last edge to iterate. If negative,
  ///                  all edges from @p begin until the end are iterated.
  inline void InitTraversal(int begin = 0, int end = -1)
  {
    _EdgeId = begin;
    _EndId  = ((end < 0 || end > _Table.NumberOfEdges()) ? _Table.NumberOfEdges() : end);
    if (_EndId > _EdgeId) {
      for (_PointId2 = 0; _PointId2 < _Table.NumberOfPoints(); ++_PointId2) {
        _Table.GetAdjacentPoints(_PointId2, _PointId1, _ListEnd);
        while (_PointId1 != _ListEnd) {
          if (begin == 0 || (*_PointId1) > _PointId2) break;
          ++_PointId1, --begin;
        }
        if (begin == 0 && _PointId1 != _ListEnd && (*_PointId1) <= _PointId2) break;
      }
    } else {
      _PointId1 = _ListEnd = NULL;
      _PointId2 = -1;
    }
  }

  /// Initialize traversal of edges in range
  inline void InitTraversal(const blocked_range<int> &re)
  {
    InitTraversal(re.begin(), re.end());
  }

  /// Get next edge
  ///
  /// @param[out] ptId1 ID of first edge point.
  /// @param[out] ptId2 ID of second edge point.
  ///
  /// @return ID of undirected edge (ptId1, ptId2) or -1 if end reached
  template <class IdType>
  int GetNextEdge(IdType &ptId1, IdType &ptId2)
  {
    if (_EdgeId >= _EndId) return -1;
    ptId1 = static_cast<IdType>(*_PointId1);
    ptId2 = static_cast<IdType>(_PointId2);
    int edgeId = _EdgeId;
    if (++_EdgeId < _EndId) {
      if (++_PointId1 == _ListEnd || (*_PointId1) > _PointId2) {
        do {
          _Table.GetAdjacentPoints(++_PointId2, _PointId1, _ListEnd);
        } while (_PointId1 == _ListEnd);
      }
    }
    return edgeId;
  }
};


} } // namespace irtk::polydata

#endif
