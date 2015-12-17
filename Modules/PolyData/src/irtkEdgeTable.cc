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

#include <irtkEdgeTable.h>

#include <vtkSmartPointer.h>
#include <vtkGenericCell.h>


namespace irtk { namespace polydata {


// -----------------------------------------------------------------------------
irtkEdgeTable::irtkEdgeTable(vtkDataSet *mesh)
{
  if (mesh) Initialize(mesh);
}

// -----------------------------------------------------------------------------
irtkEdgeTable::irtkEdgeTable(const irtkEdgeTable &other)
:
  irtkGenericSparseMatrix(other)
{
}

// -----------------------------------------------------------------------------
irtkEdgeTable &irtkEdgeTable::operator =(const irtkEdgeTable &other)
{
  irtkGenericSparseMatrix::operator =(other);
  return *this;
}

// -----------------------------------------------------------------------------
irtkEdgeTable::~irtkEdgeTable()
{
}

// -----------------------------------------------------------------------------
void irtkEdgeTable::Initialize(vtkDataSet *mesh)
{
  IRTK_START_TIMING();

  _NumberOfEdges = 0;

  const int       numPts   = static_cast<int>(mesh->GetNumberOfPoints());
  const vtkIdType numCells = mesh->GetNumberOfCells();

  typedef irtkGenericSparseMatrix<int>::Entries Entries;
  Entries *entries = new Entries[numPts];

  vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();

  Entries::const_iterator entry;
  vtkIdType cellId, edgeId, ptId1, ptId2;
  int numCellEdges, numEdgePts;
  bool new_edge;
  vtkCell *edge;

  for (cellId = 0; cellId < numCells; ++cellId) {
    mesh->GetCell(cellId, cell);
    numCellEdges = cell->GetNumberOfEdges();
    for (edgeId = 0; edgeId < numCellEdges; ++edgeId) {
      edge = cell->GetEdge(edgeId);
      if (edge->IsLinear()) {
        numEdgePts = edge->GetNumberOfPoints();
        if (numEdgePts > 1) {
          ptId1 = edge->PointIds->GetId(0);
          for (int i = 1; i < numEdgePts; ++i, ptId1 = ptId2) {
            ptId2 = edge->PointIds->GetId(i);
            new_edge = true;
            for (entry = entries[ptId1].begin(); entry != entries[ptId1].end(); ++entry) {
              if (entry->first == ptId2) {
                new_edge = false;
                break;
              }
            }
            if (new_edge) {
              // Symmetric entries such that AdjacentPoints is efficient
              ++_NumberOfEdges; // edgeId + 1 such that entries are non-zero
              entries[ptId1].push_back(make_pair(ptId2, _NumberOfEdges));
              entries[ptId2].push_back(make_pair(ptId1, _NumberOfEdges));
            }
          }
        }
      } else {
        cerr << "WARNING: irtkEdgeTable::Initialize: Only linear edges supported" << endl;
      }
    }
  }

  for (int i = 0; i < numPts; ++i) {
    sort(entries[i].begin(), entries[i].end());
  }

  irtkGenericSparseMatrix::Initialize(numPts, numPts, entries, true);
  delete[] entries;

  // Reassign edge IDs -- same order as edges are visited by irtkEdgeIterator!
  int i1, i2, j1, j2;
  const int *row = _Row;
  const int *col = _Col;
  if (_Layout == CRS) swap(row, col);

  int edgeIdPlusOne = 0;
  for (ptId2 = 0; ptId2 < numPts; ++ptId2) {
    for (i1 = col[ptId2], j1 = col[ptId2 + 1]; i1 < j1; ++i1) {
      ptId1 = row[i1];
      if (ptId1 > ptId2) break;
      i2 = col[ptId1];
      j2 = col[ptId1 + 1];
      while (i2 < j2 && row[i2] < ptId2) ++i2;
      irtkAssert(i2 < j2 && row[i2] == ptId2, "matrix is symmetric");
      // one-based IDs to have non-zero entries in sparse matrix
      _Data[i1] = _Data[i2] = ++edgeIdPlusOne;
    }
  }
  irtkAssert(edgeIdPlusOne == _NumberOfEdges, "edge ID reassigned is consistent");

  IRTK_DEBUG_TIMING(5, "initialization of edge table");
}

// -----------------------------------------------------------------------------
void irtkEdgeTable::GetAdjacentPoints(vtkPolyData *mesh, int ptId, vtkIdList *ids)
{
  unsigned short  ncells;
  int             numCellEdges, numEdgePts;
  vtkIdType      *cells, ptId1, ptId2;
  vtkCell        *cell, *edge;

  ids->Reset();
  mesh->GetPointCells(ptId, ncells, cells);

  for (unsigned short i = 0; i < ncells; ++i) {
    cell = mesh->GetCell(cells[i]);
    numCellEdges = cell->GetNumberOfEdges();
    for (int edgeId = 0; edgeId < numCellEdges; ++edgeId) {
      edge = cell->GetEdge(edgeId);
      if (edge->IsLinear()) {
        numEdgePts = edge->GetNumberOfPoints();
        if (numEdgePts > 1) {
          ptId1 = edge->PointIds->GetId(0);
          for (int i = 1; i < numEdgePts; ++i, ptId1 = ptId2) {
            ptId2 = edge->PointIds->GetId(i);
            if (ptId1 == ptId) {
              ids->InsertUniqueId(ptId2);
            } else if (ptId2 == ptId) {
              ids->InsertUniqueId(ptId1);
            }
          }
        }
      } else {
        cerr << "WARNING: irtkEdgeTable::GetAdjacentPoints: Only linear edges supported" << endl;
      }
    }
  }
}


} } // namespace irtk::polydata
