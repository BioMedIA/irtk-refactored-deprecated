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

#include <irtkCommon.h>
#include <irtkRemesher.h>

#include <vtkMath.h>
#include <vtkPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkDataArray.h>
#include <vtkCellData.h>
#include <vtkExtractEdges.h>
#include <vtkPriorityQueue.h>
#include <vtkTriangleFilter.h>
#include <vtkPointData.h>


// -----------------------------------------------------------------------------
double irtkRemesher::GetEdgeLength(double p1[], double p2[])
{
  return sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
}

// -----------------------------------------------------------------------------
void irtkRemesher::Midpoint(double p1[], double p2[], double mp[])
{
  for (int i = 0; i < 3; ++i) mp[i] = 0.5 * (p1[i] + p2[i]);
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> irtkRemesher::CleanPolyData(vtkPolyData *poly)
{
  vtkSmartPointer<vtkCleanPolyData> cleaner;
  cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
  cleaner->PointMergingOn();
  cleaner->ConvertLinesToPointsOn();
  cleaner->ConvertPolysToLinesOn();
  cleaner->ConvertStripsToPolysOn();
  cleaner->ToleranceIsAbsoluteOn();
  cleaner->SetAbsoluteTolerance(.0);
  SetVTKInput(cleaner, poly);
  cleaner->Update();
  cleaner->GetOutput()->SetVerts(NULL);
  cleaner->GetOutput()->SetLines(NULL);
  return cleaner->GetOutput();
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> irtkRemesher::ProcessInput(vtkPolyData *poly)
{
  vtkSmartPointer<vtkTriangleFilter> triangulator;
  triangulator = vtkSmartPointer<vtkTriangleFilter>::New();
  vtkSmartPointer<vtkPolyData> cleanPolyData = CleanPolyData(poly);
  SetVTKInput(triangulator, cleanPolyData);
  triangulator->Update();
  triangulator->GetOutput()->BuildCells();
  triangulator->GetOutput()->BuildLinks();
  return triangulator->GetOutput();
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> irtkRemesher::BuildEdgeStructure(vtkPolyData *poly)
{
  vtkSmartPointer<vtkExtractEdges> edges = vtkSmartPointer<vtkExtractEdges>::New();
  SetVTKInput(edges, poly);
  edges->Update();
  edges->GetOutput()->BuildCells();
  return edges->GetOutput();
}

// -----------------------------------------------------------------------------
// TODO: Need to parallalise this!
vtkSmartPointer<vtkDoubleArray> irtkRemesher::ComputeEdgeLengths(vtkPolyData *edges)
{
  double d, p0[3], p1[3];
  int    id0, id1;

  vtkSmartPointer<vtkDoubleArray> edgeLengths = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkIdList>      edge        = vtkSmartPointer<vtkIdList>::New();

  const int numEdges = edges->GetNumberOfLines();
  edgeLengths->SetNumberOfComponents(1);
  edgeLengths->SetNumberOfTuples(numEdges);
  edgeLengths->SetName("Length");

  for (int edgeId = 0; edgeId < numEdges; ++edgeId) {
    edges->GetCellPoints(edgeId, edge);
    id0 = edge->GetId(0);
    id1 = edge->GetId(1);
    edges->GetPoint(id0, p0);
    edges->GetPoint(id1, p1);
    d = GetEdgeLength(p0, p1);
    edgeLengths->InsertTuple1(edgeId, d);
  }

  return edgeLengths;
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPriorityQueue>
irtkRemesher::QueueEdgesByLength(vtkPolyData *edges, double minEdgeLength, double maxEdgeLength)
{
  int    splitNum = 0, collapseNum = 0;
  int    numEdges = edges->GetNumberOfLines();
  double priority, d;

  vtkSmartPointer<vtkPriorityQueue> edgeQueue   = vtkSmartPointer<vtkPriorityQueue>::New();
  vtkSmartPointer<vtkDoubleArray>   edgeLengths = ComputeEdgeLengths(edges);

  vtkSmartPointer<vtkIntArray> isCollapse = vtkSmartPointer<vtkIntArray>::New();
  isCollapse->SetNumberOfComponents(1);
  isCollapse->SetNumberOfTuples(numEdges);
  isCollapse->SetName("isCollapse");

  for (int edgeId = 0; edgeId < numEdges; ++edgeId) {
    d = edgeLengths->GetTuple1(edgeId);
    if (d > maxEdgeLength) {
      priority = maxEdgeLength / (d - maxEdgeLength);
      edgeQueue->Insert(priority, edgeId);
      splitNum += 1;
      isCollapse->SetTuple1(edgeId, 0);
    } else if (d < minEdgeLength) {
      priority = minEdgeLength / (minEdgeLength - d);
      edgeQueue->Insert(priority, edgeId);
      collapseNum += 1;
      isCollapse->SetTuple1(edgeId, 1);
    }
  }
  if (verbose) {
    cout << "Attempting to split " << splitNum << " edges and collapse " << collapseNum
         << " edges (" << numEdges << " total)." << endl;
  }
  edges->GetCellData()->AddArray(isCollapse);
  return edgeQueue;
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkIdList>
irtkRemesher::GetAdjacentCellIds(int pId0, int pId1, vtkPolyData *poly)
{
  vtkSmartPointer<vtkIdList> cellList0, cellList1;
  cellList0 = vtkSmartPointer<vtkIdList>::New();
  cellList1 = vtkSmartPointer<vtkIdList>::New();

  poly->GetPointCells(pId0, cellList0);
  poly->GetPointCells(pId1, cellList1);
  cellList0->IntersectWith(cellList1); // common cells that both points are in

  return cellList0;
}

// -----------------------------------------------------------------------------
void irtkRemesher::DeleteEdgesAfterCollapse(int pId0,int pId1, vtkPolyData *edges)
{
  vtkSmartPointer<vtkIdList> edgeList0, edgeList1;
  edgeList0 = vtkSmartPointer<vtkIdList>::New();
  edgeList1 = vtkSmartPointer<vtkIdList>::New();

  edges->GetPointCells(pId0, edgeList0);
  edges->GetPointCells(pId1, edgeList1);

  for (int i = 0; i < edgeList0->GetNumberOfIds(); ++i) {
    edges->DeleteCell(edgeList0->GetId(i));
  }
  for (int i = 0; i < edgeList1->GetNumberOfIds(); ++i) {
    edges->DeleteCell(edgeList1->GetId(i));
  }
}

// -----------------------------------------------------------------------------
void DeleteEdgesAfterSplit(int pId0, int pId1, int otherPId0, int otherPId1, vtkPolyData *edges)
{
  vtkSmartPointer<vtkIdList> edgeList0      = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> edgeList1      = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> otherEdgeList0 = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> otherEdgeList1 = vtkSmartPointer<vtkIdList>::New();

  edges->GetPointCells(pId0,      edgeList0);
  edges->GetPointCells(pId1,      edgeList1);
  edges->GetPointCells(otherPId0, otherEdgeList0);
  edges->GetPointCells(otherPId1, otherEdgeList1);

  for (int i = 0; i < edgeList1->GetNumberOfIds(); ++i) {
    edgeList0->InsertUniqueId(edgeList1->GetId(i));
  }
  for (int i = 0; i < otherEdgeList1->GetNumberOfIds(); ++i) {
    otherEdgeList0->InsertUniqueId(otherEdgeList1->GetId(i));
  }
  edgeList0->IntersectWith(otherEdgeList0);
  for (int i = 0; i < edgeList0->GetNumberOfIds(); ++i) {
    edges->DeleteCell(edgeList0->GetId(i));
  }
}

// -----------------------------------------------------------------------------
void irtkRemesher::CreateNewCells(int oldCellId0, int oldCellId1,int pId0, int pId1,
                                  int newPId, vtkCellArray * newCells, vtkPolyData * poly)
{
  int i;

  vtkSmartPointer<vtkIdList> oldCellPIds0 = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> oldCellPIds1 = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> newCell      = vtkSmartPointer<vtkIdList>::New();

  poly->GetCellPoints(oldCellId0, oldCellPIds0);
  poly->GetCellPoints(oldCellId1, oldCellPIds1);

  newCell->DeepCopy(oldCellPIds0);
  i = newCell->IsId(pId0);
  newCell->SetId(i, newPId);
  newCells->InsertNextCell(newCell);

  newCell->DeepCopy(oldCellPIds0);
  i = newCell->IsId(pId1);
  newCell->SetId(i, newPId);
  newCells->InsertNextCell(newCell);

  newCell->DeepCopy(oldCellPIds1);
  i = newCell->IsId(pId0);
  newCell->SetId(i, newPId);
  newCells->InsertNextCell(newCell);

  newCell->DeepCopy(oldCellPIds1);
  i = newCell->IsId(pId1);
  newCell->SetId(i, newPId);
  newCells->InsertNextCell(newCell);
}

// -----------------------------------------------------------------------------
void irtkRemesher::CombineCells(vtkPolyData *poly, vtkCellArray *newCells, vtkIntArray *toDelete)
{
  //add the new cells we created
  vtkSmartPointer<vtkIdList> newCell = vtkSmartPointer<vtkIdList>::New();
  for (int i = 0; i < poly->GetNumberOfCells(); ++i) {
    if (toDelete->GetTuple1(i) == 0) {
      poly->GetCellPoints(i,newCell);
      newCells->InsertNextCell(newCell);
    }
  }
  newCells->Squeeze();
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkIdList>
GetOtherPoints(int pId0, int pId1, int oldCellId0, int oldCellId1, vtkPolyData *poly)
{
  int otherId;

  vtkSmartPointer<vtkIdList> oldCellPIds0 = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> oldCellPIds1 = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> otherPIds    = vtkSmartPointer<vtkIdList>::New();
  otherPIds->Allocate(2);

  poly->GetCellPoints(oldCellId0, oldCellPIds0);
  poly->GetCellPoints(oldCellId1, oldCellPIds1);

  for (int i = 0; i < 3; ++i) {
    otherId = oldCellPIds0->GetId(i);
    if (otherId != pId0 && otherId != pId1) {
      otherPIds->InsertNextId(otherId);
      break;
    }
  }

  for (int i = 0; i < 3; ++i) {
    otherId = oldCellPIds1->GetId(i);
    if (otherId != pId0 && otherId != pId1){
      otherPIds->InsertNextId(otherId);
      break;
    }
  }

  return otherPIds;
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData>
irtkRemesher::ProcessEdgeQueue(vtkPriorityQueue *edgeQueue, vtkPolyData *edges, vtkPolyData *poly, vtkPolyData *newPoly)
{
  int collapseCount = 0, splitCount = 0;
  int edgeId;
  int pId0,pId1,newPId;
  int oldCellId0,oldCellId1;
  double p0 [3],p1 [3],m [3];

  // TODO: Copy/Interpolate cell data as well
  vtkSmartPointer<vtkPointData> inputPD  = poly   ->GetPointData();
  vtkSmartPointer<vtkPointData> outputPD = newPoly->GetPointData();
  outputPD->InterpolateAllocate(inputPD);
  outputPD->CopyData(inputPD, 0, poly->GetNumberOfPoints());

  vtkSmartPointer<vtkIdList>    edge       = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList>    oldCellIds = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList>    otherPIds  = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkCellArray> newCells   = vtkSmartPointer<vtkCellArray>::New();
  //TODO newPoly num cells change to number of cells to split * 4
  newCells->Allocate(newPoly->GetNumberOfCells() + 4 * edgeQueue->GetNumberOfItems());

  vtkIntArray * isCollapse = (vtkIntArray *) edges->GetCellData()->GetArray("isCollapse");
  vtkSmartPointer<vtkIntArray> cellsToDelete = vtkSmartPointer<vtkIntArray>::New();
  for (int i = 0; i < newPoly->GetNumberOfCells(); ++i) {
    cellsToDelete->InsertNextTuple1(0);
  }

  while (true) {
    edgeId = edgeQueue->Pop();
    if (edgeId < 0) break; //break loop when we've run out of edges

    //step 1 if edge has been marked as deleted skip
    if (edges->GetCellType(edgeId) == 0) continue;

    //step 2 find the pids for the edge points and the cells that border it
    edges->GetCellPoints(edgeId,edge);
    pId0 = edge->GetId(0);
    pId1 = edge->GetId(1);
    oldCellIds = GetAdjacentCellIds(pId0, pId1, newPoly);

    //if (oldCellIds->GetNumberOfIds() < 2){
    //cout << "ERROR - found an edge with less than 2 adajacent faces" << endl;
    //continue;

    //There are cases where the two edge points have more than two edges with
    //common points. We will skip this for now, because its extremely rare and
    //much more complex.
    if (oldCellIds->GetNumberOfIds() != 2) continue;

    oldCellId0 = oldCellIds->GetId(0);
    oldCellId1 = oldCellIds->GetId(1);

    //step 3 calculate edge midpoint
    newPoly->GetPoint(pId0,p0);
    newPoly->GetPoint(pId1,p1);
    Midpoint(p0,p1,m);

    //step 4 determine if we're collapsing or splitting
    if (isCollapse->GetTuple1(edgeId) == 1) { //collapse

      // Move edge points to midpoint, effectively collapsing it!
      newPoly->GetPoints()->SetPoint(pId0,m);
      newPoly->GetPoints()->SetPoint(pId1,m);

      //mark neighbouring edges as "deleted" so they aren't processed in the
      //next iteration
      DeleteEdgesAfterCollapse(pId0,pId1,edges);

      collapseCount++;
    } else { //split
        //Insert mid point as a new point and createnew cells
        newPId = newPoly->GetPoints()->InsertNextPoint(m);
        CreateNewCells(oldCellId0,oldCellId1,pId0,pId1,newPId,newCells,newPoly);
        outputPD->InterpolateEdge(inputPD, newPId, pId0, pId1, .5);

        //mark neighbouring edges as "deleted" so they aren't processed in the
        //next iteration
        otherPIds = GetOtherPoints(pId0,pId1,oldCellId0,oldCellId1,poly);
        DeleteEdgesAfterSplit(pId0,pId1,otherPIds->GetId(0),otherPIds->GetId(1),edges);
        edges->DeleteCell(edgeId);

        splitCount++;
    }

    //step 6 finally mark the old cells for deletion
    cellsToDelete->InsertTuple1(oldCellId0, 1);
    cellsToDelete->InsertTuple1(oldCellId1, 1);
  }

  // combine cells
  CombineCells(newPoly, newCells, cellsToDelete);
  newPoly->SetPolys(newCells);

  if (verbose) {
    cout << "Successfully split " << splitCount  << " edges and collapsed " << collapseCount << " edges." << endl;
    cout << "The new mesh consists of " << newPoly->GetNumberOfCells() << " cells." << endl;
  }

  return CleanPolyData(newPoly);
}

// -----------------------------------------------------------------------------
irtkRemesher::irtkRemesher()
{
  _MaxEdgeLength = numeric_limits<double>::max();
  _MinEdgeLength = -1;
}

// -----------------------------------------------------------------------------
irtkRemesher::~irtkRemesher()
{
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> irtkRemesher::Remesh(vtkPolyData *input, int *numedges)
{
  if (verbose) cout << "Remeshing" << endl;
  vtkSmartPointer<vtkPolyData> output = vtkSmartPointer<vtkPolyData>::New();
  output->DeepCopy(input);
  output->GetPointData()->Initialize();
  output->GetCellData ()->Initialize();
  vtkSmartPointer<vtkPolyData> edges = BuildEdgeStructure(output);
  vtkSmartPointer<vtkPriorityQueue> edgeQueue;
  edgeQueue = QueueEdgesByLength(edges, _MinEdgeLength, _MaxEdgeLength);
  if (numedges) *numedges = edgeQueue->GetNumberOfItems();
  return ProcessEdgeQueue(edgeQueue, edges, input, output);
}
