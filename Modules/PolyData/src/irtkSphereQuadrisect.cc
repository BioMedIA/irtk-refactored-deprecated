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
#include <irtkImage.h>

#include <map>
#include <vtkMath.h>
#include <vtkTriangle.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkPlatonicSolidSource.h>
#include <vtkPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <irtkSphereQuadrisect.h>


irtkSphereQuadrisect::irtkSphereQuadrisect(){
}

irtkSphereQuadrisect::~irtkSphereQuadrisect(){
}



typedef struct {
    int id0, id1;
} edge;

struct edgeCmp {
   bool operator()( const edge ea, const edge eb ) const {
     return ( ea.id0 < eb.id0 ) || (( ea.id0 == eb.id0 ) && ( ea.id1 < eb.id1));
   }
 };


typedef map<edge, int, edgeCmp> midpointMap;

map<edge, short>::iterator iter;


char *output_name = NULL;

//VTK_SOLID_TETRAHEDRON  0
//VTK_SOLID_CUBE         1
//VTK_SOLID_OCTAHEDRON   2
//VTK_SOLID_ICOSAHEDRON  3
//VTK_SOLID_DODECAHEDRON 4

void Midpoint(double *x, double *y, double *m){
  m[0] = 0.5 * (x[0] + y[0]);
  m[1] = 0.5 * (x[1] + y[1]);
  m[2] = 0.5 * (x[2] + y[2]);
}

void copyPoints(vtkPoints *src, vtkPoints*dest){
  int i, noOfPts;
  double p[3];

  noOfPts = src->GetNumberOfPoints();
  dest->SetNumberOfPoints(noOfPts);

  for (i = 0; i < noOfPts; ++i){
    src->GetPoint(i, p);
    dest->InsertPoint(i, p);
  }
  cerr << "Copied " << noOfPts << " points" << endl;
}

void copyFaces(vtkCellArray *src, vtkCellArray *dest){
  int i, j, noOfCells;
  vtkIdType npts = 0;
  vtkIdType *ptIds;
  noOfCells = src->GetNumberOfCells();
  dest->SetNumberOfCells(noOfCells);

  int estSize;
  int maxPts;
  maxPts = src->GetMaxCellSize();

  dest->Initialize();
//  cerr << "dest cells " << dest->GetNumberOfCells() << endl;
  estSize = dest->EstimateSize(noOfCells, maxPts);
  dest->Allocate(estSize, 0);
  dest->Reset();
//  cerr << "dest cells " << dest->GetNumberOfCells() << endl;

  dest->InitTraversal();
  src->InitTraversal();
  for (i = 0; i < noOfCells; ++i){
    src->GetNextCell(npts, ptIds);
//    cerr << "npts " << npts << endl;
    dest->InsertNextCell(npts);
    for (j = 0; j < npts; ++j){
      dest->InsertCellPoint(ptIds[j]);
//      cerr << " " << ptIds[j];
    }
//    cerr << endl;
    dest->UpdateCellCount(npts);
  }
//  cerr << "dest cells " << dest->GetNumberOfCells() << endl;
  dest->Squeeze();
  cerr << "dest cells " << dest->GetNumberOfCells() << endl;

  cerr << "Copied " << noOfCells << " cells" << endl;
}

void createIdListArray(vtkIdList** &list, int size){
  int i;
  list = new vtkIdList*[size];
  for (i = 0; i < size; ++i){
    list[i] = vtkIdList::New();
  }
}

void copyIdListArray(vtkIdList** &src, vtkIdList** &dest, int noOfVerts){
  int i, j, nPts;
  // Assume both lists initialised to same size.
  for (i = 0; i < noOfVerts; ++i){
    nPts = src[i]->GetNumberOfIds();
    for (j = 0; j < nPts; ++j){
      dest[i]->InsertUniqueId(src[i]->GetId(j));
    }
  }
}

void deleteIdListArray(vtkIdList** &list, int size){
  int i;
  for (i = 0; i < size; ++i){
    list[i]->Delete();
  }
  delete [] list;
  list = NULL;
}

void addAdjacency(vtkIdList**adj, vtkIdType a, vtkIdType b){
  int temp;
  if (a == b)
    return;
  if (a > b){
    temp = a;
    a = b;
    b = temp;
  }
  adj[a]->InsertUniqueId(b);
}

void getAdjacency(vtkCellArray *faces, vtkIdList** adj){

  int i, j, noOfFaces;
  vtkIdType ida, idb;

  vtkIdType npts = 0;
  vtkIdType *pts;

  noOfFaces = faces->GetNumberOfCells();

  faces->InitTraversal();
  for (i = 0; i < noOfFaces; ++i){

    faces->GetNextCell(npts, pts);

    ida = pts[0];
    idb = pts[npts - 1];

    if (ida < idb){
      adj[ida]->InsertUniqueId(idb);
    } else {
      adj[idb]->InsertUniqueId(ida);
    }

    for (j = 0; j < npts - 1; ++j){
      ida = pts[j];
      idb = pts[j + 1];
      if (ida < idb){
        adj[ida]->InsertUniqueId(idb);
      } else {
        adj[idb]->InsertUniqueId(ida);
      }
    }
  }
}

int getEdgeCount(vtkIdList** adj, int noOfVertices){
  int i, edgeCount;

  edgeCount = 0;
  for (i = 0; i < noOfVertices; ++i){
    edgeCount += adj[i]->GetNumberOfIds();
  }
  return edgeCount;
}

void printStuff(vtkPoints *points, vtkIdList** adj){

  int i, j, noOfPts;
  noOfPts = points->GetNumberOfPoints();
  double a[3];

  for (i = 0; i < noOfPts; ++i){
    points->GetPoint(i, a);
    cerr << i << " : " << a[0] << " " << a[1] << " " << a[2] << endl;
    for(j = 0; j < adj[i]->GetNumberOfIds(); ++j){
      cerr << " " << i << " - " << adj[i]->GetId(j) << endl;
    }
  }
  cerr << "Edge count: " << getEdgeCount(adj, noOfPts) << endl;

}





vtkSmartPointer<vtkPolyData> irtkSphereQuadrisect::Quadrisect(vtkSmartPointer<vtkPolyData> input, int levels) {
  cout << "numcells " << input->GetNumberOfCells() << endl;

  bool ok;
  int i, j, k;
  double a[3];
  double v0[3], v1[3], v2[3];
  vtkIdType npts = 0;
  vtkIdType id0, id1, id2, ida, idb, idc;
  vtkIdType *pts;
  int vertCounter;
  int noOfVertsOld = 0, noOfEdgesOld, noOfFacesOld;
  int noOfVertsNew, noOfEdgesNew, noOfFacesNew;
  midpointMap midptMap;
  edge e;
  int solidIndex = VTK_SOLID_ICOSAHEDRON;
  int maxPts, estSize;


  vtkPoints *oldPoints = vtkPoints::New();
  vtkCellArray* oldFaces = vtkCellArray::New();
  vtkIdList** oldAdj = NULL;

  vtkPoints *newPoints = vtkPoints::New();
  vtkCellArray* newFaces = vtkCellArray::New();
  vtkIdList** newAdj = NULL;

//  copyPoints(input->GetPoints(), newPoints);
  newPoints->DeepCopy(input->GetPoints());

//  copyFaces(input->GetPolys(), newFaces);
  newFaces->DeepCopy(input->GetPolys());

  noOfVertsNew = newPoints->GetNumberOfPoints();
  noOfFacesNew = newFaces->GetNumberOfCells();

  createIdListArray(newAdj, noOfVertsNew);

  // Convert the adjacency information, that is implicit in the faces,
  // into an adjacency list.
  getAdjacency(newFaces, newAdj);

  printStuff(newPoints, newAdj);

  for (k = 0; k < levels; k++){

    // Delete old adjacency if necessary
    if (oldAdj != NULL){
      deleteIdListArray(oldAdj, noOfVertsOld);
    }

    // Copy the new points into the old (current) points.
//    copyPoints(newPoints, oldPoints);
    oldPoints->DeepCopy(newPoints);

    noOfVertsOld = oldPoints->GetNumberOfPoints();

    // Same for faces.
//    copyFaces(newFaces, oldFaces);
    oldFaces->DeepCopy(newFaces);

    noOfFacesOld = oldFaces->GetNumberOfCells();

    // Same for adjacency (edge) info.
    createIdListArray(oldAdj, noOfVertsOld);
    copyIdListArray(newAdj, oldAdj, noOfVertsOld);
    noOfEdgesOld = getEdgeCount(oldAdj, noOfVertsOld);
    deleteIdListArray(newAdj, noOfVertsOld);

    // How many vertices, edges and faces after the current iteration?
    // V' = V + E
    // F' = 4F
    // E' = 2E + 3F

    cerr << "old:  v " << noOfVertsOld << " e " << noOfEdgesOld << " f " << noOfFacesOld << endl;

    noOfVertsNew = noOfVertsOld + noOfEdgesOld;
    noOfFacesNew = 4 * noOfFacesOld;
    noOfEdgesNew = 2 * noOfEdgesOld + 3 * noOfFacesOld;

    // Allocate space for new vertices edges faces
    newPoints->SetNumberOfPoints(noOfVertsNew);

    maxPts = oldFaces->GetMaxCellSize();
    newFaces->Initialize();
    estSize = newFaces->EstimateSize(noOfFacesNew, maxPts);
    newFaces->Allocate(estSize, 0);
    newFaces->Reset();

    createIdListArray(newAdj, noOfVertsNew);

    // Copy over original vertices, these will retain the same indices in the vertex array.
    for (i = 0; i < noOfVertsOld; ++i){
      oldPoints->GetPoint(i, a);
      newPoints->InsertPoint(i, a);
    }
    // Generate all the new vertices, one for each of the old edges.
    vertCounter = noOfVertsOld;
    for (i = 0; i < noOfVertsOld; ++i){
      id0 = i;
      newPoints->GetPoint(id0, v0);
      npts = oldAdj[id0]->GetNumberOfIds();
      // Loop over list of neighbours.
      for (j = 0; j < npts; ++j){
        id1 = oldAdj[i]->GetId(j);
        newPoints->GetPoint(id1, v1);
        Midpoint(v0, v1, a);
        vtkMath::Normalize(a);
        newPoints->InsertPoint(vertCounter, a);

        // We know that id0 should be < id1
        e.id0 = id0;
        e.id1 = id1;
        midptMap[e] = vertCounter;
        ++vertCounter;
      }
    }

    // Loop over old faces, subdividing each.
    oldFaces->InitTraversal();
    newFaces->InitTraversal();
    cerr << "new faces : " << newFaces->GetNumberOfCells() << endl;

    for(i = 0; i < noOfFacesOld; ++i){
      oldFaces->GetNextCell(npts, pts);

      if (npts != 3){
        cerr << "Expecting triangles only!" << endl;
        exit(1);
      }

      id0 = pts[0];
      id1 = pts[1];
      id2 = pts[2];

      newPoints->GetPoint(id0, v0);
      newPoints->GetPoint(id1, v1);
      newPoints->GetPoint(id2, v2);
      //        Subdivide each triangle in the old approximation and normalize
      //        the new points thus generated to lie on the surface of the unit
      //        sphere.
      //        Each input triangle with vertices labelled [0,1,2] as shown
      //        below will be turned into four new triangles:
      //
      //            Make new points
      //                a = (0+2)/2
      //                b = (0+1)/2
      //                c = (1+2)/2
      //             1
      //             /\   Normalize a, b, c
      //            /  \
      //          b/____\ c   Construct new triangles
      //          /\    /\        [0,b,a]
      //         /  \  /  \       [b,1,c]
      //        /____\/____\      [a,b,c]
      //       0     a      2     [a,c,2]

      // Look up the ids of the midpoints.
      e.id0 = min(id0, id2);
      e.id1 = max(id0, id2);
      ida = midptMap[e];

      e.id0 = min(id0, id1);
      e.id1 = max(id0, id1);
      idb = midptMap[e];

      e.id0 = min(id1, id2);
      e.id1 = max(id1, id2);
      idc = midptMap[e];

      // Add the new adjacencies.
      addAdjacency(newAdj, id0, ida);
      addAdjacency(newAdj, id2, ida);

      addAdjacency(newAdj, id0, idb);
      addAdjacency(newAdj, id1, idb);

      addAdjacency(newAdj, id1, idc);
      addAdjacency(newAdj, id2, idc);

      addAdjacency(newAdj, ida, idb);
      addAdjacency(newAdj, ida, idc);
      addAdjacency(newAdj, idb, idc);

      // Add the new faces.
      newFaces->InsertNextCell(3);
      newFaces->InsertCellPoint(id0);
      newFaces->InsertCellPoint(idb);
      newFaces->InsertCellPoint(ida);
      newFaces->UpdateCellCount(3);

      newFaces->InsertNextCell(3);
      newFaces->InsertCellPoint(idb);
      newFaces->InsertCellPoint(id1);
      newFaces->InsertCellPoint(idc);
      newFaces->UpdateCellCount(3);

      newFaces->InsertNextCell(3);
      newFaces->InsertCellPoint(ida);
      newFaces->InsertCellPoint(idb);
      newFaces->InsertCellPoint(idc);
      newFaces->UpdateCellCount(3);

      newFaces->InsertNextCell(3);
      newFaces->InsertCellPoint(ida);
      newFaces->InsertCellPoint(idc);
      newFaces->InsertCellPoint(id2);
      newFaces->UpdateCellCount(3);

    }

    newFaces->Squeeze();
    cerr << "new faces : " << newFaces->GetNumberOfCells() << endl;

  }

  vtkSmartPointer<vtkPolyData>  output = vtkSmartPointer<vtkPolyData> ::New();
  output->SetPoints(newPoints);
  output->SetPolys(newFaces);

  deleteIdListArray(newAdj, noOfVertsNew);
  return output;
}
