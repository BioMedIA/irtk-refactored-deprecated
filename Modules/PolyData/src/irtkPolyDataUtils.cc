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

#include <irtkPolyDataUtils.h>


#include <algorithm>

#include <irtkPath.h>
#include <irtkDefines.h>
#include <irtkParallel.h>
#include <irtkEdgeTable.h>
#include <irtkDataStatistics.h>

#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkDataSetAttributes.h>
#include <vtkCell.h>
#include <vtkTriangle.h>
#include <vtkTetra.h>

#include <vtkOBJReader.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>
#include <vtkPLYReader.h>
#include <vtkPLYWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkXMLGenericDataObjectReader.h>
#include <vtkXMLDataSetWriter.h>
#include <vtkDataSetWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>

#include <vtkCleanPolyData.h>
#include <vtkTriangleFilter.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkHull.h>
#include <vtkMaskPoints.h>
#include <vtkDelaunay3D.h>
#include <vtkMassProperties.h>

using namespace std;


namespace irtk { namespace polydata {

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
const char *DefaultExtension(vtkDataSet *dataset)
{
  if      (vtkPolyData        ::SafeDownCast(dataset)) return ".vtp";
  else if (vtkUnstructuredGrid::SafeDownCast(dataset)) return ".vtu";
  else if (vtkStructuredGrid  ::SafeDownCast(dataset)) return ".vts";
  return ".vtk";
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPointSet> ReadPointSet(const char *fname, int *ftype)
{
  vtkSmartPointer<vtkPointSet> pointset;
  const std::string ext = Extension(fname);
  if (ext == ".vtp" || ext == ".stl" || ext == ".ply" || ext == ".obj") {
    pointset = ReadPolyData(fname);
  } else if (ext.length() == 4  && ext.substr(0, 3) == ".vt" && ext != ".vtk") {
    vtkSmartPointer<vtkXMLGenericDataObjectReader> reader;
    reader = vtkSmartPointer<vtkXMLGenericDataObjectReader>::New();
    reader->SetFileName(fname);
    reader->Update();
    pointset = vtkPointSet::SafeDownCast(reader->GetOutput());
  } else {
    vtkSmartPointer<vtkGenericDataObjectReader> reader;
    reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
    reader->SetFileName(fname);
    reader->Update();
    if (ftype) *ftype = reader->GetFileType();
    pointset = vtkPointSet::SafeDownCast(reader->GetOutput());
  }
  return pointset;
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> ReadPolyData(const char *fname, int *ftype)
{
  const std::string ext = Extension(fname);
  if (ext == ".vtp") {
    vtkSmartPointer<vtkXMLPolyDataReader> reader;
    reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader->SetFileName(fname);
    reader->Update();
    return reader->GetOutput();
  } else if (ext == ".stl") {
    vtkSmartPointer<vtkSTLReader> reader;
    reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(fname);
    reader->Update();
    return reader->GetOutput();
  } else if (ext == ".ply") {
    vtkSmartPointer<vtkPLYReader> reader;
    reader = vtkSmartPointer<vtkPLYReader>::New();
    reader->SetFileName(fname);
    reader->Update();
    return reader->GetOutput();
  } else if (ext == ".obj") {
    vtkSmartPointer<vtkOBJReader> reader;
    reader = vtkSmartPointer<vtkOBJReader>::New();
    reader->SetFileName(fname);
    reader->Update();
    return reader->GetOutput();
  } else {
    vtkSmartPointer<vtkPolyDataReader> reader;
    reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(fname);
    reader->Update();
    if (ftype) *ftype = reader->GetFileType();
    return reader->GetOutput();
  }
}

// -----------------------------------------------------------------------------
/// Write dataset points to TetGen .node file
int WriteTetGenNode(std::ostream &os, vtkPolyData *polydata)
{
  vtkPointData *pd = polydata->GetPointData();
  vtkDataArray *ar;
  int nattributes = 0;
  for (int i = 0; i < pd->GetNumberOfArrays(); ++i) {
    nattributes += pd->GetArray(i)->GetNumberOfComponents();
  }
  os << polydata->GetNumberOfPoints() << " 3 " << nattributes << " 0\n";
  double p[3];
  std::streamsize precision = os.precision();
  for (vtkIdType ptId = 0; ptId < polydata->GetNumberOfPoints(); ++ptId) {
    os << (ptId + 1) << " ";
    os.precision(8); // default TetGen tolerance is 1e-8
    polydata->GetPoint(ptId, p);
    os << " " << p[0] << " " << p[1] << " " << p[2];
    os.precision(5);
    for (int i = 0; i < pd->GetNumberOfArrays(); ++i) {
      ar = pd->GetArray(i);
      for (int j = 0; j < ar->GetNumberOfComponents(); ++j) {
        os << " " << ar->GetComponent(ptId, j);
      }
    }
    os << "\n";
  }
  os.precision(precision);
  return os.fail() ? 0 : 1;
}

// -----------------------------------------------------------------------------
/// Write dataset points to TetGen .node file
int WriteTetGenNode(const char *fname, vtkPolyData *polydata)
{
  ofstream os(fname);
  if (!os.is_open()) return false;
  return WriteTetGenNode(os, polydata);
}

// -----------------------------------------------------------------------------
/// Write polygonal dataset to TetGen .poly file
int WriteTetGenPoly(const char *fname, vtkPolyData *polydata)
{
  ofstream os(fname);
  if (!os.is_open()) return 0;
  os << "# part 1: nodes\n";
  WriteTetGenNode(os, polydata);
  vtkIdType npts, *pts;
  vtkCellArray *verts  = polydata->GetVerts();
  vtkCellArray *lines  = polydata->GetLines();
  vtkCellArray *polys  = polydata->GetPolys();
  vtkCellArray *strips = polydata->GetStrips();
  int nfacets = 0;
  if (verts ->GetNumberOfCells() > 0) nfacets += 1;
  if (lines ->GetNumberOfCells() > 0) nfacets += 1;
  if (polys ->GetNumberOfCells() > 0) nfacets += 1;
  if (strips->GetNumberOfCells() > 0) nfacets += 1;
  os << "\n# part 2: facets\n";
  os << nfacets << " 0\n";
  if (verts->GetNumberOfCells() > 0) {
    os << "# verts\n";
    os << verts->GetNumberOfCells() << "\n";
    verts->InitTraversal();
    while (verts->GetNextCell(npts, pts)) {
      os << npts << " ";
      for (vtkIdType i = 0; i < npts; ++i) os << " " << (pts[i] + 1);
      os << "\n";
    }
  }
  if (lines->GetNumberOfCells() > 0) {
    os << "# lines\n";
    os << lines->GetNumberOfCells() << "\n";
    lines->InitTraversal();
    while (lines->GetNextCell(npts, pts)) {
      os << npts << " ";
      for (vtkIdType i = 0; i < npts; ++i) os << " " << (pts[i] + 1);
      os << "\n";
    }
  }
  if (polys->GetNumberOfCells() > 0) {
    os << "# polys\n";
    os << polys->GetNumberOfCells() << "\n";
    polys->InitTraversal();
    while (polys->GetNextCell(npts, pts)) {
      os << npts << " ";
      for (vtkIdType i = 0; i < npts; ++i) os << " " << (pts[i] + 1);
      os << "\n";
    }
  }
  if (strips->GetNumberOfCells() > 0) {
    os << "# strips\n";
    os << strips->GetNumberOfCells() << "\n";
    strips->InitTraversal();
    while (strips->GetNextCell(npts, pts)) {
      os << npts << " ";
      for (vtkIdType i = 0; i < npts; ++i) os << " " << (pts[i] + 1);
      os << "\n";
    }
  }
  os << "\n# part 3: hole list\n";
  os << "0\n";
  os << "\n# part 4: region list\n";
  os << "0\n";
  return os.fail() ? 0 : 1;
}

// -----------------------------------------------------------------------------
/// Write polygonal dataset to TetGen .smesh file
int WriteTetGenSMesh(const char *fname, vtkPolyData *polydata)
{
  ofstream os(fname);
  if (!os.is_open()) return 0;
  os << "# part 1: nodes\n";
  WriteTetGenNode(os, polydata);
  vtkIdType npts, *pts;
  vtkCellArray *verts  = polydata->GetVerts();
  vtkCellArray *lines  = polydata->GetLines();
  vtkCellArray *polys  = polydata->GetPolys();
  vtkCellArray *strips = polydata->GetStrips();
  int nfacets = 0;
  if (verts ->GetNumberOfCells() > 0) nfacets += verts ->GetNumberOfCells();
  if (lines ->GetNumberOfCells() > 0) nfacets += lines ->GetNumberOfCells();
  if (polys ->GetNumberOfCells() > 0) nfacets += polys ->GetNumberOfCells();
  if (strips->GetNumberOfCells() > 0) nfacets += strips->GetNumberOfCells();
  os << "\n# part 2: facets\n";
  os << nfacets << " 0\n";
  if (verts->GetNumberOfCells() > 0) {
    verts->InitTraversal();
    while (verts->GetNextCell(npts, pts)) {
      os << npts << " ";
      for (vtkIdType i = 0; i < npts; ++i) os << " " << (pts[i] + 1);
      os << "\n";
    }
  }
  if (lines->GetNumberOfCells() > 0) {
    lines->InitTraversal();
    while (lines->GetNextCell(npts, pts)) {
      os << npts << " ";
      for (vtkIdType i = 0; i < npts; ++i) os << " " << (pts[i] + 1);
      os << "\n";
    }
  }
  if (polys->GetNumberOfCells() > 0) {
    polys->InitTraversal();
    while (polys->GetNextCell(npts, pts)) {
      os << npts << " ";
      for (vtkIdType i = 0; i < npts; ++i) os << " " << (pts[i] + 1);
      os << "\n";
    }
  }
  if (strips->GetNumberOfCells() > 0) {
    strips->InitTraversal();
    while (strips->GetNextCell(npts, pts)) {
      os << npts;
      for (vtkIdType i = 0; i < npts; ++i) os << " " << (pts[i] + 1);
      os << "\n";
    }
  }
  os << "\n# part 3: hole list\n";
  os << "0\n";
  os << "\n# part 4: region list\n";
  os << "0\n";
  return os.fail() ? 0 : 1;
}

// -----------------------------------------------------------------------------
bool WritePointSet(const char *fname, vtkPointSet *pointset, bool compress, bool ascii)
{
  vtkPolyData *polydata = vtkPolyData::SafeDownCast(pointset);
  if (polydata) return WritePolyData(fname, polydata, compress, ascii);
  const std::string ext = Extension(fname);
  int success = 0;
  if (ext.length() == 4 && ext.substr(0, 3) == ".vt" && ext != ".vtk") {
    vtkSmartPointer<vtkXMLDataSetWriter> writer;
    writer = vtkSmartPointer<vtkXMLDataSetWriter>::New();
    SetVTKInput(writer, pointset);
    writer->SetFileName(fname);
    if (compress) writer->SetCompressorTypeToZLib();
    else          writer->SetCompressorTypeToNone();
    success = writer->Write();
  } else {
    vtkSmartPointer<vtkDataSetWriter> writer;
    writer = vtkSmartPointer<vtkDataSetWriter>::New();
    SetVTKInput(writer, pointset);
    writer->SetFileName(fname);
    if (ascii) writer->SetFileTypeToASCII();
    else       writer->SetFileTypeToBinary();
    success = writer->Write();
  }
  return (success == 1);
}

// -----------------------------------------------------------------------------
bool WritePolyData(const char *fname, vtkPolyData *polydata, bool compress, bool ascii)
{
  const std::string ext = Extension(fname);
  int success = 0;
  if (ext == ".vtp") {
    vtkSmartPointer<vtkXMLPolyDataWriter> writer;
    writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    SetVTKInput(writer, polydata);
    writer->SetFileName(fname);
    if (compress) writer->SetCompressorTypeToZLib();
    else          writer->SetCompressorTypeToNone();
    success = writer->Write();
  } else if (ext == ".stl") {
    vtkSmartPointer<vtkSTLWriter> writer;
    writer = vtkSmartPointer<vtkSTLWriter>::New();
    SetVTKInput(writer, polydata);
    if (ascii) writer->SetFileTypeToASCII();
    else       writer->SetFileTypeToBinary();
    writer->SetFileName(fname);
    success = writer->Write();
  } else if (ext == ".ply") {
    vtkSmartPointer<vtkPLYWriter> writer;
    writer = vtkSmartPointer<vtkPLYWriter>::New();
    SetVTKInput(writer, polydata);
    if (ascii) writer->SetFileTypeToASCII();
    else       writer->SetFileTypeToBinary();
    writer->SetFileName(fname);
    success = writer->Write();
  } else if (ext == ".node") {
    success = WriteTetGenNode(fname, polydata);
  } else if (ext == ".poly") {
    success = WriteTetGenPoly(fname, polydata);
  } else if (ext == ".smesh") {
    success = WriteTetGenSMesh(fname, polydata);
  } else {
    vtkSmartPointer<vtkPolyDataWriter> writer;
    writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    SetVTKInput(writer, polydata);
    writer->SetFileName(fname);
    if (ascii) writer->SetFileTypeToASCII();
    else       writer->SetFileTypeToBinary();
    success = writer->Write();
  }
  return (success == 1);
}

// =============================================================================
// Point/cell data
// =============================================================================

// -----------------------------------------------------------------------------
int PolyDataAttributeType(const char *type)
{
  std::string ltype(type);
  transform(ltype.begin(), ltype.end(), ltype.begin(), ::tolower);
  if (ltype == "scalars")     return vtkDataSetAttributes::SCALARS;
  if (ltype == "vectors")     return vtkDataSetAttributes::VECTORS;
  if (ltype == "normals")     return vtkDataSetAttributes::NORMALS;
  if (ltype == "tcoords")     return vtkDataSetAttributes::TCOORDS;
  if (ltype == "tensors")     return vtkDataSetAttributes::TENSORS;
  if (ltype == "globalids")   return vtkDataSetAttributes::GLOBALIDS;
  if (ltype == "pedigreeids") return vtkDataSetAttributes::PEDIGREEIDS;
  if (ltype == "edgeflag")    return vtkDataSetAttributes::EDGEFLAG;
  return -1;
}

// -----------------------------------------------------------------------------
vtkDataArray *GetArrayByCaseInsensitiveName(vtkDataSetAttributes *data, const char *name, int *loc)
{
  std::string lname(name), lower_name;
  std::transform(lname.begin(), lname.end(), lname.begin(), ::tolower);
  for (int i = 0; i < data->GetNumberOfArrays(); ++i) {
    const char *array_name = data->GetArrayName(i);
    if (array_name) {
      lower_name = array_name;
      transform(lower_name.begin(), lower_name.end(), lower_name.begin(), ::tolower);
      if (lower_name == lname) {
        if (loc) *loc = i;
        return data->GetArray(i);
      }
    }
  }
  if (loc) *loc = -1;
  return NULL;
}

// -----------------------------------------------------------------------------
int DeepCopyArrayUsingCaseInsensitiveName(vtkDataSetAttributes *dst, vtkDataSetAttributes *src, const char *name, bool deep)
{
  int loc = -1;
  vtkDataArray *src_array = GetArrayByCaseInsensitiveName(src, name);
  if (src_array) {
    vtkDataArray *dst_array = GetArrayByCaseInsensitiveName(dst, name, &loc);
    if (dst_array) {
      dst_array->DeepCopy(src_array);
    } else {
      vtkSmartPointer<vtkDataArray> copy = src_array->NewInstance();
      copy->DeepCopy(src_array);
      loc = dst->AddArray(copy);
    }
  }
  return loc;
}

// =============================================================================
// Cells
// =============================================================================

// -----------------------------------------------------------------------------
// Attention: vtkGenericCell is not derived from the respective cell types
double ComputeArea(vtkCell *cell)
{
  if (cell->GetCellDimension() < 2) return .0;
  vtkPoints *points = cell->GetPoints();
  switch (cell->GetCellType()) {
    case VTK_TRIANGLE: {
      double p1[3], p2[3], p3[3];
      points->GetPoint(0, p1);
      points->GetPoint(1, p2);
      points->GetPoint(2, p3);
      return vtkTriangle::TriangleArea(p1, p2, p3);
    }
    default: return numeric_limits<double>::quiet_NaN();
  }
}

// -----------------------------------------------------------------------------
// Attention: vtkGenericCell is not derived from the respective cell types
double ComputeVolume(vtkCell *cell)
{
  if (cell->GetCellDimension() < 3) return .0;
  vtkPoints *points = cell->GetPoints();
  switch (cell->GetCellType()) {
    case VTK_TETRA: {
      double p1[3], p2[3], p3[3], p4[3];
      points->GetPoint(0, p1);
      points->GetPoint(1, p2);
      points->GetPoint(2, p3);
      points->GetPoint(3, p4);
      return vtkTetra::ComputeVolume(p1, p2, p3, p4);
    }
    default: return numeric_limits<double>::quiet_NaN();
  }
}

// =============================================================================
// Surface meshes
// =============================================================================

namespace EdgeLengthUtils {

// -----------------------------------------------------------------------------
/// Compute lengths of all edges
struct ComputeEdgeLength
{
  vtkPoints           *_Points;
  const irtkEdgeTable *_EdgeTable;
  double              *_EdgeLength;

  void operator ()(const blocked_range<int> &re) const
  {
    int    ptId1, ptId2, edgeId;
    double p1[3], p2[3];

    irtkEdgeIterator it(*_EdgeTable);
    for (it.InitTraversal(re); (edgeId = it.GetNextEdge(ptId1, ptId2)) != -1;) {
      _Points->GetPoint(ptId1, p1);
      _Points->GetPoint(ptId2, p2);
      _EdgeLength[edgeId] = sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
    }
  }
};

// -----------------------------------------------------------------------------
/// Determine minimum/maximum edge length
struct MinMaxEdgeLength
{
  vtkPoints           *_Points;
  const irtkEdgeTable *_EdgeTable;
  double               _Min;
  double               _Max;

  MinMaxEdgeLength() : _Min(numeric_limits<double>::infinity()), _Max(-_Min) {}

  MinMaxEdgeLength(const MinMaxEdgeLength &other, split)
  :
    _Points(other._Points),
    _EdgeTable(other._EdgeTable),
    _Min(other._Min),
    _Max(other._Max)
  {}

  void join(const MinMaxEdgeLength &other)
  {
    if (other._Min < _Min) _Min = other._Min;
    if (other._Max < _Max) _Max = other._Max;
  }

  void operator ()(const blocked_range<int> &re)
  {
    int    ptId1, ptId2;
    double p1[3], p2[3], d;

    irtkEdgeIterator it(*_EdgeTable);
    for (it.InitTraversal(re); it.GetNextEdge(ptId1, ptId2) != -1;) {
      _Points->GetPoint(ptId1, p1);
      _Points->GetPoint(ptId2, p2);
      d = sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
      if (d < _Min) _Min = d;
      if (d > _Max) _Max = d;
    }
  }
};

// -----------------------------------------------------------------------------
/// Calculate sum of edge lengths
struct SumEdgeLengths
{
  vtkPoints           *_Points;
  const irtkEdgeTable *_EdgeTable;
  double               _Sum;

  SumEdgeLengths() : _Sum(.0) {}

  SumEdgeLengths(const SumEdgeLengths &other, split)
  :
    _Points(other._Points),
    _EdgeTable(other._EdgeTable),
    _Sum(.0)
  {}

  void join(const SumEdgeLengths &other)
  {
    _Sum += other._Sum;
  }

  void operator ()(const blocked_range<int> &re)
  {
    int    ptId1, ptId2;
    double p1[3], p2[3];

    irtkEdgeIterator it(*_EdgeTable);
    for (it.InitTraversal(re); it.GetNextEdge(ptId1, ptId2) != -1;) {
      _Points->GetPoint(ptId1, p1);
      _Points->GetPoint(ptId2, p2);
      _Sum += sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
    }
  }
};


} // namespace EdgeLengthUtils

// -----------------------------------------------------------------------------
double AverageEdgeLength(vtkSmartPointer<vtkPoints> points, const irtkEdgeTable &edgeTable)
{
  if (edgeTable.NumberOfEdges() == 0) return .0;
  EdgeLengthUtils::SumEdgeLengths eval;
  eval._Points    = points;
  eval._EdgeTable = &edgeTable;
  parallel_reduce(blocked_range<int>(0, edgeTable.NumberOfEdges()), eval);
  return eval._Sum / edgeTable.NumberOfEdges();
}

// -----------------------------------------------------------------------------
double AverageEdgeLength(vtkSmartPointer<vtkPointSet> pointset)
{
  irtkEdgeTable edgeTable(pointset);
  return AverageEdgeLength(pointset->GetPoints(), edgeTable);
}

// -----------------------------------------------------------------------------
double RobustAverageEdgeLength(vtkSmartPointer<vtkPoints> points, const irtkEdgeTable &edgeTable)
{
  const int n = edgeTable.NumberOfEdges();
  if (n == 0) return .0;
  double *edgeLength = new double[n];
  EdgeLengthUtils::ComputeEdgeLength eval;
  eval._Points     = points;
  eval._EdgeTable  = &edgeTable;
  eval._EdgeLength = edgeLength;
  parallel_for(blocked_range<int>(0, n), eval);
  double mean = irtk::data::statistic::RobustMean::Calculate(5, n, edgeLength);
  delete[] edgeLength;
  return mean;
}

// -----------------------------------------------------------------------------
double RobustAverageEdgeLength(vtkSmartPointer<vtkPointSet> pointset)
{
  irtkEdgeTable edgeTable(pointset);
  return RobustAverageEdgeLength(pointset->GetPoints(), edgeTable);
}

// -----------------------------------------------------------------------------
double MedianEdgeLength(vtkSmartPointer<vtkPoints> points, const irtkEdgeTable &edgeTable)
{
  const int n = edgeTable.NumberOfEdges();
  if (n == 0) return .0;
  double *edgeLength = new double[n];
  EdgeLengthUtils::ComputeEdgeLength eval;
  eval._Points     = points;
  eval._EdgeTable  = &edgeTable;
  eval._EdgeLength = edgeLength;
  parallel_for(blocked_range<int>(0, n), eval);
  sort(edgeLength, edgeLength + n);
  double median = edgeLength[n / 2];
  delete[] edgeLength;
  return median;
}

// -----------------------------------------------------------------------------
double MedianEdgeLength(vtkSmartPointer<vtkPointSet> pointset)
{
  irtkEdgeTable edgeTable(pointset);
  return MedianEdgeLength(pointset->GetPoints(), edgeTable);
}

// -----------------------------------------------------------------------------
void GetMinMaxEdgeLength(vtkSmartPointer<vtkPoints> points, const irtkEdgeTable &edgeTable, double &min, double &max)
{
  min = max = .0;
  if (edgeTable.NumberOfEdges() == 0) return;
  EdgeLengthUtils::MinMaxEdgeLength eval;
  eval._Points    = points;
  eval._EdgeTable = &edgeTable;
  parallel_reduce(blocked_range<int>(0, edgeTable.NumberOfEdges()), eval);
  min = eval._Min;
  max = eval._Max;
}

// -----------------------------------------------------------------------------
void GetMinMaxEdgeLength(vtkSmartPointer<vtkPointSet> pointset, double &min, double &max)
{
  irtkEdgeTable edgeTable(pointset);
  GetMinMaxEdgeLength(pointset->GetPoints(), edgeTable, min, max);
}

// -----------------------------------------------------------------------------
double MinEdgeLength(vtkSmartPointer<vtkPoints> points, const irtkEdgeTable &edgeTable)
{
  double min, max;
  GetMinMaxEdgeLength(points, edgeTable, min, max);
  return min;
}

// -----------------------------------------------------------------------------
double MinEdgeLength(vtkSmartPointer<vtkPointSet> pointset)
{
  irtkEdgeTable edgeTable(pointset);
  return MinEdgeLength(pointset->GetPoints(), edgeTable);
}

// -----------------------------------------------------------------------------
double MaxEdgeLength(vtkSmartPointer<vtkPoints> points, const irtkEdgeTable &edgeTable)
{
  double min, max;
  GetMinMaxEdgeLength(points, edgeTable, min, max);
  return max;
}

// -----------------------------------------------------------------------------
double MaxEdgeLength(vtkSmartPointer<vtkPointSet> pointset)
{
  irtkEdgeTable edgeTable(pointset);
  return MaxEdgeLength(pointset->GetPoints(), edgeTable);
}

// -----------------------------------------------------------------------------
void EdgeLengthNormalDistribution(vtkSmartPointer<vtkPoints> points, const irtkEdgeTable &edgeTable, double &mean, double &sigma)
{
  mean = sigma = .0;
  if (edgeTable.NumberOfEdges() == 0) return;

  int    ptId1, ptId2, n = 0;
  double p1[3], p2[3], d, delta;

  irtkEdgeIterator it(edgeTable);
  for (it.InitTraversal(); it.GetNextEdge(ptId1, ptId2) != -1;) {
    points->GetPoint(ptId1, p1);
    points->GetPoint(ptId2, p2);
    d = sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
    ++n;
    delta = d - mean;
    mean  += delta / n;
    sigma += delta * (d - mean);
  }

  if (n > 1) sigma /= n - 1;
  sigma = sqrt(sigma);
}

// -----------------------------------------------------------------------------
void EdgeLengthNormalDistribution(vtkSmartPointer<vtkPointSet> pointset, double &mean, double &sigma)
{
  irtkEdgeTable edgeTable(pointset);
  EdgeLengthNormalDistribution(pointset->GetPoints(), edgeTable, mean, sigma);
}

// -----------------------------------------------------------------------------
double GetVolume(vtkSmartPointer<vtkPolyData> surface)
{
  vtkSmartPointer<vtkMassProperties> mp = vtkMassProperties::New();
  SetVTKInput(mp, surface);
  return mp->GetVolume();
}

// -----------------------------------------------------------------------------
bool IsSurfaceMesh(vtkDataSet *dataset)
{
  vtkPolyData *poly = vtkPolyData::SafeDownCast(dataset);
  return poly && poly->GetPolys() && poly->GetPolys()->GetNumberOfCells() > 0 &&
         (!poly->GetLines() || poly->GetLines()->GetNumberOfCells() == 0) &&
         (!poly->GetVerts() || poly->GetVerts()->GetNumberOfCells() == 0);
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> DataSetSurface(vtkSmartPointer<vtkDataSet> dataset, bool passPtIds, bool passCellIds)
{
  vtkSmartPointer<vtkDataSetSurfaceFilter> surface;
  surface = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
  SetVTKInput(surface, dataset);
  surface->SetPassThroughPointIds(passPtIds);
  surface->SetPassThroughCellIds(passCellIds);
  surface->Update();
  return surface->GetOutput();
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> Triangulate(vtkSmartPointer<vtkPolyData> mesh)
{
  vtkSmartPointer<vtkTriangleFilter> filter;
  filter = vtkSmartPointer<vtkTriangleFilter>::New();
  SetVTKInput(filter, mesh);
  filter->PassVertsOff();
  filter->PassLinesOff();
  vtkSmartPointer<vtkCleanPolyData> cleaner;
  cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
  cleaner->PointMergingOff();
  cleaner->ConvertLinesToPointsOn();
  cleaner->ConvertPolysToLinesOn();
  cleaner->ConvertStripsToPolysOn();
  SetVTKConnection(cleaner, filter);
  cleaner->Update();
  cleaner->GetOutput()->SetLines(NULL);
  cleaner->GetOutput()->SetVerts(NULL);
  return cleaner->GetOutput();
}

// -----------------------------------------------------------------------------
// Implements two different filter pipelines to extract the convex hull
vtkSmartPointer<vtkPolyData> ConvexHull(vtkSmartPointer<vtkPointSet> pointset, bool delaunay)
{
  const int levels = 3;

  if (delaunay) {

    // as12312: This pipeline is based on my polydatacortex tool

    // Spatially stratify points to prevent "Unable to factor linear system"
    // warning of vtkDelaunay3D filter due to numerical imprecisions
    vtkSmartPointer<vtkMaskPoints> stratify;
    stratify = vtkSmartPointer<vtkMaskPoints>::New();
    stratify->RandomModeOn();
    stratify->SetRandomModeType(2);
    stratify->SetMaximumNumberOfPoints(.75 * pointset->GetNumberOfPoints());
    SetVTKInput(stratify, pointset);

    // Get convex hull of largest component
    vtkSmartPointer<vtkHull> hull;
    hull = vtkSmartPointer<vtkHull>::New();
    hull->AddRecursiveSpherePlanes(levels);
    SetVTKConnection(hull, stratify);

    // Compute Delaunay triangulation
    vtkSmartPointer<vtkDelaunay3D> delaunay;
    delaunay = vtkSmartPointer<vtkDelaunay3D>::New();
    SetVTKConnection(delaunay, hull);

    // Construct surface mesh
    vtkSmartPointer<vtkDataSetSurfaceFilter> mesher;
    mesher = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    SetVTKConnection(mesher, delaunay);

    mesher->Update();
    return mesher->GetOutput();

  } else {

    // as12312: Copied the following from Paul Aljabar's polydatacurvature tool

    // Get convex hull of largest component
    vtkSmartPointer<vtkHull> hull;
    hull = vtkSmartPointer<vtkHull>::New();
    hull->AddRecursiveSpherePlanes(levels);
    SetVTKInput(hull, pointset);

    // Clean polydata
    vtkSmartPointer<vtkCleanPolyData> cleaner;
    cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
    cleaner->ConvertLinesToPointsOff();
    cleaner->ConvertPolysToLinesOff();
    cleaner->ConvertStripsToPolysOff();
    cleaner->PointMergingOn();
    cleaner->ToleranceIsAbsoluteOn();
    cleaner->SetAbsoluteTolerance(0);
    SetVTKConnection(cleaner, hull);

    // Triangulate convex hull
    cleaner->Update();
    return Triangulate(cleaner->GetOutput());
  }
}

// -----------------------------------------------------------------------------
bool IsTriangularMesh(vtkDataSet *input)
{
  if (vtkPointSet::SafeDownCast(input) == NULL) return false;
  for (vtkIdType cellId = 0; cellId < input->GetNumberOfCells(); ++cellId) {
    int type = input->GetCellType(cellId);
    if (type != VTK_EMPTY_CELL && type != VTK_TRIANGLE) return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
bool IsTetrahedralMesh(vtkDataSet *input)
{
  if (vtkPointSet::SafeDownCast(input) == NULL) return false;
  for (vtkIdType cellId = 0; cellId < input->GetNumberOfCells(); ++cellId) {
    int type = input->GetCellType(cellId);
    if (type != VTK_EMPTY_CELL && type != VTK_TETRA) return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPointSet> Tetrahedralize(vtkSmartPointer<vtkPointSet> input)
{
  vtkSmartPointer<vtkPointSet> mesh;

  if (IsTetrahedralMesh(input)) {
    mesh = input->NewInstance();
    mesh->ShallowCopy(input);
  } else {
    // TODO: Use TetGen library to tetrahedralize interior of input PLC
    cerr << "irtkPolyDataUtils::Tetrahedralize: Not implemented, use TetGen command instead" << endl;
    exit(1);
  }

  return mesh;
}


} } // namespace irtk::polydata
