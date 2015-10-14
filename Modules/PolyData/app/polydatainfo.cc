// =============================================================================
// Project: Image Registration Toolkit (IRTK)
// Package: PolyData
//
// Copyright (c) 2015 Imperial College London
// Copyright (c) 2015 Andreas Schuh
// =============================================================================

#ifdef HAS_VTK
#include <irtkCommon.h>

// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << "usage: " << name << " <input> [options]" << endl;
  cout << endl;
  cout << "Print information about polygonal dataset." << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -edgelength              Report edge length statistics. (default: off)" << endl;
  cout << "  -self-intersections      Check for self-intersections. (default: off)" << endl;
  cout << "  -output-surface <file>   Write surface mesh to specified file. (default: none)" << endl;
  PrintStandardOptions(cout);
}

// =============================================================================
// Includes
// =============================================================================

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkDataArray.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkGenericCell.h>
#include <vtkOctreePointLocator.h>
#include <vtkMath.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>

#include <vtkDataSetSurfaceFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkIntersectionPolyDataFilter.h>

#include <irtkEdgeTable.h>
#include <irtkTriangle.h>
#include <irtkPolyhedron.h>
#include <irtkPolyDataUtils.h>

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
/// Convert input dataset to polygonal data, i.e., boundary surface
vtkSmartPointer<vtkPolyData> ConvertToPolyData(vtkDataSet *dataset)
{
  vtkSmartPointer<vtkPolyData> polydata = vtkPolyData::SafeDownCast(dataset);
  if (polydata) return polydata;
  vtkSmartPointer<vtkDataSetSurfaceFilter> surface;
  surface = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
  SetVTKInput(surface, dataset);
  surface->PassThroughCellIdsOff();
  surface->PassThroughPointIdsOff();
  surface->Update();
  return surface->GetOutput();
}

// -----------------------------------------------------------------------------
/// Clean polydata
vtkSmartPointer<vtkPolyData> CleanPolyData(vtkSmartPointer<vtkPolyData> polydata)
{
  vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
  cleaner->SetAbsoluteTolerance(.0);
  cleaner->PointMergingOn();
  cleaner->ConvertLinesToPointsOn();
  cleaner->ConvertPolysToLinesOn();
  cleaner->ConvertStripsToPolysOn();
  SetVTKInput(cleaner, polydata);
  cleaner->Update();
  return cleaner->GetOutput();
}

// -----------------------------------------------------------------------------
/// Get number of edges
int NumberOfEdges(vtkDataSet *polydata)
{
  irtkEdgeTable edgeTable(polydata);
  return edgeTable.NumberOfEdges();
}

// -----------------------------------------------------------------------------
/// Determine number of connected components
int NumberOfComponents(vtkPolyData *polydata)
{
  vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivity;
  connectivity = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
  connectivity->SetExtractionModeToAllRegions();
  SetVTKInput(connectivity, polydata);
  connectivity->Update();
  return connectivity->GetNumberOfExtractedRegions();
}

// -----------------------------------------------------------------------------
/// Get center point and radius of bounding sphere
inline double GetBoundingSphereRadius(vtkGenericCell *cell, double c[3])
{
  double p[3];
  // Get center of bounding sphere
  c[0] = .0, c[1] = .0, c[2] = .0;
  double *weights = new double[cell->GetNumberOfPoints()];
  int subId = cell->GetParametricCenter(p);
  cell->EvaluateLocation(subId, p, c, weights);
  delete[] weights;
  // Get radius of bounding sphere
  double r = .0;
  vtkPoints *points = cell->GetPoints();
  for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i) {
    points->GetPoint(i, p);
    r = max(r, vtkMath::Distance2BetweenPoints(c, p));
  }
  return sqrt(r);
}

// -----------------------------------------------------------------------------
/// Auxiliary functor which counts the number of self-intersections of a triangulated surface mesh
struct CountTriangleTriangleIntersections
{
  vtkPolyData             *_DataSet;
  vtkAbstractPointLocator *_PointLocator;
  vtkDataArray            *_Mask;
  int                      _NumberOfIntersections;

  CountTriangleTriangleIntersections() : _NumberOfIntersections(0) {}

  CountTriangleTriangleIntersections(const CountTriangleTriangleIntersections &other, split)
  :
    _DataSet(other._DataSet),
    _PointLocator(other._PointLocator),
    _Mask(other._Mask),
    _NumberOfIntersections(0)
  {}

  void join(const CountTriangleTriangleIntersections &other)
  {
    _NumberOfIntersections += other._NumberOfIntersections;
  }

  void operator ()(const blocked_range<vtkIdType> &re)
  {
    vtkSmartPointer<vtkGenericCell> cell1   = vtkSmartPointer<vtkGenericCell>::New();
    vtkSmartPointer<vtkGenericCell> cell2   = vtkSmartPointer<vtkGenericCell>::New();
    vtkSmartPointer<vtkIdList>      ptIds   = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList>      cellIds = vtkSmartPointer<vtkIdList>::New();

    unsigned short ncells;
    vtkIdType *cells, cellId2;
    double r1, c1[3], v1[3], v2[3], v3[3], u1[3], u2[3], u3[3];

    for (vtkIdType cellId1 = re.begin(); cellId1 != re.end(); ++cellId1) {
      _DataSet->GetCell(cellId1, cell1);
      if (cell1->GetNumberOfPoints() != 3) continue;
      _DataSet->GetPoint(cell1->GetPointId(0), v1);
      _DataSet->GetPoint(cell1->GetPointId(1), v2);
      _DataSet->GetPoint(cell1->GetPointId(2), v3);
      r1 = GetBoundingSphereRadius(cell1, c1);
      _PointLocator->FindPointsWithinRadius(3.0 * r1, c1, ptIds);
      cellIds->Reset();
      for (vtkIdType i = 0; i < ptIds->GetNumberOfIds(); ++i) {
        _DataSet->GetPointCells(ptIds->GetId(i), ncells, cells);
        for (unsigned short j = 0; j < ncells; ++j) {
          if (cells[j] != cellId1) cellIds->InsertUniqueId(cells[j]);
        }
      }
      for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {
        cellId2 = cellIds->GetId(i);
        _DataSet->GetCell(cellId2, cell2);
        if (cell2->GetNumberOfPoints() != 3) continue;
        if (cell1->GetPointIds()->IsId(cell2->GetPointId(0)) >= 0) continue;
        if (cell1->GetPointIds()->IsId(cell2->GetPointId(1)) >= 0) continue;
        if (cell1->GetPointIds()->IsId(cell2->GetPointId(2)) >= 0) continue;
        _DataSet->GetPoint(cell2->GetPointId(0), u1);
        _DataSet->GetPoint(cell2->GetPointId(1), u2);
        _DataSet->GetPoint(cell2->GetPointId(2), u3);
        if (irtkTriangle::TriangleTriangleIntersection(v1, v2, v3, u1, u2, u3) != 0) {
          _Mask->SetComponent(cellId1, 0, 1);
          _Mask->SetComponent(cellId2, 0, 1);
          ++_NumberOfIntersections;
        }
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Count number of self-intersections of triangulated surface mesh
int NumberOfTriangleTriangleIntersections(vtkPolyData *polydata)
{
  vtkSmartPointer<vtkDataArray> mask = vtkSmartPointer<vtkUnsignedCharArray>::New();
  mask->SetName("TriangleTriangleIntersection");
  mask->SetNumberOfComponents(1);
  mask->SetNumberOfTuples(polydata->GetNumberOfCells());
  for (vtkIdType cellId = 0; cellId < polydata->GetNumberOfCells(); ++cellId) {
    mask->SetComponent(cellId, 0, 0);
  }
  polydata->GetCellData()->AddArray(mask);
  vtkSmartPointer<vtkOctreePointLocator> octree = vtkSmartPointer<vtkOctreePointLocator>::New();
  octree->SetDataSet(polydata);
  octree->BuildLocator();
  CountTriangleTriangleIntersections count;
  count._DataSet      = polydata;
  count._PointLocator = octree;
  count._Mask         = mask;
  parallel_reduce(blocked_range<vtkIdType>(0, polydata->GetNumberOfCells()), count);
  return count._NumberOfIntersections;
}

// -----------------------------------------------------------------------------
/// Auxiliary functor which counts the number of vertices inside the polyhedron
struct CountVerticesInsidePolyhedron
{
  irtkPolyhedron *_Polyhedron;
  vtkDataArray   *_Mask;
  int             _Num;

  CountVerticesInsidePolyhedron() : _Num(0) {}
  CountVerticesInsidePolyhedron(const CountVerticesInsidePolyhedron &other, split)
  :
    _Polyhedron(other._Polyhedron), _Mask(other._Mask), _Num(0)
  {}

  void join(const CountVerticesInsidePolyhedron &other)
  {
    _Num += other._Num;
  }

  void operator ()(const blocked_range<int> &re)
  {
    double p[3];
    for (int ptId = re.begin(); ptId != re.end(); ++ptId) {
      _Polyhedron->GetPoint(ptId, p);
      if (_Polyhedron->IsInside(p)) {
        _Mask->SetComponent(ptId, 0, 1);
        ++_Num;
      } else {
        _Mask->SetComponent(ptId, 0, 0);
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Point-inside-polyhedron test
int NumberOfVerticesInsidePolyhedron(vtkPolyData *polydata)
{
  vtkSmartPointer<vtkDataArray> mask = vtkSmartPointer<vtkUnsignedCharArray>::New();
  mask->SetName("VertexInsidePolyhedron");
  mask->SetNumberOfComponents(1);
  mask->SetNumberOfTuples(polydata->GetNumberOfPoints());
  polydata->GetPointData()->AddArray(mask);
  irtkPolyhedron polyhedron(polydata);
  CountVerticesInsidePolyhedron count;
  count._Polyhedron = &polyhedron;
  count._Mask       = mask;
  parallel_reduce(blocked_range<int>(0, polyhedron.NumberOfPoints()), count);
  return count._Num;
}

// -----------------------------------------------------------------------------
/// Compute volume of cells
struct ComputeCellVolumes
{
  vtkPointSet  *_PointSet;
  vtkDataArray *_Volume;

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    double volume;
    vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
    for (vtkIdType cellId = re.begin(); cellId != re.end(); ++cellId) {
      _PointSet->GetCell(cellId, cell);
      volume = ComputeVolume(cell);
      _Volume->SetComponent(cellId, 0, volume);
    }
  }

  static void Run(vtkPointSet *pointset)
  {
    vtkSmartPointer<vtkDataArray> volume;
    volume = vtkSmartPointer<vtkFloatArray>::New();
    volume->SetName("Volume");
    volume->SetNumberOfComponents(1);
    volume->SetNumberOfTuples(pointset->GetNumberOfCells());
    pointset->GetCellData()->AddArray(volume);
    ComputeCellVolumes body;
    body._PointSet = pointset;
    body._Volume   = volume;
    parallel_for(blocked_range<vtkIdType>(0, pointset->GetNumberOfCells()), body);
  }
};

// -----------------------------------------------------------------------------
const char *DataObjectTypeString(int type)
{
  switch (type) {
    case VTK_IMAGE_DATA:        return "Image";
    case VTK_POLY_DATA:         return "Polydata";
    case VTK_UNSTRUCTURED_GRID: return "Unstructured grid";
    case VTK_STRUCTURED_GRID:   return "Structured grid";
    case VTK_STRUCTURED_POINTS: return "Structured points";
    default: return "unknown";
  }
}

// -----------------------------------------------------------------------------
const char *CellTypeString(int type)
{
  switch (type) {
    case VTK_EMPTY_CELL:      return "empty cells";
    case VTK_VERTEX:          return "vertices";
    case VTK_POLY_VERTEX:     return "poly-vertices";
    case VTK_LINE:            return "lines";
    case VTK_TRIANGLE:        return "triangles";
    case VTK_TRIANGLE_STRIP:  return "triangle strips";
    case VTK_POLYGON:         return "polygons";
    case VTK_QUAD:            return "quads";
    case VTK_TETRA:           return "tetrahedra";
    case VTK_VOXEL:           return "voxels";
    case VTK_HEXAHEDRON:      return "hexadrons";
    case VTK_HEXAGONAL_PRISM: return "hexagonal prisms";
    case VTK_PYRAMID:         return "pyramids";
    default:                  return "unknown cells";
  }
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  EXPECTS_POSARGS(1);

  // Parse options
  bool check_intersections = false;
  bool compute_volumes     = false;
  bool report_edge_lengths = false;
  const char *output_pointset_name = NULL;
  const char *output_surface_name  = NULL;

  for (ALL_OPTIONS) {
    if (OPTION("-self-intersections")) check_intersections = true;
    else if (OPTION("-edgelength")) report_edge_lengths = true;
    else if (OPTION("-vol") || OPTION("-volumes")) compute_volumes = true;
    else if (OPTION("-o") || OPTION("-out") || OPTION("-output")) {
      output_pointset_name = ARGUMENT;
    }
    else if (OPTION("-s") || OPTION("-output-surface")) {
      output_surface_name = ARGUMENT;
    }
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (verbose) {
    check_intersections = true;
    report_edge_lengths = true;
  }

  // Read input dataset
  vtkSmartPointer<vtkPointSet> pointset = ReadPointSet(POSARG(1));
  vtkSmartPointer<vtkPolyData> polydata = ConvertToPolyData(pointset);
  polydata->BuildLinks();

  // Cell types
  cout << DataObjectTypeString(pointset->GetDataObjectType()) << ":" << endl;
  cout << "  No. of points:             " << pointset->GetNumberOfPoints() << endl;
  if (verbose) {
    cout << "  No. of edges:              " << NumberOfEdges(pointset) << endl;
  }
  cout << "  No. of cells:              " << pointset->GetNumberOfCells() << endl;
  if (pointset->GetNumberOfCells() > 0) {
    map<int, int> cellTypeHist;
    map<int, int>::iterator it;
    for (vtkIdType cellId = 0; cellId < pointset->GetNumberOfCells(); ++cellId) {
      int type = pointset->GetCellType(cellId);
      it = cellTypeHist.find(type);
      if (it == cellTypeHist.end()) cellTypeHist[type] = 1;
      else ++it->second;
    }
    for (it = cellTypeHist.begin(); it != cellTypeHist.end(); ++it) {
      string s = CellTypeString(it->first);
      cout << "    No. of " << setw(18) << left << (s + ": ") << it->second << endl;
    }
    cout << "  Maximum cell size:         " << pointset->GetMaxCellSize() << endl;
  }

  // Attributes
  if (pointset->GetPointData()->GetNumberOfArrays() > 0) {
  	cout << "\nPoint data arrays:" << endl;
    for (int i = 0; i < pointset->GetPointData()->GetNumberOfArrays(); ++i) {
      vtkDataArray *scalars = pointset->GetPointData()->GetArray(i);
      printf(" %2d %-24s (dim: %2d)\n", i, scalars->GetName(), scalars->GetNumberOfComponents());
    }
  }
  if (pointset->GetCellData()->GetNumberOfArrays() > 0) {
    cout << "\nCell data arrays:" << endl;
    for (int i = 0; i < pointset->GetCellData()->GetNumberOfArrays(); ++i) {
      vtkDataArray *scalars = pointset->GetCellData()->GetArray(i);
      printf(" %2d %-24s (dim: %2d)\n", i, scalars->GetName(), scalars->GetNumberOfComponents());
    }
  }

  if (pointset->GetNumberOfCells() > 0) {

    // Compute cell volumes
    if (compute_volumes) ComputeCellVolumes::Run(pointset);

    // Discard degenerate cells and unused points
    polydata = CleanPolyData(polydata);
    polydata->SetVerts(NULL);
    polydata->SetLines(NULL);
    polydata = CleanPolyData(polydata);
    polydata->BuildLinks();

    // Compute edge table
    irtkEdgeTable edgeTable(polydata);

    // Euler characteristic / Genus
    const int    noOfVerts = polydata->GetNumberOfPoints();
    const int    noOfFaces = polydata->GetNumberOfCells();
    const int    noOfComps = NumberOfComponents(polydata);
    const int    noOfEdges = edgeTable.NumberOfEdges();
    const int    euler     = noOfVerts - noOfEdges + noOfFaces;
    const double genus     = 0.5 * (2 * noOfComps - euler);

    cout << endl;
    cout << "Surface mesh:" << endl;
    cout << "  V " << noOfVerts << endl;
    cout << "  E " << noOfEdges << endl;
    cout << "  F " << noOfFaces << endl;
    cout << "  C " << noOfComps << endl;
    cout << endl;
    cout << "  Euler characteristic / Genus (V - E + F = 2C - 2g)" << endl;
    cout << "    Euler: " << euler << endl;
    cout << "    Genus: " << genus << endl;

    if (report_edge_lengths) {
      double min_length, max_length, mean, sigma;
      EdgeLengthNormalDistribution(polydata->GetPoints(), edgeTable, mean, sigma);
      GetMinMaxEdgeLength(polydata->GetPoints(), edgeTable, min_length, max_length);
      cout << endl;
      cout << "  Average edge length: " << mean << endl;
      cout << "  Edge length StDev:   " << sigma << endl;
      cout << "  Minimum edge length: " << min_length << endl;
      cout << "  Maximum edge length: " << max_length << endl;
    }

    // Self-intersections
    if (check_intersections) {
      cout << endl;
      cout << "  No. of triangle/triangle intersections: ", cout.flush();
      cout << NumberOfTriangleTriangleIntersections(polydata) << endl;
  //    cout << "  No. of vertices inside polyhedron: ", cout.flush();
  //    cout << NumberOfVerticesInsidePolyhedron(polydata) << endl;
    }
  }

  // Write output point set (possibly with additional attributes)
  if (output_pointset_name) WritePointSet(output_pointset_name, pointset);

  // Write output surface mesh (possibly with additional attributes)
  if (output_surface_name && polydata->GetNumberOfCells() > 0) {
    WritePointSet(output_surface_name, polydata);
  }

  return 0;
}

#else // HAS_VTK
int main(int, char *argv[])
{
  cerr << argv[0] << " needs to be compiled with the VTK" << endl;
  exit(1);
}
#endif // HAS_VTK
