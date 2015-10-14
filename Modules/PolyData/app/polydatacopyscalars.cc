// =============================================================================
// Project: Image Registration Toolkit (IRTK)
// Package: PolyData
//
// Copyright (c) 2015 Andreas Schuh
// Copyright (c) 2015 Imperial College London
// =============================================================================

#include <irtkCommon.h>

// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << "usage: " << name << " <source> <target> [<output>] [options]" << endl;
  cout << endl;
  cout << "Copies point and/or cell data from a source data set to a target data set." << endl;
  cout << "If the points sets have differing number of points or cells, respectively," << endl;
  cout << "zero entries are either added to the target arrays or only the first n" << endl;
  cout << "tuples of the source arrays are copied." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  source   Point set from which to copy the point/cell data." << endl;
  cout << "  target   Point set to which to add the point/cell data." << endl;
  cout << "  output   Output point set with copied point/cell data." << endl;
  cout << "           If not specified, the target file is overwritten." << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -points [<target_name>]                   Add points of source as point data of target. (default: off)" << endl;
  cout << "  -pointdata <name|index> [<target_name>]   Name or index of point data to copy. (default: all)" << endl;
  cout << "  -celldata  <name|index> [<target_name>]   Name or index of cell  data to copy. (default: all)" << endl;
  cout << "  -case-insensitive                         Lookup source arrays by case insensitive name. (default: on)" << endl;
  cout << "  -case-sensitive                           Lookup source arrays by case sensitive name. (default: off)" << endl;
  PrintCommonOptions(cout);
}

// =============================================================================
// Includes
// =============================================================================

#ifdef HAS_VTK

#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  vtkDataArray *array;

  // Positional arguments
  REQUIRES_POSARGS(2);

  if (NUM_POSARGS > 3) {
    PrintHelp(EXECNAME);
    exit(1);
  }

  const char *source_name = POSARG(1);
  const char *target_name = POSARG(2);
  const char *output_name = (NUM_POSARGS == 3 ? POSARG(3) : target_name);

  // Optional arguments
  const char *points_name = NULL;
  vector<pair<const char *, const char *> > pdname;
  vector<pair<int,          const char *> > pdindex;
  vector<pair<const char *, const char *> > cdname;
  vector<pair<int,          const char *> > cdindex;
  bool case_sensitive = false;

  for (ALL_OPTIONS) {
    if (OPTION("-points")) {
      if (HAS_ARGUMENT) points_name = ARGUMENT;
      else              points_name = "Coords";
    }
    else if (OPTION("-pointdata")) {
      int index;
      const char *arg  = ARGUMENT;
      const char *name = (HAS_ARGUMENT ? ARGUMENT : NULL);
      if (FromString(arg, index)) {
        pdindex.push_back(make_pair(index, name));
      } else {
        pdname.push_back(make_pair(arg, name));
      }
    }
    else if (OPTION("-celldata")) {
      int index;
      const char *arg  = ARGUMENT;
      const char *name = (HAS_ARGUMENT ? ARGUMENT : NULL);
      if (FromString(arg, index)) {
        cdindex.push_back(make_pair(index, name));
      } else {
        cdname.push_back(make_pair(arg, name));
      }
    }
    else if (OPTION("-case-sensitive"))   case_sensitive = true;
    else if (OPTION("-case-insensitive")) case_sensitive = false;
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // Read input data sets
  vtkSmartPointer<vtkPointSet> source = ReadPointSet(source_name);
  vtkSmartPointer<vtkPointSet> target = ReadPointSet(target_name);

  // Add point data
  vtkPointData *sourcePD = source->GetPointData();
  vtkPointData *targetPD = target->GetPointData();
  if (pdname.empty() && pdindex.empty()) {
    pdindex.resize(sourcePD->GetNumberOfArrays());
    for (int i = 0; i < sourcePD->GetNumberOfArrays(); ++i) {
      pdindex[i].first = i;
    }
  }
  for (size_t i = 0; i < pdname.size(); ++i) {
    if (case_sensitive) {
      array = sourcePD->GetArray(pdname[i].first);
    } else {
      array = GetArrayByCaseInsensitiveName(sourcePD, pdname[i].first);
    }
    if (array == NULL) {
      cerr << "Source has no point data array named " << pdname[i].first << endl;
      exit(1);
    }
    vtkSmartPointer<vtkDataArray> copy = array->NewInstance();
    copy->SetName(array->GetName());
    copy->SetNumberOfComponents(array->GetNumberOfComponents());
    copy->SetNumberOfTuples(target->GetNumberOfPoints());
    vtkIdType end = min(array->GetNumberOfTuples(), copy->GetNumberOfTuples());
    for (vtkIdType ptId = 0; ptId < end; ++ptId) {
      copy->SetTuple(ptId, array->GetTuple(ptId));
    }
    for (vtkIdType ptId = end; ptId < copy->GetNumberOfTuples(); ++ptId) {
      for (int j = 0; j < copy->GetNumberOfComponents(); ++j) {
        copy->SetComponent(ptId, j, .0);
      }
    }
    if (pdname[i].second) copy->SetName(pdname[i].second);
    targetPD->AddArray(copy);
  }
  for (size_t i = 0; i < pdindex.size(); ++i) {
    if (pdindex[i].first < 0 || pdindex[i].first > sourcePD->GetNumberOfArrays()) {
      cerr << "Source has no point data array with index " << pdindex[i].first << endl;
      exit(1);
    }
    array = sourcePD->GetArray(pdindex[i].first);
    vtkSmartPointer<vtkDataArray> copy = array->NewInstance();
    copy->SetName(array->GetName());
    copy->SetNumberOfComponents(array->GetNumberOfComponents());
    copy->SetNumberOfTuples(target->GetNumberOfPoints());
    vtkIdType end = min(array->GetNumberOfTuples(), copy->GetNumberOfTuples());
    for (vtkIdType ptId = 0; ptId < end; ++ptId) {
      copy->SetTuple(ptId, array->GetTuple(ptId));
    }
    for (vtkIdType ptId = end; ptId < copy->GetNumberOfTuples(); ++ptId) {
      for (int j = 0; j < copy->GetNumberOfComponents(); ++j) {
        copy->SetComponent(ptId, j, .0);
      }
    }
    if (pdindex[i].second) copy->SetName(pdindex[i].second);
    targetPD->AddArray(copy);
  }

  // Add cell data
  vtkCellData *sourceCD = source->GetCellData();
  vtkCellData *targetCD = target->GetCellData();
  if (cdname.empty() && cdindex.empty()) {
    cdindex.resize(sourceCD->GetNumberOfArrays());
    for (int i = 0; i < sourceCD->GetNumberOfArrays(); ++i) {
      cdindex[i].first = i;
    }
  }
  for (size_t i = 0; i < cdname.size(); ++i) {
    if (case_sensitive) {
      array = sourceCD->GetArray(cdname[i].first);
    } else {
      array = GetArrayByCaseInsensitiveName(sourceCD, cdname[i].first);
    }
    if (array == NULL) {
      cerr << "Source has no cell data array named " << cdname[i].first << endl;
      exit(1);
    }
    vtkSmartPointer<vtkDataArray> copy = array->NewInstance();
    copy->SetName(array->GetName());
    copy->SetNumberOfComponents(array->GetNumberOfComponents());
    copy->SetNumberOfTuples(target->GetNumberOfCells());
    vtkIdType end = min(array->GetNumberOfTuples(), copy->GetNumberOfTuples());
    for (vtkIdType ptId = 0; ptId < end; ++ptId) {
      copy->SetTuple(ptId, array->GetTuple(ptId));
    }
    for (vtkIdType ptId = end; ptId < copy->GetNumberOfTuples(); ++ptId) {
      for (int j = 0; j < copy->GetNumberOfComponents(); ++j) {
        copy->SetComponent(ptId, j, .0);
      }
    }
    if (cdname[i].second) copy->SetName(cdname[i].second);
    targetCD->AddArray(copy);
  }
  for (size_t i = 0; i < cdindex.size(); ++i) {
    if (cdindex[i].first < 0 || cdindex[i].first > sourceCD->GetNumberOfArrays()) {
      cerr << "Source has no cell data array with index " << cdindex[i].first << endl;
      exit(1);
    }
    array = sourceCD->GetArray(cdindex[i].first);
    vtkSmartPointer<vtkDataArray> copy = array->NewInstance();
    copy->SetName(array->GetName());
    copy->SetNumberOfComponents(array->GetNumberOfComponents());
    copy->SetNumberOfTuples(target->GetNumberOfCells());
    vtkIdType end = min(array->GetNumberOfTuples(), copy->GetNumberOfTuples());
    for (vtkIdType ptId = 0; ptId < end; ++ptId) {
      copy->SetTuple(ptId, array->GetTuple(ptId));
    }
    for (vtkIdType ptId = end; ptId < copy->GetNumberOfTuples(); ++ptId) {
      for (int j = 0; j < copy->GetNumberOfComponents(); ++j) {
        copy->SetComponent(ptId, j, .0);
      }
    }
    if (cdindex[i].second) copy->SetName(cdindex[i].second);
    targetCD->AddArray(copy);
  }

  // Add source point coordinates
  if (points_name) {
    vtkSmartPointer<vtkDataArray> coords = vtkSmartPointer<vtkFloatArray>::New();
    coords->SetName(points_name);
    coords->SetNumberOfComponents(3);
    coords->SetNumberOfTuples(target->GetNumberOfPoints());
    vtkIdType end = min(source->GetNumberOfPoints(), target->GetNumberOfPoints());
    for (vtkIdType ptId = 0; ptId < end; ++ptId) {
      coords->SetTuple(ptId, source->GetPoint(ptId));
    }
    const double zero[3] = {.0, .0, .0};
    for (vtkIdType ptId = end; ptId < target->GetNumberOfPoints(); ++ptId) {
      coords->SetTuple(ptId, zero);
    }
    targetPD->AddArray(coords);
  }

  // Write resulting data set
  WritePointSet(output_name, target);
}

#else // HAS_VTK
int main(int argc, char *argv[])
{
  PrintHelp(EXECNAME);
  cout << endl;
  cerr << "Error: This program cannot be used as it was not compiled with the VTK library!" << endl;
  return 1;
}
#endif
