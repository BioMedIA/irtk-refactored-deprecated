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
  cout << "usage: " << name << " <input> <output> [options]" << endl;
  cout << endl;
  cout << "Convert dataset from one (file) format to another." << endl;
  cout << endl;
  cout << "The current implementation can only convert between different" << endl;
  cout << "polygonal data file formats based on the file name extension." << endl;
  cout << "Besides the common formats supported by VTK, it can also write" << endl;
  cout << "a Piecewise Linear Complex (PLC) B-Rep description in the TetGen" << endl;
  cout << "input formats .poly and .smesh." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  input   Input  dataset file (.vtk, .vtp, .obj, .stl, .ply)." << endl;
  cout << "  output  Output dataset file (.vtk, .vtp, .stl, .ply, .node, .poly, .smesh)." << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -nocelldata    Do not write cell  data to output file. (default: off)" << endl;
  cout << "  -nopointdata   Do not write point data to output file. (default: off)" << endl;
  cout << "  -ascii         Write legacy VTK files encoded in ASCII. (default: off)" << endl;
  cout << "  -binary        Write legacy VTK files in binary form. (default: on)" << endl;
  cout << "  -compress      Compress XML VTK files. (default: on)" << endl;
  cout << "  -nocompress    Do not compress XML VTK files. (default: off)" << endl;
  PrintStandardOptions(cout);
}

// =============================================================================
// Includes
// =============================================================================
#ifdef HAS_VTK

#include <vtkPointData.h>
#include <vtkCellData.h>

#endif // HAS_VTK
// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
#ifdef HAS_VTK

  // Positional arguments
  EXPECTS_POSARGS(2);

  const char *input_name  = POSARG(1);
  const char *output_name = POSARG(2);

  // Optional arguments
  bool pointdata = true;
  bool celldata  = true;
  bool ascii     = false;
  bool compress  = true;

  for (ALL_OPTIONS) {
    if      (OPTION("-nopointdata")) pointdata = false;
    else if (OPTION("-nocelldata"))  celldata  = false;
    else if (OPTION("-ascii"))       ascii = true;
    else if (OPTION("-binary"))      ascii = false;
    else if (OPTION("-compress"))    compress = true;
    else if (OPTION("-nocompress"))  compress = false;
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // Read data set
  vtkSmartPointer<vtkPointSet> dataset = ReadPointSet(input_name);

  // Reset point/cell data
  if (!pointdata) dataset->GetPointData()->Initialize();
  if (!celldata ) dataset->GetCellData ()->Initialize();

  // Write data set
  WritePointSet(output_name, dataset, compress, ascii);
  return 0;

#else // HAS_VTK

  PrintHelp(EXECNAME);
  cout << endl;
  cerr << "Error: This program cannot be used as it was not compiled with the VTK library!" << endl;
  return 1;

#endif // HAS_VTK
}

