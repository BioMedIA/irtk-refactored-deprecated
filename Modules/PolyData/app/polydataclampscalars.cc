// =============================================================================
// Project: Image Registration Toolkit (IRTK)
// Package: PolyData
//
// Copyright (c) 2015 Imperial College London
// Copyright (c) 2015 Andreas Schuh
// Copyright (c) Paul Aljabar
// =============================================================================

#include <irtkCommon.h>

// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << "Usage: " << name << " <in> <out> <lo> <hi> [options]" << endl;
  cout << endl;
  cout << "Clamp scalar point data of input point set to given value range." << endl;
  cout << "By default, the <lo> and <hi> values are the percentiles. Use the" << endl;
  cout << "-minmax option to use these directly as lower and upper thresholds." << endl;
  cout << endl;
  cout << "Options: " << endl;
  cout << "  -name         Case insensitive name of point data array to clamp. (default: SCALARS)" << endl;
  cout << "  -minmax       Use input thresholds directly. (default: off)" << endl;
  cout << "  -percentile   Compute percentile thresholds. (default: on)" << endl;
  PrintStandardOptions(cout);
}

// =============================================================================
// Includes
// =============================================================================

#ifdef HAS_VTK
#  include <vtkDataArray.h>
#  include <vtkPointData.h>
#  include <irtkPolyDataUtils.h>
#endif

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  EXPECTS_POSARGS(4);

  const char *input_name  = POSARG(1);
  const char *output_name = POSARG(2);

  double minVal, maxVal;
  if (!FromString(POSARG(3), minVal)) {
    cerr << "Invalid <lo> argument: " << POSARG(3) << endl;
    exit(1);
  }
  if (!FromString(POSARG(4), maxVal)) {
    cerr << "Invalid <hi> argument: " << POSARG(4) << endl;
    exit(1);
  }

  const char *scalar_name = NULL;  // Name of scalar point data array
  bool        minmax      = false; // false: percentiles, true: thresholds

  for (ALL_OPTIONS) {
    if      (OPTION("-q"))           verbose     = 2;
    else if (OPTION("-minmax"))      minmax      = true;
    else if (OPTION("-percentiles")) minmax      = false;
    else if (OPTION("-name"))        scalar_name = ARGUMENT;
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

#if !defined(HAS_VTK)
  cerr << EXECNAME << ": Must be compiled and linked with the VTK library" << endl;
  return 1;
#else // HAS_VTK
  vtkSmartPointer<vtkPointSet> pointset = ReadPointSet(input_name);
  const vtkIdType noOfPoints = pointset->GetNumberOfPoints();

  vtkDataArray *scalars = NULL;
  if (scalar_name) {
    scalars = GetArrayByCaseInsensitiveName(pointset->GetPointData(), scalar_name);
    if (!scalars) {
      cerr << "Input point set has no point data array named " << scalar_name << " (case insensitive)" << endl;
      exit(1);
    }
  } else {
    scalars = pointset->GetPointData()->GetScalars();
    if (!scalars) {
      cerr << "Input point set has no SCALARS point data" << endl;
      exit(1);
    }
  }

  if (!minmax) {
    vector<double> data(noOfPoints);
    for (vtkIdType i = 0; i < noOfPoints; ++i) data[i] = scalars->GetTuple1(i);
    sort(data.begin(), data.end());
    minVal = data[int(round(minVal * (noOfPoints - 1) / 100.0))];
    maxVal = data[int(round(maxVal * (noOfPoints - 1) / 100.0))];
  }
  if (verbose > 1) {
    cout << minVal << "," << maxVal << endl;
  } else if (verbose) {
    cout << "No of points " << noOfPoints << endl;
    cout << "Clamping to min and max : " << minVal << ", " << maxVal << endl;
  }

  for (vtkIdType i = 0; i < noOfPoints; ++i) {
    scalars->SetTuple1(i, clamp(scalars->GetTuple1(i), minVal, maxVal));
  }

  return WritePointSet(output_name, pointset) ? 0 : 1;
#endif // HAS_VTK
}

