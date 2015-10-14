// =============================================================================
// Project: Image Registration Toolkit (IRTK)
// Package: PolyData
//
// Copyright (c) 2015 Imperial College London
// Copyright (c) 2015 Andreas Schuh
// =============================================================================

#include <irtkCommon.h>

// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
/// Print help screen
void PrintHelp(const char *name)
{
  cout << "usage: " << name << " <input> <mat> <output> [options]" << endl;
  cout << endl;
  cout << "Sets/Adds variables stored in a MAT file as point data of a VTK dataset." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  input   File name of input dataset." << endl;
  cout << "  mat     File name of MAT file." << endl;
  cout << "  output  File name of output dataset." << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << "  -var <name> [<type>]   Set/Add variable with given name as point data" << endl;
  cout << "                         of specified type. If no type is specified, the" << endl;
  cout << "                         array is added as non-attribute point data array." << endl;
  cout << "                         Type can be 'scalars', 'vectors', 'normals', 'tcoords'," << endl;
  cout << "                         'tensors', 'globalids', or 'pedigreeids' (case insensitive)." << endl;
}

// =============================================================================
// Includes
// =============================================================================

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkDataSetAttributes.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkCharArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkShortArray.h>
#include <vtkUnsignedShortArray.h>
#include <vtkIntArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkLongArray.h>
#include <vtkUnsignedLongArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>

#include <vtkPolyDataReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkOBJReader.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>
#include <vtkPLYReader.h>
#include <vtkPLYWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyDataWriter.h>

#include <mclmcrrt.h>

// =============================================================================
// Types
// =============================================================================

// -----------------------------------------------------------------------------
/// Enumeration of explicit requested output array type
enum OutputType
{
  InputType = 0,
  OutputChar,
  OutputUChar,
  OutputShort,
  OutputUShort,
  OutputInt,
  OutputUInt,
  OutputLong,
  OutputULong,
  OutputFloat,
  OutputDouble
};

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
/// Convert mxArray to vtkDataArray
template <typename TData, class TArray = vtkFloatArray>
vtkSmartPointer<vtkDataArray> ToTypedDataArray(const mxArray *pa, bool transpose)
{
  const mwSize ndim = mxGetNumberOfDimensions(pa);
  if (ndim > 2) {
    cerr << "Error: MAT variable has unsupported number of dimensions: " << ndim << endl;
    exit(1);
  }
  vtkSmartPointer<TArray> data = vtkSmartPointer<TArray>::New();
  const TData *real = reinterpret_cast<TData *>(mxGetData(pa));
  const size_t m = mxGetM(pa);
  const size_t n = mxGetN(pa);
  mwIndex i, subs[2];
  if (transpose) {
    data->SetNumberOfComponents(m);
    data->SetNumberOfTuples    (n);
    for (subs[0] = 0; subs[0] < m; ++subs[0]) {
      for (subs[1] = 0; subs[1] < n; ++subs[1]) {
        i = mxCalcSingleSubscript(pa, ndim, subs);
        data->SetComponent(subs[1], subs[0], real[i]);
      }
    }
  } else {
    data->SetNumberOfComponents(n);
    data->SetNumberOfTuples    (m);
    for (subs[0] = 0; subs[0] < m; ++subs[0]) {
      for (subs[1] = 0; subs[1] < n; ++subs[1]) {
        i = mxCalcSingleSubscript(pa, ndim, subs);
        data->SetComponent(subs[0], subs[1], real[i]);
      }
    }
  }
  return data;
}

// -----------------------------------------------------------------------------
/// Convert mxArray to vtkDataArray of type specified as template argument
template <class TArray>
vtkSmartPointer<vtkDataArray> ToCastedDataArray(const mxArray *pa, bool transpose = false)
{
  const mxClassID classID = mxGetClassID(pa);
  switch (classID) {
    case mxINT8_CLASS:   return ToTypedDataArray<char,           TArray>(pa, transpose);
    case mxUINT8_CLASS:  return ToTypedDataArray<unsigned char,  TArray>(pa, transpose);
    case mxINT16_CLASS:  return ToTypedDataArray<short,          TArray>(pa, transpose);
    case mxUINT16_CLASS: return ToTypedDataArray<unsigned short, TArray>(pa, transpose);
    case mxINT32_CLASS:  return ToTypedDataArray<int,            TArray>(pa, transpose);
    case mxUINT32_CLASS: return ToTypedDataArray<unsigned int,   TArray>(pa, transpose);
    case mxINT64_CLASS:  return ToTypedDataArray<long,           TArray>(pa, transpose);
    case mxUINT64_CLASS: return ToTypedDataArray<unsigned long,  TArray>(pa, transpose);
    case mxSINGLE_CLASS: return ToTypedDataArray<float,          TArray>(pa, transpose);
    case mxDOUBLE_CLASS: return ToTypedDataArray<double,         TArray>(pa, transpose);
    default:
      cerr << "Error: MAT variable has unsupported type: " << classID << endl;
      exit(1);
  }
  return NULL;
}

// -----------------------------------------------------------------------------
/// Convert mxArray to vtkDataArray of corresponding type
vtkSmartPointer<vtkDataArray> ToDataArray(const mxArray *pa, bool transpose = false)
{
  const mxClassID classID = mxGetClassID(pa);
  switch (classID) {
    case mxINT8_CLASS:   return ToTypedDataArray<char,           vtkCharArray         >(pa, transpose);
    case mxUINT8_CLASS:  return ToTypedDataArray<unsigned char,  vtkUnsignedCharArray >(pa, transpose);
    case mxINT16_CLASS:  return ToTypedDataArray<short,          vtkShortArray        >(pa, transpose);
    case mxUINT16_CLASS: return ToTypedDataArray<unsigned short, vtkUnsignedShortArray>(pa, transpose);
    case mxINT32_CLASS:  return ToTypedDataArray<int,            vtkIntArray          >(pa, transpose);
    case mxUINT32_CLASS: return ToTypedDataArray<unsigned int,   vtkUnsignedIntArray  >(pa, transpose);
    case mxINT64_CLASS:  return ToTypedDataArray<long,           vtkLongArray         >(pa, transpose);
    case mxUINT64_CLASS: return ToTypedDataArray<unsigned long,  vtkUnsignedLongArray >(pa, transpose);
    case mxSINGLE_CLASS: return ToTypedDataArray<float,          vtkFloatArray        >(pa, transpose);
    case mxDOUBLE_CLASS: return ToTypedDataArray<double,         vtkDoubleArray       >(pa, transpose);
    default:
      cerr << "Error: MAT variable has unsupported type: " << classID << endl;
      exit(1);
  }
  return NULL;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  REQUIRES_POSARGS(3);

  const char *input_name  = POSARG(1);
  const char *mat_name    = POSARG(2);
  const char *output_name = POSARG(3);

  // Optional arguments
  vector<string> var_name;
  vector<int>    var_type;
  bool           ascii    = false;
  bool           compress = true;
  OutputType     out_type = InputType;

  for (ALL_OPTIONS) {
    if      (OPTION("-var")) {
      var_name.push_back(ARGUMENT);
      if (HAS_ARGUMENT) var_type.push_back(PolyDataAttributeType(ARGUMENT));
      else              var_type.push_back(-1);
    }
    else if (OPTION("-char"))   out_type = OutputChar;
    else if (OPTION("-uchar"))  out_type = OutputUChar;
    else if (OPTION("-short"))  out_type = OutputShort;
    else if (OPTION("-ushort")) out_type = OutputUShort;
    else if (OPTION("-int"))    out_type = OutputInt;
    else if (OPTION("-uint"))   out_type = OutputUInt;
    else if (OPTION("-long"))   out_type = OutputLong;
    else if (OPTION("-float"))  out_type = OutputFloat;
    else if (OPTION("-double")) out_type = OutputDouble;
    else if (OPTION("-ascii")  || OPTION("-nobinary")) ascii = true;
    else if (OPTION("-binary") || OPTION("-noascii") ) ascii = false;
    else if (OPTION("-compress"))   compress = true;
    else if (OPTION("-nocompress")) compress = false;
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // Initialize MCR
  mclmcrInitialize();

  // Read input dataset
  vtkSmartPointer<vtkPolyData> dataset = ReadPolyData(input_name);
  if (verbose > 1) cout << "Reading dataset from " << input_name << endl;
  if (dataset == NULL || dataset->GetNumberOfPoints() == 0) {
    cerr << "Error: Failed to read input dataset or dataset has no points!" << endl;
    exit(1);
  }

  // Open MAT file
  if (verbose > 1) cout << "Reading MAT variables from " << mat_name << endl;
  MATFile *mat_file = matOpen(mat_name, "r");
  if (mat_file == NULL) {
    cerr << "Error: Failed to open MAT file " << mat_name << endl;
    exit(1);
  }

  // List directory of MAT file if no specific variable chosen by user
  if (var_name.empty()) {
    int         ndir;
    const char **dir = const_cast<const char **>(reinterpret_cast<char **>(matGetDir(mat_file, &ndir)));
    if (dir == NULL) {
      cerr << "Error: Failed to read directory of MAT file " << mat_name << endl;
      exit(1);
    }
    for (int i = 0; i < ndir; ++i) {
      var_name.push_back(dir[i]);
      var_type.push_back(-1);
    }
    mxFree(dir);
    if (matClose(mat_file) != 0) {
      cerr << "Error: Failed to close MAT file" << endl;
      exit(1);
    }
    mat_file = matOpen(mat_name, "r");
    if (mat_file == NULL) {
      cerr << "Error: Failed to reopen MAT file " << mat_name << endl;
      exit(1);
    }
  }

  // Set/Add variables in MAT file as point data
  const char                   *name;
  mxArray                      *var;
  vtkSmartPointer<vtkDataArray> array;
  vtkPointData                 *pd = dataset->GetPointData();
  int                           pos;

  for (size_t i = 0; i < var_name.size(); ++i) {
    name  = var_name[i].c_str();
    var   = matGetVariable(mat_file, name);
    if (var == NULL) {
      cerr << "Error: Unknown MAT variable " << name << endl;
      exit(1);
    }
    const mwSize ndim = mxGetNumberOfDimensions(var);
    if (ndim > 2) {
      cerr << "Error: MAT variable has unsupported number of dimensions: " << ndim << endl;
      exit(1);
    }
    bool transpose = false;
    if (mxGetM(var) != dataset->GetNumberOfPoints()) {
      if (mxGetN(var) != dataset->GetNumberOfPoints()) {
        cerr << "Error: MAT variable " << name << " has invalid number of rows (cols)" << endl;
        exit(1);
      }
      transpose = true;
    }
    switch (out_type) {
      case OutputChar:   array = ToCastedDataArray<vtkCharArray         >(var, transpose); break;
      case OutputUChar:  array = ToCastedDataArray<vtkUnsignedCharArray >(var, transpose); break;
      case OutputShort:  array = ToCastedDataArray<vtkShortArray        >(var, transpose); break;
      case OutputUShort: array = ToCastedDataArray<vtkUnsignedShortArray>(var, transpose); break;
      case OutputInt:    array = ToCastedDataArray<vtkIntArray          >(var, transpose); break;
      case OutputUInt:   array = ToCastedDataArray<vtkUnsignedIntArray  >(var, transpose); break;
      case OutputLong:   array = ToCastedDataArray<vtkLongArray         >(var, transpose); break;
      case OutputULong:  array = ToCastedDataArray<vtkUnsignedLongArray >(var, transpose); break;
      case OutputFloat:  array = ToCastedDataArray<vtkFloatArray        >(var, transpose); break;
      case OutputDouble: array = ToCastedDataArray<vtkDoubleArray       >(var, transpose); break;
      default: array = ToDataArray(var, transpose);
    }
    mxDestroyArray(var);
    array->SetName(name);
    switch (var_type[i]) {
      case vtkDataSetAttributes::SCALARS:     pos = pd->SetScalars    (array); break;
      case vtkDataSetAttributes::VECTORS:     pos = pd->SetVectors    (array); break;
      case vtkDataSetAttributes::NORMALS:     pos = pd->SetNormals    (array); break;
      case vtkDataSetAttributes::TCOORDS:     pos = pd->SetTCoords    (array); break;
      case vtkDataSetAttributes::TENSORS:     pos = pd->SetTensors    (array); break;
      case vtkDataSetAttributes::GLOBALIDS:   pos = pd->SetGlobalIds  (array); break;
      case vtkDataSetAttributes::PEDIGREEIDS: pos = pd->SetPedigreeIds(array); break;
      default: pos = pd->AddArray(array);
    }
    if (verbose) cout << pos << ": " << var_name[i] << endl;
  }

  // Write modified dataset
  WritePolyData(output_name, dataset, compress, ascii);
  if (verbose > 1) cout << "Wrote modified dataset to " << output_name << endl;

  // Close MAT file
  if (matClose(mat_file) != 0) {
    cerr << "Error: Failed to close MAT file" << endl;
    exit(1);
  }
  return 0;
}
