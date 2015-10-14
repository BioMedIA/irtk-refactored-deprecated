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

/// @file  irtkDataOp.h
/// @brief Defines base class and I/O functions for arbitrary 1D data sequences
///
/// Functions to manipulate the data are defined in irtkDataFunctions.h.
/// Statistics of the data sequence such as mean and variance or percentiles
/// can be computed using the operators found in irtkDataStatistics.h.
/// The data operators are used in particular by the calculate.cc tool for
/// which they were originally developed. They were added to the linear
/// algebra library because they are useful to compute common statistics or
/// perform basic mathematical operations on a data sequence such as an image
/// or attributes of a VTK point set.
///
/// @sa irtkDataFunctions.h
/// @sa irtkDataStatistics.h
#ifndef _IRTKDATAOP_H
#define _IRTKDATAOP_H

#include <irtkImage.h>

#ifdef HAS_VTK
#include <vtkSmartPointer.h>
#include <vtkDataSet.h>
#include <vtkDataArray.h>
#endif


namespace irtk { namespace data {


// =============================================================================
// Base class of data operations
// =============================================================================

// -----------------------------------------------------------------------------
/// Base class of all data operations
class Op
{
public:

  /// Destructor
  virtual ~Op() {}

  /// Process given data
  virtual void Process(int, double *, bool * = NULL) = 0;

#ifdef HAS_VTK
  /// Process given vtkDataArray
  virtual void Process(vtkDataArray *data, bool *mask = NULL)
  {
    const int n = static_cast<int>(data->GetNumberOfTuples() * data->GetNumberOfComponents());
    if (data->GetDataType() == VTK_DOUBLE) {
      this->Process(n, reinterpret_cast<double *>(data->GetVoidPointer(0)), mask);
    } else {
      double *d = new double[n];
      double *tuple = d;
      for (vtkIdType i = 0; i < data->GetNumberOfTuples(); ++i) {
        data->GetTuple(i, tuple);
        tuple += data->GetNumberOfComponents();
      }
      this->Process(n, d, mask);
      for (vtkIdType i = 0; i < data->GetNumberOfTuples(); ++i) {
        data->SetTuple(i, tuple);
        tuple += data->GetNumberOfComponents();
      }
      delete[] d;
    }
  }
#endif
};

// =============================================================================
// I/O functions
// =============================================================================

// -----------------------------------------------------------------------------
/// Enumeration of supported input data file types
enum DataFileType
{
  UnknownDataFile,
  IMAGE,
  LEGACY_VTK,
  XML_VTK
};

/// Get (or guess) type of input file
DataFileType FileType(const char *name);

// -----------------------------------------------------------------------------
/// Read data sequence from any supported input file type
#ifdef HAS_VTK
int Read(const char *name, double *&data, int *dtype = NULL, irtkImageAttributes *attr = NULL,
         vtkSmartPointer<vtkDataSet> *dataset= NULL, const char *scalar_name = NULL);
#else
int Read(const char *name, double *&data, int *dtype = NULL, irtkImageAttributes *attr = NULL);
#endif

// -----------------------------------------------------------------------------
/// Write data sequence
class Write : public Op
{
  /// Name of output file
  irtkPublicAttributeMacro(string, FileName);

#ifdef HAS_VTK
  /// VTK input dataset whose scalar data was modified
  irtkPublicAttributeMacro(vtkSmartPointer<vtkDataSet>, DataSet);

  /// Name of input/output point data array
  irtkPublicAttributeMacro(string, ArrayName);
#endif

  /// Attributes of input image whose data was modified
  irtkPublicAttributeMacro(irtkImageAttributes, Attributes);

  /// Output data type
  irtkPublicAttributeMacro(int, DataType);

public:

  /// Constructor
#ifdef HAS_VTK
  Write(const char *fname, int dtype = IRTK_VOXEL_DOUBLE,
        irtkImageAttributes attr = irtkImageAttributes(),
        vtkDataSet *dataset      = NULL,
        const char *scalars_name = NULL)
  :
    _FileName(fname),
    _DataSet(dataset),
    _Attributes(attr),
    _DataType(dtype)
  {
    if (scalars_name) _ArrayName = scalars_name;
  }
#else
  Write(const char *fname, int dtype = IRTK_VOXEL_DOUBLE,
        irtkImageAttributes attr = irtkImageAttributes())
  :
    _FileName(fname),
    _Attributes(attr),
    _DataType(dtype)
  {}
#endif

  /// Process given data
  virtual void Process(int n, double *data, bool * = NULL);

};

// =============================================================================
// Auxiliary macros for subclass implementation
// =============================================================================

// -----------------------------------------------------------------------------
// Add Calculate function that takes a vtkDataArray as argument
// and computes a single return value
//
// Subclass must implement:
//   template <class T> double Calculate(int n, const T *data, const bool *mask)
#ifdef HAS_VTK
  #define irtkCalculateVtkDataArray1() \
    static double Calculate(vtkDataArray *data, const bool *mask = NULL) \
    { \
      const void *p = data->GetVoidPointer(0); \
      const int   n = static_cast<int>(data->GetNumberOfTuples() * \
                                       data->GetNumberOfComponents()); \
      switch (data->GetDataType()) { \
        case VTK_SHORT:  return Calculate(n, reinterpret_cast<const short  *>(p), mask); \
        case VTK_INT:    return Calculate(n, reinterpret_cast<const int    *>(p), mask); \
        case VTK_FLOAT:  return Calculate(n, reinterpret_cast<const float  *>(p), mask); \
        case VTK_DOUBLE: return Calculate(n, reinterpret_cast<const double *>(p), mask); \
        default: \
          cerr << "Unsupported vtkDataArray type: " << data->GetDataType() << endl; \
          exit(1); \
      } \
      return numeric_limits<double>::quiet_NaN(); \
    }
#else
  #define irtkCalculateVtkDataArray1()
#endif

// -----------------------------------------------------------------------------
// Add Calculate function that takes a vtkDataArray as argument
// and computes two return values
//
// Subclass must implement:
//   template <class T>
//   void Calculate(double &, double &, int n, const T *data, const bool *mask)
#ifdef HAS_VTK
  #define irtkCalculateVtkDataArray2() \
    static void Calculate(double &v1, double &v2, vtkDataArray *data, const bool *mask = NULL) \
    { \
      const void *p = data->GetVoidPointer(0); \
      const int   n = static_cast<int>(data->GetNumberOfTuples() * \
                                       data->GetNumberOfComponents()); \
      switch (data->GetDataType()) { \
        case VTK_SHORT:  Calculate(v1, v2, n, reinterpret_cast<const short  *>(p), mask); break; \
        case VTK_INT:    Calculate(v1, v2, n, reinterpret_cast<const int    *>(p), mask); break; \
        case VTK_FLOAT:  Calculate(v1, v2, n, reinterpret_cast<const float  *>(p), mask); break; \
        case VTK_DOUBLE: Calculate(v1, v2, n, reinterpret_cast<const double *>(p), mask); break; \
        default: \
          cerr << "Unsupported vtkDataArray type: " << data->GetDataType() << endl; \
          exit(1); \
      } \
    }
#else
  #define irtkCalculateVtkDataArray2()
#endif


} } // namespace irtk::data

#endif
