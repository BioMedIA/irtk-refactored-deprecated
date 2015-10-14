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

#ifndef _IRTKVTK_H
#define _IRTKVTK_H

#include <vtkSmartPointer.h>
#include <vtkCharArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkShortArray.h>
#include <vtkUnsignedShortArray.h>
#include <vtkIntArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>


// VTK magic number
#define VTK_MAGIC1        "# vtk DataFile Version 1.0"
#define VTK_MAGIC2        "# vtk DataFile Version 2.0"
#define VTK_MAGIC3        "# vtk DataFile Version 3.0"

// VTK data types
#define VTK_DATA_CHAR          "char"
#define VTK_DATA_U_CHAR        "unsigned_char"
#define VTK_DATA_SHORT         "short"
#define VTK_DATA_U_SHORT       "unsigned_short"
#define VTK_DATA_FLOAT         "float"

// -----------------------------------------------------------------------------
/// Get VTK data type from IRTK voxel type
inline int ToVTKDataType(int type)
{
  switch (type) {
    case IRTK_VOXEL_CHAR:           return VTK_CHAR;
    case IRTK_VOXEL_UNSIGNED_CHAR:  return VTK_UNSIGNED_CHAR;
    case IRTK_VOXEL_SHORT:          return VTK_SHORT;
    case IRTK_VOXEL_UNSIGNED_SHORT: return VTK_UNSIGNED_SHORT;
    case IRTK_VOXEL_INT:            return VTK_INT;
    case IRTK_VOXEL_UNSIGNED_INT:   return VTK_UNSIGNED_INT;
    case IRTK_VOXEL_FLOAT:          return VTK_FLOAT;
    case IRTK_VOXEL_DOUBLE:         return VTK_DOUBLE;
    default:                        return VTK_VOID;
  }
}

// -----------------------------------------------------------------------------
/// Get IRTK voxel type from VTK data type
inline int FromVTKDataType(int type)
{
  switch (type) {
    case VTK_CHAR:           return IRTK_VOXEL_CHAR;
    case VTK_UNSIGNED_CHAR:  return IRTK_VOXEL_UNSIGNED_CHAR;
    case VTK_SHORT:          return IRTK_VOXEL_SHORT;
    case VTK_UNSIGNED_SHORT: return IRTK_VOXEL_UNSIGNED_SHORT;
    case VTK_INT:            return IRTK_VOXEL_INT;
    case VTK_UNSIGNED_INT:   return IRTK_VOXEL_UNSIGNED_INT;
    case VTK_FLOAT:          return IRTK_VOXEL_FLOAT;
    case VTK_DOUBLE:         return IRTK_VOXEL_DOUBLE;
    default:                 return IRTK_VOXEL_UNKNOWN;
  }
}

// -----------------------------------------------------------------------------
/// Instantiate new VTK data array of given type
inline vtkSmartPointer<vtkDataArray> NewVTKDataArray(int vtkType)
{
  switch (vtkType) {
    case VTK_CHAR:           return vtkSmartPointer<vtkCharArray>::New();
    case VTK_UNSIGNED_CHAR:  return vtkSmartPointer<vtkUnsignedCharArray>::New();
    case VTK_SHORT:          return vtkSmartPointer<vtkShortArray>::New();
    case VTK_UNSIGNED_SHORT: return vtkSmartPointer<vtkUnsignedShortArray>::New();
    case VTK_INT:            return vtkSmartPointer<vtkIntArray>::New();
    case VTK_UNSIGNED_INT:   return vtkSmartPointer<vtkUnsignedIntArray>::New();
    case VTK_FLOAT:          return vtkSmartPointer<vtkFloatArray>::New();
    case VTK_DOUBLE:         return vtkSmartPointer<vtkDoubleArray>::New();
    default:
      cerr << "Invalid VTK data type: " << vtkType << endl;
      exit(1);
  }
}


#endif
