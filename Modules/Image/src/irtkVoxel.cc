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

#include <irtkVoxel.h>

// -----------------------------------------------------------------------------
int DataTypeSize(int type)
{
  switch (type){
    case IRTK_VOXEL_CHAR:           return sizeof(char);
    case IRTK_VOXEL_UNSIGNED_CHAR:  return sizeof(unsigned char);
    case IRTK_VOXEL_SHORT:          return sizeof(short);
    case IRTK_VOXEL_UNSIGNED_SHORT: return sizeof(unsigned short);
    case IRTK_VOXEL_INT:            return sizeof(int);
    case IRTK_VOXEL_UNSIGNED_INT:   return sizeof(unsigned int);
    case IRTK_VOXEL_FLOAT:          return sizeof(float);
    case IRTK_VOXEL_DOUBLE:         return sizeof(double);
    case IRTK_VOXEL_FLOAT1:         return sizeof(float1);
    case IRTK_VOXEL_FLOAT2:         return sizeof(float2);
    case IRTK_VOXEL_FLOAT3:         return sizeof(float3);
    case IRTK_VOXEL_FLOAT4:         return sizeof(float4);
    case IRTK_VOXEL_FLOAT2x2:       return sizeof(float2x2);
    case IRTK_VOXEL_FLOAT3x3:       return sizeof(float3x3);
    case IRTK_VOXEL_FLOAT3x4:       return sizeof(float3x4);
    case IRTK_VOXEL_FLOAT4x4:       return sizeof(float4x4);
    case IRTK_VOXEL_DOUBLE1:        return sizeof(double1);
    case IRTK_VOXEL_DOUBLE2:        return sizeof(double2);
    case IRTK_VOXEL_DOUBLE3:        return sizeof(double3);
    case IRTK_VOXEL_DOUBLE4:        return sizeof(double4);
    case IRTK_VOXEL_DOUBLE2x2:      return sizeof(double2x2);
    case IRTK_VOXEL_DOUBLE3x3:      return sizeof(double3x3);
    case IRTK_VOXEL_DOUBLE3x4:      return sizeof(double3x4);
    case IRTK_VOXEL_DOUBLE4x4:      return sizeof(double4x4);
    default:                        return 0;
  }
}

// -----------------------------------------------------------------------------
std::string DataTypeName(int type)
{
  switch (type){
    case IRTK_VOXEL_CHAR:           return "char";
    case IRTK_VOXEL_UNSIGNED_CHAR:  return "uchar";
    case IRTK_VOXEL_SHORT:          return "short";
    case IRTK_VOXEL_UNSIGNED_SHORT: return "ushort";
    case IRTK_VOXEL_INT:            return "int";
    case IRTK_VOXEL_UNSIGNED_INT:   return "uint";
    case IRTK_VOXEL_FLOAT:          return "float";
    case IRTK_VOXEL_DOUBLE:         return "double";
    case IRTK_VOXEL_RGB:            return "RGB";
    case IRTK_VOXEL_FLOAT1:         return "float1";
    case IRTK_VOXEL_FLOAT2:         return "float2";
    case IRTK_VOXEL_FLOAT3:         return "float3";
    case IRTK_VOXEL_FLOAT4:         return "float4";
    case IRTK_VOXEL_FLOAT1x1:       return "float1x1";
    case IRTK_VOXEL_FLOAT2x2:       return "float2x2";
    case IRTK_VOXEL_FLOAT3x3:       return "float3x3";
    case IRTK_VOXEL_FLOAT3x4:       return "float3x4";
    case IRTK_VOXEL_FLOAT4x4:       return "float4x4";
    case IRTK_VOXEL_DOUBLE1:        return "double1";
    case IRTK_VOXEL_DOUBLE2:        return "double2";
    case IRTK_VOXEL_DOUBLE3:        return "double3";
    case IRTK_VOXEL_DOUBLE4:        return "double4";
    case IRTK_VOXEL_DOUBLE1x1:      return "double1x1";
    case IRTK_VOXEL_DOUBLE2x2:      return "double2x2";
    case IRTK_VOXEL_DOUBLE3x3:      return "double3x3";
    case IRTK_VOXEL_DOUBLE3x4:      return "double3x4";
    case IRTK_VOXEL_DOUBLE4x4:      return "double4x4";
    default:                        return "unknown";
  }
}

// -----------------------------------------------------------------------------
int ToDataType(const char *name)
{
  int type = 0;
  while (type <= IRKT_VOXEL_LAST && DataTypeName(type) != name) ++type;
  if (DataTypeName(type) == "unknown") return IRTK_VOXEL_UNKNOWN;
  return type;
}

// -----------------------------------------------------------------------------
int ToDataType(const string &name)
{
  return ToDataType(name.c_str());
}
