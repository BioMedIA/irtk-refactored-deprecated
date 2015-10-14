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

#ifndef IRTKVOXEL_H_
#define IRTKVOXEL_H_

#include <irtkVector.h>


// =============================================================================
// Voxel type enumeration
// =============================================================================

#define IRTK_VOXEL_UNKNOWN             0
#define IRTK_VOXEL_CHAR                1
#define IRTK_VOXEL_UNSIGNED_CHAR       2
#define IRTK_VOXEL_SHORT               3
#define IRTK_VOXEL_UNSIGNED_SHORT      4
#define IRTK_VOXEL_INT                 5
#define IRTK_VOXEL_UNSIGNED_INT        6
#define IRTK_VOXEL_FLOAT               7
#define IRTK_VOXEL_DOUBLE              8
#define IRTK_VOXEL_RGB                 9
#define IRTK_VOXEL_FLOAT1             10 // unused
#define IRTK_VOXEL_FLOAT2             11
#define IRTK_VOXEL_FLOAT3             12
#define IRTK_VOXEL_FLOAT4             13
#define IRTK_VOXEL_DOUBLE1            20 // unused
#define IRTK_VOXEL_DOUBLE2            21
#define IRTK_VOXEL_DOUBLE3            22
#define IRTK_VOXEL_DOUBLE4            23
#define IRTK_VOXEL_FLOAT1x1           30 // unused
#define IRTK_VOXEL_FLOAT2x2           31
#define IRTK_VOXEL_FLOAT3x3           32
#define IRTK_VOXEL_FLOAT3x4           33
#define IRTK_VOXEL_FLOAT4x4           34
#define IRTK_VOXEL_DOUBLE1x1          40 // unused
#define IRTK_VOXEL_DOUBLE2x2          41
#define IRTK_VOXEL_DOUBLE3x3          42
#define IRTK_VOXEL_DOUBLE3x4          43
#define IRTK_VOXEL_DOUBLE4x4          44
#define IRKT_VOXEL_LAST IRTK_VOXEL_DOUBLE4x4

#define IRTK_VOXEL_BINARY     IRTK_VOXEL_UNSIGNED_CHAR
#define IRTK_VOXEL_BYTE       IRTK_VOXEL_UNSIGNED_CHAR
#define IRTK_VOXEL_GREY       IRTK_VOXEL_SHORT
#if USE_FLOAT_BY_DEFAULT
#  define IRTK_VOXEL_REAL     IRTK_VOXEL_FLOAT
#  define IRTK_VOXEL_REAL1    IRTK_VOXEL_FLOAT1
#  define IRTK_VOXEL_REAL2    IRTK_VOXEL_FLOAT2
#  define IRTK_VOXEL_REAL3    IRTK_VOXEL_FLOAT3
#  define IRTK_VOXEL_REAL4    IRTK_VOXEL_FLOAT4
#  define IRTK_VOXEL_REAL1x1  IRTK_VOXEL_FLOAT1x1
#  define IRTK_VOXEL_REAL2x2  IRTK_VOXEL_FLOAT2x2
#  define IRTK_VOXEL_REAL3x3  IRTK_VOXEL_FLOAT3x3
#  define IRTK_VOXEL_REAL3x4  IRTK_VOXEL_FLOAT3x4
#  define IRTK_VOXEL_REAL4x4  IRTK_VOXEL_FLOAT4x4
#else
#  define IRTK_VOXEL_REAL     IRTK_VOXEL_DOUBLE
#  define IRTK_VOXEL_REAL1    IRTK_VOXEL_DOUBLE1
#  define IRTK_VOXEL_REAL2    IRTK_VOXEL_DOUBLE2
#  define IRTK_VOXEL_REAL3    IRTK_VOXEL_DOUBLE3
#  define IRTK_VOXEL_REAL4    IRTK_VOXEL_DOUBLE4
#  define IRTK_VOXEL_REAL1x1  IRTK_VOXEL_DOUBLE1x1
#  define IRTK_VOXEL_REAL2x2  IRTK_VOXEL_DOUBLE2x2
#  define IRTK_VOXEL_REAL3x3  IRTK_VOXEL_DOUBLE3x3
#  define IRTK_VOXEL_REAL3x4  IRTK_VOXEL_DOUBLE3x4
#  define IRTK_VOXEL_REAL4x4  IRTK_VOXEL_DOUBLE4x4
#endif

// =============================================================================
// Voxel types
// =============================================================================

typedef unsigned char irtkBinaryPixel;
typedef unsigned char irtkBytePixel;
typedef short         irtkGreyPixel;
typedef realt         irtkRealPixel;

// TODO: Remove, and use float3, float4, double3, double4 instead
typedef irtkVector3D<float>    irtkFloat3;
typedef irtkVector4D<float>    irtkFloat4;
typedef irtkVector3D<double>   irtkDouble3;
typedef irtkVector4D<double>   irtkDouble4;

// =============================================================================
// Voxel type limits
// =============================================================================

#define MIN_GREY std::numeric_limits<irtkGreyPixel>::min()
#define MAX_GREY std::numeric_limits<irtkGreyPixel>::max()

// -----------------------------------------------------------------------------
template <class T>
struct voxel_limits
{
  /// Minimum value that can be represented by this voxel type as double
  /// \note For vector types, this corresponds to the minium value that can
  ///       be represented by each component of the vector.
  static double min() throw();
  /// Maximum value that can be represented by this voxel type as double
  /// \note For vector types, this corresponds to the maxium value that can
  ///       be represented by each component of the vector.
  static double max() throw();
  /// Minimum value that can be represented by this voxel type
  static T min_value() throw();
  /// Maximum value that can be represented by this voxel type
  static T max_value() throw();
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<char>
{
  static char   min_value() throw() { return static_cast<char  >(0x80); }
  static char   max_value() throw() { return static_cast<char  >(0x7f); }
	static double min()       throw() { return static_cast<double>(min_value()); }
	static double max()       throw() { return static_cast<double>(max_value()); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<unsigned char>
{
	static unsigned char min_value() throw() { return static_cast<unsigned char>(0u); }
	static unsigned char max_value() throw() { return static_cast<unsigned char>(0xffu); }
  static double        min()       throw() { return static_cast<double       >(min_value()); }
	static double        max()       throw() { return static_cast<double       >(max_value()); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<short>
{
	static short  min_value() throw() { return static_cast<short >(0x8000); }
	static short  max_value() throw() { return static_cast<short >(0x7fff); }
  static double min()       throw() { return static_cast<double>(min_value()); }
	static double max()       throw() { return static_cast<double>(max_value()); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<unsigned short>
{
	static unsigned short min_value() throw() { return static_cast<unsigned short>(0u); }
	static unsigned short max_value() throw() { return static_cast<unsigned short>(0xffffu); }
  static double         min()       throw() { return static_cast<double        >(min_value()); }
	static double         max()       throw() { return static_cast<double        >(max_value()); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<int> {
	static int    min_value() throw() { return static_cast<int   >(~(~0u >> 1)); }
	static int    max_value() throw() { return static_cast<int   >( ~0u >> 1); }
  static double min()       throw() { return static_cast<double>(min_value()); }
	static double max()       throw() { return static_cast<double>(max_value()); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<unsigned int>
{
	static unsigned int min_value() throw() { return static_cast<unsigned int>(0u); }
	static unsigned int max_value() throw() { return static_cast<unsigned int>( ~0u >> 1); }
  static double       min()       throw() { return static_cast<double      >(min_value()); }
	static double       max()       throw() { return static_cast<double      >(max_value()); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<float>
{
	static float  min_value() throw() { return static_cast<float >(-1.0e+38f); }
	static float  max_value() throw() { return static_cast<float >( 1.0e+38f); }
  static double min()       throw() { return static_cast<double>(min_value()); }
	static double max()       throw() { return static_cast<double>(max_value()); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<double>
{
	static double min_value() throw() { return static_cast<double>(-1.0e+299); }
	static double max_value() throw() { return static_cast<double>( 1.0e+299); }
  static double min()       throw() { return static_cast<double>(min_value()); }
	static double max()       throw() { return static_cast<double>(max_value()); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<float1>
{
	static float1 min_value() throw() { return make_float1(voxel_limits<float>::min_value()); }
	static float1 max_value() throw() { return make_float1(voxel_limits<float>::max_value()); }
  static double min() throw() { return voxel_limits<float>::min(); }
	static double max() throw() { return voxel_limits<float>::max(); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<float2>
{
	static float2 min_value() throw() { return make_float2(voxel_limits<float>::min_value()); }
	static float2 max_value() throw() { return make_float2(voxel_limits<float>::max_value()); }
  static double min() throw() { return voxel_limits<float>::min(); }
	static double max() throw() { return voxel_limits<float>::max(); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<float3>
{
	static float3 min_value() throw() { return make_float3(voxel_limits<float>::min_value()); }
	static float3 max_value() throw() { return make_float3(voxel_limits<float>::max_value()); }
  static double min() throw() { return voxel_limits<float>::min(); }
	static double max() throw() { return voxel_limits<float>::max(); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<float4>
{
	static float4 min_value() throw() { return make_float4(voxel_limits<float>::min_value()); }
	static float4 max_value() throw() { return make_float4(voxel_limits<float>::max_value()); }
  static double min() throw() { return voxel_limits<float>::min(); }
	static double max() throw() { return voxel_limits<float>::max(); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<float3x3>
{
	static float3x3 min_value() throw() { return make_float3x3(voxel_limits<float>::min_value()); }
	static float3x3 max_value() throw() { return make_float3x3(voxel_limits<float>::max_value()); }
  static double min() throw() { return voxel_limits<float>::min(); }
	static double max() throw() { return voxel_limits<float>::max(); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<double1>
{
	static double1 min_value() throw() { return make_double1(voxel_limits<double>::min_value()); }
	static double1 max_value() throw() { return make_double1(voxel_limits<double>::max_value()); }
  static double min() throw() { return voxel_limits<double>::min(); }
	static double max() throw() { return voxel_limits<double>::max(); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<double2>
{
	static double2 min_value() throw() { return make_double2(voxel_limits<double>::min_value()); }
	static double2 max_value() throw() { return make_double2(voxel_limits<double>::max_value()); }
  static double min() throw() { return voxel_limits<double>::min(); }
	static double max() throw() { return voxel_limits<double>::max(); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<double3>
{
	static double3 min_value() throw() { return make_double3(voxel_limits<double>::min_value()); }
	static double3 max_value() throw() { return make_double3(voxel_limits<double>::max_value()); }
  static double min() throw() { return voxel_limits<double>::min(); }
	static double max() throw() { return voxel_limits<double>::max(); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<double4>
{
	static double4 min_value() throw() { return make_double4(voxel_limits<double>::min_value()); }
	static double4 max_value() throw() { return make_double4(voxel_limits<double>::max_value()); }
  static double min() throw() { return voxel_limits<double>::min(); }
	static double max() throw() { return voxel_limits<double>::max(); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<double3x3>
{
	static double3x3 min_value() throw() { return make_double3x3(voxel_limits<double>::min_value()); }
	static double3x3 max_value() throw() { return make_double3x3(voxel_limits<double>::max_value()); }
  static double min() throw() { return voxel_limits<double>::min(); }
	static double max() throw() { return voxel_limits<double>::max(); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<irtkVector3D<float> >
{
	static irtkVector3D<float> min_value() throw()
  { return irtkVector3D<float>(voxel_limits<float>::min_value()); }
	static irtkVector3D<float> max_value() throw()
  { return irtkVector3D<float>(voxel_limits<float>::max_value()); }
  static double min() throw() { return voxel_limits<float>::min(); }
	static double max() throw() { return voxel_limits<float>::max(); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<irtkVector3D<double> >
{
	static irtkVector3D<double> min_value() throw()
  { return irtkVector3D<double>(voxel_limits<double>::min_value()); }
	static irtkVector3D<double> max_value() throw()
  { return irtkVector3D<double>(voxel_limits<double>::max_value()); }
  static double min() throw() { return voxel_limits<double>::min(); }
	static double max() throw() { return voxel_limits<double>::max(); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<irtkVector4D<float> >
{
	static irtkVector4D<float> min_value() throw()
  { return irtkVector4D<float>(voxel_limits<float>::min_value()); }
	static irtkVector4D<float> max_value() throw()
  { return irtkVector4D<float>(voxel_limits<float>::max_value()); }
  static double min() throw() { return voxel_limits<float>::min(); }
	static double max() throw() { return voxel_limits<float>::max(); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<irtkVector4D<double> >
{
	static irtkVector4D<double> min_value() throw()
  { return irtkVector4D<double>(voxel_limits<double>::min_value()); }
	static irtkVector4D<double> max_value() throw()
  { return irtkVector4D<double>(voxel_limits<double>::max_value()); }
  static double min() throw() { return voxel_limits<double>::min(); }
	static double max() throw() { return voxel_limits<double>::max(); }
};

// -----------------------------------------------------------------------------
// Variable length vector type not allowed as actual voxel type of an image
// instance. Only used as voxel type by base class methods and general
// interpolators. Treat voxel type as if it was a scalar type here.
template <> struct voxel_limits<irtkVector>
{
	static irtkVector min_value() throw() { return irtkVector(min()); }
	static irtkVector max_value() throw() { return irtkVector(max()); }
  static double min()           throw() { return static_cast<double>(-1.0e+299); }
	static double max()           throw() { return static_cast<double>( 1.0e+299); }
};

// =============================================================================
// Voxel type information
// =============================================================================

// -----------------------------------------------------------------------------
extern int         DataTypeSize(int);
extern std::string DataTypeName(int);
extern int         ToDataType(const char *);
extern int         ToDataType(const std::string &);

// -----------------------------------------------------------------------------
template <class T>
struct voxel_info : public voxel_limits<T>
{
  /// Scalar type compatible with this voxel type
  typedef T             ScalarType;
  /// Floating point type compatible with this voxel type
  typedef irtkRealPixel RealType;
  /// Number of (vector) elements stored by this voxel
  static int vector_size() throw();
  /// Enumeration value corresponding to voxel type of (vector) elements
  static int element_type() throw();
  /// Enumeration value corresponding to this voxel type
  static int type() throw();
  /// Minimum value that can be represented by this voxel type as double
  /// \note For vector types, this corresponds to the minium value that can
  ///       be represented by each component of the vector.
  using voxel_limits<T>::min;
  /// Maximum value that can be represented by this voxel type as double
  /// \note For vector types, this corresponds to the maxium value that can
  ///       be represented by each component of the vector.
  using voxel_limits<T>::max;
  /// Minimum value that can be represented by this voxel type
  using voxel_limits<T>::min_value;
  /// Maximum value that can be represented by this voxel type
  using voxel_limits<T>::max_value;
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<char>
{
  typedef char            ScalarType;
  typedef irtkRealPixel   RealType;
  static int vector_size()  throw() { return 1; }
  static int element_type() throw() { return IRTK_VOXEL_CHAR; }
  static int type()         throw() { return IRTK_VOXEL_CHAR; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<unsigned char>
{
  typedef unsigned char   ScalarType;
  typedef irtkRealPixel   RealType;
  static int vector_size()  throw() { return 1; }
  static int element_type() throw() { return IRTK_VOXEL_UNSIGNED_CHAR; }
  static int type()         throw() { return IRTK_VOXEL_UNSIGNED_CHAR; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<short>
{
  typedef short           ScalarType;
  typedef irtkRealPixel   RealType;
  static int vector_size()  throw() { return 1; }
  static int element_type() throw() { return IRTK_VOXEL_SHORT; }
  static int type()         throw() { return IRTK_VOXEL_SHORT; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<unsigned short>
{
  typedef unsigned short   ScalarType;
  typedef irtkRealPixel    RealType;
  static int vector_size()  throw() { return 1; }
  static int element_type() throw() { return IRTK_VOXEL_UNSIGNED_SHORT; }
  static int type()         throw() { return IRTK_VOXEL_UNSIGNED_SHORT; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<int>
{
  typedef int             ScalarType;
  typedef irtkRealPixel   RealType;
  static int vector_size()  throw() { return 1; }
  static int element_type() throw() { return IRTK_VOXEL_INT; }
  static int type()         throw() { return IRTK_VOXEL_INT; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<unsigned int>
{
  typedef unsigned int    ScalarType;
  typedef irtkRealPixel   RealType;
  static int vector_size()  throw() { return 1; }
  static int element_type() throw() { return IRTK_VOXEL_UNSIGNED_INT; }
  static int type()         throw() { return IRTK_VOXEL_UNSIGNED_INT; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<float>
{
  typedef float1   Type1;
  typedef float2   Type2;
  typedef float3   Type3;
  typedef float4   Type4;
  typedef float2x2 Type2x2;
  typedef float3x3 Type3x3;
  typedef float3x4 Type3x4;
  typedef float4x4 Type4x4;
  typedef float    ScalarType;
  typedef float    RealType;
  static int vector_size()  throw() { return 1; }
  static int element_type() throw() { return IRTK_VOXEL_FLOAT; }
  static int type()         throw() { return IRTK_VOXEL_FLOAT; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<float1>
{
  typedef float  ScalarType;
  typedef float1 RealType;
  static int vector_size()  throw() { return 1; }
  static int element_type() throw() { return IRTK_VOXEL_FLOAT; }
  static int type()         throw() { return IRTK_VOXEL_FLOAT1; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<float2>
{
  typedef float  ScalarType;
  typedef float2 RealType;
  static int vector_size()  throw() { return 2; }
  static int element_type() throw() { return IRTK_VOXEL_FLOAT; }
  static int type()         throw() { return IRTK_VOXEL_FLOAT2; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<float3>
{
  typedef float  ScalarType;
  typedef float3 RealType;
  static int vector_size()  throw() { return 3; }
  static int element_type() throw() { return IRTK_VOXEL_FLOAT; }
  static int type()         throw() { return IRTK_VOXEL_FLOAT3; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<float4>
{
  typedef float  ScalarType;
  typedef float4 RealType;
  static int vector_size()  throw() { return 4; }
  static int element_type() throw() { return IRTK_VOXEL_FLOAT; }
  static int type()         throw() { return IRTK_VOXEL_FLOAT4; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<float3x3>
{
  typedef float    ScalarType;
  typedef float3x3 RealType;
  static int vector_size()  throw() { return 9; }
  static int element_type() throw() { return IRTK_VOXEL_FLOAT; }
  static int type()         throw() { return IRTK_VOXEL_FLOAT3x3; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<double>
{
  typedef double1   Type1;
  typedef double2   Type2;
  typedef double3   Type3;
  typedef double4   Type4;
  typedef double2x2 Type2x2;
  typedef double3x3 Type3x3;
  typedef double3x4 Type3x4;
  typedef double4x4 Type4x4;
  typedef double    ScalarType;
  typedef double    RealType;
  static int vector_size()  throw() { return 1; }
  static int element_type() throw() { return IRTK_VOXEL_DOUBLE; }
  static int type()         throw() { return IRTK_VOXEL_DOUBLE; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<double1>
{
  typedef double  ScalarType;
  typedef double1 RealType;
  static int vector_size()  throw() { return 1; }
  static int element_type() throw() { return IRTK_VOXEL_DOUBLE; }
  static int type()         throw() { return IRTK_VOXEL_DOUBLE1; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<double2>
{
  typedef double  ScalarType;
  typedef double2 RealType;
  static int vector_size()  throw() { return 2; }
  static int element_type() throw() { return IRTK_VOXEL_DOUBLE; }
  static int type()         throw() { return IRTK_VOXEL_DOUBLE2; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<double3>
{
  typedef double  ScalarType;
  typedef double3 RealType;
  static int vector_size()  throw() { return 3; }
  static int element_type() throw() { return IRTK_VOXEL_DOUBLE; }
  static int type()         throw() { return IRTK_VOXEL_DOUBLE3; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<double4>
{
  typedef double  ScalarType;
  typedef double4 RealType;
  static int vector_size()  throw() { return 4; }
  static int element_type() throw() { return IRTK_VOXEL_DOUBLE; }
  static int type()         throw() { return IRTK_VOXEL_DOUBLE4; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<double3x3>
{
  typedef double    ScalarType;
  typedef double3x3 RealType;
  static int vector_size()  throw() { return 9; }
  static int element_type() throw() { return IRTK_VOXEL_DOUBLE; }
  static int type()         throw() { return IRTK_VOXEL_DOUBLE3x3; }
};

// -----------------------------------------------------------------------------
// Variable length vector type not allowed as actual voxel type of an image
// instance. Only used as voxel type by base class methods and general
// interpolators. Treat voxel type as if it was a scalar type here.
template <> struct voxel_info<irtkVector>
{
  typedef double ScalarType;
  typedef double RealType;
  static int vector_size()  throw() { return 1; }
  static int element_type() throw() { return IRTK_VOXEL_DOUBLE; }
  static int type()         throw() { return IRTK_VOXEL_DOUBLE; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<irtkVector3D<float> >
{
  typedef float                 ScalarType;
  typedef irtkVector3D<float>   RealType;
  static int vector_size()  throw() { return 3; }
  static int element_type() throw() { return IRTK_VOXEL_FLOAT; }
  static int type()         throw() { return IRTK_VOXEL_FLOAT3; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<irtkVector3D<double> >
{
  typedef double                 ScalarType;
  typedef irtkVector3D<double>   RealType;
  static int vector_size()  throw() { return 3; }
  static int element_type() throw() { return IRTK_VOXEL_DOUBLE; }
  static int type()         throw() { return IRTK_VOXEL_DOUBLE3; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<irtkVector4D<float> >
{
  typedef float                 ScalarType;
  typedef irtkVector4D<float>   RealType;
  static int vector_size()  throw() { return 4; }
  static int element_type() throw() { return IRTK_VOXEL_FLOAT; }
  static int type()         throw() { return IRTK_VOXEL_FLOAT4; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<irtkVector4D<double> >
{
  typedef double                 ScalarType;
  typedef irtkVector4D<double>   RealType;
  static int vector_size()  throw() { return 4; }
  static int element_type() throw() { return IRTK_VOXEL_DOUBLE; }
  static int type()         throw() { return IRTK_VOXEL_DOUBLE4; }
};


#endif
