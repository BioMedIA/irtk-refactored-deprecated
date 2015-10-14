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

#ifndef IRTKVOXELCAST_H_
#define IRTKVOXELCAST_H_

#include <irtkVector.h>


// Overloading the C++ conversion operators for vector types would be possible
// as well, however, this can lead to ambiguitis in mathematical expressions
// when a vector can be constructed from a scalar and converted to a scalar at
// the same time. Then it is no longer clear whether the arithmetic operation
// should be performed between scalars, vectors, or a mix of both. Therefore,
// only support the needed voxel type conversions as template specializations
// of our own voxel_cast template function.

// -----------------------------------------------------------------------------
/// Auxiliary template class for partial voxel_cast specialization
template <class TIn, class TOut>
struct VoxelCaster
{
  /// By default, use static_cast to convert voxel value from one type to another.
  /// If the value is out of the range which can be represented by the output
  /// type, the minium/maximum value of the output type is returned instead.
  static TOut Convert(const TIn &value)
  {
    if      (static_cast<double>(value) < voxel_limits<TOut>::min()) return voxel_limits<TOut>::min_value();
    else if (static_cast<double>(value) > voxel_limits<TOut>::max()) return voxel_limits<TOut>::max_value();
    else                                                             return static_cast <TOut>(value);
  }
};

// -----------------------------------------------------------------------------
template <class T>
struct VoxelCaster<T, T>
{
  static T Convert(const T &value)
  {
    return value;
  }
};

// -----------------------------------------------------------------------------
template <class TIn>
struct VoxelCaster<TIn, float1>
{
  static float1 Convert(const TIn &value)
  {
    return make_float1(value);
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<float1, TOut>
{
  static TOut Convert(const float1 &value)
  {
    return VoxelCaster<float, TOut>::Convert(value.x);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float1, float1>
{
  static float1 Convert(const float1 &value)
  {
    return value;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float1, double1>
{
  static double1 Convert(const float1 &value)
  {
    return make_double1(value);
  }
};

// -----------------------------------------------------------------------------
template <class TIn>
struct VoxelCaster<TIn, float2>
{
  static float2 Convert(const TIn &value)
  {
    return make_float2(value);
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<float2, TOut>
{
  static TOut Convert(const float2 &)
  {
    cerr << "Cannot cast 2D vector to a scalar!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float2, float2>
{
  static float2 Convert(const float2 &value)
  {
    return value;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float2, double2>
{
  static double2 Convert(const float2 &value)
  {
    return make_double2(value);
  }
};

// -----------------------------------------------------------------------------
template <class TIn>
struct VoxelCaster<TIn, float3>
{
  static float3 Convert(const TIn &value)
  {
    return make_float3(value);
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<float3, TOut>
{
  static TOut Convert(const float3 &)
  {
    cerr << "Cannot cast 3D vector to a scalar!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float3, float3>
{
  static float3 Convert(const float3 &value)
  {
    return value;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float3, double3>
{
  static double3 Convert(const float3 &value)
  {
    return make_double3(value);
  }
};

// -----------------------------------------------------------------------------
template <class TIn>
struct VoxelCaster<TIn, float4>
{
  static float4 Convert(const TIn &value)
  {
    return make_float4(value);
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<float4, TOut>
{
  static TOut Convert(const float4 &)
  {
    cerr << "Cannot cast 4D vector to a scalar!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float4, float4>
{
  static float4 Convert(const float4 &value)
  {
    return value;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float4, double4>
{
  static double4 Convert(const float4 &value)
  {
    return make_double4(value);
  }
};

// -----------------------------------------------------------------------------
template <class TIn>
struct VoxelCaster<TIn, double1>
{
  static double1 Convert(const TIn &value)
  {
    return make_double1(value);
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<double1, TOut>
{
  static TOut Convert(const double1 &value)
  {
    return VoxelCaster<double, TOut>::Convert(value.x);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double1, float1>
{
  static float1 Convert(const double1 &value)
  {
    return make_float1(value);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double1, double1>
{
  static double1 Convert(const double1 &value)
  {
    return value;
  }
};

// -----------------------------------------------------------------------------
template <class TIn>
struct VoxelCaster<TIn, double2>
{
  static double2 Convert(const TIn &value)
  {
    return make_double2(value);
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<double2, TOut>
{
  static TOut Convert(const double2 &)
  {
    cerr << "Cannot cast 2D vector to a scalar!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double2, float2>
{
  static float2 Convert(const double2 &value)
  {
    return make_float2(value);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double2, double2>
{
  static double2 Convert(const double2 &value)
  {
    return value;
  }
};

// -----------------------------------------------------------------------------
template <class TIn>
struct VoxelCaster<TIn, double3>
{
  static double3 Convert(const TIn &value)
  {
    return make_double3(value);
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<double3, TOut>
{
  static TOut Convert(const double3 &)
  {
    cerr << "Cannot cast 3D vector to a scalar!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double3, float3>
{
  static float3 Convert(const double3 &value)
  {
    return make_float3(value);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double3, double3>
{
  static double3 Convert(const double3 &value)
  {
    return value;
  }
};

// -----------------------------------------------------------------------------
template <class TIn>
struct VoxelCaster<TIn, double4>
{
  static double4 Convert(const TIn &value)
  {
    return make_double4(value);
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<double4, TOut>
{
  static TOut Convert(const double4 &)
  {
    cerr << "Cannot cast 4D vector to a scalar!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double4, float4>
{
  static float4 Convert(const double4 &value)
  {
    return make_float4(value);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double4, double4>
{
  static double4 Convert(const double4 &value)
  {
    return value;
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<irtkVector, TOut>
{
  static TOut Convert(const irtkVector &value)
  {
    if (value.Rows() == 1) return VoxelCaster<double, TOut>::Convert(value(0));
    cerr << "Can only cast vector with exactly one element to a scalar!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <class TIn>
struct VoxelCaster<TIn, irtkVector>
{
  static irtkVector Convert(const TIn &value)
  {
    return irtkVector(1, VoxelCaster<TIn, double>::Convert(value));
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<irtkVector, irtkVector>
{
  static irtkVector Convert(const irtkVector &value)
  {
    return value;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<irtkVector, float1>
{
  static float1 Convert(const irtkVector &value)
  {
    if (value.Rows() == 1) return make_float1(value(0));
    cerr << "Can only cast vector with exactly one element to a 1D vector!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float1, irtkVector>
{
  static irtkVector Convert(const float1 &value)
  {
    irtkVector v(1);
    v(0) = value.x;
    return v;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<irtkVector, float2>
{
  static float2 Convert(const irtkVector &value)
  {
    if (value.Rows() == 2) return make_float2(value(0), value(1));
    cerr << "Can only cast vector with exactly two elements to a 2D vector!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float2, irtkVector>
{
  static irtkVector Convert(const float2 &value)
  {
    irtkVector v(2);
    v(0) = value.x;
    v(1) = value.y;
    return v;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<irtkVector, float3>
{
  static float3 Convert(const irtkVector &value)
  {
    if (value.Rows() == 3) return make_float3(value(0), value(1), value(2));
    cerr << "Can only cast vector with exactly three elements to a 3D vector!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float3, irtkVector>
{
  static irtkVector Convert(const float3 &value)
  {
    irtkVector v(3);
    v(0) = value.x;
    v(1) = value.y;
    v(2) = value.z;
    return v;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<irtkVector, float4>
{
  static float4 Convert(const irtkVector &value)
  {
    if (value.Rows() == 4) return make_float4(value(0), value(1), value(2), value(3));
    cerr << "Can only cast vector with exactly four elements to a 4D vector!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float4, irtkVector>
{
  static irtkVector Convert(const float4 &value)
  {
    irtkVector v(4);
    v(0) = value.x;
    v(1) = value.y;
    v(2) = value.z;
    v(3) = value.w;
    return v;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<irtkVector, double1>
{
  static double1 Convert(const irtkVector &value)
  {
    if (value.Rows() == 1) return make_double1(value(0));
    cerr << "Can only cast vector with exactly one element to a 1D vector!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double1, irtkVector>
{
  static irtkVector Convert(const double1 &value)
  {
    irtkVector v(1);
    v(0) = value.x;
    return v;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<irtkVector, double2>
{
  static double2 Convert(const irtkVector &value)
  {
    if (value.Rows() == 2) return make_double2(value(0), value(1));
    cerr << "Can only cast vector with exactly two elements to a 2D vector!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double2, irtkVector>
{
  static irtkVector Convert(const double2 &value)
  {
    irtkVector v(2);
    v(0) = value.x;
    v(1) = value.y;
    return v;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<irtkVector, double3>
{
  static double3 Convert(const irtkVector &value)
  {
    if (value.Rows() == 3) return make_double3(value(0), value(1), value(2));
    cerr << "Can only cast vector with exactly three elements to a 3D vector!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double3, irtkVector>
{
  static irtkVector Convert(const double3 &value)
  {
    irtkVector v(3);
    v(0) = value.x;
    v(1) = value.y;
    v(2) = value.z;
    return v;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<irtkVector, double4>
{
  static double4 Convert(const irtkVector &value)
  {
    if (value.Rows() == 4) return make_double4(value(0), value(1), value(2), value(3));
    cerr << "Can only cast vector with exactly four elements to a 4D vector!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double4, irtkVector>
{
  static irtkVector Convert(const double4 &value)
  {
    irtkVector v(4);
    v(0) = value.x;
    v(1) = value.y;
    v(2) = value.z;
    v(3) = value.w;
    return v;
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<irtkVector, irtkVector3D<TOut> >
{
  static irtkVector3D<TOut> Convert(const irtkVector &value)
  {
    if (value.Rows() == 3) {
      return irtkVector3D<TOut>(VoxelCaster<TOut, double>::Convert(value(0)),
                                VoxelCaster<TOut, double>::Convert(value(1)),
                                VoxelCaster<TOut, double>::Convert(value(2)));
    }
    cerr << "Can only cast vector with exactly three elements to a 3D vector!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<irtkVector, irtkVector4D<TOut> >
{
  static irtkVector4D<TOut> Convert(const irtkVector &value)
  {
    if (value.Rows() == 4) {
      return irtkVector4D<TOut>(VoxelCaster<TOut, double>::Convert(value(0)),
                                VoxelCaster<TOut, double>::Convert(value(1)),
                                VoxelCaster<TOut, double>::Convert(value(2)),
                                VoxelCaster<TOut, double>::Convert(value(3)));
    }
    cerr << "Can only cast vector with exactly four elements to a 4D vector!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <class TOut, class TIn>
struct VoxelCaster<irtkVector3D<TIn>, TOut>
{
  static TOut Convert(const irtkVector3D<TIn> &)
  {
    cerr << "Cannot cast 3D vector to a scalar!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <class TOut, class TIn>
struct VoxelCaster<irtkVector3D<TIn>, irtkVector3D<TOut> >
{
  static irtkVector3D<TOut> Convert(const irtkVector3D<TIn> &value)
  {
    return irtkVector3D<TOut>(VoxelCaster<TOut, TIn>::Convert(value._x),
                              VoxelCaster<TOut, TIn>::Convert(value._y),
                              VoxelCaster<TOut, TIn>::Convert(value._z));
  }
};

// -----------------------------------------------------------------------------
template <class T>
struct VoxelCaster<irtkVector3D<T>, irtkVector3D<T> >
{
  static irtkVector3D<T> Convert(const irtkVector3D<T> &value)
  {
    return value;
  }
};

// -----------------------------------------------------------------------------
template <class TIn>
struct VoxelCaster<irtkVector3D<TIn>, irtkVector>
{
  static irtkVector Convert(const irtkVector3D<TIn> &value)
  {
    irtkVector v(3);
    v.Put(value);
    return v;
  }
};

// -----------------------------------------------------------------------------
template <class TOut, class TIn>
struct VoxelCaster<irtkVector4D<TIn>, TOut>
{
  static TOut Convert(const irtkVector4D<TIn> &)
  {
    cerr << "Cannot cast 4D vector to a scalar!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <class TOut, class TIn>
struct VoxelCaster<irtkVector4D<TIn>, irtkVector4D<TOut> >
{
  static irtkVector4D<TOut> Convert(const irtkVector4D<TIn> &value)
  {
    return irtkVector4D<TOut>(VoxelCaster<TOut, TIn>::Convert(value._x),
                              VoxelCaster<TOut, TIn>::Convert(value._y),
                              VoxelCaster<TOut, TIn>::Convert(value._z),
                              VoxelCaster<TOut, TIn>::Convert(value._t));
  }
};

// -----------------------------------------------------------------------------
template <class T>
struct VoxelCaster<irtkVector4D<T>, irtkVector4D<T> >
{
  static irtkVector4D<T> Convert(const irtkVector4D<T> &value)
  {
    return value;
  }
};

// -----------------------------------------------------------------------------
template <class TIn>
struct VoxelCaster<irtkVector4D<TIn>, irtkVector>
{
  static irtkVector Convert(const irtkVector4D<TIn> &value)
  {
    irtkVector v(4);
    v.Put(value);
    return v;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double3x3, float3x3>
{
  static float3x3 Convert(const double3x3 &value)
  {
    return make_float3x3(value);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float3x3, double3x3>
{
  static double3x3 Convert(const float3x3 &value)
  {
    return make_double3x3(value);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float3x3, float3x3>
{
  static float3x3 Convert(const float3x3 &value)
  {
    return value;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double3x3, double3x3>
{
  static double3x3 Convert(const double3x3 &value)
  {
    return value;
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<float3x3, TOut>
{
  static TOut Convert(const float3x3 &)
  {
    cerr << "Cannot cast 3x3 matrix to a scalar!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<double3x3, TOut>
{
  static TOut Convert(const double3x3 &)
  {
    cerr << "Cannot cast 3x3 matrix to a scalar!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float3x3, irtkVector>
{
  static irtkVector Convert(const float3x3 &value)
  {
    irtkVector v(9);
    v(0) = value.a.x;
    v(1) = value.a.y;
    v(2) = value.a.z;
    v(3) = value.b.x;
    v(4) = value.b.y;
    v(5) = value.b.z;
    v(6) = value.c.x;
    v(7) = value.c.y;
    v(8) = value.c.z;
    return v;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double3x3, irtkVector>
{
  static irtkVector Convert(const double3x3 &value)
  {
    irtkVector v(9);
    v(0) = value.a.x;
    v(1) = value.a.y;
    v(2) = value.a.z;
    v(3) = value.b.x;
    v(4) = value.b.y;
    v(5) = value.b.z;
    v(6) = value.c.x;
    v(7) = value.c.y;
    v(8) = value.c.z;
    return v;
  }
};

// -----------------------------------------------------------------------------
template <class TIn>
struct VoxelCaster<TIn, float3x3>
{
  static float3x3 Convert(const TIn &value)
  {
    return make_float3x3(value);
  }
};

// -----------------------------------------------------------------------------
template <class TIn>
struct VoxelCaster<TIn, double3x3>
{
  static double3x3 Convert(const TIn &value)
  {
    return make_double3x3(value);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<irtkVector, float3x3>
{
  static float3x3 Convert(const irtkVector &v)
  {
    float3x3 m;
    if (v.Rows() == 9) {
      m.a.x = v(0);
      m.a.y = v(1);
      m.a.z = v(2);
      m.b.x = v(3);
      m.b.y = v(4);
      m.b.z = v(5);
      m.c.x = v(6);
      m.c.y = v(7);
      m.c.z = v(8);
    } else if (v.Rows() == 6) {
      m.a.x = v(0);
      m.a.y = v(1);
      m.a.z = v(2);
      m.b.x = m.a.y;
      m.b.y = v(3);
      m.b.z = v(4);
      m.c.x = m.a.z;
      m.c.y = m.b.z;
      m.c.z = v(5);
    } else {
      cerr << "Can only cast vector of size 6 or 9 to a 3x3 matrix!" << endl;
      cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
      exit(1);
    }
    return m;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<irtkVector, double3x3>
{
  static double3x3 Convert(const irtkVector &v)
  {
    double3x3 m;
    if (v.Rows() == 9) {
      m.a.x = v(0);
      m.a.y = v(1);
      m.a.z = v(2);
      m.b.x = v(3);
      m.b.y = v(4);
      m.b.z = v(5);
      m.c.x = v(6);
      m.c.y = v(7);
      m.c.z = v(8);
    } else if (v.Rows() == 6) {
      m.a.x = v(0);
      m.a.y = v(1);
      m.a.z = v(2);
      m.b.x = m.a.y;
      m.b.y = v(3);
      m.b.z = v(4);
      m.c.x = m.a.z;
      m.c.y = m.b.z;
      m.c.z = v(5);
    } else {
      cerr << "Can only cast vector of size 6 or 9 to a 3x3 matrix!" << endl;
      cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
      exit(1);
    }
    return m;
  }
};

// -----------------------------------------------------------------------------
template <class TOut, class TIn>
TOut voxel_cast(const TIn &value)
{
  return VoxelCaster<TIn, TOut>::Convert(value);
}


#endif
