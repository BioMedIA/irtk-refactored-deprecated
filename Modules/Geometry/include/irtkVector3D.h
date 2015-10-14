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

#ifndef _IRTKVECTOR3D_H
#define _IRTKVECTOR3D_H

#include <irtkCommon.h>
#include <irtkPoint.h>


/// Represents a 3D vector
///
/// Must be a primitive type which can be treated as an array of three
/// values of type T such that sizeof(irtkVector<T>) == 3 * sizeof(T).
/// Thus, this primitive vector type may not have any other data members
/// besides the three vector components. This is required especially
/// when irtkVector3D is used as voxel type of an image and further an
/// externally allocated continuous block of 3 * sizeof(T) bytes used
/// internally by the image instance which only reinterprets the memory
/// as consecutive irtkVector3D<T> instances, i.e.,
///
/// \code
/// const int X   = 256;
/// const int Y   = 256;
/// const int Z   = 128;
/// const int num = X * Y * Z;
/// double *data = new double[3 * num];
/// irtkGenericImage<irtkVector3D<double> > image(X, Y, Z, data);
/// \endcode

template <typename T>
struct irtkVector3D
{
  typedef T ComponentType;

  // ---------------------------------------------------------------------------
  // Attributes

  T _x; ///< The x component
  T _y; ///< The y component
  T _z; ///< The z component

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  irtkVector3D();

  /// Construct from scalar
  irtkVector3D(T);

  /// Construct from vector components
  irtkVector3D(T, T, T);

  /// Construct from 3D point
  irtkVector3D(const irtkPoint &);

  /// Assignment operator
  irtkVector3D &operator =(const irtkVector3D &);

  /// Assignment operator
  irtkVector3D &operator =(const irtkPoint &);

  /// Copy constructor
  template <typename T2> irtkVector3D(const irtkVector3D<T2> &);

  // ---------------------------------------------------------------------------
  // Accessors

  /// Number of vector components
  static int Rows() { return 3; }

  /// Set/get vector component at index 0: _x, 1: _y, or 2: _z
  T &operator ()(int);

  /// Get vector component at index 0: _x, 1: _y, or 2: _z
  T operator ()(int) const;

  // ---------------------------------------------------------------------------
  // Vector/integral-valued scalar operators

  /// Assign integral valued scalar
  irtkVector3D &operator =(int);

  /// Add integral valued scalar
  irtkVector3D &operator +=(int);

  /// Subtract integral valued scalar
  irtkVector3D &operator -=(int);

  /// Multiply by integral valued scalar
  irtkVector3D &operator *=(int);

  /// Divide by integral valued scalar
  irtkVector3D &operator /=(int);

  /// Add integral valued scalar to vector
  irtkVector3D operator +(int) const;

  /// Subtract integral valued scalar to vector
  irtkVector3D operator -(int) const;

  /// Multiply vector by integral valued scalar
  irtkVector3D operator *(int) const;

  /// Divide vector by integral valued scalar
  irtkVector3D operator /(int) const;

  // ---------------------------------------------------------------------------
  // Vector/real-valued scalar operators

  /// Assign real valued scalar
  irtkVector3D &operator =(double);

  /// Add real valued scalar
  irtkVector3D &operator +=(double);

  /// Subtract real valued scalar
  irtkVector3D &operator -=(double);

  /// Multiply by real valued scalar
  irtkVector3D &operator *=(double);

  /// Divide by real valued scalar
  irtkVector3D &operator /=(double);

  /// Add real valued scalar to vector
  irtkVector3D operator +(double) const;

  /// Subtract real valued scalar to vector
  irtkVector3D operator -(double) const;

  /// Multiply vector by real valued scalar
  irtkVector3D operator *(double) const;

  /// Divide vector by real valued scalar
  irtkVector3D operator /(double) const;

  // ---------------------------------------------------------------------------
  // Vector/vector operators

  /// Unary negation operator
  irtkVector3D operator -() const;

  /// Assignment from other vector
  template <typename T2> irtkVector3D &operator =(const irtkVector3D<T2> &);

  /// Addition of other vector
  template <typename T2> irtkVector3D &operator +=(const irtkVector3D<T2> &);

  /// Subtraction of other vector
  template <typename T2> irtkVector3D &operator -=(const irtkVector3D<T2> &);

  /// Element-wise multiplication with other vector
  template <typename T2> irtkVector3D &operator *=(const irtkVector3D<T2> &);

  /// Element-wise division by other vector
  template <typename T2> irtkVector3D &operator /=(const irtkVector3D<T2> &);

  /// Addition of two vectors
  template <typename T2> irtkVector3D operator +(const irtkVector3D<T2> &) const;

  /// Subtraction of two vectors
  template <typename T2> irtkVector3D operator -(const irtkVector3D<T2> &) const;

  /// Element-wise multiplication of two vectors
  template <typename T2> irtkVector3D operator *(const irtkVector3D<T2> &) const;

  /// Element-wise division of two vectors
  template <typename T2> irtkVector3D operator /(const irtkVector3D<T2> &) const;

  // ---------------------------------------------------------------------------
  // Vector/integral-valued scalar comparison

  /// Element-wise equality comparison with inegral-valued scalar
  bool operator ==(int) const;

  /// Element-wise inequality comparison with inegral-valued scalar
  bool operator !=(int) const;

  /// Element-wise less than comparison to inegral-valued scalar
  bool operator <(int) const;

  /// Element-wise greater than comparison to inegral-valued scalar
  bool operator >(int) const;

  /// Element-wise less or equal than comparison to inegral-valued scalar
  bool operator <=(int) const;

  /// Element-wise greater or equal than comparison to inegral-valued scalar
  bool operator >=(int) const;

  // ---------------------------------------------------------------------------
  // Vector/real-valued scalar comparison

  /// Element-wise equality comparison with real-valued scalar
  bool operator ==(double) const;

  /// Element-wise inequality comparison with real-valued scalar
  bool operator !=(double) const;

  /// Element-wise less than comparison to real-valued scalar
  bool operator <(double) const;

  /// Element-wise greater than comparison to real-valued scalar
  bool operator >(double) const;

  /// Element-wise less or equal than comparison to real-valued scalar
  bool operator <=(double) const;

  /// Element-wise greater or equal than comparison to real-valued scalar
  bool operator >=(double) const;

  // ---------------------------------------------------------------------------
  // Vector/vector comparison

  /** Operator for testing equality of two vectors. */
  template <typename T2> bool operator ==(const irtkVector3D<T2> &) const;

  /** Operator for testing non-equality of two vector. */
  template <typename T2> bool operator !=(const irtkVector3D<T2> &) const;

  /** Operator for comparing sizes of vectors. */
  template <typename T2> bool operator <(const irtkVector3D<T2> &) const;

  /** Operator for comparing sizes of vectors. */
  template <typename T2> bool operator >(const irtkVector3D<T2> &) const;

  /** Operator for comparing sizes of vectors. */
  template <typename T2> bool operator <=(const irtkVector3D<T2> &) const;

  /** Operator for comparing sizes of vectors. */
  template <typename T2> bool operator >=(const irtkVector3D<T2> &) const;

  // ---------------------------------------------------------------------------
  // Other vector functions

  /// Compute length of vector
  double Length() const;

  /// Normalize vector to length one
  void Normalize();

  /// Cross-product with other vector
  irtkVector3D CrossProduct(const irtkVector3D &) const;

  /// Dot-product with other vector
  double DotProduct(const irtkVector3D &) const;

  /// Cross-product of two vectors
  ///
  /// \deprecated Use v1.CrossProduct(v2) instead which is less to type because
  ///             it does not require the irtkVector<T>:: type prefix and further
  ///             puts the name of the operation in between the arguments.
  static irtkVector3D CrossProduct(const irtkVector3D &, const irtkVector3D &);

  /// Dot-product of two vectors
  ///
  /// \deprecated Use v1.DotProduct(v2) instead which is less to type because
  ///             it does not require the irtkVector<T>:: type prefix and further
  ///             puts the name of the operation in between the arguments.
  static double DotProduct(const irtkVector3D &, const irtkVector3D &);

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T>::irtkVector3D()
{
  _x = _y = _z = T();
}

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T>::irtkVector3D(T s)
{
  _x = _y = _z = s;
}

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T>::irtkVector3D(T x, T y, T z)
{
  _x = x;
  _y = y;
  _z = z;
}

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T>::irtkVector3D(const irtkPoint &p)
{
  _x = p._x;
  _y = p._y;
  _z = p._z;
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline irtkVector3D<T1>::irtkVector3D(const irtkVector3D<T2> &v)
{
  _x = v._x;
  _y = v._y;
  _z = v._z;
}

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T> &irtkVector3D<T>::operator =(const irtkVector3D &v)
{
  _x = v._x;
  _y = v._y;
  _z = v._z;
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T> &irtkVector3D<T>::operator =(const irtkPoint &p)
{
  _x = p._x;
  _y = p._y;
  _z = p._z;
  return *this;
}

// =============================================================================
// Accessors
// =============================================================================

// -----------------------------------------------------------------------------
template <typename T>
inline T& irtkVector3D<T>::operator ()(int i)
{
  switch (i) {
    case 0: return _x;
    case 1: return _y;
    case 2: return _z;
    default:
      cerr << "irtkVector3D::operator(): Invalid index " << i << endl;
      exit(1);
  }
}

// -----------------------------------------------------------------------------
template <typename T>
inline T irtkVector3D<T>::operator ()(int i) const
{
  return const_cast<irtkVector3D<T> *>(this)->operator()(i);
}

// =============================================================================
// 3D vector/integral-valued scalar operators
// =============================================================================

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T>& irtkVector3D<T>::operator =(int s)
{
  _x = s;
  _y = s;
  _z = s;
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T>& irtkVector3D<T>::operator +=(int s)
{
  _x += s;
  _y += s;
  _z += s;
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T>& irtkVector3D<T>::operator -=(int s)
{
  _x -= s;
  _y -= s;
  _z -= s;
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T>& irtkVector3D<T>::operator *=(int s)
{
  _x *= s;
  _y *= s;
  _z *= s;
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T>& irtkVector3D<T>::operator /=(int s)
{
  _x /= s;
  _y /= s;
  _z /= s;
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T> irtkVector3D<T>::operator +(int s) const
{
  irtkVector3D<T> r(*this);
  r += s;
  return r;
}

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T> irtkVector3D<T>::operator -(int s) const
{
  irtkVector3D<T> r(*this);
  r -= s;
  return r;
}

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T> irtkVector3D<T>::operator *(int s) const
{
  irtkVector3D<T> r(*this);
  r *= s;
  return r;
}

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T> irtkVector3D<T>::operator /(int s) const
{
  irtkVector3D<T> r(*this);
  r /= s;
  return r;
}

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T> operator +(int s, const irtkVector3D<T> &v)
{
  return v + s;
}

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T> operator -(int s, const irtkVector3D<T> &v)
{
  return v - s;
}

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T> operator *(int s, const irtkVector3D<T> &v)
{
  return v * s;
}

// =============================================================================
// 3D vector/real-valued scalar operators
// =============================================================================

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T>& irtkVector3D<T>::operator =(double s)
{
  _x = s;
  _y = s;
  _z = s;
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T>& irtkVector3D<T>::operator +=(double s)
{
  _x += s;
  _y += s;
  _z += s;
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T>& irtkVector3D<T>::operator -=(double s)
{
  _x -= s;
  _y -= s;
  _z -= s;
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T>& irtkVector3D<T>::operator *=(double s)
{
  _x *= s;
  _y *= s;
  _z *= s;
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T>& irtkVector3D<T>::operator /=(double s)
{
  _x /= s;
  _y /= s;
  _z /= s;
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T> irtkVector3D<T>::operator +(double s) const
{
  irtkVector3D<T> r(*this);
  r += s;
  return r;
}

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T> irtkVector3D<T>::operator -(double s) const
{
  irtkVector3D<T> r(*this);
  r -= s;
  return r;
}

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T> irtkVector3D<T>::operator *(double s) const
{
  irtkVector3D<T> r(*this);
  r *= s;
  return r;
}

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T> irtkVector3D<T>::operator /(double s) const
{
  irtkVector3D<T> r(*this);
  r /= s;
  return r;
}

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T> operator +(double s, const irtkVector3D<T> &v)
{
  return v + s;
}

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T> operator -(double s, const irtkVector3D<T> &v)
{
  return v - s;
}

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T> operator *(double s, const irtkVector3D<T> &v)
{
  return v * s;
}

// =============================================================================
// 3D vector/vector operators
// =============================================================================

// -----------------------------------------------------------------------------
template <typename T>
inline irtkVector3D<T> irtkVector3D<T>::operator -() const
{
  return irtkVector3D<T>(-_x, -_y, -_z);
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline irtkVector3D<T1> &irtkVector3D<T1>::operator =(const irtkVector3D<T2> &v)
{
  _x = v._x;
  _y = v._y;
  _z = v._z;
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline irtkVector3D<T1> &irtkVector3D<T1>::operator +=(const irtkVector3D<T2>& v)
{
  _x += v._x;
  _y += v._y;
  _z += v._z;
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline irtkVector3D<T1> &irtkVector3D<T1>::operator -=(const irtkVector3D<T2> &v)
{
  _x -= v._x;
  _y -= v._y;
  _z -= v._z;
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline irtkVector3D<T1> &irtkVector3D<T1>::operator *=(const irtkVector3D<T2> &v)
{
  _x *= v._x;
  _y *= v._y;
  _z *= v._z;
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline irtkVector3D<T1> &irtkVector3D<T1>::operator /=(const irtkVector3D<T2> &v)
{
  _x /= v._x;
  _y /= v._y;
  _z /= v._z;
  return *this;
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline irtkVector3D<T1> irtkVector3D<T1>::operator +(const irtkVector3D<T2> &v) const
{
  irtkVector3D<T1> r(*this);
  r += v;
  return r;
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline irtkVector3D<T1> irtkVector3D<T1>::operator -(const irtkVector3D<T2> &v) const
{
  irtkVector3D<T1> r(*this);
  r -= v;
  return r;
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline irtkVector3D<T1> irtkVector3D<T1>::operator *(const irtkVector3D<T2> &v) const
{
  irtkVector3D<T1> r(*this);
  r *= v;
  return r;
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline irtkVector3D<T1> irtkVector3D<T1>::operator /(const irtkVector3D<T2> &v) const
{
  irtkVector3D<T1> r(*this);
  r /= v;
  return r;
}

// =============================================================================
// 3D Vector/integral-valued scalar comparison
// =============================================================================

// -----------------------------------------------------------------------------
template <typename T>
inline bool irtkVector3D<T>::operator ==(int s) const
{
  return (_x == s && _y == s && _z == s);
}

// -----------------------------------------------------------------------------
template <typename T>
inline bool irtkVector3D<T>::operator !=(int s) const
{
  return !(*this == s);
}

// -----------------------------------------------------------------------------
template <typename T>
inline bool irtkVector3D<T>::operator <(int s) const
{
  return (_x < s && _y < s && _z < s);
}

// -----------------------------------------------------------------------------
template <typename T>
inline bool irtkVector3D<T>::operator >(int s) const
{
  return (_x > s && _y > s && _z > s);
}

// -----------------------------------------------------------------------------
template <typename T>
inline bool irtkVector3D<T>::operator <=(int s) const
{
  return (_x <= s && _y <= s && _z <= s);
}

// -----------------------------------------------------------------------------
template <typename T>
inline bool irtkVector3D<T>::operator >=(int s) const
{
  return (_x >= s && _y >= s && _z >= s);
}

// =============================================================================
// 3D Vector/real-valued scalar comparison
// =============================================================================

// -----------------------------------------------------------------------------
template <typename T>
inline bool irtkVector3D<T>::operator ==(double s) const
{
  return (_x == s && _y == s && _z == s);
}

// -----------------------------------------------------------------------------
template <typename T>
inline bool irtkVector3D<T>::operator !=(double s) const
{
  return !(*this == s);
}

// -----------------------------------------------------------------------------
template <typename T>
inline bool irtkVector3D<T>::operator <(double s) const
{
  return (_x < s && _y < s && _z < s);
}

// -----------------------------------------------------------------------------
template <typename T>
inline bool irtkVector3D<T>::operator >(double s) const
{
  return (_x > s && _y > s && _z > s);
}

// -----------------------------------------------------------------------------
template <typename T>
inline bool irtkVector3D<T>::operator <=(double s) const
{
  return (_x <= s && _y <= s && _z <= s);
}

// -----------------------------------------------------------------------------
template <typename T>
inline bool irtkVector3D<T>::operator >=(double s) const
{
  return (_x >= s && _y >= s && _z >= s);
}

// =============================================================================
// 3D Vector/vector comparison
// =============================================================================

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline bool irtkVector3D<T1>::operator ==(const irtkVector3D<T2> &v) const
{
  return ((_z == v._z) && (_y == v._y) && (_x == v._x));
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline bool irtkVector3D<T1>::operator !=(const irtkVector3D<T2> &v) const
{
  return ((_z != v._z) || (_y != v._y) || (_x != v._x));
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline bool irtkVector3D<T1>::operator <(const irtkVector3D<T2> &v) const
{
  return ((_z < v._z) ||
          ((_z == v._z) && (_y < v._y)) ||
          ((_z == v._z) && (_y == v._y) && (_x < v._x)));
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline bool irtkVector3D<T1>::operator >(const irtkVector3D<T2> &v) const
{
  return ((_z > v._z) ||
          ((_z == v._z) && (_y > v._y)) ||
          ((_z == v._z) && (_y == v._y) && (_x > v._x)));
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline bool irtkVector3D<T1>::operator <=(const irtkVector3D<T2> &v) const
{
  return ((*this < v) || (*this == v));
}

// -----------------------------------------------------------------------------
template <typename T1> template <typename T2>
inline bool irtkVector3D<T1>::operator >=(const irtkVector3D<T2> &v) const
{
  return ((*this > v) || (*this == v));
}

// =============================================================================
// 3D vector functions
// =============================================================================

// -----------------------------------------------------------------------------
template <typename T>
inline void irtkVector3D<T>::Normalize()
{
  double length = sqrt(static_cast<double>(_x*_x + _y*_y + _z*_z));
  if (length != .0) (*this) /= length;
}

// -----------------------------------------------------------------------------
template <typename T>
inline double irtkVector3D<T>::Length() const
{
  return sqrt(static_cast<double>(_x*_x + _y*_y + _z*_z));
}

// -----------------------------------------------------------------------------
template<typename T>
inline irtkVector3D<T>
irtkVector3D<T>::CrossProduct(const irtkVector3D<T> &v) const
{
  return irtkVector3D<T>((_y * v._z - _z * v._y),
                         (_z * v._x - _x * v._z),
                         (_x * v._y - _y * v._x));
}

// -----------------------------------------------------------------------------
template<typename T>
inline irtkVector3D<T> irtkVector3D<T>::CrossProduct(const irtkVector3D<T> &v1,
                                                     const irtkVector3D<T> &v2)
{
  return v1.CrossProduct(v2);
}

// -----------------------------------------------------------------------------
template <typename T>
inline double irtkVector3D<T>::DotProduct(const irtkVector3D<T> &v) const
{
  return (_x * v._x + _y * v._y + _z * v._z);
}

// -----------------------------------------------------------------------------
template<typename T>
inline double irtkVector3D<T>::DotProduct(const irtkVector3D<T> &v1,
                                          const irtkVector3D<T> &v2)
{
  return v1.DotProduct(v2);
}

// -----------------------------------------------------------------------------
template <typename T>
irtkVector3D<T> pow(const irtkVector3D<T> &v, int e)
{
  return irtkVector3D<T>(pow(v._x, e), pow(v._y, e), pow(v._z, e));
}

// -----------------------------------------------------------------------------
template <typename T>
irtkVector3D<T> pow(const irtkVector3D<T> &v, double e)
{
  return irtkVector3D<T>(pow(v._x, e), pow(v._y, e), pow(v._z, e));
}

// -----------------------------------------------------------------------------
template <typename T>
irtkVector3D<T> sqrt(const irtkVector3D<T> &v)
{
  return irtkVector3D<T>(sqrt(v._x), sqrt(v._y), sqrt(v._z));
}

// =============================================================================
// Indexed element access
// =============================================================================

// -----------------------------------------------------------------------------
template <class T>
inline T get(const irtkVector3D<T> &v, int n)
{
  switch (n) {
    case 0: return v._x;
    case 1: return v._y;
    case 2: return v._z;
    default:
      cerr << "Invalid 3D vector element index: " << n << endl;
      cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
      exit(1);
  }
}

// -----------------------------------------------------------------------------
template <class T>
inline T put(irtkVector3D<T> &v, int n, const T &value)
{
  switch (n) {
    case 0: v._x = value;
    case 1: v._y = value;
    case 2: v._z = value;
    default:
      cerr << "Invalid 3D vector element index: " << n << endl;
      cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
      exit(1);
  }
}


#endif // __IRTKVECTOR3D_H
