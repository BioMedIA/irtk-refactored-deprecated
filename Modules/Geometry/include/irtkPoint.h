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

#ifndef _IRTKPOINT_H

#define _IRTKPOINT_H

#include <irtkObject.h>


// Forward declaration of irtkVector in particular, as its header
// includes irtkVector3D.h which in turn includes this file again.
// Therefore, only include irtkVector.h after irtkPoint is declared.
class irtkVector;
class irtkMatrix;


/**

  Point class.

*/

class irtkPoint : public irtkObject
{
  irtkObjectMacro(irtkPoint);

public:

  /// x coordinate of Point
  double _x;

  /// y coordinate of Point
  double _y;

  /// z coordinate of Point
  double _z;

  //
  // Constructors and destructor
  //

  /// Constructor
  irtkPoint();

  /// Constructor with three coordinates
  irtkPoint(double, double, double);

  /// Constructor with three coordinates
  irtkPoint(double [3]);

  /// Constructor with Point
  irtkPoint(const irtkPoint &);

  /// Constructor with Vector
  irtkPoint(const irtkVector&);

  /// Default destructor
  virtual ~irtkPoint();

  //
  // Operators for Point
  //

  /// Copy operator for point
  irtkPoint& operator =(const irtkPoint&);

  /// Substraction operator for point
  irtkPoint& operator-=(const irtkPoint&);

  /// Addition operator for point
  irtkPoint& operator+=(const irtkPoint&);

  /// Multiplication operator for point
  irtkPoint& operator*=(const irtkPoint&);

  /// Division operator for point
  irtkPoint& operator/=(const irtkPoint&);

  /// Return result of point substraction
  irtkPoint  operator- (const irtkPoint&) const;

  /// Return result of point addition
  irtkPoint  operator+ (const irtkPoint&) const;

  /// Return result of point multiplication
  irtkPoint  operator* (const irtkPoint&) const;

  /// Return result of point division
  irtkPoint  operator/ (const irtkPoint&) const;

  //
  // Operators for comparison
  //

  /// Comparison operator ==
  int    operator==(const irtkPoint&) const;

  /// Comparison operator != (if USE_STL is defined, negate == operator)
  int    operator!=(const irtkPoint&) const;

  /// Comparison operator <
  int    operator<(const irtkPoint&) const;

  /// Comparison operator >
  int    operator>(const irtkPoint&) const;

  //
  // Operators for double
  //

  /// Substraction of double
  irtkPoint& operator-=(double);

  /// Addition of double
  irtkPoint& operator+=(double);

  /// Multiplication with double
  irtkPoint& operator*=(double);

  /// Division by double
  irtkPoint& operator/=(double);

  // Return result of substraction of double
  irtkPoint  operator- (double) const;

  // Return result of addition of double
  irtkPoint  operator+ (double) const;

  // Return result of multiplication with double
  irtkPoint  operator* (double) const;

  // Return result of division by double
  irtkPoint  operator/ (double) const;

  //
  // Operators for Vector
  //

  /// Copy operator for vectors
  irtkPoint& operator =(const irtkVector&);

  /// Substraction operator for vectors
  irtkPoint& operator-=(const irtkVector&);

  /// Addition operator for vectors
  irtkPoint& operator+=(const irtkVector&);

  /// Multiplication operator for vectors (componentwise)
  irtkPoint& operator*=(const irtkVector&);

  /// Division operator for vectors (componentwise)
  irtkPoint& operator/=(const irtkVector&);

  // Return result of vector substraction
  irtkPoint  operator- (const irtkVector&) const;

  // Return result of vector addition
  irtkPoint  operator+ (const irtkVector&) const;

  // Return result of vector multiplication
  irtkPoint  operator* (const irtkVector&) const;

  // Return result of vector division
  irtkPoint  operator/ (const irtkVector&) const;

  //
  // Operators for Matrix
  //

  /// Point multiplication operator for matrices
  irtkPoint& operator*=(const irtkMatrix&);

  /// Return result from Matrix multiplication
  irtkPoint  operator* (const irtkMatrix&) const;

  //
  // Distance methods
  //

  /// Squared distance from origin
  double  SquaredDistance() const;

  /// Squared distance from point
  double  SquaredDistance(const irtkPoint&) const;

  /// Distance from origin
  double  Distance() const;

  /// Distance from point
  double  Distance(const irtkPoint&) const;

  //
  // I/O methods
  //

  // Interface to output stream
  friend ostream& operator<< (ostream&, const irtkPoint&);

  /// Interface to input stream
  friend istream& operator>> (istream&, irtkPoint&);
};

// =============================================================================
// Inline definitions
// =============================================================================

#include <irtkVector.h>
#include <irtkMatrix.h>


inline irtkPoint::irtkPoint()
{
  _x = 0;
  _y = 0;
  _z = 0;
}

inline irtkPoint::irtkPoint(double x, double y, double z)
{
  _x = x;
  _y = y;
  _z = z;
}

inline irtkPoint::irtkPoint(double p[3])
{
  _x = p[0];
  _y = p[1];
  _z = p[2];
}

inline irtkPoint::irtkPoint(const irtkPoint& p) : irtkObject(p)
{
  _x = p._x;
  _y = p._y;
  _z = p._z;
}

inline irtkPoint::irtkPoint(const irtkVector& v)
{
  if ((v.Rows() < 0) || (v.Rows() > 3)) {
    cerr << "irtkPoint::irtkPoint(const irtkVector&) Illegal dimension: " << v.Rows() << endl;
    exit(1);
  } else {
    if (v.Rows() == 1) {
      _x = v(0);
      _y = 0;
      _z = 0;
    }
    if (v.Rows() == 2) {
      _x = v(0);
      _y = v(1);
      _z = 0;
    }
    if (v.Rows() == 3) {
      _x = v(0);
      _y = v(1);
      _z = v(2);
    }
  }
}

inline irtkPoint::~irtkPoint()
{
}

inline irtkPoint& irtkPoint::operator =(const irtkPoint& p)
{
  _x = p._x;
  _y = p._y;
  _z = p._z;
  return *this;
}

inline irtkPoint& irtkPoint::operator+=(const irtkPoint& p)
{
  _x += p._x;
  _y += p._y;
  _z += p._z;
  return *this;
}

inline irtkPoint& irtkPoint::operator-=(const irtkPoint& p)
{
  _x -= p._x;
  _y -= p._y;
  _z -= p._z;
  return *this;
}

inline irtkPoint& irtkPoint::operator*=(const irtkPoint& p)
{
  _x *= p._x;
  _y *= p._y;
  _z *= p._z;
  return *this;
}

inline irtkPoint& irtkPoint::operator/=(const irtkPoint& p)
{
  _x /= p._x;
  _y /= p._y;
  _z /= p._z;
  return *this;
}

inline irtkPoint irtkPoint::operator+ (const irtkPoint& p) const
{
  irtkPoint tmp;

  tmp._x = _x + p._x;
  tmp._y = _y + p._y;
  tmp._z = _z + p._z;
  return tmp;
}

inline irtkPoint irtkPoint::operator- (const irtkPoint& p) const
{
  irtkPoint tmp;

  tmp._x = _x - p._x;
  tmp._y = _y - p._y;
  tmp._z = _z - p._z;
  return tmp;
}

inline irtkPoint irtkPoint::operator* (const irtkPoint& p) const
{
  irtkPoint tmp;

  tmp._x = _x * p._x;
  tmp._y = _y * p._y;
  tmp._z = _z * p._z;
  return tmp;
}

inline irtkPoint irtkPoint::operator/ (const irtkPoint& p) const
{
  irtkPoint tmp;

  tmp._x = _x / p._x;
  tmp._y = _y / p._y;
  tmp._z = _z / p._z;
  return tmp;
}

inline int irtkPoint::operator==(irtkPoint const &p) const
{
  return ((_x == p._x) && (_y == p._y) && (_z == p._z));
}

inline int irtkPoint::operator!=(irtkPoint const &p) const
{
  return ((_x != p._x) || (_y != p._y) || (_z != p._z));
}

inline int irtkPoint::operator<(irtkPoint const& p) const
{
  return ((_x < p._x) && (_y < p._y) && (_z < p._z));
}

inline int irtkPoint::operator>(irtkPoint const& p) const
{
  return ((_x > p._x) && (_y > p._y) && (_z > p._z));
}

inline irtkPoint& irtkPoint::operator+=(double x)
{
  _x += x;
  _y += x;
  _z += x;
  return *this;
}

inline irtkPoint& irtkPoint::operator-=(double x)
{
  _x -= x;
  _y -= x;
  _z -= x;
  return *this;
}

inline irtkPoint& irtkPoint::operator/=(double x)
{
  _x /= x;
  _y /= x;
  _z /= x;
  return *this;
}

inline irtkPoint& irtkPoint::operator*=(double x)
{
  _x *= x;
  _y *= x;
  _z *= x;
  return *this;
}

inline irtkPoint irtkPoint::operator+ (double x) const
{
  irtkPoint p;

  p._x = _x + x;
  p._y = _y + x;
  p._z = _z + x;
  return p;
}

inline irtkPoint irtkPoint::operator- (double x) const
{
  irtkPoint p;

  p._x = _x - x;
  p._y = _y - x;
  p._z = _z - x;
  return p;
}

inline irtkPoint irtkPoint::operator* (double x) const
{
  irtkPoint p;

  p._x = _x * x;
  p._y = _y * x;
  p._z = _z * x;
  return p;
}

inline irtkPoint irtkPoint::operator/ (double x) const
{
  irtkPoint p;

  p._x = _x / x;
  p._y = _y / x;
  p._z = _z / x;
  return p;
}

inline irtkPoint& irtkPoint::operator+=(const irtkVector& v)
{
  if (v.Rows() != 3) {
    cerr << "irtkPoint::operator+=(const irtkVector& v): Size mismatch" << endl;
    exit(1);
  }
  _x += v(0);
  _y += v(1);
  _z += v(2);
  return *this;
}

inline irtkPoint& irtkPoint::operator-=(const irtkVector& v)
{
  if (v.Rows() != 3) {
    cerr << "irtkPoint::operator-=(const irtkVector& v): Size mismatch" << endl;
    exit(1);
  }
  _x -= v(0);
  _y -= v(1);
  _z -= v(2);
  return *this;
}

inline irtkPoint& irtkPoint::operator*=(const irtkVector& v)
{
  if (v.Rows() != 3) {
    cerr << "irtkPoint::operator*=(const irtkVector& v): Size mismatch" << endl;
    exit(1);
  }
  _x *= v(0);
  _y *= v(1);
  _z *= v(2);
  return *this;
}

inline irtkPoint& irtkPoint::operator/=(const irtkVector& v)
{
  if (v.Rows() != 3) {
    cerr << "irtkPoint::operator/=(const irtkVector& v): Size mismatch" << endl;
    exit(1);
  }
  _x /= v(0);
  _y /= v(1);
  _z /= v(2);
  return *this;
}

inline irtkPoint irtkPoint::operator+ (const irtkVector& v) const
{
  irtkPoint tmp;

  if (v.Rows() != 3) {
    cerr << "irtkPoint::operator+(const irtkVector& v): Size mismatch" << endl;
    exit(1);
  }
  tmp._x = _x + v(0);
  tmp._y = _y + v(1);
  tmp._z = _z + v(2);
  return tmp;
}

inline irtkPoint irtkPoint::operator- (const irtkVector& v) const
{
  irtkPoint tmp;

  if (v.Rows() != 3) {
    cerr << "irtkPoint::operator-(const irtkVector& v): Size mismatch" << endl;
    exit(1);
  }
  tmp._x = _x - v(0);
  tmp._y = _y - v(1);
  tmp._z = _z - v(2);
  return tmp;
}

inline irtkPoint irtkPoint::operator* (const irtkVector& v) const
{
  irtkPoint tmp;

  if (v.Rows() != 3) {
    cerr << "irtkPoint::operator*(const irtkVector& v): Size mismatch" << endl;
    exit(1);
  }
  tmp._x = _x * v(0);
  tmp._y = _y * v(1);
  tmp._z = _z * v(2);
  return tmp;
}

inline irtkPoint irtkPoint::operator/ (const irtkVector& v) const
{
  irtkPoint tmp;

  if (v.Rows() != 3) {
    cerr << "irtkPoint::operator/(const irtkVector& v): Size mismatch" << endl;
    exit(1);
  }
  tmp._x = _x / v(0);
  tmp._y = _y / v(1);
  tmp._z = _z / v(2);
  return tmp;
}

inline irtkPoint irtkPoint::operator* (const irtkMatrix& m) const
{
  irtkPoint tmp;

  if ((m.Rows() != 4) && (m.Cols() != 4)) {
    cerr << "irtkPoint::operator*(const irtkMatrix& m): Size mismatch" << endl;
    exit(1);
  }
  tmp._x = _x * m(0, 0) + _y * m(0, 1) + _z * m(0, 2) + m(0, 3);
  tmp._y = _x * m(1, 0) + _y * m(1, 1) + _z * m(1, 2) + m(1, 3);
  tmp._z = _x * m(2, 0) + _y * m(2, 1) + _z * m(2, 2) + m(2, 3);

  return tmp;
}

inline irtkPoint& irtkPoint::operator*=(const irtkMatrix& m)
{
  irtkPoint tmp;

  if ((m.Rows() != 4) && (m.Cols() != 4)) {
    cerr << "irtkPoint::operator*(const irtkMatrix& m): Size mismatch" << endl;
    exit(1);
  }
  tmp._x = _x * m(0, 0) + _y * m(0, 1) + _z * m(0, 2) + m(0, 3);
  tmp._y = _x * m(1, 0) + _y * m(1, 1) + _z * m(1, 2) + m(1, 3);
  tmp._z = _x * m(2, 0) + _y * m(2, 1) + _z * m(2, 2) + m(2, 3);
  *this  = tmp;

  return *this;
}

inline double irtkPoint::SquaredDistance() const
{
  return _x*_x + _y*_y  + _z *_z;
}

inline double irtkPoint::SquaredDistance(const irtkPoint &p) const
{
  const double dx = _x - p._x;
  const double dy = _y - p._y;
  const double dz = _z - p._z;
  return dx*dx + dy*dy + dz*dz;
}

inline double irtkPoint::Distance() const
{
  return sqrt(SquaredDistance());
}

inline double irtkPoint::Distance(const irtkPoint &p) const
{
  return sqrt(SquaredDistance(p));
}

#endif
