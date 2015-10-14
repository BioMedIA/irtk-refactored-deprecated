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

#ifndef _IRTKVECTOR_H

#define _IRTKVECTOR_H

#include <irtkCommon.h>

#ifdef HAVE_EIGEN
#  include <Eigen/Core>
#endif

#ifdef HAVE_MATLAB
#  include <mclmcrrt.h>
#endif

// Forward declaration of specialized vector types included after
// declaration of irtkVector for definition of inline functions
template <class T> struct irtkVector3D;
template <class T> struct irtkVector4D;

/**

  Vector class.

*/

class irtkVector : public irtkObject
{
  irtkObjectMacro(irtkVector);

  // ---------------------------------------------------------------------------
  // Data members

protected:

  /// Number of rows
  int _rows;

  /// Vector elements
  double *_vector;

public:

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  irtkVector();

  /// Constructor for given row dimensions
  explicit irtkVector(int);

  /// Constructor for given row dimensions
  irtkVector(int, double);

  /// Constructor for given row dimensions
  irtkVector(int, double *);

  /// Copy constructor
  irtkVector(const irtkVector &);

  /// Construct vector from 3D vector
  template <class T> irtkVector(const irtkVector3D<T> &);

  /// Construct vector from 4D vector
  template <class T> irtkVector(const irtkVector4D<T> &);

  /// Destructor
  ~irtkVector();

  /// Intialize matrix with number of rows
  void Initialize(int);

  /// Intialize matrix with number of rows
  void Initialize(int, double);

  /// Intialize matrix with number of rows
  void Initialize(int, double *);

  /// Change size of vector, preserving existing rows
  void Resize(int, double = .0);

  /// Free vector
  void Clear();

  /// Initialize from 3D vector
  template <class T> irtkVector &Put(const irtkVector3D<T> &);

  /// Initialize from 4D vector
  template <class T> irtkVector &Put(const irtkVector4D<T> &);

  /// Whether vector is non-empty, i.e., initialized and number of rows greater zero
  operator bool() const;

  // ---------------------------------------------------------------------------
  // Vector access functions

  /// Returns number of rows
  int Rows() const;

  /// Puts vector value
  void Put(int, double);

  /// Gets vector value
  const double &Get(int) const;

  // ---------------------------------------------------------------------------
  // Element access

  /// Get pointer to linear memory which stores vector elements
  double *RawPointer(int r = 0);

  /// Get pointer to linear memory which stores vector elements
  const double *RawPointer(int r = 0) const;

  /// Puts vector value
  double &operator()(int);

  /// Gets vector value
  const double &operator()(int) const;

  // ---------------------------------------------------------------------------
  // Vector/scalar operations

  /// Assignment of double
  irtkVector &operator =(double);

  /// Subtraction of a double
  irtkVector &operator-=(double);

  /// Addition of a double
  irtkVector &operator+=(double);

  /// Multiplication with a double
  irtkVector &operator*=(double);

  /// Division by a double
  irtkVector &operator/=(double);

  /// Return result of subtraction of a double
  irtkVector operator- (double) const;

  /// Return result of addition of a double
  irtkVector operator+ (double) const;

  /// Return result of multiplication with a double
  irtkVector operator* (double) const;

  /// Return result of division by a double
  irtkVector operator/ (double) const;

  // ---------------------------------------------------------------------------
  // Element-wise vector/vector operations

  /// Unary negation operator
  irtkVector operator -() const;

  /// Vector copy operator
  irtkVector &operator =(const irtkVector &);

  /// Vector subtraction operator
  irtkVector &operator-=(const irtkVector &);

  /// Vector addition operator
  irtkVector &operator+=(const irtkVector &);

  /// Vector componentwise multiplication operator (no scalar nor cross product)
  irtkVector &operator*=(const irtkVector &);

  /// Vector componentwise division operator
  irtkVector &operator/=(const irtkVector &);

  /// Return result for vector subtraction
  irtkVector operator- (const irtkVector &) const;

  /// Return result for vector addition
  irtkVector operator+ (const irtkVector &) const;

  /// Return result for componentwise vector multiplication (no scalar nor cross product)
  irtkVector operator* (const irtkVector &) const;

  /// Return result for componentwise vector division
  irtkVector operator/ (const irtkVector &) const;

  // ---------------------------------------------------------------------------
  // Comparison

  /// Comparison operator ==
  bool operator==(const irtkVector &) const;

#ifndef USE_STL
  /// Comparison operator != (if USE_STL is defined, negate == operator)
  bool operator!=(const irtkVector &) const;
#endif

  /// Comparison operator <
  bool operator<(const irtkVector &) const;

  // ---------------------------------------------------------------------------
  // Vector products

  /// Scalar/dot product
  double ScalarProduct(const irtkVector &) const;

  /// Scalar/dot product
  double DotProduct(const irtkVector &) const;

  /// Vector/cross product
  irtkVector CrossProduct(const irtkVector &) const;

  // ---------------------------------------------------------------------------
  // Other vector functions

  /// Returns norm of a vector
  double Norm() const;

  /// Normalize vector
  irtkVector &Normalize();

  /// Replace each vector element by its inverse
  irtkVector &Inverse();

  /// Permute vector elements
  void PermuteRows(vector<int>);

  /// Permute vector elements (cf. PermuteRows)
  void Permute(const vector<int> &);

  // ---------------------------------------------------------------------------
  // I/O

  /// Print vector
  void Print(irtkIndent = 0) const;

  /// Read vector from file
  void Read(const char *);

  /// Write vector to file
  void Write(const char *) const;

  /// Interface to output stream
  friend ostream &operator<< (ostream &, const irtkVector &);

  /// Interface to input stream
  friend istream &operator>> (istream &, irtkVector &);

  /// Interface to output stream
  friend irtkCofstream &operator<< (irtkCofstream &, const irtkVector &);

  /// Interface to input stream
  friend irtkCifstream &operator>> (irtkCifstream &, irtkVector &);

#ifdef HAVE_EIGEN

  /// Conversion to Eigen library vector
  void Vector2Eigen(Eigen::VectorXd &) const;

  /// Conversion from Eigen library vector
  void Eigen2Vector(Eigen::VectorXd &);

#endif

#ifdef HAVE_MATLAB
  /// Create mxArray from vector
  /// \returns Array created using mxCreateDoubleMatrix.
  ///          Must be deleted by caller using mxDestroyArray.
  mxArray *MxArray() const;

  /// Write vector to MAT-file
  bool WriteMAT(const char *, const char * = "A") const;
#endif

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

#include <irtkVector3D.h>
#include <irtkVector4D.h>

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
inline irtkVector::irtkVector()
:
  _rows  (0),
  _vector(NULL)
{
}

// -----------------------------------------------------------------------------
inline irtkVector::irtkVector(int rows)
:
  _rows  (0),
  _vector(NULL)
{
  Initialize(rows);
}

// -----------------------------------------------------------------------------
inline irtkVector::irtkVector(int rows, double s)
:
  _rows  (0),
  _vector(NULL)
{
  Initialize(rows, s);
}

// -----------------------------------------------------------------------------
inline irtkVector::irtkVector(int rows, double *v)
:
  _rows  (0),
  _vector(NULL)
{
  Initialize(rows, v);
}

// -----------------------------------------------------------------------------
inline irtkVector::irtkVector(const irtkVector& v)
:
  irtkObject(v),
  _rows  (0),
  _vector(NULL)
{
  Initialize(v._rows, v._vector);
}

// -----------------------------------------------------------------------------
template <class T>
inline irtkVector::irtkVector(const irtkVector3D<T> &v)
{
  _rows      = 3;
  _vector    = new double[3];
  _vector[0] = static_cast<double>(v._x);
  _vector[1] = static_cast<double>(v._y);
  _vector[2] = static_cast<double>(v._z);
}

// -----------------------------------------------------------------------------
template <class T>
inline irtkVector::irtkVector(const irtkVector4D<T> &v)
{
  _rows      = 4;
  _vector    = new double[4];
  _vector[0] = static_cast<double>(v._x);
  _vector[1] = static_cast<double>(v._y);
  _vector[2] = static_cast<double>(v._z);
  _vector[3] = static_cast<double>(v._t);
}

// -----------------------------------------------------------------------------
inline irtkVector::~irtkVector()
{
  delete[] _vector;
}

// -----------------------------------------------------------------------------
inline void irtkVector::Initialize(int rows)
{
  if (_rows != rows) {
    delete[] _vector;
    _rows   = rows;
    _vector = (_rows > 0) ? new double[_rows] : NULL;
  }
  memset(_vector, 0, _rows * sizeof(double));
}

// -----------------------------------------------------------------------------
inline void irtkVector::Initialize(int rows, double s)
{
  if (_rows != rows) {
    delete[] _vector;
    _rows   = rows;
    _vector = (_rows > 0) ? new double[_rows] : NULL;
  }
  for (int i = 0; i < _rows; i++) _vector[i] = s;
}

// -----------------------------------------------------------------------------
inline void irtkVector::Initialize(int rows, double *v)
{
  if (_rows != rows) {
    delete[] _vector;
    _rows   = rows;
    _vector = (_rows > 0) ? new double[_rows] : NULL;
  }
  if (v) memcpy(_vector, v, _rows * sizeof(double));
}

// -----------------------------------------------------------------------------
inline void irtkVector::Clear()
{
  delete[] _vector;
  _vector = NULL;
  _rows   = 0;
}

// -----------------------------------------------------------------------------
inline void irtkVector::Resize(int n, double value)
{
  if (n <= 0) {
    Clear();
  } else if (_rows != n) {
    double *vector = new double[n];
    const int m = min(n, _rows);
    for (int i = 0; i < m; ++i) vector[i] = _vector[i];
    for (int i = m; i < n; ++i) vector[i] = value;
    delete[] _vector;
    _vector = vector;
    _rows   = n;
  }
}

// -----------------------------------------------------------------------------
inline irtkVector::operator bool() const
{
  return _rows != 0;
}

// =============================================================================
// Access operators
// =============================================================================

// -----------------------------------------------------------------------------
inline int irtkVector::Rows() const
{
  return _rows;
}

// -----------------------------------------------------------------------------
inline double *irtkVector::RawPointer(int r)
{
  return &_vector[r];
}

// -----------------------------------------------------------------------------
inline const double *irtkVector::RawPointer(int r) const
{
  return &_vector[r];
}

// -----------------------------------------------------------------------------
inline void irtkVector::Put(int rows, double vector)
{
  _vector[rows] = vector;
}

// -----------------------------------------------------------------------------
template <class T>
inline irtkVector &irtkVector::Put(const irtkVector3D<T> &v)
{
  if (_rows != 3) {
    delete[] _vector;
    _rows      = 3;
    _vector    = new double[3];
  }
  _vector[0] = static_cast<double>(v._x);
  _vector[1] = static_cast<double>(v._y);
  _vector[2] = static_cast<double>(v._z);
  return *this;
}

// -----------------------------------------------------------------------------
template <class T>
inline irtkVector &irtkVector::Put(const irtkVector4D<T> &v)
{
  if (_rows != 4) {
    delete[] _vector;
    _rows      = 4;
    _vector    = new double[4];
  }
  _vector[0] = static_cast<double>(v._x);
  _vector[1] = static_cast<double>(v._y);
  _vector[2] = static_cast<double>(v._z);
  _vector[3] = static_cast<double>(v._t);
  return *this;
}

// -----------------------------------------------------------------------------
inline const double &irtkVector::Get(int rows) const
{
  return _vector[rows];
}

// -----------------------------------------------------------------------------
inline double &irtkVector::operator()(int rows)
{
  return _vector[rows];
}

// -----------------------------------------------------------------------------
inline const double &irtkVector::operator()(int rows) const
{
  return _vector[rows];
}

// =============================================================================
// Vector/scalar operations
// =============================================================================

// -----------------------------------------------------------------------------
inline irtkVector &irtkVector::operator =(double x)
{
  for (int i = 0; i < _rows; ++i) _vector[i] = x;
  return *this;
}

// -----------------------------------------------------------------------------
inline irtkVector &irtkVector::operator-=(double x)
{
  for (int i = 0; i < _rows; ++i) _vector[i] -= x;
  return *this;
}

// -----------------------------------------------------------------------------
inline irtkVector &irtkVector::operator+=(double x)
{
  for (int i = 0; i < _rows; ++i) _vector[i] += x;
  return *this;
}

// -----------------------------------------------------------------------------
inline irtkVector &irtkVector::operator*=(double x)
{
  for (int i = 0; i < _rows; ++i) _vector[i] *= x;
  return *this;
}

// -----------------------------------------------------------------------------
inline irtkVector &irtkVector::operator/=(double x)
{
  for (int i = 0; i < _rows; ++i) _vector[i] /= x;
  return *this;
}

// -----------------------------------------------------------------------------
inline irtkVector irtkVector::operator- (double x) const
{
  return (irtkVector(*this) -= x);
}

// -----------------------------------------------------------------------------
inline irtkVector irtkVector::operator+ (double x) const
{
  return (irtkVector(*this) += x);
}

// -----------------------------------------------------------------------------
inline irtkVector irtkVector::operator* (double x) const
{
  return (irtkVector(*this) *= x);
}

// -----------------------------------------------------------------------------
inline irtkVector irtkVector::operator/ (double x) const
{
  return (irtkVector(*this) /= x);
}

// -----------------------------------------------------------------------------
inline irtkVector operator- (double x, const irtkVector &v)
{
  return v - x;
}

// -----------------------------------------------------------------------------
inline irtkVector operator+ (double x, const irtkVector &v)
{
  return v + x;
}

// -----------------------------------------------------------------------------
inline irtkVector operator* (double x, const irtkVector &v)
{
  return v * x;
}

// =============================================================================
// Element-wise vector/vector operations
// =============================================================================

// -----------------------------------------------------------------------------
inline irtkVector irtkVector::operator -() const
{
  irtkVector negative(_rows);
  for (int i = 0; i < _rows; ++i) negative._vector[i] = -_vector[i];
  return negative;
}

// -----------------------------------------------------------------------------
inline irtkVector &irtkVector::operator =(const irtkVector &v)
{
  Initialize(v._rows, v._vector);
  return *this;
}

// -----------------------------------------------------------------------------
inline irtkVector &irtkVector::operator-=(const irtkVector &v)
{
  if (_rows != v._rows) {
    cerr << "irtkVector::operator-=: Size mismatch" << endl;
    exit(1);
  }
  for (int i = 0; i < _rows; ++i) _vector[i] -= v._vector[i];
  return *this;
}

// -----------------------------------------------------------------------------
inline irtkVector &irtkVector::operator+=(const irtkVector &v)
{
  if (_rows != v._rows) {
    cerr << "irtkVector::operator+=: Size mismatch" << endl;
    exit(1);
  }
  for (int i = 0; i < _rows; ++i) _vector[i] += v._vector[i];
  return *this;
}

// -----------------------------------------------------------------------------
inline irtkVector &irtkVector::operator*=(const irtkVector &v)
{
  if (_rows != v._rows) {
    cerr << "irtkVector::operator*=: Size mismatch" << endl;
    exit(1);
  }
  for (int i = 0; i < _rows; ++i) _vector[i] *= v._vector[i];
  return *this;
}

// -----------------------------------------------------------------------------
inline irtkVector &irtkVector::operator/=(const irtkVector &v)
{
  if (_rows != v._rows) {
    cerr << "irtkVector::operator/=: Size mismatch" << endl;
    exit(1);
  }
  for (int i = 0; i < _rows; ++i) _vector[i] /= v._vector[i];
  return *this;
}

// -----------------------------------------------------------------------------
inline irtkVector irtkVector::operator- (const irtkVector &v) const
{
  return (irtkVector(*this) -= v);
}

// -----------------------------------------------------------------------------
inline irtkVector irtkVector::operator+ (const irtkVector &v) const
{
  return (irtkVector(*this) += v);
}

// -----------------------------------------------------------------------------
inline irtkVector irtkVector::operator* (const irtkVector &v) const
{
  return (irtkVector(*this) *= v);
}

// -----------------------------------------------------------------------------
inline irtkVector irtkVector::operator/ (const irtkVector& v) const
{
  return (irtkVector(*this) /= v);
}

// =============================================================================
// Comparison
// =============================================================================

// -----------------------------------------------------------------------------
inline bool irtkVector::operator==(const irtkVector &v) const
{
  if (_rows != v._rows) return false;
  for (int i = 0; i < _rows; ++i) {
    if (_vector[i] != v._vector[i]) return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
#ifndef USE_STL
inline bool irtkVector::operator!=(const irtkVector &v) const
{
  if (_rows != v._rows) return true;
  for (int i = 0; i < _rows; ++i) {
    if (_vector[i] != v._vector[i]) return true;
  }
  return false;
}
#endif

// -----------------------------------------------------------------------------
inline bool irtkVector::operator<(const irtkVector &v) const
{
  if (_rows > v._rows) return false;
  for (int i = 0; i < _rows; ++i) {
    if (_vector[i] >= v._vector[i]) return false;
  }
  return true;
}

// =============================================================================
// Vector products
// =============================================================================

// -----------------------------------------------------------------------------
inline double irtkVector::ScalarProduct(const irtkVector &v) const
{
  if (_rows != v._rows) {
    cerr << "irtkVector::ScalarProduct: Size mismatch" << endl;
    exit(1);
  }
  double s = .0;
  for (int i = 0; i < _rows; ++i) {
    s += _vector[i] * v._vector[i];
  }
  return s;
}

// -----------------------------------------------------------------------------
inline double ScalarProduct(const irtkVector &a, const irtkVector &b)
{
  return a.ScalarProduct(b);
}

// -----------------------------------------------------------------------------
inline double irtkVector::DotProduct(const irtkVector &v) const
{
  return ScalarProduct(v);
}

// -----------------------------------------------------------------------------
inline double DotProduct(const irtkVector &a, const irtkVector &b)
{
  return a.DotProduct(b);
}

// -----------------------------------------------------------------------------
inline irtkVector irtkVector::CrossProduct(const irtkVector &v) const
{
  if (_rows != v._rows) {
    cerr << "irtkVector::CrossProduct: Size mismatch" << endl;
    exit(1);
  }
  int        a, b;
  irtkVector c(_rows, (double*)NULL); // allocate without initialization
  for (int i = 0; i < _rows; ++i) {
    a = (i+1) % _rows;
    b = (i+2) % _rows;
    c._vector[i] = _vector[a] * v._vector[b] - _vector[b] * v._vector[a];
  }
  return c;
}

// -----------------------------------------------------------------------------
inline irtkVector CrossProduct(const irtkVector &a, const irtkVector &b)
{
  return a.CrossProduct(b);
}

// =============================================================================
// Functions
// =============================================================================

// -----------------------------------------------------------------------------
inline double irtkVector::Norm() const
{
  double norm = .0;
  for (int i = 0; i < _rows; i++) norm += _vector[i] * _vector[i];
  return sqrt(norm);
}

// -----------------------------------------------------------------------------
inline irtkVector &irtkVector::Normalize()
{
  double norm = Norm();
  if (norm != .0) (*this) /= norm;
  return *this;
}

// -----------------------------------------------------------------------------
inline irtkVector &irtkVector::Inverse()
{
  for (int i = 0; i < _rows; i++) {
    if (_vector[i] != .0) _vector[i] = 1.0 / _vector[i];
  }
  return *this;
}

// -----------------------------------------------------------------------------
inline void irtkVector::Permute(const vector<int> &idx)
{
  PermuteRows(idx);
}

// =============================================================================
// Conversion (vnl)
// =============================================================================
#ifdef USE_VXL

// -----------------------------------------------------------------------------
template <class T>
inline void irtkVector::Vector2Vnl(vnl_diag_matrix<T>* m) const
{
  unsigned i;

  for (i = 0; i < (unsigned) _rows; i++) {
    (*m)(i) = (T) _vector[i];
  }
}

// -----------------------------------------------------------------------------
template <class T>
inline void irtkVector::Vnl2Vector(vnl_diag_matrix<T>* m)
{
  unsigned i;

  for (i = 0; i < (unsigned) _rows; i++) {
    _vector[i] = (float) (*m)(i);
  }
}


#endif

#endif
