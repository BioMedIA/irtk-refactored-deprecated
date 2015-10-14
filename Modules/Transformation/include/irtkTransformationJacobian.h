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

#ifndef _IRTKTRANSFORMATIONJACOBIAN_H

#define _IRTKTRANSFORMATIONJACOBIAN_H

#include <map>


/**
 * Sparse matrix for the transformation Jacobian of derivatives w.r.t the parameters.
 *
 * This matrix type only stores the non-zero columns of the Jacobian matrix of
 * the transformation, which contains the derivatives of the transformation
 * w.r.t the transformation parameters. The full Jacobian matrix has dimension
 * 3xN, where N is the number of transformation parameters and the number of
 * rows correspond to the deformation in each spatial dimension (T_x, T_y, T_z).
 */
class irtkTransformationJacobian : public irtkObject
{
  irtkObjectMacro(irtkTransformationJacobian);

public:
  typedef irtkVector3D<double>               ColumnType;
  typedef std::map<int, ColumnType>          SparseMatrixType;
  typedef SparseMatrixType::iterator         ColumnIterator;
  typedef SparseMatrixType::const_iterator   ConstColumnIterator;
  typedef irtkMatrix                         DenseMatrixType;

protected:

  /// Non-zero columns of transformation Jacobian
  SparseMatrixType _Columns;

public:

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Constructor
  irtkTransformationJacobian();

  /// Destructor  
  ~irtkTransformationJacobian();

  /// Remove all non-zero columns
  void Clear();

  // ---------------------------------------------------------------------------
  // Element access

  /// Get number of non-zero columns
  int NumberOfNonZeroColumns() const;

  /// Get iterator to first non-zero column
  ColumnIterator Begin();

  /// Get iterator to first non-zero column
  ConstColumnIterator Begin() const;

  /// Get iterator to position after last non-zero column
  ColumnIterator End();

  /// Get iterator to position after last non-zero column
  ConstColumnIterator End() const;

  /// Get iterator to i-th non-zero column
  ColumnIterator GetNonZeroColumn(int);

  /// Get iterator to i-th non-zero column
  ConstColumnIterator GetNonZeroColumn(int) const;

  /// Get i-th non-zero column vector
  ColumnType &operator [](int);

  /// Get i-th non-zero column vector
  const ColumnType &operator [](int) const;

  /// Get i-th column vector
  ColumnType &operator ()(int);

  /// Get column index of the i-th non-zero column.
  int ColumnIndex(int) const;

  /// Get i-th non-zero column vector.
  ColumnType &ColumnVector(int);

  /// Get i-th non-zero column vector.
  const ColumnType &ColumnVector(int) const;

  /// Get i-th column vector, inserts new zero column if necessary.
  ColumnType &Column(int);

  /// Get iterator to i-th column.
  ColumnIterator Find(int);

  /// Get iterator to i-th column.
  ConstColumnIterator Find(int) const;

  // ---------------------------------------------------------------------------
  // Operators
  
  /// Add transformation Jacobian to this Jacobian matrix
  irtkTransformationJacobian &operator +=(const irtkTransformationJacobian &);

  /// Multiply this transformation Jacobian by a scalar
  irtkTransformationJacobian &operator *=(const double);

  /// Pre-multiply (!) this transformation Jacobian with the given 3x3 matrix
  irtkTransformationJacobian &operator *=(const irtkMatrix &);

  // ---------------------------------------------------------------------------
  // Mathematical operations

  /// Add scaled transformation Jacobian to this Jacobian matrix
  irtkTransformationJacobian &add(const irtkTransformationJacobian &, double);

};


////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
inline irtkTransformationJacobian::irtkTransformationJacobian()
{
}

// -----------------------------------------------------------------------------
inline irtkTransformationJacobian::~irtkTransformationJacobian()
{
}

// =============================================================================
// Element access
// =============================================================================

// -----------------------------------------------------------------------------
inline void irtkTransformationJacobian::Clear()
{
  return _Columns.clear();
}

// -----------------------------------------------------------------------------
inline int irtkTransformationJacobian::NumberOfNonZeroColumns() const
{
  return static_cast<int>(_Columns.size());
}

// -----------------------------------------------------------------------------
inline irtkTransformationJacobian::ColumnIterator irtkTransformationJacobian::Begin()
{
  return _Columns.begin();
}

// -----------------------------------------------------------------------------
inline irtkTransformationJacobian::ConstColumnIterator irtkTransformationJacobian::Begin() const
{
  return _Columns.begin();
}

// -----------------------------------------------------------------------------
inline irtkTransformationJacobian::ColumnIterator irtkTransformationJacobian::End()
{
  return _Columns.end();
}

// -----------------------------------------------------------------------------
inline irtkTransformationJacobian::ConstColumnIterator irtkTransformationJacobian::End() const
{
  return _Columns.end();
}

// -----------------------------------------------------------------------------
inline irtkTransformationJacobian::ColumnIterator irtkTransformationJacobian::GetNonZeroColumn(int i)
{
  ColumnIterator it = _Columns.begin();
  for (int j = 0; j < i; j++) {
    ++it;
    if (it == _Columns.end()) {
      cerr << "irtkTransformationJacobian::GetNonZeroColumn: Index is out of bounds: " << i << endl;
      exit(1);
    }
  }
  return it;
}

// -----------------------------------------------------------------------------
inline irtkTransformationJacobian::ConstColumnIterator irtkTransformationJacobian::GetNonZeroColumn(int i) const
{
  ConstColumnIterator it = _Columns.begin();
  for (int j = 0; j < i; j++) {
    ++it;
    if (it == _Columns.end()) {
      cerr << "irtkTransformationJacobian::GetNonZeroColumn: Index is out of bounds: " << i << endl;
      exit(1);
    }
  }
  return it;
}

// -----------------------------------------------------------------------------
inline int irtkTransformationJacobian::ColumnIndex(int i) const
{
  return GetNonZeroColumn(i)->first;
}

// -----------------------------------------------------------------------------
inline irtkTransformationJacobian::ColumnType &irtkTransformationJacobian::ColumnVector(int i)
{
  return GetNonZeroColumn(i)->second;
}

// -----------------------------------------------------------------------------
inline const irtkTransformationJacobian::ColumnType &irtkTransformationJacobian::ColumnVector(int i) const
{
  return GetNonZeroColumn(i)->second;
}

// -----------------------------------------------------------------------------
inline irtkTransformationJacobian::ColumnType &irtkTransformationJacobian::Column(int c)
{
  return _Columns[c];
}

// -----------------------------------------------------------------------------
inline irtkTransformationJacobian::ColumnType &irtkTransformationJacobian::operator [](int i)
{
  return GetNonZeroColumn(i)->second;
}

// -----------------------------------------------------------------------------
inline const irtkTransformationJacobian::ColumnType &irtkTransformationJacobian::operator [](int i) const
{
  return GetNonZeroColumn(i)->second;
}

// -----------------------------------------------------------------------------
inline irtkTransformationJacobian::ColumnType &irtkTransformationJacobian::operator ()(int i)
{
  return Column(i);
}

// -----------------------------------------------------------------------------
inline irtkTransformationJacobian::ColumnIterator irtkTransformationJacobian::Find(int c)
{
  return _Columns.find(c);
}

// -----------------------------------------------------------------------------
inline irtkTransformationJacobian::ConstColumnIterator irtkTransformationJacobian::Find(int c) const
{
  return _Columns.find(c);
}

// =============================================================================
// Unary operators
// =============================================================================

// -----------------------------------------------------------------------------
inline irtkTransformationJacobian &irtkTransformationJacobian::operator +=(const irtkTransformationJacobian &b)
{
  for (ConstColumnIterator it = b.Begin(); it != b.End(); ++it) {
    _Columns[it->first] += it->second;
  }
  return *this;
}

// -----------------------------------------------------------------------------
inline irtkTransformationJacobian &irtkTransformationJacobian::operator *=(const double s)
{
  for (ColumnIterator it = Begin(); it != End(); ++it) it->second *= s;
  return *this;
}

// -----------------------------------------------------------------------------
inline irtkTransformationJacobian &irtkTransformationJacobian::operator *=(const irtkMatrix &a)
{
  // Note: Read this operator as: a * (*this)!
  ColumnType v;
  for (ColumnIterator it = Begin(); it != End(); ++it) {
    v._x = a(0, 0) * it->second._x + a(0, 1) * it->second._y + a(0, 2) * it->second._z;
    v._y = a(1, 0) * it->second._x + a(1, 1) * it->second._y + a(1, 2) * it->second._z;
    v._z = a(2, 0) * it->second._x + a(2, 1) * it->second._y + a(2, 2) * it->second._z;
    it->second = v;
  }
  return *this;
}

// =============================================================================
// Binary operators
// =============================================================================

// -----------------------------------------------------------------------------
/// Calculate column-by-column sum of transformation Jacobian
inline irtkTransformationJacobian operator +(irtkTransformationJacobian &a, irtkTransformationJacobian &b)
{
  irtkTransformationJacobian c = a;
  c += b;
  return c;
}

// -----------------------------------------------------------------------------
/// Multiply transformation Jacobian and scalar
inline irtkTransformationJacobian operator *(irtkTransformationJacobian &a, double s)
{
  irtkTransformationJacobian b = a;
  b *= s;
  return b;
}

// -----------------------------------------------------------------------------
/// Multiply transformation Jacobian and scalar
inline irtkTransformationJacobian operator *(double s, irtkTransformationJacobian &a)
{
  return a * s;
}

// -----------------------------------------------------------------------------
/// Calculate product of 3x3 matrix and transformation Jacobian
inline irtkTransformationJacobian operator *(irtkMatrix &a, irtkTransformationJacobian &b)
{
  irtkTransformationJacobian c = b;
  c *= a;
  return c;
}

// =============================================================================
// Mathematical operations
// =============================================================================

// -----------------------------------------------------------------------------
inline irtkTransformationJacobian &irtkTransformationJacobian::add(const irtkTransformationJacobian &b, double s)
{
  for (ConstColumnIterator it = b.Begin(); it != b.End(); ++it) {
    _Columns[it->first] += it->second * s;
  }
  return *this;
}

// =============================================================================
// Debugging
// =============================================================================

inline bool has_nan(const irtkTransformationJacobian &a)
{
  for (irtkTransformationJacobian::ConstColumnIterator it = a.Begin(); it != a.End(); ++it) {
    if (IsNaN(it->second._x) || IsNaN(it->second._y) || IsNaN(it->second._z)) {
      cerr << "irtkTransformationJacobian::has_nan: Found NaN in column " << it->first << endl;
      return true;
    }
  }
  return false;
}


#endif
