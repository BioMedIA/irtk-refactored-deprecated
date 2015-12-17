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

#ifndef _IRTKMATRIX_H
#define _IRTKMATRIX_H

#include <irtkObject.h>
#include <irtkIndent.h>
#include <irtkMath.h>

#ifdef HAVE_EIGEN
#  include <Eigen/Core>
#  include <Eigen/LU>
#  include <Eigen/SVD>
#  include <Eigen/Eigenvalues>
#endif

#ifdef HAVE_VNL
#  include <vnl/vnl_matrix.h>
#endif

#ifdef HAVE_MATLAB
#  include <irtkMatlab.h>
#endif

// Forward declaration to reduce cyclic header dependencies
class irtkPointSet;
class irtkVector;
class irtkCifstream;
class irtkCofstream;

/**
 * Dense matrix
 *
 * \sa irtkSparseMatrix
 */

class irtkMatrix : public irtkObject
{
  irtkObjectMacro(irtkMatrix);

  // ---------------------------------------------------------------------------
  // Types
public:

  typedef double ElementType;

  // ---------------------------------------------------------------------------
  // Attributes
protected:

  /// Number of rows
  int _rows;

  /// Number of colums
  int _cols;

  /// Matrix values
  double **_matrix;

  /// Whether this instance owns the memory of the matrix elements
  bool _owner;

  // ---------------------------------------------------------------------------
  // Construction/destruction
public:

  /// Default constructor
  irtkMatrix();

  /// Constructor for given number of rows and columns
  irtkMatrix(int, int = -1, double * = NULL);

  /// Copy constructor
  irtkMatrix(const irtkMatrix &);

  /// Destructor
  ~irtkMatrix();

  /// Assign scalar value to all elements
  irtkMatrix& operator =(double);

  /// Assignment operator
  irtkMatrix& operator =(const irtkMatrix &);

  /// Initialize matrix with number of rows and columns
  void Initialize(int, int = -1, double * = NULL);

  /// Resize matrix preserving existing entries
  void Resize(int, int = -1);

  /// Free memory
  void Clear();

  /// Set elements to zero
  void Zero();

  // ---------------------------------------------------------------------------
  // Indexing

  /// Get number of elements
  int NumberOfElements() const;

  /// Get matrix size in each dimension
  pair<int, int> Size() const;

  /// Get number of rows
  int Rows() const;

  /// Get number of columns
  int Cols() const;

  /// Get linear index of element given its row and column indices
  int Index(int, int) const;

  /// Get row index of element given its linear index
  int RowIndex(int) const;

  /// Get column index of element given its linear index
  int ColIndex(int) const;

  /// Get row and column index of element given its linear index
  void SubIndex(int, int &, int &) const;

  /// Get row and column index of element given its linear index
  pair<int, int> SubIndex(int) const;

  // ---------------------------------------------------------------------------
  // Element access

  /// Get pointer to linear memory which stores matrix elements in row-major order
  double *RawPointer(int r = 0, int c = 0);

  /// Get pointer to linear memory which stores matrix elements in row-major order
  const double *RawPointer(int r = 0, int c = 0) const;

  /// Get pointer to linear memory which stores matrix elements in row-major order
  /// \deprecated Use RawPointer instead.
  double *GetPointerToElements(int r = 0, int c = 0);

  /// Get pointer to linear memory which stores matrix elements in row-major order
  /// \deprecated Use RawPointer instead.
  const double *GetPointerToElements(int r = 0, int c = 0) const;

  /// Get reference to element with specified linear index
  double &operator ()(int);

  /// Get const reference to element with specified linear index
  const double &operator ()(int) const;

  /// Get reference to element in specified row and column
  double &operator ()(int, int);

  /// Get const reference to element in specified row and column
  const double &operator ()(int, int) const;

  /// Set value of element with specified linear index
  void Put(int, double);

  /// Get value of element with specified linear index
  double Get(int) const;

  /// Set value of element in specified row and column
  void Put(int, int, double);

  /// Get value of element in specified row and column
  double Get(int, int) const;

  /// Get submatrix
  irtkMatrix operator ()(int, int, int, int) const;

  /// Set submatrix
  void operator ()(irtkMatrix &, int, int);

  // ---------------------------------------------------------------------------
  // Scalar matrix operations

  /// Scale row by a scalar
  irtkMatrix &ScaleRow(int, double);

  /// Scale column by a scalar
  irtkMatrix &ScaleCol(int, double);

  /// Scalar assignment operator
  irtkMatrix &operator =(const double&);

  /// Subtraction of a double
  irtkMatrix &operator -=(const double&);

  /// Addition of a double
  irtkMatrix &operator +=(const double&);

  /// Multiplication with a double
  irtkMatrix &operator *=(const double&);

  /// Division by a double
  irtkMatrix &operator /=(const double&);

  /// Return result of subtraction of a double
  irtkMatrix operator -(const double&) const;

  /// Return result of addition of a double
  irtkMatrix operator +(const double&) const;

  /// Return result of multiplication with a double
  irtkMatrix operator *(const double&) const;

  /// Return result of division by a double
  irtkMatrix operator /(const double&) const;

  // ---------------------------------------------------------------------------
  // Matrix-vector operations

  /// Right-multiply matrix with vector
  irtkVector operator* (const irtkVector&) const;

  // ---------------------------------------------------------------------------
  // Matrix-matrix operations

  /// Matrix subtraction operator
  irtkMatrix& operator -=(const irtkMatrix&);

  /// Matrix addition operator
  irtkMatrix& operator +=(const irtkMatrix&);

  /// Matrix multiplication operator
  irtkMatrix& operator *=(const irtkMatrix&);

  /// Return result of matrix subtraction
  irtkMatrix operator -(const irtkMatrix&) const;

  /// Return result of matrix addition
  irtkMatrix operator +(const irtkMatrix&) const;

  /// Return result of matrix multiplication
  irtkMatrix operator *(const irtkMatrix&) const;

  // ---------------------------------------------------------------------------
  // Comparison operations

  /// Matrix equality
  bool operator ==(const irtkMatrix &) const;

  /// Matrix inequality
  bool operator !=(const irtkMatrix &) const;

  /// Find elements with specific value
  /// \returns Sorted linear indices of found elements
  vector<int> operator ==(double) const;

  /// Find elements with value not equal to specified value
  /// \returns Sorted linear indices of found elements
  vector<int> operator !=(double) const;

  /// Find elements with value less than specified value
  /// \returns Sorted linear indices of found elements
  vector<int> operator <(double) const;

  /// Find elements with value less or equal than specified value
  /// \returns Sorted linear indices of found elements
  vector<int> operator <=(double) const;

  /// Find elements with value greater than specified value
  /// \returns Sorted linear indices of found elements
  vector<int> operator >(double) const;

  /// Find elements with value greater or equal than specified value
  /// \returns Sorted linear indices of found elements
  vector<int> operator >=(double) const;

  // ---------------------------------------------------------------------------
  // Matrix functions

  /// Minimum value in specified row
  double RowMin(int) const;

  /// Maximum value in specified row
  double RowMax(int) const;

  /// Minimum and maximum value in specified row
  void RowRange(int, double &, double &) const;

  /// Minimum value in specified column
  double ColMin(int) const;

  /// Maximum value in specified column
  double ColMax(int) const;

  /// Minimum and maximum value in specified column
  void ColRange(int, double &, double &) const;

  /// Matrix exponential via Pade approximation
  /// (cf. Golub and Van Loan, Matrix Computations, Algorithm 11.3-1)
  irtkMatrix Exp() const;

  /// Matrix logarithm
  irtkMatrix Log() const;

  /// Matrix square root
  irtkMatrix Sqrt() const;

  /// Calculate norm of matrix
  double Norm() const;

  /// Calculate trace of matrix
  double Trace() const;

  // The infinity norm is the maximum of the absolute value row sums.
  double InfinityNorm() const;

  /// Calculate determinant of matrix
  double Det() const;

  /// Calculate determinant of a 3x3 matrix
  double Det3x3() const;

  /// Set to identity matrix
  irtkMatrix &Ident();

  /// Returns true if the matrix is an identity matrix.
  bool IsIdentity() const;

  /// Whether matrix is square
  bool IsSquare() const;

  /// Whether matrix is symmetric
  bool IsSymmetric() const;

  /// Whether matrix is diagonalizable
  bool IsDiagonalizable() const;

  /// Make square matrix symmetric by adding its transpose and divide by 2
  void MakeSymmetric();

  /// Invert matrix
  irtkMatrix &Invert();

  /// Get inverse matrix
  irtkMatrix Inverse() const;

  /// Matrix inversion operator
  irtkMatrix operator !();

  /// Adjugate matrix and return determinant
  irtkMatrix &Adjugate(double &);

  /// Transpose matrix
  irtkMatrix &Transpose();

  /// Matrix transpose operator
  irtkMatrix operator ~();

  /// Permute rows
  irtkMatrix &PermuteRows(vector<int>);

  /// Permute columns
  irtkMatrix &PermuteCols(vector<int>);

  /// Calculate LU decomposition of square matrix
  void LU(irtkMatrix &, irtkMatrix &, double &) const;

  /// Calculate singular value decomposition
  void SVD(irtkMatrix &, irtkVector &, irtkMatrix &) const;

  /// Calculate eigenvalues and eigenvectors of matrix
  ///
  /// This function chooses the appropriate algorithm
  /// according to matrix properties.
  ///
  /// \returns Whether eigenvalue decomposition exists.
  ///          If \c false, squared singular values are returned.
  bool Eigenvalues(irtkMatrix &, irtkVector &, irtkMatrix &) const;

  /// Calculate eigendecomposition of symmetric matrix
  void SymmetricEigen(irtkMatrix &, irtkVector &) const;

  /// Calculate least square fit via SVD
  void LeastSquaresFit(const irtkVector &, irtkVector &) const;

  // ---------------------------------------------------------------------------
  // I/O

  /// Interface to output stream
  friend ostream& operator<< (ostream&, const irtkMatrix&);

  /// Interface to input stream
  friend istream& operator>> (istream&, irtkMatrix&);

  /// Interface to output stream
  friend irtkCofstream& operator<< (irtkCofstream&, const irtkMatrix&);

  /// Interface to input stream
  friend irtkCifstream& operator>> (irtkCifstream&, irtkMatrix&);

  /// Print matrix
  void Print(irtkIndent = 0) const;

  /// Read matrix from file
  void Read(const char *);

  /// Write matrix to file
  void Write(const char *) const;

  /// Import matrix from text file (requires no. of expected rows and cols)
  void Import(const char *, int, int);

#ifdef HAVE_EIGEN

  /// Conversion to Eigen library matrix
  void Matrix2Eigen(Eigen::MatrixXd &) const;

  /// Conversion from Eigen library matrix
  void Eigen2Matrix(Eigen::MatrixXd &);

#endif

#ifdef HAVE_VNL

  /// Conversion to VNL matrix
  template <class T> void Matrix2Vnl(vnl_matrix<T> *m) const;

  /// Conversion from VNL matrix
  template <class T> void Vnl2Matrix(vnl_matrix<T> *m);

#endif

#ifdef HAVE_MATLAB

  /// Create mxArray from dense matrix
  /// \returns Array created using mxCreateDoubleMatrix.
  ///          Must be deleted by caller using mxDestroyArray.
  mxArray *MxArray() const;

  /// Write dense matrix to MAT-file
  bool WriteMAT(const char *, const char * = "A") const;

#endif

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

#include <irtkPoint.h>

// =============================================================================
// Indexing
// =============================================================================

// -----------------------------------------------------------------------------
inline int irtkMatrix::NumberOfElements() const
{
  return _rows * _cols;
}

// -----------------------------------------------------------------------------
inline pair<int, int> irtkMatrix::Size() const
{
  return make_pair(_rows, _cols);
}

// -----------------------------------------------------------------------------
inline int irtkMatrix::Rows() const
{
  return _rows;
}

// -----------------------------------------------------------------------------
inline int irtkMatrix::Cols() const
{
  return _cols;
}

// -----------------------------------------------------------------------------
inline int irtkMatrix::Index(int r, int c) const
{
  return c * _cols + r;
}

// -----------------------------------------------------------------------------
inline int irtkMatrix::RowIndex(int i) const
{
  return i % _cols;
}

// -----------------------------------------------------------------------------
inline int irtkMatrix::ColIndex(int i) const
{
  return i / _cols;
}

// -----------------------------------------------------------------------------
inline void irtkMatrix::SubIndex(int i, int &r, int &c) const
{
  r = RowIndex(i);
  c = ColIndex(i);
}

// -----------------------------------------------------------------------------
inline pair<int, int> irtkMatrix::SubIndex(int i) const
{
  return make_pair(RowIndex(i), ColIndex(i));
}

// =============================================================================
// Element access
// =============================================================================

// -----------------------------------------------------------------------------
inline double *irtkMatrix::RawPointer(int r, int c)
{
  return &_matrix[c][r];
}

// -----------------------------------------------------------------------------
inline const double *irtkMatrix::RawPointer(int r, int c) const
{
  return &_matrix[c][r];
}

// -----------------------------------------------------------------------------
inline double *irtkMatrix::GetPointerToElements(int r, int c)
{
  return &_matrix[c][r];
}

// -----------------------------------------------------------------------------
inline const double *irtkMatrix::GetPointerToElements(int r, int c) const
{
  return &_matrix[c][r];
}

// -----------------------------------------------------------------------------
inline double &irtkMatrix::operator()(int i)
{
#ifndef NO_BOUNDS
  if (0 <= i && i < Rows()) {
#endif
    return _matrix[0][i];
#ifndef NO_BOUNDS
  } else {
    cerr << "irtkMatrix::operator(int): index out of range ";
    exit(1);
  }
#endif
}

// -----------------------------------------------------------------------------
inline const double &irtkMatrix::operator()(int i) const
{
  return const_cast<const double &>(const_cast<irtkMatrix *>(this)->operator()(i));
}

// -----------------------------------------------------------------------------
inline double &irtkMatrix::operator()(int r, int c)
{
#ifndef NO_BOUNDS
  if (0 <= r && r < Rows() && 0 <= c && c < Cols()) {
#endif
    return _matrix[c][r];
#ifndef NO_BOUNDS
  } else {
    cerr << "irtkMatrix::operator(int, int): index out of range ";
    exit(1);
  }
#endif
}

// -----------------------------------------------------------------------------
inline const double &irtkMatrix::operator()(int r, int c) const
{
  return const_cast<const double &>(const_cast<irtkMatrix *>(this)->operator()(r, c));
}

// -----------------------------------------------------------------------------
inline void irtkMatrix::Put(int i, double v)
{
  this->operator()(i) = v;
}

// -----------------------------------------------------------------------------
inline double irtkMatrix::Get(int i) const
{
  return this->operator()(i);
}

// -----------------------------------------------------------------------------
inline void irtkMatrix::Put(int r, int c, double v)
{
  this->operator()(r, c) = v;
}

// -----------------------------------------------------------------------------
inline double irtkMatrix::Get(int r, int c) const
{
  return this->operator()(r, c);
}

// =============================================================================
// Scalar matrix operations
// =============================================================================

// -----------------------------------------------------------------------------
inline irtkMatrix &irtkMatrix::ScaleRow(int r, double s)
{
  for (int c = 0; c < _cols; ++c) _matrix[c][r] *= s;
  return *this;
}

// -----------------------------------------------------------------------------
inline irtkMatrix &irtkMatrix::ScaleCol(int c, double s)
{
  for (int r = 0; r < _rows; ++r) _matrix[c][r] *= s;
  return *this;
}

// -----------------------------------------------------------------------------
inline irtkMatrix &irtkMatrix::operator =(const double &v)
{
  const int n = this->NumberOfElements();
  double   *p = this->RawPointer();
  for (int i = 0; i < n; ++i, ++p) (*p) = v;
  return *this;
}

// -----------------------------------------------------------------------------
inline irtkMatrix &irtkMatrix::operator -=(const double &x)
{
  const int n = this->NumberOfElements();
  double   *p = this->RawPointer();
  for (int i = 0; i < n; ++i, ++p) (*p) -= x;
  return *this;
}

// -----------------------------------------------------------------------------
inline irtkMatrix &irtkMatrix::operator +=(const double &x)
{
  const int n = this->NumberOfElements();
  double   *p = this->RawPointer();
  for (int i = 0; i < n; ++i, ++p) (*p) += x;
  return *this;
}

// -----------------------------------------------------------------------------
inline irtkMatrix &irtkMatrix::operator *=(const double &x)
{
  const int n = this->NumberOfElements();
  double   *p = this->RawPointer();
  for (int i = 0; i < n; ++i, ++p) (*p) *= x;
  return *this;
}

// -----------------------------------------------------------------------------
inline irtkMatrix &irtkMatrix::operator /=(const double &x)
{
  const int n = this->NumberOfElements();
  double   *p = this->RawPointer();
  for (int i = 0; i < n; ++i, ++p) (*p) /= x;
  return *this;
}

// -----------------------------------------------------------------------------
inline irtkMatrix irtkMatrix::operator -(const double &x) const
{
  return irtkMatrix(*this) -= x;
}

// -----------------------------------------------------------------------------
inline irtkMatrix irtkMatrix::operator +(const double &x) const
{
  return irtkMatrix(*this) += x;
}

// -----------------------------------------------------------------------------
inline irtkMatrix irtkMatrix::operator *(const double &x) const
{
  return irtkMatrix(*this) *= x;
}

// -----------------------------------------------------------------------------
inline irtkMatrix irtkMatrix::operator /(const double &x) const
{
  return irtkMatrix(*this) /= x;
}

// =============================================================================
// Comparison operations
// =============================================================================

// -----------------------------------------------------------------------------
inline bool irtkMatrix::operator ==(const irtkMatrix &m) const
{
  if ((m._rows != _rows) || (m._cols != _cols)) return false;
  const int n = this->NumberOfElements();
  const double *ptr1 = m   . RawPointer();
  const double *ptr2 = this->RawPointer();
  for (int i = 0; i < n; ++i, ++ptr1, ++ptr2) {
    if ((*ptr2) != (*ptr1)) return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
inline bool irtkMatrix::operator !=(const irtkMatrix &m) const
{
  return !(*this == m);
}

// -----------------------------------------------------------------------------
inline vector<int> irtkMatrix::operator ==(double x) const
{
  vector<int> idx;
  const int     n = this->NumberOfElements();
  const double *p = this->RawPointer();
  for (int i = 0; i < n; ++i, ++p) {
    if (*p == x) idx.push_back(i);
  }
  return idx;
}

// -----------------------------------------------------------------------------
inline vector<int> irtkMatrix::operator !=(double x) const
{
  vector<int> idx;
  const int     n = this->NumberOfElements();
  const double *p = this->RawPointer();
  for (int i = 0; i < n; ++i, ++p) {
    if (*p != x) idx.push_back(i);
  }
  return idx;
}

// -----------------------------------------------------------------------------
inline vector<int> irtkMatrix::operator <(double x) const
{
  vector<int> idx;
  const int     n = this->NumberOfElements();
  const double *p = this->RawPointer();
  for (int i = 0; i < n; ++i, ++p) {
    if (*p < x) idx.push_back(i);
  }
  return idx;
}

// -----------------------------------------------------------------------------
inline vector<int> irtkMatrix::operator <=(double x) const
{
  vector<int> idx;
  const int     n = this->NumberOfElements();
  const double *p = this->RawPointer();
  for (int i = 0; i < n; ++i, ++p) {
    if (*p <= x) idx.push_back(i);
  }
  return idx;
}

// -----------------------------------------------------------------------------
inline vector<int> irtkMatrix::operator >(double x) const
{
  vector<int> idx;
  const int     n = this->NumberOfElements();
  const double *p = this->RawPointer();
  for (int i = 0; i < n; ++i, ++p) {
    if (*p > x) idx.push_back(i);
  }
  return idx;
}

// -----------------------------------------------------------------------------
inline vector<int> irtkMatrix::operator >=(double x) const
{
  vector<int> idx;
  const int     n = this->NumberOfElements();
  const double *p = this->RawPointer();
  for (int i = 0; i < n; ++i, ++p) {
    if (*p >= x) idx.push_back(i);
  }
  return idx;
}

// =============================================================================
// Matrix functions
// =============================================================================

// -----------------------------------------------------------------------------
inline double irtkMatrix::RowMin(int r) const
{
  if (_rows == 0) return .0;
  double v = _matrix[0][r];
  for (int c = 1; c < _cols; ++c) v = min(v, _matrix[c][r]);
  return v;
}

// -----------------------------------------------------------------------------
inline double irtkMatrix::RowMax(int r) const
{
  if (_rows == 0) return .0;
  double v = _matrix[0][r];
  for (int c = 1; c < _cols; ++c) v = max(v, _matrix[c][r]);
    return v;
}

// -----------------------------------------------------------------------------
inline void irtkMatrix::RowRange(int r, double &v1, double &v2) const
{
  if (_rows > 0) {
    v1 = v2 = _matrix[0][r];
    for (int c = 1; c < _cols; ++c) {
      const double &v = _matrix[c][r];
      v1 = min(v1, v), v2 = max(v2, v);
    }
  } else {
    v1 = v2 = .0;
  }
}

// -----------------------------------------------------------------------------
inline double irtkMatrix::ColMin(int c) const
{
  if (_cols == 0) return .0;
  double v = _matrix[c][0];
  for (int r = 1; r < _rows; ++r) v = min(v, _matrix[c][r]);
  return v;
}

// -----------------------------------------------------------------------------
inline double irtkMatrix::ColMax(int c) const
{
  if (_cols == 0) return .0;
  double v = _matrix[c][0];
  for (int r = 1; r < _rows; ++r) v = min(v, _matrix[c][r]);
  return v;
}

// -----------------------------------------------------------------------------
inline void irtkMatrix::ColRange(int c, double &v1, double &v2) const
{
  if (_cols > 0) {
    v1 = v2 = _matrix[c][0];
    for (int r = 1; r < _rows; ++r) {
      const double &v = _matrix[c][r];
      v1 = min(v1, v), v2 = max(v2, v);
    }
  } else {
    v1 = v2 = .0;
  }
}

// -----------------------------------------------------------------------------
inline double irtkMatrix::Trace() const
{
    int i, j;
    double trace = 0;

    if(_rows == _cols){

        // The trace of a matrix
        for (j = 0; j < _cols; j++) {
            for (i = 0; i < _rows; i++) {
                trace += _matrix[j][i];
            }
        }
        return trace;
    }else{
        cerr << "irtkMatrix::Trace() matrix number of col != row" << endl;
        return 0;
    }
}

// -----------------------------------------------------------------------------
inline double irtkMatrix::Norm() const
{
  int i, j;
  double norm = 0;

  // The norm of a matrix M is defined as trace(M * M~)
  for (j = 0; j < _cols; j++) {
    for (i = 0; i < _rows; i++) {
      norm += _matrix[j][i]*_matrix[j][i];
    }
  }
  return sqrt(norm);
}

// -----------------------------------------------------------------------------
inline double irtkMatrix::InfinityNorm() const
{
  int i, j;
  double normInf = -1.0 * DBL_MAX;
  double sum;

  for (i = 0; i < _rows; ++i) {
    sum = 0;
    for (j = 0; j < _cols; ++j) {
      sum += abs(_matrix[j][i]);
    }
    if (sum > normInf)
      normInf = sum;
  }
  return normInf;
}

// -----------------------------------------------------------------------------
inline double irtkMatrix::Det3x3() const
{
  return _matrix[0][0] * _matrix[1][1] * _matrix[2][2]
       + _matrix[1][0] * _matrix[2][1] * _matrix[0][2]
       + _matrix[2][0] * _matrix[0][1] * _matrix[1][2]
       - _matrix[2][0] * _matrix[1][1] * _matrix[0][2]
       - _matrix[0][0] * _matrix[2][1] * _matrix[1][2]
       - _matrix[1][0] * _matrix[0][1] * _matrix[2][2];
}

#ifdef HAVE_VNL

// -----------------------------------------------------------------------------
template <class T>
void irtkMatrix::Matrix2Vnl(vnl_matrix<T> *m) const
{
  unsigned r, c;

  for (c = 0; c < (unsigned) _cols; c++) {
    for (r = 0; r < (unsigned) _rows; r++) {
      (*m)(r,c) = (T) _matrix[c][r];
    }
  }
}

// -----------------------------------------------------------------------------
template <class T>
void irtkMatrix::Vnl2Matrix(vnl_matrix<T> *m)
{
  unsigned r, c;

  for (c = 0; c < (unsigned) _cols; c++) {
    for (r = 0; r < (unsigned) _rows; r++) {
      _matrix[c][r] = (float) (*m)(r,c);
    }
  }
}

#endif

// -----------------------------------------------------------------------------
inline irtkMatrix irtkMatrix::Inverse() const
{
  return irtkMatrix(*this).Invert();
}

// -----------------------------------------------------------------------------
inline irtkMatrix irtkMatrix::operator !()
{
  return Inverse();
}

// -----------------------------------------------------------------------------
inline irtkMatrix irtkMatrix::operator ~()
{
  return irtkMatrix(*this).Transpose();
}

// =============================================================================
// Backwards compatibility
// =============================================================================

// -----------------------------------------------------------------------------
inline irtkMatrix expm(const irtkMatrix &m)
{
  return m.Exp();
}

// -----------------------------------------------------------------------------
inline irtkMatrix logm(const irtkMatrix &m)
{
  return m.Log();
}

// -----------------------------------------------------------------------------
inline irtkMatrix sqrtm(const irtkMatrix &m)
{
  return m.Sqrt();
}

////////////////////////////////////////////////////////////////////////////////
// Conversion to CUDA vector types
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
// Construct float3x3 matrix from upper left 3x3 sub-matrix of irtkMatrix
IRTKCU_HOST_API inline float3x3 make_float3x3(const irtkMatrix &m)
{
  float3x3 f;
  f.a = make_float3(m(0, 0), m(0, 1), m(0, 2));
  f.b = make_float3(m(1, 0), m(1, 1), m(1, 2));
  f.c = make_float3(m(2, 0), m(2, 1), m(2, 2));
  return f;
}

// -----------------------------------------------------------------------------
// Construct float3x4 matrix from upper left 3x4 sub-matrix of irtkMatrix
IRTKCU_HOST_API inline float3x4 make_float3x4(const irtkMatrix &m)
{
  float3x4 f;
  f.a = make_float4(m(0, 0), m(0, 1), m(0, 2), m(0, 3));
  f.b = make_float4(m(1, 0), m(1, 1), m(1, 2), m(1, 3));
  f.c = make_float4(m(2, 0), m(2, 1), m(2, 2), m(2, 3));
  return f;
}

// -----------------------------------------------------------------------------
// Construct float4x4 matrix from upper left 4x4 sub-matrix of irtkMatrix
IRTKCU_HOST_API inline float4x4 make_float4x4(const irtkMatrix &m)
{
  float4x4 d;
  d.a = make_float4(m(0, 0), m(0, 1), m(0, 2), m(0, 3));
  d.b = make_float4(m(1, 0), m(1, 1), m(1, 2), m(1, 3));
  d.c = make_float4(m(2, 0), m(2, 1), m(2, 2), m(2, 3));
  d.d = make_float4(m(3, 0), m(3, 1), m(3, 2), m(3, 3));
  return d;
}

// -----------------------------------------------------------------------------
// Construct double3x3 matrix from upper left 3x3 sub-matrix of irtkMatrix
IRTKCU_HOST_API inline double3x3 make_double3x3(const irtkMatrix &m)
{
  double3x3 d;
  d.a = make_double3(m(0, 0), m(0, 1), m(0, 2));
  d.b = make_double3(m(1, 0), m(1, 1), m(1, 2));
  d.c = make_double3(m(2, 0), m(2, 1), m(2, 2));
  return d;
}

// -----------------------------------------------------------------------------
// Construct double3x4 matrix from upper left 3x4 sub-matrix of irtkMatrix
IRTKCU_HOST_API inline double3x4 make_double3x4(const irtkMatrix &m)
{
  double3x4 d;
  d.a = make_double4(m(0, 0), m(0, 1), m(0, 2), m(0, 3));
  d.b = make_double4(m(1, 0), m(1, 1), m(1, 2), m(1, 3));
  d.c = make_double4(m(2, 0), m(2, 1), m(2, 2), m(2, 3));
  return d;
}

// -----------------------------------------------------------------------------
// Construct double4x4 matrix from upper left 4x4 sub-matrix of irtkMatrix
IRTKCU_HOST_API inline double4x4 make_double4x4(const irtkMatrix &m)
{
  double4x4 d;
  d.a = make_double4(m(0, 0), m(0, 1), m(0, 2), m(0, 3));
  d.b = make_double4(m(1, 0), m(1, 1), m(1, 2), m(1, 3));
  d.c = make_double4(m(2, 0), m(2, 1), m(2, 2), m(2, 3));
  d.d = make_double4(m(3, 0), m(3, 1), m(3, 2), m(3, 3));
  return d;
}

////////////////////////////////////////////////////////////////////////////////
// Means on Lie groups of transformation matrices
////////////////////////////////////////////////////////////////////////////////

/// Compute log-Euclidean mean of transformation matrices
///
/// \f[
///   mu = \exp\left( \sum w_i \log\left( matrices_i \right) \right)
/// \f]
///
/// \param[in] n        Number of matrices.
/// \param[in] matrices Transformation matrices.
/// \param[in] weights  Weights of transformation matrices.
///                     Uniform weighting if \c NULL.
///
/// \return Mean transformation matrix.
irtkMatrix LogEuclideanMean(int n, const irtkMatrix *matrices, const double *weights = NULL);

/// Compute exponential barycenter / bi-invariant mean of transformation matrices
///
/// This function finds the exponential barycenter of a set of rigid/affine
/// transformation matrices using a fixed point iteration (Gauss-Newton;
/// Barycentric fixed point iteration on Lie groups).
///
/// \f[
///   mu_{t+1} = mu_t \exp\left( \sum w_i \log\left( mu_t^{-1} matrices_i \right) \right)
/// \f]
///
/// \see Arsigny and Pennec, Exponential Barycenters of the Canonical Cartan
///      Connection and Invariant Means on Lie Groups,
///      Matrix Information Geometry (2012)
///
/// \param[in] n        Number of matrices.
/// \param[in] matrices Transformation matrices.
/// \param[in] weights  Weights of transformation matrices.
///                     Uniform weighting if \c NULL.
/// \param[in] niter    Maximum number of fixed point iterations.
/// \param[in] tol      Tolerance of residual infinity norm.
/// \param[in] mu0      Start value of fixed point iteration.
///                     First matrix if \c NULL.
///
/// \return Mean transformation matrix.
irtkMatrix BiInvariantMean(int n, const irtkMatrix *matrices,
                           const double *weights = NULL,
                           int niter = 20, double tol = 1e-12,
                           const irtkMatrix *mu0 = NULL);

// -----------------------------------------------------------------------------
/// \deprecated Use BiInvariantMean instead
inline irtkMatrix FrechetMean(const irtkMatrix *matrices, const double *weights, int n,
                              int niter = 20, double tol = 1e-12,
                              const irtkMatrix *mu0 = NULL)
{
  return BiInvariantMean(n, matrices, weights, niter, tol, mu0);
}

// -----------------------------------------------------------------------------
/// \deprecated Use BiInvariantMean instead
inline irtkMatrix FrechetMean(const irtkMatrix *matrices, int n,
                              int niter = 20, double tol = 1e-12,
                              const irtkMatrix *mu0 = NULL)
{
  return BiInvariantMean(n, matrices, NULL, niter, tol, mu0);
}

////////////////////////////////////////////////////////////////////////////////
// Affine transformation matrices
////////////////////////////////////////////////////////////////////////////////

/// Muliply point by homogeneous transformation matrix
inline void Transform(const irtkMatrix &m, double  x,  double  y,  double  z,
                                           double &mx, double &my, double &mz)
{
  mx = m(0, 0) * x + m(0, 1) * y + m(0, 2) * z + m(0, 3);
  my = m(1, 0) * x + m(1, 1) * y + m(1, 2) * z + m(1, 3);
  mz = m(2, 0) * x + m(2, 1) * y + m(2, 2) * z + m(2, 3);
}

/// Muliply point by homogeneous transformation matrix
inline void Transform(const irtkMatrix &m, double &x, double &y, double &z)
{
  double mx, my, mz;
  Transform(m, x, y, z, mx, my, mz);
  x = mx, y = my, z = mz;
}

/// Muliply point by homogeneous transformation matrix
inline irtkPoint Transform(const irtkMatrix &m, const irtkPoint &p)
{
  irtkPoint p2;
  Transform(m, p._x, p._y, p._z, p2._x, p2._y, p2._z);
  return p2;
}

/// Construct 4x4 homogeneous coordinate transformation matrix from rigid
/// transformation parameters. The output transformation matrix is the
/// composition of a rotation followed by a translation, i.e.,
///
///   T = Translate * Rotate, where Rotate = (Rz Ry Rx)^T
///
/// \note This overloaded function is mainly used by the linear transformation
///       classes such as in particular irtkRigidTransformation as these
///       store the cosine and sine values of the rotation angles as members.
///
/// \param[in]  tx    Translation along x axis.
/// \param[in]  ty    Translation along y axis.
/// \param[in]  tz    Translation along z axis.
/// \param[in]  cosrx Cosine of rotation around x axis in radians.
/// \param[in]  cosry Cosine of rotation around y axis in radians.
/// \param[in]  cosrz Cosine of rotation around z axis in radians.
/// \param[in]  sinrx Sine   of rotation around x axis in radians.
/// \param[in]  sinry Sine   of rotation around y axis in radians.
/// \param[in]  sinrz Sine   of rotation around z axis in radians.
/// \param[out] m     Homogeneous transformation matrix.
void RigidParametersToMatrix(double tx,     double ty,     double tz,
                             double cosrx,  double cosry,  double cosrz,
                             double sinrx,  double sinry,  double sinrz,
                             irtkMatrix &m);

/// Construct 4x4 homogeneous coordinate transformation matrix from rigid
/// transformation parameters. The output transformation matrix is the
/// composition of a rotation followed by a translation, i.e.,
///
///   T = Translate * Rotate
///
/// \param[in]  tx Translation along x axis.
/// \param[in]  ty Translation along y axis.
/// \param[in]  tz Translation along z axis.
/// \param[in]  rx Rotation around x axis in radians.
/// \param[in]  ry Rotation around y axis in radians.
/// \param[in]  rz Rotation around z axis in radians.
/// \param[out] m  Homogeneous transformation matrix.
void RigidParametersToMatrix(double tx,  double ty,  double tz,
                             double rx,  double ry,  double rz, irtkMatrix &m);

/// Construct 4x4 homogeneous coordinate transformation matrix from rigid
/// transformation parameters. The output transformation matrix is the
/// composition of a rotation followed by a translation, i.e.,
///
///   T = Translate * Rotate, where Rotate = (Rz Ry Rx)^T
///
/// \param[in]  tx Translation along x axis.
/// \param[in]  ty Translation along y axis.
/// \param[in]  tz Translation along z axis.
/// \param[in]  rx Rotation around x axis in radians.
/// \param[in]  ry Rotation around y axis in radians.
/// \param[in]  rz Rotation around z axis in radians.
///
/// \returns Homogeneous transformation matrix.
inline irtkMatrix RigidParametersToMatrix(double tx,  double ty,  double tz,
                                          double rx,  double ry,  double rz)
{
  irtkMatrix m;
  RigidParametersToMatrix(tx, ty, tz, rx, ry, rz, m);
  return m;
}

/// Extract Euler angles from 4x4 homogeneous coordinate transformation matrix.
/// The input transformation matrix is assumed to be the composition of a
/// rotation followed by a translation, i.e.,
///
///   T = Translate * Rotate, where Rotate = (Rz Ry Rx)^T
///
/// \param[in]  m  Homogeneous coordinate transformation matrix.
/// \param[out] rx Rotation around x axis in radians.
/// \param[out] ry Rotation around y axis in radians.
/// \param[out] rz Rotation around z axis in radians.
void MatrixToEulerAngles(const irtkMatrix &m,
                         double &rx,  double &ry,  double &rz);

/// Extract rigid transformation parameters from 4x4 homogeneous coordinate
/// transformation matrix. The input transformation matrix is assumed to be
/// the composition of a rotation followed by a translation, i.e.,
///
///   T = Translate * Rotate, where Rotate = (Rz Ry Rx)^T
///
/// \param[in]  m  Homogeneous coordinate transformation matrix.
/// \param[out] tx Translation along x axis.
/// \param[out] ty Translation along y axis.
/// \param[out] tz Translation along z axis.
/// \param[out] rx Rotation around x axis in radians.
/// \param[out] ry Rotation around y axis in radians.
/// \param[out] rz Rotation around z axis in radians.
void MatrixToRigidParameters(const irtkMatrix &m,
                             double &tx,  double &ty,  double &tz,
                             double &rx,  double &ry,  double &rz);

/// Construct 4x4 homogeneous coordinate transformation matrix from affine
/// transformation parameters. The output transformation matrix is the
/// composition of a shearing, followed by a scaling, followed by a
/// rotation, followed by a translation, i.e.,
///
///   T = Translate * Rotate * Scale * Shear, where Rotate = (Rz Ry Rx)^T
///
/// \param[in]  tx  Translation along x axis.
/// \param[in]  ty  Translation along y axis.
/// \param[in]  tz  Translation along z axis.
/// \param[in]  rx  Rotation around x axis in radians.
/// \param[in]  ry  Rotation around y axis in radians.
/// \param[in]  rz  Rotation around z axis in radians.
/// \param[in]  sx  Scaling of x axis (factor, not percentage).
/// \param[in]  sy  Scaling of y axis (factor, not percentage).
/// \param[in]  sz  Scaling of z axis (factor, not percentage).
/// \param[in]  sxy Skew between x and y axes in radians.
/// \param[in]  sxz Skew between x and z axes in radians.
/// \param[in]  syz Skew between y and z axes in radians.
/// \param[out] m  Homogeneous transformation matrix.
void AffineParametersToMatrix(double tx,  double ty,  double tz,
                              double rx,  double ry,  double rz,
                              double sx,  double sy,  double sz,
                              double sxy, double sxz, double syz, irtkMatrix &m);

/// Construct 4x4 homogeneous coordinate transformation matrix from affine
/// transformation parameters. The output transformation matrix is the
/// composition a scaling, followed by a rotation, followed by a translation, i.e.,
///
///   T = Translate * Rotate * Scale, where Rotate = (Rz Ry Rx)^T
///
/// \param[in]  tx  Translation along x axis.
/// \param[in]  ty  Translation along y axis.
/// \param[in]  tz  Translation along z axis.
/// \param[in]  rx  Rotation around x axis in radians.
/// \param[in]  ry  Rotation around y axis in radians.
/// \param[in]  rz  Rotation around z axis in radians.
/// \param[in]  sx  Scaling of x axis (factor, not percentage).
/// \param[in]  sy  Scaling of y axis (factor, not percentage).
/// \param[in]  sz  Scaling of z axis (factor, not percentage).
/// \param[in]  sxy Skew between x and y axes in radians.
/// \param[in]  sxz Skew between x and z axes in radians.
/// \param[in]  syz Skew between y and z axes in radians.
///
/// \returns Homogeneous transformation matrix.
inline irtkMatrix AffineParametersToMatrix(double tx,  double ty,  double tz,
                                           double rx,  double ry,  double rz,
                                           double sx,  double sy,  double sz,
                                           double sxy, double sxz, double syz)
{
  irtkMatrix m(4, 4);
  AffineParametersToMatrix(tx, ty, tz, rx, ry, rz, sx, sy, sz, sxy, sxz, syz, m);
  return m;
}

/// Construct 4x4 homogeneous coordinate transformation matrix from affine
/// transformation parameters. The output transformation matrix is the
/// composition a scaling, followed by a rotation, followed by a translation, i.e.,
///
///   T = Translate * Rotate * Scale, where Rotate = (Rz Ry Rx)^T
///
/// \param[in]  tx  Translation along x axis.
/// \param[in]  ty  Translation along y axis.
/// \param[in]  tz  Translation along z axis.
/// \param[in]  rx  Rotation around x axis in radians.
/// \param[in]  ry  Rotation around y axis in radians.
/// \param[in]  rz  Rotation around z axis in radians.
/// \param[in]  sx  Scaling of x axis (factor, not percentage).
/// \param[in]  sy  Scaling of y axis (factor, not percentage).
/// \param[in]  sz  Scaling of z axis (factor, not percentage).
/// \param[out] m  Homogeneous transformation matrix.
inline void AffineParametersToMatrix(double tx,  double ty,  double tz,
                                     double rx,  double ry,  double rz,
                                     double sx,  double sy,  double sz, irtkMatrix &m)
{
  AffineParametersToMatrix(tx, ty, tz, rx, ry, rz, sx, sy, sz, .0, .0, .0, m);
}

/// Construct 4x4 homogeneous coordinate transformation matrix from affine
/// transformation parameters. The output transformation matrix is the
/// composition a scaling, followed by a rotation, followed by a translation, i.e.,
///
///   T = Translate * Rotate * Scale, where Rotate = (Rz Ry Rx)^T
///
/// \param[in]  tx  Translation along x axis.
/// \param[in]  ty  Translation along y axis.
/// \param[in]  tz  Translation along z axis.
/// \param[in]  rx  Rotation around x axis in radians.
/// \param[in]  ry  Rotation around y axis in radians.
/// \param[in]  rz  Rotation around z axis in radians.
/// \param[in]  sx  Scaling of x axis (factor, not percentage).
/// \param[in]  sy  Scaling of y axis (factor, not percentage).
/// \param[in]  sz  Scaling of z axis (factor, not percentage).
///
/// \returns Homogeneous transformation matrix.
inline irtkMatrix AffineParametersToMatrix(double tx,  double ty,  double tz,
                                           double rx,  double ry,  double rz,
                                           double sx,  double sy,  double sz)
{
  irtkMatrix m(4, 4);
  AffineParametersToMatrix(tx, ty, tz, rx, ry, rz, sx, sy, sz, .0, .0, .0, m);
  return m;
}

/// Extract affine transformation parameters from 4x4 homogeneous coordinate
/// transformation matrix. The input transformation matrix is assumed to be
/// the composition of a shearing, followed by a scaling, followed by a
/// rotation, followed by a translation, i.e.,
///
///   T = Translate * Rotate * Scale * Shear, where Rotate = (Rz Ry Rx)^T
///
///  @sa Graphicx Gems II, with a 12 DoF model without perspective distortion
///      https://github.com/erich666/GraphicsGems/blob/master/gemsii/unmatrix.c
///
/// \param[in]  m   Homogeneous coordinate transformation matrix.
/// \param[out] tx  Translation along x axis.
/// \param[out] ty  Translation along y axis.
/// \param[out] tz  Translation along z axis.
/// \param[out] rx  Rotation around x axis in radians.
/// \param[out] ry  Rotation around y axis in radians.
/// \param[out] rz  Rotation around z axis in radians.
/// \param[out] sx  Scaling of x axis (factor, not percentage).
/// \param[out] sy  Scaling of y axis (factor, not percentage).
/// \param[out] sz  Scaling of z axis (factor, not percentage).
/// \param[out] sxy Skew between x and y axes in radians.
/// \param[out] sxz Skew between x and z axes in radians.
/// \param[out] syz Skew between y and z axes in radians.
void MatrixToAffineParameters(const irtkMatrix &m,
                              double &tx,  double &ty,  double &tz,
                              double &rx,  double &ry,  double &rz,
                              double &sx,  double &sy,  double &sz,
                              double &sxy, double &sxz, double &syz);

/// Find affine transformation matrix which minimizes the mean squared distance
/// between two given sets of corresponding points (e.g., landmarks, fiducial markers)
///
/// \param[in] target (Transformed) Target point set.
/// \param[in] source Fixed source point set.
/// \param[in] weight Weight of corresponding point pair. All correspondences
///                   are weighted equally if vector is empty.
///
/// \returns Homogeneous transformation matrix of the target points.
irtkMatrix ApproximateAffineMatrix(const irtkPointSet &target,
                                   const irtkPointSet &source,
                                   const irtkVector   &weight);

/// Find affine transformation matrix which minimizes the mean squared distance
/// between two given sets of corresponding points (e.g., landmarks, fiducial markers)
///
/// \param[in] target (Transformed) Target point set.
/// \param[in] source Fixed source point set.
///
/// \returns Homogeneous transformation matrix of the target points.
irtkMatrix ApproximateAffineMatrix(const irtkPointSet &target,
                                   const irtkPointSet &source);

/// Orthonormalize upper 3x3 matrix using stabilized Gram-Schmidt
void OrthoNormalize3x3(irtkMatrix &);


#endif
