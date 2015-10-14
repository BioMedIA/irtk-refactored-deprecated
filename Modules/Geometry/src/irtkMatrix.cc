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

#include <irtkGeometry.h>

#if defined(HAVE_EIGEN)
#  include <Eigen/Dense>
#elif defined(HAVE_VNL)
#  include <vnl/algo/vnl_qr.h>
#endif

#define EIGEN_TOL 1.0e-5
//#define JACOBI


// =============================================================================
// Auxiliary macros
// =============================================================================

// -----------------------------------------------------------------------------
#define IRTK_CHECK_MATRIX_SIZE(op, m) \
  while (_rows != m._rows || _cols != m._cols) { \
    cerr << "irtkMatrix::" op ": Right-hand matrix must have same size" << endl; \
    exit(1); \
  }

// -----------------------------------------------------------------------------
#define IRTK_CHECK_MATRIX_IS_SQUARE(op) \
   while (_rows != _cols) { \
    cerr << "irtkMatrix::" op ": Matrix must be square" << endl; \
    exit(1); \
  }

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkMatrix::irtkMatrix()
:
  _rows  (0),
  _cols  (0),
  _matrix(NULL),
  _owner (false)
{
}

// -----------------------------------------------------------------------------
irtkMatrix::irtkMatrix(int rows, int cols, double *data)
:
  _rows  (rows),
  _cols  (cols < 0 ? rows : cols),
  _matrix(NULL),
  _owner (false)
{
  if (_rows * _cols > 0) {
    Allocate(_matrix, _rows, _cols, data);
    _owner = (this->RawPointer() != data);
    if (_owner) {
      memset(this->RawPointer(), 0, _rows * _cols * sizeof(double));
    }
  }
}

// -----------------------------------------------------------------------------
irtkMatrix::irtkMatrix(const irtkMatrix& m)
:
  irtkObject(m),
  _rows  (m._rows),
  _cols  (m._cols),
  _matrix(NULL),
  _owner (false)
{
  if (_rows * _cols > 0) {
    Allocate(_matrix, _rows, _cols);
    _owner = true;
    memcpy(this->RawPointer(), m.RawPointer(), _rows * _cols * sizeof(double));
  }
}

// -----------------------------------------------------------------------------
irtkMatrix::~irtkMatrix()
{
  Clear();
}

// -----------------------------------------------------------------------------
irtkMatrix& irtkMatrix::operator =(double s)
{
  for (int c = 0; c < _cols; ++c)
  for (int r = 0; r < _rows; ++r) {
    _matrix[c][r] = s;
  }
  return *this;
}

// -----------------------------------------------------------------------------
irtkMatrix& irtkMatrix::operator =(const irtkMatrix& m)
{
  if (_rows != m._rows || _cols != m._cols) {
    _rows = m._rows;
    _cols = m._cols;
    if (_owner) Deallocate(_matrix);
    Allocate(_matrix, _rows, _cols);
    _owner = true;
  }
  if (_matrix && this->RawPointer() != m.RawPointer()) {
    memcpy(this->RawPointer(), m.RawPointer(), _rows * _cols * sizeof(double));
  }
  return *this;
}

// -----------------------------------------------------------------------------
void irtkMatrix::Initialize(int rows, int cols, double *data)
{
  if (cols < 0) cols = rows;
  if (_rows != rows || _cols != cols || (data && this->RawPointer() != data)) {
    if (_owner) Deallocate(_matrix);
    _rows = rows;
    _cols = cols;
    Allocate(_matrix, _rows, _cols, data);
    _owner = (this->RawPointer() != data);
  }
  if (_owner) {
    memset(this->RawPointer(), 0, _rows * _cols * sizeof(double));
  }
}

// -----------------------------------------------------------------------------
void irtkMatrix::Resize(int rows, int cols)
{
  if (rows <= 0 || cols <= 0) {
    Clear();
  } else if (_rows != rows || _cols != cols) {
    double **matrix = Allocate<double>(rows, cols);
    const int m = min(rows, _rows);
    const int n = min(cols, _cols);
    for (int c = 0; c < n; ++c) {
      for (int r = 0; r < m; ++r) {
        matrix[c][r] = _matrix[c][r];
      }
      for (int r = m; r < rows; ++r) {
        matrix[c][r] = .0;
      }
    }
    for (int c = n; c < cols; ++c) {
      for (int r = 0; r < rows; ++r) {
        matrix[c][r] = .0;
      }
    }
    if (_owner) Deallocate(_matrix);
    _matrix = matrix;
    _owner  = true;
    _rows   = rows;
    _cols   = cols;
  }
}

// -----------------------------------------------------------------------------
void irtkMatrix::Clear()
{
  if (_owner) Deallocate(_matrix);
  _rows = 0;
  _cols = 0;
}

// -----------------------------------------------------------------------------
void irtkMatrix::Zero()
{
  memset(this->RawPointer(), 0, _rows * _cols * sizeof(double));
}

// =============================================================================
// Element access
// =============================================================================

// -----------------------------------------------------------------------------
void irtkMatrix::operator()(irtkMatrix &m, int rows, int cols)
{
  int i, j;

  if ((m._rows + rows > _rows) || (m._cols + cols > _cols)) {
    cerr << "irtkMatrix::operator(): Invalid range" << endl;
    exit(1);
  }
  for (j = cols; j < m._cols + cols; j++) {
    for (i = rows; i < m._rows + rows; i++) {
      _matrix[j][i] = m._matrix[j - cols][i - rows];
    }
  }
}

// -----------------------------------------------------------------------------
irtkMatrix irtkMatrix::operator()(int rows1, int cols1, int rows2, int cols2) const
{
  int i, j;

  if ((rows1 < 0) || (rows2 > _rows) || (cols1 < 0) || (cols2 > _cols)) {
    cerr << "irtkMatrix::operator(): Invalid range" << endl;
    exit(1);
  }

  // Create new matrix
  irtkMatrix m(rows2-rows1, cols2-cols1);

  for (j = cols1; j < cols2; j++) {
    for (i = rows1; i < rows2; i++) {
      m._matrix[j-cols1][i-rows1] = _matrix[j][i];
    }
  }

  return m;
}

// =============================================================================
// Matrix-vector operations
// =============================================================================

// -----------------------------------------------------------------------------
irtkVector irtkMatrix::operator *(const irtkVector &v) const
{
  if (_cols != v.Rows()) {
    cerr << "irtkMatrix::operator*(irtkVector): Size mismatch" << endl;
    exit(1);
  }

  irtkVector p(_rows);
  double sum;
  for (int r = 0; r < _rows; ++r) {
    sum = .0;
    for (int c = 0; c < _cols; ++c) {
      sum += _matrix[c][r] * v(c);
    }
    p(r) = sum;
  }

  return p;
}

// =============================================================================
// Matrix-matrix operations
// =============================================================================

// -----------------------------------------------------------------------------
irtkMatrix &irtkMatrix::operator -=(const irtkMatrix &m)
{
  IRTK_CHECK_MATRIX_SIZE("operator -=", m);
  const int n = this->NumberOfElements();
  const double *ptr1 = m   . RawPointer();
  double       *ptr2 = this->RawPointer();
  for (int i = 0; i < n; ++i, ++ptr1, ++ptr2) (*ptr2) -= (*ptr1);
  return *this;
}

// -----------------------------------------------------------------------------
irtkMatrix &irtkMatrix::operator +=(const irtkMatrix &m)
{
  IRTK_CHECK_MATRIX_SIZE("operator +=", m);
  const int n = this->NumberOfElements();
  const double *ptr1 = m   . RawPointer();
  double       *ptr2 = this->RawPointer();
  for (int i = 0; i < n; ++i, ++ptr1, ++ptr2) (*ptr2) += (*ptr1);
  return *this;
}

// -----------------------------------------------------------------------------
irtkMatrix& irtkMatrix::operator *=(const irtkMatrix& m)
{
  (*this) = (*this) * m;
  return *this;
}

// -----------------------------------------------------------------------------
irtkMatrix irtkMatrix::operator -(const irtkMatrix& m) const
{
  IRTK_CHECK_MATRIX_SIZE("operator -", m);
  irtkMatrix m2(*this);
  m2 -= m;
  return m2;
}

// -----------------------------------------------------------------------------
irtkMatrix irtkMatrix::operator +(const irtkMatrix& m) const
{
  IRTK_CHECK_MATRIX_SIZE("operator +", m);
  irtkMatrix m2(*this);
  m2 += m;
  return m2;
}

// -----------------------------------------------------------------------------
irtkMatrix irtkMatrix::operator *(const irtkMatrix& m) const
{
  if (_cols != m.Rows()) {
    cerr << "irtkMatrix::operator *: Matrix size mismatch" << endl;
    exit(1);
  }
  irtkMatrix m2(_rows, m._cols);
  for (int i = 0; i < _rows; ++i) {
    for (int j = 0; j < m._cols; ++j) {
      m2._matrix[j][i] = 0.0;
      for (int k = 0; k < _cols; ++k) {
        m2._matrix[j][i] += _matrix[k][i] * m._matrix[j][k];
      }
    }
  }
  return m2;
}

// =============================================================================
// Matrix functions
// =============================================================================

// -----------------------------------------------------------------------------
irtkMatrix irtkMatrix::Exp() const
{
  // Matrix exponential via Pade approximation.
  // See Golub and Van Loan, Matrix Computations, Algorithm 11.3-1.

  // Number of iterations (determines accuracy).
  const int q = 6;

  // First this matrix is scaled by a power of 2 so its norm is < 1/2.
  // j is the index of the power required.
  double norm = InfinityNorm();
  int e = (int) ceil(log(norm) / log(2.0));
  int j = max(0, 1 + e);

  irtkMatrix A = (*this) / (pow(2.0, j));
  irtkMatrix D(_rows, _cols);
  irtkMatrix N(_rows, _cols);
  irtkMatrix X(_rows, _cols);

  D.Ident();
  N.Ident();
  X.Ident();

  double c             = 1.0;
  int    minusOnePower = 1;
  for (int k = 1; k <= q; ++k) {
    c = c * (q - k + 1) / ((double) k * (2*q - k + 1));
    X = A * X;
    N = N +  X * c;
    minusOnePower *= -1;
    D = D + X * minusOnePower * c;
  }

  D.Invert();
  X = D * N;

  // Squaring steps
  for (int k = 1; k <= j; ++k) X = X * X;
  return X;
}

// -----------------------------------------------------------------------------
irtkMatrix irtkMatrix::Log() const
{
  IRTK_CHECK_MATRIX_IS_SQUARE("Log");

  irtkVector eigval;
  irtkMatrix eigvec, tmp;
  Eigenvalues(eigvec, eigval, tmp);
  for (int i = 0; i < eigval.Rows(); ++i) {
    if (eigval(i) <= 0) {
      cerr << "irtkMatrix::Log: Found non-positive eigenvalue: e(" << (i+1) << ") = " << eigval(i) << endl;
      this->Print(2);
      exit(1);
    }
  }

  const int    maxit = 100;
  const double tol   = 0.00001;
  int          i, k, n;

  irtkMatrix A(*this);
  irtkMatrix I(_rows, _cols);
  irtkMatrix Z(_rows, _cols);
  irtkMatrix X(_rows, _cols);
  irtkMatrix D(_rows, _cols);

  I.Ident();

  D = A - I;
  k = 0, n = 0;
  while (D.InfinityNorm() > 0.5 && n < maxit) {
    A = A.Sqrt(), ++k;
    D = A - I,    ++n;
  }

  A = X = Z = I - A;
  i = 1, n = 0;
  while (Z.InfinityNorm() > tol && n < maxit) {
    Z = Z * A,       ++i;
    X = X + (Z / i), ++n;
  }

  X = X * -pow(2.0, k);
  return X;
}

// -----------------------------------------------------------------------------
irtkMatrix irtkMatrix::Sqrt() const
{
  const int    maxit = 100;
  const double tol   = 0.0001;

  irtkMatrix X   (*this);
  irtkMatrix Y   (_rows, _cols);
  irtkMatrix D   (_rows, _cols);
  irtkMatrix invX(_rows, _cols);
  irtkMatrix invY(_rows, _cols);

  Y.Ident();
  D = X * X - (*this);

  for (int i = 0; i < maxit && D.InfinityNorm() > tol; ++i) {
    invX = X.Inverse();
    invY = Y.Inverse();
    X = (X + invY) * 0.5;
    Y = (Y + invX) * 0.5;
    D = X * X - (*this);
  }

  return X;
}

// -----------------------------------------------------------------------------
double irtkMatrix::Det() const
{
  IRTK_CHECK_MATRIX_IS_SQUARE("Det");

  double d, sign;
  int i;
  irtkMatrix lu, perm;

  this->LU(lu, perm, sign);

  d = sign;
  for (i = 0; i < _rows; i++) {
    d *= lu(i, i);
  }

  return d;
}

// -----------------------------------------------------------------------------
irtkMatrix &irtkMatrix::Invert()
{
  IRTK_CHECK_MATRIX_IS_SQUARE("Invert");

  double sign;
  int i, j, p, q;
  irtkMatrix lu, perm, A, inv;
  irtkVector b, x, y;

  if (this->Det() == 0) {
    cerr << "irtkMatrix::Invert: Zero determinant\n";
    exit(1);
  }

  this->LU(lu, perm, sign);

  A.Initialize(_rows, _rows);
  b.Initialize(_rows);
  x.Initialize(_rows);
  y.Initialize(_rows);

  for (j = 0; j < _rows; j++) {
    for (i = 0; i < _rows; i++) {
      b(i) = 0.0;
    }
    b(j) = 1.0;

    // Forward substitution
    for (p = 0; p < _rows; p++)
    {
      y(p) = b(p);
      for (q = 0; q < p; q++)
      {
        y(p) -= lu(p, q) * y(q);
      }
    }

    // Back substitution
    for (p = _rows - 1; p >= 0; p--)
    {
      x(p) = y(p);
      for (q = p + 1; q < _rows; q++)
      {
        x(p) -= lu(p, q) * x(q);
      }
      x(p) /= lu(p, p);
    }

    // Write column j
    for (i = 0; i < _rows; i++) {
      A(i, j) = x(i);
    }
  }

  // Multiply with permutation matrix to get actual inverse
  inv = A * perm;

  // Copy into *this and return
  return (*this = inv);
}

// -----------------------------------------------------------------------------
irtkMatrix &irtkMatrix::Adjugate(double &d)
{
  IRTK_CHECK_MATRIX_IS_SQUARE("Adjugate");

  double sign;
  int i, j, p, q;
  irtkMatrix lu, perm, A, inv;
  irtkVector b, x, y;

  d = this->Det();

  if (d == 0) {
    cerr << "irtkMatrix::Invert: Zero determinant\n";
    exit(1);
  }

  this->LU(lu, perm, sign);

  A.Initialize(_rows, _rows);
  b.Initialize(_rows);
  x.Initialize(_rows);
  y.Initialize(_rows);

  for (j = 0; j < _rows; j++) {
    for (i = 0; i < _rows; i++) {
      b(i) = 0.0;
    }
    b(j) = 1.0;

    // Forward substitution
    for (p = 0; p < _rows; p++)
    {
      y(p) = b(p);
      for (q = 0; q < p; q++)
      {
        y(p) -= lu(p, q) * y(q);
      }
    }

    // Back substitution
    for (p = _rows - 1; p >= 0; p--)
    {
      x(p) = y(p);
      for (q = p + 1; q < _rows; q++)
      {
        x(p) -= lu(p, q) * x(q);
      }
      x(p) /= lu(p, p);
    }

    // Write column j
    for (i = 0; i < _rows; i++) {
      A(i, j) = x(i) * d;
    }
  }

  // Multiply with permutation matrix to get actual inverse
  inv = A * perm;

  // Copy into *this and return
  return (*this = inv);
}

// -----------------------------------------------------------------------------
irtkMatrix &irtkMatrix::Transpose()
{
  int i, j;

  if (_rows == _cols) {
    double tmp;
    for (i = 0; i < _rows; i++) {
      for (j = 0; j < i; j++) {
        tmp=_matrix[i][j];
        _matrix[i][j] = _matrix[j][i];
        _matrix[j][i] = tmp;
      }
    }
  } else {
    irtkMatrix tmp(_rows, _cols);
    for (j = 0; j < _cols; j++) {
      for (i = 0; i < _rows; i++) {
        tmp._matrix[j][i] = _matrix[j][i];
      }
    }
    Deallocate(_matrix);
    Allocate(_matrix, tmp._cols, tmp._rows);
    _rows = tmp._cols;
    _cols = tmp._rows;
    for (j = 0; j < _cols; j++) {
      for (i = 0; i < _rows; i++) {
        _matrix[j][i] = tmp._matrix[i][j];
      }
    }
  }

  return *this;
}

// -----------------------------------------------------------------------------
irtkMatrix &irtkMatrix::PermuteRows(std::vector<int> idx)
{
  irtkAssert(idx.size() <= static_cast<size_t>(_cols), "valid permutation");
  for (int r1 = 0; r1 < static_cast<int>(idx.size()); ++r1) {
    int r2 = idx[r1];
    if (r2 == r1) continue;
    for (int c = 0; c < _cols; ++c) swap(_matrix[c][r1], _matrix[c][r2]);
    for (int r = r1 + 1; r < _rows; ++r) {
      if (idx[r] == r1) {
        swap(idx[r], idx[r1]);
        break;
      }
    }
  }
  return *this;
}

// -----------------------------------------------------------------------------
irtkMatrix &irtkMatrix::PermuteCols(std::vector<int> idx)
{
  irtkAssert(idx.size() <= static_cast<size_t>(_cols), "valid permutation");
  for (int c1 = 0; c1 < static_cast<int>(idx.size()); ++c1) {
    int c2 = idx[c1];
    if (c2 == c1) continue;
    for (int r = 0; r < _rows; ++r) swap(_matrix[c1][r], _matrix[c2][r]);
    for (int c = c1 + 1; c < _cols; ++c) {
      if (idx[c] == c1) {
        swap(idx[c], idx[c1]);
        break;
      }
    }
  }
  return *this;
}

// -----------------------------------------------------------------------------
void irtkMatrix::LU(irtkMatrix &lu, irtkMatrix &perm, double &sign) const
{
  IRTK_CHECK_MATRIX_IS_SQUARE("LU");

#ifdef HAVE_EIGEN
  Eigen::MatrixXd m(_rows, _rows);
  Eigen::PartialPivLU<Eigen::MatrixXd> ppiv_lu;

  Matrix2Eigen(m);

  ppiv_lu.compute(m);

  Eigen::MatrixXd eigen_lu    = ppiv_lu.matrixLU();
  Eigen::MatrixXd eigen_perm  = ppiv_lu.permutationP();

  lu = irtkMatrix(eigen_lu.rows(), eigen_lu.cols());
  lu.Eigen2Matrix(eigen_lu);
  perm = irtkMatrix(eigen_perm.rows(), eigen_perm.cols());
  perm.Eigen2Matrix(eigen_perm);

  // Trick to determine the sign (either 1 or -1) of the LU decomposition
  // since the PartialPivLU class interface does not provide a way to get it
  sign = ppiv_lu.determinant() / eigen_lu.diagonal().prod();
#else
  cerr << "irtkMatrix::LU: IRTK must be compiled with Eigen library to have this functionality" << endl;
#endif
}

// -----------------------------------------------------------------------------
void irtkMatrix::SVD(irtkMatrix &u, irtkVector &w, irtkMatrix &v) const
{
#ifdef HAVE_EIGEN
  Eigen::MatrixXd m(_rows, _cols);
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(_rows, _cols, Eigen::ComputeThinU | Eigen::ComputeThinV);

  Matrix2Eigen(m);

  svd.compute(m);

  Eigen::MatrixXd eigen_u = svd.matrixU();
  Eigen::MatrixXd eigen_v = svd.matrixV();
  Eigen::VectorXd eigen_w = svd.singularValues();

  u = irtkMatrix(eigen_u.rows(), eigen_u.cols());
  u.Eigen2Matrix(eigen_u);
  v = irtkMatrix(eigen_v.rows(), eigen_v.cols());
  v.Eigen2Matrix(eigen_v);
  w = irtkVector(eigen_w.size());
  w.Eigen2Vector(eigen_w);
#else
  cerr << "irtkMatrix::SVD: IRTK must be compiled with Eigen library to have this functionality" << endl;
#endif
}

// -----------------------------------------------------------------------------
bool irtkMatrix::Eigenvalues(irtkMatrix &E1, irtkVector &e, irtkMatrix &E2) const
{
  bool ok = true;

  if (this->IsSymmetric()) {

    this->SymmetricEigen(E1, e);
    E2 = E1.Inverse();

  } else {

    // TODO: Consider use of Eigen::EigenSolver instead.
    //       http://eigen.tuxfamily.org/dox/classEigen_1_1EigenSolver.html

    this->SVD(E1, e, E2);
    for (int i = 0; i < e.Rows(); ++i) e(i) *= e(i);

    // Singular value decomposition differs from eigenvalue decomposition
    // if matrix is not square and diagonalizable
    ok = this->IsDiagonalizable();

  }

  return ok;
}

// -----------------------------------------------------------------------------
void irtkMatrix::SymmetricEigen(irtkMatrix &E, irtkVector &e) const
{
  IRTK_CHECK_MATRIX_IS_SQUARE("SymmetricEigen");

#ifdef HAVE_EIGEN
  Eigen::MatrixXd m(_rows, _rows);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_solver;

  Matrix2Eigen(m);

  eig_solver.compute(m);

  Eigen::MatrixXd eigen_evec = eig_solver.eigenvectors();
  Eigen::VectorXd eigen_eval = eig_solver.eigenvalues();

  E = irtkMatrix(eigen_evec.rows(), eigen_evec.cols());
  E.Eigen2Matrix(eigen_evec);
  e = irtkVector(eigen_eval.size());
  e.Eigen2Vector(eigen_eval);
#else
  cerr << "irtkMatrix::SymmetricEigen: IRTK must be compiled with Eigen library to have this functionality" << endl;
#endif
}

// -----------------------------------------------------------------------------
void irtkMatrix::LeastSquaresFit(const irtkVector &y, irtkVector &x) const
{
  int i, j;
  float wmax, thresh;
  irtkMatrix u, v, m_inv;
  irtkVector w;

  // nmatrix should be rows
  // ma should be cols

  if (y.Rows() != _rows) {
    cerr << "irtkMatrix::LeastSquaresFit: Y has wrong dimensions" << endl;
    exit(1);
  }

  if (x.Rows() != _cols) {
    cerr << "irtkMatrix::LeastSquaresFit: X has wrong dimensions" << endl;
    exit(1);
  }

  // Calculate least squares fit via SVD
  wmax = 0.0;

  this->SVD(u, w, v);

  for (j = 0; j < _cols; j++) if (w(j) > wmax) wmax = w(j);
  thresh = EIGEN_TOL * wmax;
  for (j = 0; j < _cols; j++) if (w(j) < thresh) w(j) = 0.0;

  // Compute u_new = {w}^-1 * u^T
  u.Transpose();
  for (i = 0; i < _cols; i++) {
    for (j = 0; j < _rows; j++) {
      if (w(i) != 0.0) {
        u(i, j) /= w(i);
      } else {
        u(i, j) = 0.0;
      }
    }
  }

  // Compute x = V * u_new * y
  x = (v * u) * y;
}

// -----------------------------------------------------------------------------
irtkMatrix &irtkMatrix::Ident()
{
  for (int c = 0; c < _cols; ++c) {
    for (int r = 0; r < _rows; ++r) {
      _matrix[c][r] = static_cast<double>(r == c);
    }
  }
  return *this;
}

// -----------------------------------------------------------------------------
bool irtkMatrix::IsIdentity() const
{
  if (_rows == 0 || _cols == 0 || _rows != _cols) return false;

  for (int c = 0; c < _cols; ++c) {
    for (int r = 0; r < _rows; ++r) {
      if (_matrix[c][r] != static_cast<double>(r == c)) return false;
    }
  }

  return true;
}

// -----------------------------------------------------------------------------
bool irtkMatrix::IsSquare() const
{
  return _rows == _cols;
}

// -----------------------------------------------------------------------------
bool irtkMatrix::IsSymmetric() const
{
  if (_rows == 0 || _cols == 0 || _rows != _cols) return false;
  for (int c = 0; c < _cols; ++c)
  for (int r = 0; r < _rows; ++r) {
    if (!fequal(_matrix[c][r], _matrix[r][c], 1e-6)) return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
bool irtkMatrix::IsDiagonalizable() const
{
  if (_rows == 0 || _cols == 0 || _rows != _cols) return false;
  irtkMatrix AT  = (*this); AT.Transpose();
  irtkMatrix ATA = AT * (*this);
  irtkMatrix AAT = (*this) * AT;
  for (int c = 0; c < _cols; ++c)
  for (int r = 0; r < _rows; ++r) {
    if (!fequal(AAT(r, c), ATA(r, c), 1e-6)) return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
struct MakeMatrixSymmetric
{
  irtkMatrix &_Matrix;

  MakeMatrixSymmetric(irtkMatrix &m) : _Matrix(m) {}

  void operator ()(const blocked_range2d<int> &re) const
  {
    for (int c = re.cols().begin(); c != re.cols().end(); ++c) {
    for (int r = re.rows().begin(); r != re.rows().end(); ++r) {
      _Matrix(r, c) = _Matrix(c, r) = 0.5 * (_Matrix(r, c) + _Matrix(c, r));
    }
    }
  }
};

void irtkMatrix::MakeSymmetric()
{
  if (_rows != _cols) {
    cerr << "irtkMatrix::MakeSymmetric: Matrix must be square" << endl;
    exit(1);
  }
  parallel_for(blocked_range2d<int>(0, _rows, 0, _cols), MakeMatrixSymmetric(*this));
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
ostream& operator<< (ostream& os, const irtkMatrix &m)
{
  int index, i, j;

  // Write header in ascii
  os << "irtkMatrix " << m._rows << " x " << m._cols << endl;

  // Allocate temporary memory
  double *data  = new double [m._rows*m._cols];

  // Convert data
  index = 0;
  for (j = 0; j < m._cols; j++) {
    for (i = 0; i < m._rows; i++) {
      data[index] = m._matrix[j][i];
      index++;
    }
  }

#ifndef WORDS_BIGENDIAN
  swap64((char *)data, (char *)data, m._rows*m._cols);
#endif

  // Write binary data
  os.write((char *)data, m._rows*m._cols*sizeof(double));

  // Free temporary memory
  delete [] data;

  return os;
}

// -----------------------------------------------------------------------------
istream& operator>> (istream& is, irtkMatrix &m)
{
  int index, i, j, cols, rows;
  char buffer[255];

  // Read header
  is >> buffer;
  if (strcmp(buffer, "irtkMatrix") != 0) {
    cerr << "irtkMatrix: Can't read file " << buffer << endl;
    exit(1);
  }

  // Read size
  is >> rows;
  is >> buffer;
  is >> cols;

  // Allocate matrix
  m = irtkMatrix(rows, cols);

  // Read header, skip comments
  is.get(buffer, 255);
  is.clear();
  is.seekg(1, ios::cur);

  // Allocate temporary memory
  double *data  = new double [m._rows*m._cols];

  // Read binary data
  is.read((char *)data, m._rows*m._cols*sizeof(double));

#ifndef WORDS_BIGENDIAN
  swap64((char *)data, (char *)data, m._rows*m._cols);
#endif

  // Convert data
  index = 0;
  for (j = 0; j < m._cols; j++) {
    for (i = 0; i < m._rows; i++) {
      m._matrix[j][i] = data[index];
      index++;
    }
  }

  // Free temporary memory
  delete []data;

  return is;
}

// -----------------------------------------------------------------------------
irtkCofstream& operator<< (irtkCofstream& to, const irtkMatrix &m)
{
  to.WriteAsChar("irtkMatrix", 11);

  to.WriteAsInt(&m._rows, 1);
  to.WriteAsInt(&m._cols, 1);

  to.WriteAsDouble(m.RawPointer(), m._rows * m._cols);

  return to;
}

// -----------------------------------------------------------------------------
irtkCifstream& operator>> (irtkCifstream& from, irtkMatrix &m)
{
  char keyword[11];
  from.ReadAsChar(keyword, 11);
  if (strncmp(keyword, "irtkMatrix", 11) != 0) {
    keyword[10] = '\0'; // ensure it is null terminated
    cerr << "irtkMatrix: Can't read file " << keyword << endl;
    exit(1);
  }

  int rows = 0, cols = 0;
  from.ReadAsInt(&rows, 1);
  from.ReadAsInt(&cols, 1);

  m.Initialize(rows, cols);

  from.ReadAsDouble(m.RawPointer(), rows * cols);

  return from;
}

// -----------------------------------------------------------------------------
void irtkMatrix::Print(irtkIndent indent) const
{
  int i, j;

  cout << indent << "irtkMatrix " << _rows << " x " << _cols << endl;
  ++indent;
  cout.setf(ios::right);
  cout.setf(ios::fixed);
  cout.precision(4);
  for (i = 0; i < _rows; i++) {
    cout << indent;
    for (j = 0; j < _cols; j++) {
      cout << setw(15) << _matrix[j][i] << " ";
    }
    cout << endl;
  }
  cout.precision(6);
  cout.unsetf(ios::right);
  cout.unsetf(ios::fixed);
}

// -----------------------------------------------------------------------------
void irtkMatrix::Read(const char *filename)
{
  // Open file stream
  ifstream from(filename, ios::in | ios::binary);

  // Check whether file opened ok
  if (!from) {
    cerr << "irtkMatrix::Read: Can't open file " << filename << endl;
    exit(1);
  }

  // Read matrix
  from >> *this;
}

// -----------------------------------------------------------------------------
void irtkMatrix::Write(const char *filename) const
{
  // Open file stream
  ofstream to(filename, ios::out | ios::binary);

  // Check whether file opened ok
  if (!to) {
    cerr << "irtkMatrix::Write: Can't open file " << filename << endl;
    exit(1);
  }

  // Write matrix
  to << *this;
}

// -----------------------------------------------------------------------------
void irtkMatrix::Import(const char *filename, int rows, int cols)
{
  int i, j;

  // Open file stream
  ifstream from(filename);

  // Check whether file opened ok
  if (!from) {
    cerr << "irtkMatrix::Read: Can't open file " << filename << endl;
    exit(1);
  }

  // Initialize matrix
  this->Initialize(rows, cols);

  // Read matrix
  for (i = 0; i < _rows; i++) {
    for (j = 0; j < _cols; j++) {
      from >> _matrix[j][i];
    }
  }
}

#ifdef HAVE_EIGEN
// -----------------------------------------------------------------------------
void irtkMatrix::Matrix2Eigen(Eigen::MatrixXd &m) const
{
   int i, j;

   for (i = 0; i < _rows; i++) {
     for (j = 0; j < _cols; j++) {
       m(i, j) = _matrix[j][i];
     }
   }
}

// -----------------------------------------------------------------------------
void irtkMatrix::Eigen2Matrix(Eigen::MatrixXd &m)
{
  int i, j;

  for (i = 0; i < _rows; i++) {
    for (j = 0; j < _cols; j++) {
      _matrix[j][i] = m(i, j);
    }
  }
}
#endif

#ifdef HAVE_MATLAB
// -----------------------------------------------------------------------------
mxArray *irtkMatrix::MxArray() const
{
  irtkMatlab::Initialize();
  mxArray *m = mxCreateDoubleMatrix(_rows, _cols, mxREAL);
  memcpy(mxGetPr(m), RawPointer(), NumberOfElements() * sizeof(double));
  return m;
}

// -----------------------------------------------------------------------------
bool irtkMatrix::WriteMAT(const char *fname, const char *varname) const
{
  irtkMatlab::Initialize();
  MATFile *fp = matOpen(fname, "w");
  if (fp == NULL) return false;
  mxArray *m = MxArray();
  if (matPutVariable(fp, varname, m) != 0) {
    mxDestroyArray(m);
    matClose(fp);
    return false;
  }
  mxDestroyArray(m);
  return (matClose(fp) == 0);
}
#endif

////////////////////////////////////////////////////////////////////////////////
// Means
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
irtkMatrix LogEuclideanMean(int n, const irtkMatrix *matrices, const double *weights)
{
  irtkAssert(n > 0, "at least one matrix must be given");
  if (n == 1) return matrices[0];

  double totalWeight;
  if (weights) {
    totalWeight = .0;
    for (int i = 0; i < n; ++i) {
      totalWeight += weights[i];
    }
  } else {
    totalWeight = static_cast<double>(n);
  }
  if (totalWeight <= .0) {
    cerr << "itkMatrix::LogEuclideanMean: Sum of weights must be positive" << endl;
    exit(1);
  }

  irtkMatrix sumLogs(matrices[0].Rows(), matrices[0].Cols());
  for (int i = 0; i < n; ++i) {
    sumLogs += matrices[i].Log() * weights[i];
  }
  sumLogs /= totalWeight;

  return sumLogs.Exp();
}

// -----------------------------------------------------------------------------
irtkMatrix BiInvariantMean(int               n,
                           const irtkMatrix *matrices,
                           const double     *weights,
                           int               niter,
                           double            tolerance,
                           const irtkMatrix *mu0)
{
  irtkAssert(n > 0, "at least one matrix must be given");
  if (n == 1) return matrices[0];

  int    i, iter = 0;
  double totalWeight, normLogDeltaMu;
  irtkMatrix mu, muInv, deltaMu, deltaM, sumLogs;

  // Normalize weights
  if (weights) {
    totalWeight = .0;
    for (i = 0; i < n; ++i) {
      totalWeight += weights[i];
    }
  } else {
    totalWeight = static_cast<double>(n);
  }
  if (totalWeight <= .0) {
    cerr << "itkMatrix::BiInvariantMean: Sum of weights must be positive" << endl;
    exit(1);
  }

  // Check that determinants of all matrices have the same sign
  double sign = sgn(matrices[0].Det());
  for (i = 1; i < n; ++i) {
    if (sgn(matrices[i].Det()) != sign) {
      cerr << "irtkMatrix::BiInvariantMean: Sign of determinant is not the same for all matrices" << endl;
      exit(1);
    }
  }

  // Barycentric fixed point iteration
  mu = (mu0 ? *mu0 : matrices[0]);

  do {

    muInv = mu.Inverse();

    sumLogs.Initialize(mu.Rows(), mu.Cols()); // Reset sum to zero
    if (weights) {
      for (i = 0; i < n; ++i) {
        if (weights[i] == .0) continue;
        deltaM = muInv * matrices[i];
        sumLogs += deltaM.Log() * weights[i];
      }
    } else {
      for (i = 0; i < n; ++i) {
        deltaM = muInv * matrices[i];
        sumLogs += deltaM.Log();
      }
    }
    sumLogs /= totalWeight;
    deltaMu = sumLogs.Exp();

    mu = mu * deltaMu;

    normLogDeltaMu = deltaMu.Log().InfinityNorm();

  } while (normLogDeltaMu > tolerance && ++iter < niter);

  return mu;
}

////////////////////////////////////////////////////////////////////////////////
// Affine transformation matrices
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
// R = (Rz Ry Rx)^T
void RigidParametersToMatrix(double tx,    double ty,    double tz,
                             double cosrx, double cosry, double cosrz,
                             double sinrx, double sinry, double sinrz, irtkMatrix &m)
{
  m.Initialize(4, 4);

  m(0, 0) = cosry * cosrz;
  m(0, 1) = cosry * sinrz;
  m(0, 2) = -sinry;
  m(0, 3) = tx;

  m(1, 0) = (sinrx * sinry * cosrz - cosrx * sinrz);
  m(1, 1) = (sinrx * sinry * sinrz + cosrx * cosrz);
  m(1, 2) = sinrx * cosry;
  m(1, 3) = ty;

  m(2, 0) = (cosrx * sinry * cosrz + sinrx * sinrz);
  m(2, 1) = (cosrx * sinry * sinrz - sinrx * cosrz);
  m(2, 2) = cosrx * cosry;
  m(2, 3) = tz;

  m(3, 0) = 0.0;
  m(3, 1) = 0.0;
  m(3, 2) = 0.0;
  m(3, 3) = 1.0;
}

// -----------------------------------------------------------------------------
void RigidParametersToMatrix(double tx,  double ty,  double tz,
                             double rx,  double ry,  double rz, irtkMatrix &m)
{
  const double cosrx = cos(rx);
  const double cosry = cos(ry);
  const double cosrz = cos(rz);
  const double sinrx = sin(rx);
  const double sinry = sin(ry);
  const double sinrz = sin(rz);

  RigidParametersToMatrix(tx, ty, tz, cosrx, cosry, cosrz, sinrx, sinry, sinrz, m);
}

// -----------------------------------------------------------------------------
void MatrixToEulerAngles(const irtkMatrix &m, double &rx,  double &ry,  double &rz)
{
  const double TOL = 0.000001;
  double tmp;

  tmp = asin(-1 * m(0, 2));

  // asin returns values for tmp in range -pi/2 to +pi/2, i.e. cos(tmp) >=
  // 0 so the division by cos(tmp) in the first part of the if clause was
  // not needed.
  if (fabs(cos(tmp)) > TOL) {
    rx = atan2(m(1,2), m(2,2));
    ry = tmp;
    rz = atan2(m(0,1), m(0,0));
  } else {
    //m(0,2) is close to +1 or -1
    rx = atan2(-1.0*m(0,2)*m(1,0), -1.0*m(0,2)*m(2,0));
    ry = tmp;
    rz = 0;
  }
}

// -----------------------------------------------------------------------------
void MatrixToRigidParameters(const irtkMatrix &m,
                             double &tx,  double &ty,  double &tz,
                             double &rx,  double &ry,  double &rz)
{
  const double TOL = 0.000001;
  double tmp;

  tx = m(0, 3);
  ty = m(1, 3);
  tz = m(2, 3);

  tmp = asin(-1 * m(0, 2));

  // asin returns values for tmp in range -pi/2 to +pi/2, i.e. cos(tmp) >=
  // 0 so the division by cos(tmp) in the first part of the if clause was
  // not needed.
  if (fabs(cos(tmp)) > TOL) {
    rx = atan2(m(1,2), m(2,2));
    ry = tmp;
    rz = atan2(m(0,1), m(0,0));
  } else {
    //m(0,2) is close to +1 or -1
    rx = atan2(-1.0*m(0,2)*m(1,0), -1.0*m(0,2)*m(2,0));
    ry = tmp;
    rz = 0;
  }
}

// -----------------------------------------------------------------------------
void AffineParametersToMatrix(double tx,  double ty,  double tz,
                              double rx,  double ry,  double rz,
                              double sx,  double sy,  double sz,
                              double sxy, double sxz, double syz, irtkMatrix &m)
{
  irtkMatrix tmp(4, 4);

  // Construct rigid transformation matrix
  RigidParametersToMatrix(tx, ty, tz, rx, ry, rz, m);

  // Pre-multiply with shearing transformation
  tmp.Ident();
  tmp(0, 1) = tan(sxy);
  tmp(0, 2) = tan(sxz);
  tmp(1, 2) = tan(syz);
  m *= tmp;

  // Pre-multiply with scaling transformation
  tmp.Ident();
  tmp(0, 0) = sx;
  tmp(1, 1) = sy;
  tmp(2, 2) = sz;
  m *= tmp;
}

// -----------------------------------------------------------------------------
void MatrixToAffineParameters(const irtkMatrix &m,
                              double &tx,  double &ty,  double &tz,
                              double &rx,  double &ry,  double &rz,
                              double &sx,  double &sy,  double &sz,
                              double &sxy, double &sxz, double &syz)
{
  const double TOL = 0.000001;
  double tansxy, tansxz, tansyz;

  if (fabs(m(3, 3) - 1.0) > TOL) {
    cerr << "MatrixToAffineParameters: Value at m(3, 3) must equal 1." << endl;
    exit(1);
  }
  if (fabs(m.Det()) < TOL) {
    cerr << "MatrixToAffineParameters: Matrix is singular or very close to singular!." << endl;
    exit(1);
  }

  // First Part Of Graphics Gems Code Ignored Because It Relates To
  // Perspective Transformation.
  if (fabs(m(3, 0)) > TOL ||
      fabs(m(3, 1)) > TOL ||
      fabs(m(3, 2)) > TOL ) {
    cerr << "MatrixToAffineParameters: Matrix contains perspective distortion." << endl;
    exit(1);
  }

  irtkMatrix copy(m);

  // Get scale and shear by manipulating the columns of the upper left 3x3
  // sub-matrix.
  irtkVector col_0, col_1, col_2;
  col_0.Initialize(3);
  col_1.Initialize(3);
  col_2.Initialize(3);
  for (int i = 0; i < 3; ++i) {
    col_0(i) = copy(i, 0);
    col_1(i) = copy(i, 1);
    col_2(i) = copy(i, 2);
  }

  // Compute X scale factor and normalize first col.
  sx = col_0.Norm();
  col_0 /= sx;

  // Compute XY shear factor and make 2nd col orthogonal to 1st.
  tansxy = col_0.ScalarProduct(col_1);
  col_1 = col_1 - col_0 * tansxy;

  // Actually, tansxy and col_1 are still to large by a factor of sy.
  // Now, compute Y scale and normalize 2nd col and rescale tansxy.
  sy = col_1.Norm();
  col_1  /= sy;
  tansxy /= sy;

  // Compute XZ and YZ shears, orthogonalize 3rd col
  tansxz = col_0.ScalarProduct(col_2);
  col_2 = col_2 - col_0 * tansxz;

  tansyz = col_1.ScalarProduct(col_2);
  col_2 = col_2 - col_1 * tansyz;

  // Actually, tansxz, tansyz and col_2 are still to large by a factor of
  // sz.  Next, get Z scale, normalize 3rd col and scale tansxz and tansyz.
  sz = col_2.Norm();
  col_2  /= sz;
  tansxz /= sz;
  tansyz /= sz;

  // At this point, the columns are orthonormal.  Check for a coordinate
  // system flip.  If the determinant is -1, then negate the matrix and the
  // scaling factors.
  irtkVector col_1_x_col_2;
  col_1_x_col_2.Initialize(3);
  col_1_x_col_2 = col_1.CrossProduct(col_2);

  if (col_0.ScalarProduct(col_1_x_col_2) < 0) {
    sx *= -1;
    sy *= -1;
    sz *= -1;
    col_0 *= -1;
    col_1 *= -1;
    col_2 *= -1;
  }

  // Retrieve the shear angles in degrees.
  sxy = atan(tansxy);
  sxz = atan(tansxz);
  syz = atan(tansyz);

  // Now get the rigid transformation parameters.
  // Put the rotation matrix components into the upper left 3x3 submatrix.
  for (int i = 0; i < 3; ++i) {
    copy(i, 0) = col_0(i);
    copy(i, 1) = col_1(i);
    copy(i, 2) = col_2(i);
  }

  MatrixToRigidParameters(copy, tx, ty, tz, rx, ry, rz);
}

// -----------------------------------------------------------------------------
// Based on Insight Journal code: http://hdl.handle.net/10380/3299
//
// Solve Qa = C
//
//          | sum( q(1,i)*q(1,i) ) sum( q(1,i)*q(2,i) ) sum( q(1,i)*q(3,i) ) sum( q(1,i) ) |
// Q[4x4] = | sum( q(2,i)*q(1,i) ) sum( q(2,i)*q(2,i) ) sum( q(2,i)*q(3,i) ) sum( q(2,i) ) |
//          | sum( q(3,i)*q(1,i) ) sum( q(3,i)*q(2,i) ) sum( q(3,i)*q(3,i) ) sum( q(3,i) ) |
//          | sum( q(1,i) )        sum( q(2,i) )        sum( q(3,i) )        sum( q(4,i) ) |
//
//          | a(1,1) a(2,1) a(3,1) |
// a[4x3] = | a(1,2) a(2,2) a(3,2) |
//          | a(1,3) a(2,3) a(3,3) |
//          | a(1,4) a(2,4) a(3,4) |
//
//          | sum( q(1,i)*p(1,i) )  sum( q(1,i)*p(2,i) ) sum( q(1,i)*p(3,i) )  |
// C[4x3] = | sum( q(2,i)*p(1,i) )  sum( q(2,i)*p(2,i) ) sum( q(2,i)*p(3,i) )  |
//          | sum( q(3,i)*p(1,i) )  sum( q(3,i)*p(2,i) ) sum( q(3,i)*p(3,i) )  |
//          | sum( p(1,i) )         sum( p(2,i) )        sum( p(3,i) )         |
//
irtkMatrix ApproximateAffineMatrix(const irtkPointSet &target, const irtkPointSet &source, const irtkVector &weight)
{
  // Check arguments
  const int no = source.Size();
  if (target.Size() != no) {
    cerr << "ApproximateAffineMatrix: Size mismatch between number of landark pairs!" << endl;
    exit(1);
  }
  if (no < 4) {
    cerr << "ApproximateAffineMatrix: Must have at least four points" << endl;
    exit(1);
  }

  double w, wnorm = .0;
  if (weight.Rows() != 0) {
    if (weight.Rows() != no) {
      cerr << "ApproximateAffineMatrix: Size mismatch between number of landmark pairs and weights" << endl;
      exit(1);
    }
    for (int i = 0; i < no; ++i) {
      wnorm += weight(i) * weight(i);
    }
    wnorm = sqrt(wnorm);
  }

  // Output matrix
  irtkMatrix A(4, 4);
  A.Ident();

  // ---------------------------------------------------------------------------
  // Implementation using Eigen library
#if defined(HAVE_EIGEN)

  typedef Eigen::Matrix<double, 4, Eigen::Dynamic> Matrix4xN;
  typedef Eigen::Matrix<double, 3, Eigen::Dynamic> Matrix3xN;
  typedef Eigen::Matrix<double, 4, 4>              Matrix4x4;
  typedef Eigen::Matrix<double, 4, 3>              Matrix4x3;
  typedef Eigen::Matrix<double, 4, 1>              Vector4;
  typedef Eigen::FullPivHouseholderQR<Matrix4x4>   DenseSolver;

  // Convert point sets to matrices and apply landmark weights
  Matrix4xN q(4, no);
  Matrix3xN p(3, no);

  for (int i = 0; i < no; ++i) {
    // Normalized landmark weight
    w = (weight.Rows() != 0 ? weight(i) / wnorm : 1.0);
    // Target point
    q(0, i) = w * target(i)._x;
    q(1, i) = w * target(i)._y;
    q(2, i) = w * target(i)._z;
    q(3, i) = w;
    // Source point
    p(0, i) = w * source(i)._x;
    p(1, i) = w * source(i)._y;
    p(2, i) = w * source(i)._z;
  }

  // Solve Qa = C
  Matrix4x4 Q;
  Matrix4x3 C;
  Vector4   a;

  Q = Matrix4x4::Zero();
  C = Matrix4x3::Zero();
  for (int i = 0; i < no; ++i) {
    Q += q.col(i) * q.col(i).transpose();
    C += q.col(i) * p.col(i).transpose();
  }

  DenseSolver solver(Q);
  for (int i = 0; i < 3; ++i) {
    a = solver.solve(C.col(i));
    A(i, 0) = a(0), A(i, 1) = a(1), A(i, 2) = a(2), A(i, 3) = a(3);
  }

  // ---------------------------------------------------------------------------
  // Implementation using VNL/VXL library
#elif defined(HAVE_VNL)

  // Convert point sets to matrices and apply landmark weights
  vnl_matrix<double> q(4, no, 0);
  vnl_matrix<double> p(3, no, 0);

  for (int i = 0; i < no; ++i) {
    // Normalized landmark weight
    w = (weight.Rows() != 0 ? weight(i) / wnorm : 1.0);
    // Target point
    q(0, i) = w * target(i)._x;
    q(1, i) = w * target(i)._y;
    q(2, i) = w * target(i)._z;
    q(3, i) = w;
    // Source point
    p(0, i) = w * source(i)._x;
    p(1, i) = w * source(i)._y;
    p(2, i) = w * source(i)._z;
  }

  // Solve Qa = C
  vnl_matrix<double> Q(4, 4, 0);
  vnl_matrix<double> C(4, 3, 0);
  for (int i = 0; i < no; ++i) {
    vnl_matrix<double> qTemp(4, 1);
    vnl_matrix<double> pTemp(1, 3);
    for (int k = 0; k < 4; ++k) qTemp(k, 0) = q.get(k, i);
    for (int k = 0; k < 3; ++k) pTemp(0, k) = p.get(k, i);
    Q = Q + qTemp * qTemp.transpose();
    C = C + qTemp * pTemp;
  }

  vnl_matrix<double> matT = vnl_qr<double>(Q).solve(C);

  // Set output matrix
  for (int r = 0; r < 3; ++r)
  for (int c = 0; c < 4; ++c) {
    A(r, c) = matT(c, r);
  }

#else
  cerr << "ApproximateAffineMatrix requires either Eigen3 (recommended) or VNL/VXL library" << endl;
  exit(1);
#endif

  return A;
}

// -----------------------------------------------------------------------------
irtkMatrix ApproximateAffineMatrix(const irtkPointSet &target, const irtkPointSet &source)
{
  return ApproximateAffineMatrix(target, source, irtkVector());
}

// -----------------------------------------------------------------------------
void OrthoNormalize3x3(irtkMatrix &mat)
{
  for (int i = 0; i < 3; i++) {
    irtkVector3D<double> vi(mat(0, i), mat(1, i), mat(2, i));
    // normalize column i
    vi.Normalize();
    // make other vectors linearly independent of column i
    for (int j = i + 1; j < 3; j++) {
      irtkVector3D<double> vj(mat(0, j), mat(1, j), mat(2, j));
      const double s = irtkVector3D<double>::DotProduct(vj, vi);
      mat(0, j) = vj._x - s * vi._x;
      mat(1, j) = vj._y - s * vi._y;
      mat(2, j) = vj._z - s * vi._z;
    }
  }
}
