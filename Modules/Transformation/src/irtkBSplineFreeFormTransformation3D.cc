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

#include <irtkTransformation.h>
#include <irtkImageToInterpolationCoefficients.h>


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkBSplineFreeFormTransformation3D::irtkBSplineFreeFormTransformation3D()
:
  irtkFreeFormTransformation3D(_FFD, &_FFD2D)
{
}

// -----------------------------------------------------------------------------
irtkBSplineFreeFormTransformation3D
::irtkBSplineFreeFormTransformation3D(double x1, double y1, double z1,
                                      double x2, double y2, double z2,
                                      double dx, double dy, double dz,
                                      double *xaxis, double *yaxis, double *zaxis)
:
  irtkFreeFormTransformation3D(_FFD, &_FFD2D)
{
  Initialize(DefaultAttributes(x1, y1, z1, x2, y2, z2, dx, dy, dz, xaxis, yaxis, zaxis));
}

// -----------------------------------------------------------------------------
irtkBSplineFreeFormTransformation3D
::irtkBSplineFreeFormTransformation3D(const irtkImageAttributes &attr, double dx, double dy, double dz)
:
  irtkFreeFormTransformation3D(_FFD, &_FFD2D)
{
  Initialize(attr, dx, dy, dz);
}

// -----------------------------------------------------------------------------
irtkBSplineFreeFormTransformation3D
::irtkBSplineFreeFormTransformation3D(const irtkBaseImage &target, double dx, double dy, double dz)
:
  irtkFreeFormTransformation3D(_FFD, &_FFD2D)
{
  Initialize(target.Attributes(), dx, dy, dz);
}

// -----------------------------------------------------------------------------
irtkBSplineFreeFormTransformation3D
::irtkBSplineFreeFormTransformation3D(const irtkGenericImage<double> &image, bool disp)
:
  irtkFreeFormTransformation3D(_FFD, &_FFD2D)
{
  Initialize(image, disp);
}

// -----------------------------------------------------------------------------
irtkBSplineFreeFormTransformation3D
::irtkBSplineFreeFormTransformation3D(const irtkBSplineFreeFormTransformation3D &ffd)
:
  irtkFreeFormTransformation3D(ffd, _FFD, &_FFD2D)
{
  if (_NumberOfDOFs > 0) InitializeInterpolator();
}

// -----------------------------------------------------------------------------
irtkBSplineFreeFormTransformation3D::~irtkBSplineFreeFormTransformation3D()
{
}

// =============================================================================
// Approximation/Interpolation
// =============================================================================

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformation3D
::ApproximateDOFs(const double *wx, const double *wy, const double *wz, const double *,
                  const double *dx, const double *dy, const double *dz, int no)
{
  int    i, j, k, ci, cj, ck, A, B, C;
  double x, y, z, w[3], basis, sum;

  // Allocate memory
  Vector ***data = CAllocate<Vector>(_x, _y, _z);
  double ***norm = CAllocate<double>(_x, _y, _z);

  // Initial loop: Calculate change of control points
  for (int idx = 0; idx < no; ++idx) {
    if (dx[idx] == .0 && dy[idx] == .0 && dz[idx] == .0) continue;

    x = wx[idx], y = wy[idx], z = wz[idx];
    this->WorldToLattice(x, y, z);

    i = static_cast<int>(floor(x));
    j = static_cast<int>(floor(y));

    A = Kernel::VariableToIndex(x - i);
    B = Kernel::VariableToIndex(y - j);
    --i, --j;

    // 2D
    if (_z == 1) {

      sum = .0;
      for (int b = 0; b <= 3; ++b) {
        w[1] = Kernel::LookupTable[B][b];
        for (int a = 0; a <= 3; ++a) {
          w[0] = Kernel::LookupTable[A][a] * w[1];
          sum += w[0] * w[0];
        }
      }

      for (int b = 0; b <= 3; ++b) {
        cj = j + b;
        if (cj < 0 || cj >= _y) continue;
        w[1] = Kernel::LookupTable[B][b];
        for (int a = 0; a <= 3; ++a) {
          ci = i + a;
          if (ci < 0 || ci >= _x) continue;
          w[0]  = Kernel::LookupTable[A][a] * w[1];
          basis = w[0] * w[0];
          norm[0][cj][ci] += basis;
          basis *= w[0] / sum;
          data[0][cj][ci]._x += basis * dx[idx];
          data[0][cj][ci]._y += basis * dy[idx];
          data[0][cj][ci]._z += basis * dz[idx];
        }
      }

    // 3D
    } else {

      k = static_cast<int>(floor(z));
      C = Kernel::VariableToIndex(z - k);
      --k;

      sum = .0;
      for (int c = 0; c <= 3; ++c) {
        w[2] = Kernel::LookupTable[C][c];
        for (int b = 0; b <= 3; ++b) {
          w[1] = Kernel::LookupTable[B][b] * w[2];
          for (int a = 0; a <= 3; ++a) {
            w[0] = Kernel::LookupTable[A][a] * w[1];
            sum += w[0] * w[0];
          }
        }
      }

      for (int c = 0; c <= 3; ++c) {
        ck = k + c;
        if (ck < 0 || ck >= _z) continue;
        w[2] = Kernel::LookupTable[C][c];
        for (int b = 0; b <= 3; ++b) {
          cj = j + b;
          if (cj < 0 || cj >= _y) continue;
          w[1] = Kernel::LookupTable[B][b] * w[2];
          for (int a = 0; a <= 3; ++a) {
            ci = i + a;
            if (ci < 0 || ci >= _x) continue;
            w[0]  = Kernel::LookupTable[A][a] * w[1];
            basis = w[0] * w[0];
            norm[ck][cj][ci] += basis;
            basis *= w[0] / sum;
            data[ck][cj][ci]._x += basis * dx[idx];
            data[ck][cj][ci]._y += basis * dy[idx];
            data[ck][cj][ci]._z += basis * dz[idx];
          }
        }
      }

    }
  }

  // Final loop: Calculate new control points
  const Vector zero(.0);

  Vector *in  = data[0][0];
  double *div = norm[0][0];
  Vector *out = _CPImage.Data();

  for (int cp = 0; cp < NumberOfCPs(); ++cp, ++out, ++in, ++div) {
    (*out) = ((*div) ? ((*in) / (*div)) : zero);
  }

  // Deallocate memory
  Deallocate(data);
  Deallocate(norm);

  this->Changed(true);
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformation3D
::ApproximateDOFsGradient(const double *wx, const double *wy, const double *wz, const double *,
                          const double *dx, const double *dy, const double *dz,
                          int no, double *gradient, double weight) const
{
  int    i, j, k, ci, cj, ck, A, B, C;
  double x, y, z, w[3];

  // Allocate memory
  Vector ***data = CAllocate<Vector>(_x, _y, _z);

  // Initial loop: Calculate change of control points
  for (int idx = 0; idx < no; ++idx) {
    if (dx[idx] == .0 && dy[idx] == .0 && dz[idx] == .0) continue;

    x = wx[idx], y = wy[idx], z = wz[idx];
    this->WorldToLattice(x, y, z);

    i = static_cast<int>(floor(x));
    j = static_cast<int>(floor(y));
    A = Kernel::VariableToIndex(x - i);
    B = Kernel::VariableToIndex(y - j);
    --i, --j;

    // 2D
    if (_z == 1) {

      for (int b = 0; b <= 3; ++b) {
        cj = j + b;
        if (cj < 0 || cj >= _y) continue;
        w[1] = Kernel::LookupTable[B][b];
        for (int a = 0; a <= 3; ++a) {
          ci = i + a;
          if (ci < 0 || ci >= _x) continue;
          w[0] = Kernel::LookupTable[A][a] * w[1];
          data[0][cj][ci]._x += w[0] * dx[idx];
          data[0][cj][ci]._y += w[0] * dy[idx];
          data[0][cj][ci]._z += w[0] * dz[idx];
        }
      }

    // 3D
    } else {

      k = static_cast<int>(floor(z));
      C = Kernel::VariableToIndex(z - k);
      --k;

      for (int c = 0; c <= 3; ++c) {
        ck = k + c;
        if (ck < 0 || ck >= _z) continue;
        w[2] = Kernel::LookupTable[C][c];
        for (int b = 0; b <= 3; ++b) {
          cj = j + b;
          if (cj < 0 || cj >= _y) continue;
          w[1] = Kernel::LookupTable[B][b] * w[2];
          for (int a = 0; a <= 3; ++a) {
            ci = i + a;
            if (ci < 0 || ci >= _x) continue;
            w[0] = Kernel::LookupTable[A][a] * w[1];
            data[ck][cj][ci]._x += w[0] * dx[idx];
            data[ck][cj][ci]._y += w[0] * dy[idx];
            data[ck][cj][ci]._z += w[0] * dz[idx];
          }
        }
      }

    }
  }

  // Final loop
  int           xdof, ydof, zdof;
  const Vector *grad = data[0][0];

  for (int cp = 0; cp < NumberOfCPs(); ++cp, ++grad) {
    if (!IsActive(cp)) continue;
    this->IndexToDOFs(cp, xdof, ydof, zdof);
    gradient[xdof] += weight * grad->_x;
    gradient[ydof] += weight * grad->_y;
    gradient[zdof] += weight * grad->_z;
  }

  // Deallocate memory
  Deallocate(data);
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformation3D
::Interpolate(const double *dx, const double *dy, const double *dz)
{
  if (dz) {
    for (int idx = 0; idx < NumberOfCPs(); ++idx) {
      _CPImage(idx) = Vector(dx[idx], dy[idx], dz[idx]);
    }
  } else {
    for (int idx = 0; idx < NumberOfCPs(); ++idx) {
      _CPImage(idx) = Vector(dx[idx], dy[idx], .0);
    }
  }
  if (_z == 1) ConvertToSplineCoefficients(3, _CPImage, 0, 0);
  else         ConvertToSplineCoefficients(3, _CPImage,    0);
  this->Changed(true);
}

// =============================================================================
// Lattice
// =============================================================================

// -----------------------------------------------------------------------------
int irtkBSplineFreeFormTransformation3D::GetXAfterSubdivision() const
{
  return (_x == 1 ? _x : 2 * _x - 1);
}

// -----------------------------------------------------------------------------
int irtkBSplineFreeFormTransformation3D::GetYAfterSubdivision() const
{
  return (_y == 1 ? _y : 2 * _y - 1);
}

// -----------------------------------------------------------------------------
int irtkBSplineFreeFormTransformation3D::GetZAfterSubdivision() const
{
  return (_z == 1 ? _z : 2 * _z - 1);
}

// -----------------------------------------------------------------------------
int irtkBSplineFreeFormTransformation3D::GetTAfterSubdivision() const
{
  return _t;
}

// -----------------------------------------------------------------------------
double irtkBSplineFreeFormTransformation3D::GetXSpacingAfterSubdivision() const
{
  return (_x == 1 ? _dx : 0.5 * _dx);
}

// -----------------------------------------------------------------------------
double irtkBSplineFreeFormTransformation3D::GetYSpacingAfterSubdivision() const
{
  return (_y == 1 ? _dy : 0.5 * _dy);
}

// -----------------------------------------------------------------------------
double irtkBSplineFreeFormTransformation3D::GetZSpacingAfterSubdivision() const
{
  return (_z == 1 ? _dz : 0.5 * _dz);
}

// -----------------------------------------------------------------------------
double irtkBSplineFreeFormTransformation3D::GetTSpacingAfterSubdivision() const
{
  return _dt;
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformation3D
::Subdivide(bool subdivide_x, bool subdivide_y, bool subdivide_z, bool)
{
  if (!subdivide_x && !subdivide_y && !subdivide_z) return;

  if (_x == 1) subdivide_x = false;
  if (_y == 1) subdivide_y = false;
  if (_z == 1) subdivide_z = false;

  // Weights for subdivision
  const double w1[2][3] = {{0.125, 0.75, 0.125}, {0.0, 0.5, 0.5}};
  const double w2[2][3] = {{0.0,   1.0,  0.0},   {0.0, 0.0, 0.0}};

  const double (*wx)[3] = (subdivide_x ? w1 : w2);
  const double (*wy)[3] = (subdivide_y ? w1 : w2);
  const double (*wz)[3] = (subdivide_z ? w1 : w2);

  // Size of new control point grid
  const int X = (subdivide_x ? (2 * _x - 1) : _x);
  const int Y = (subdivide_y ? (2 * _y - 1) : _y);
  const int Z = (subdivide_z ? (2 * _z - 1) : _z);

  // Allocate memory for new control points
  CPValue ***data = Allocate<CPValue>(X, Y, Z);

  // Limits for inner loops
  const int I2 = (subdivide_x ? 2 : 1); // s.t. i+i2-1 == i if no subdivision
  const int J2 = (subdivide_y ? 2 : 1);
  const int K2 = (subdivide_z ? 2 : 1);

  // Compute control point values on subdivided grid
  int si, sj, sk, I1, J1, K1;
  for (int k = 0; k < _z; ++k) {
    if (subdivide_z) sk = 2 * k, K1 = int(k < _z-1);
    else             sk =     k, K1 = 0;
    for (int j = 0; j < _y; ++j) {
      if (subdivide_y) sj = 2 * j, J1 = int(j < _y-1);
      else             sj =     j, J1 = 0;
      for (int i = 0; i < _x; ++i) {
        if (subdivide_x) si = 2 * i, I1 = int(i < _x-1);
        else             si =     i, I1 = 0;
        for (int k1 = 0; k1 <= K1; ++k1)
        for (int j1 = 0; j1 <= J1; ++j1)
        for (int i1 = 0; i1 <= I1; ++i1) {
          for (int k2 = (subdivide_z ? 0 : K2); k2 <= K2; ++k2)
          for (int j2 = (subdivide_y ? 0 : J2); j2 <= J2; ++j2)
          for (int i2 = (subdivide_x ? 0 : I2); i2 <= I2; ++i2) {
            // If a dimension is not subdivided,
            // 1) si + i1    == i
            // 2) i + i2 - 1 == i
            // 3) w[i1][i2]  == 1.0
            data[sk+k1][sj+j1][si+i1] += wx[i1][i2] * wy[j1][j2] * wz[k1][k2]
                                       * _CPValue->Get(i+i2-1, j+j2-1, k+k2-1);
          }
        }
      }
    }
  }

  // Initialize subdivided free-form transformation
  irtkImageAttributes attr = this->Attributes();
  attr._x = X;
  attr._y = Y;
  attr._z = Z;
  if (subdivide_x) attr._dx *= 0.5;
  if (subdivide_y) attr._dy *= 0.5;
  if (subdivide_z) attr._dz *= 0.5;
  this->Initialize(attr);

  // Copy subdivided control point data
  memcpy(_CPImage.Data(), data[0][0], X * Y * Z * sizeof(CPValue));

  // Deallocate temporary memory
  Deallocate(data);
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformation3D
::BoundingBox(int cp, double &x1, double &y1, double &z1,
                      double &x2, double &y2, double &z2, double fraction) const
{
  int i, j, k;
  IndexToLattice(cp, i, j, k);

  fraction *= 2;
  x1 = i - fraction;
  y1 = j - fraction;
  z1 = k - fraction;
  x2 = i + fraction;
  y2 = j + fraction;
  z2 = k + fraction;

  this->LatticeToWorld(x1, y1, z1);
  this->LatticeToWorld(x2, y2, z2);

  if (x1 > x2) swap(x1, x2);
  if (y1 > y2) swap(y1, y2);
  if (z1 > z2) swap(z1, z2);
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
template <class CPImage>
void Evaluate(const CPImage *coeff, double &dx, double &dy, double &dz, int i, int j)
{
  typedef irtkBSplineFreeFormTransformation3D::Kernel Kernel;

  const double *w = Kernel::LatticeWeights;

  typename CPImage::VoxelType d;
  int                         jb;

  --i, --j;
  for (int b = 0; b < 3; ++b) {
    jb = j + b;
    for (int a = 0; a < 3; ++a) {
      d += (w[a] * w[b]) * coeff->Get(i + a, jb);
    }
  }

  dx = d._x, dy = d._y, dz = d._z;
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformation3D
::Evaluate(double &dx, double &dy, double &dz, int i, int j) const
{
  if (_FFD.IsInside(i, j)) ::Evaluate(&_CPImage, dx, dy, dz, i, j);
  else                     ::Evaluate( _CPValue, dx, dy, dz, i, j);
}

// -----------------------------------------------------------------------------
template <class CPImage>
void Evaluate(const CPImage *coeff, double &dx, double &dy, double &dz, int i, int j, int k)
{
  typedef irtkBSplineFreeFormTransformation3D::Kernel Kernel;

  const double *w = Kernel::LatticeWeights;

  typename CPImage::VoxelType d;
  int                         jb, kc;

  --i, --j, --k;
  for (int c = 0; c < 3; ++c) {
    kc = k + c;
    for (int b = 0; b < 3; ++b) {
      jb = j + b;
      for (int a = 0; a < 3; ++a) {
        d += (w[a] * w[b] * w[c]) * coeff->Get(i + a, jb, kc);
      }
    }
  }

  dx = d._x, dy = d._y, dz = d._z;
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformation3D
::Evaluate(double &dx, double &dy, double &dz, int i, int j, int k) const
{
  if (_FFD.IsInside(i, j, k)) ::Evaluate(&_CPImage, dx, dy, dz, i, j, k);
  else                        ::Evaluate( _CPValue, dx, dy, dz, i, j, k);
}

// -----------------------------------------------------------------------------
template <class CPImage>
void EvaluateJacobian(const CPImage *coeff, irtkMatrix &jac, int i, int j)
{
  typedef irtkBSplineFreeFormTransformation3D::Kernel Kernel;

  const double *w[2] = {
    Kernel::LatticeWeights,
    Kernel::LatticeWeights_I
  };

  typename CPImage::VoxelType dx, dy;
  int                         ia, jb;

  --i, --j;
  for (int b = 0; b < 3; ++b) {
    jb = j + b;
    for (int a = 0; a < 3; ++a) {
      ia = i + a;
      dx += (w[1][a] * w[0][b]) * coeff->Get(ia, jb);
      dy += (w[0][a] * w[1][b]) * coeff->Get(ia, jb);
    }
  }

  jac.Initialize(3, 3);
  jac(0, 0) = dx._x; jac(0, 1) = dy._x;
  jac(1, 0) = dx._y; jac(1, 1) = dy._y;
  jac(2, 0) = dx._z; jac(2, 1) = dy._z;
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformation3D
::EvaluateJacobian(irtkMatrix &jac, int i, int j) const
{
  if (_FFD.IsInside(i, j)) ::EvaluateJacobian(&_CPImage, jac, i, j);
  else                     ::EvaluateJacobian( _CPValue, jac, i, j);
}

// -----------------------------------------------------------------------------
template <class CPImage>
void EvaluateJacobian(const CPImage *coeff, irtkMatrix &jac, int i, int j, int k)
{
  typedef irtkBSplineFreeFormTransformation3D::Kernel Kernel;

  const double *w[2] = {
    Kernel::LatticeWeights,
    Kernel::LatticeWeights_I
  };

  typename CPImage::VoxelType dx, dy, dz;
  int                         ia, jb, kc;

  --i, --j, --k;
  for (int c = 0; c < 3; ++c) {
    kc = k + c;
    for (int b = 0; b < 3; ++b) {
      jb = j + b;
      for (int a = 0; a < 3; ++a) {
        ia = i + a;
        dx += (w[1][a] * w[0][b] * w[0][c]) * coeff->Get(ia, jb, kc);
        dy += (w[0][a] * w[1][b] * w[0][c]) * coeff->Get(ia, jb, kc);
        dz += (w[0][a] * w[0][b] * w[1][c]) * coeff->Get(ia, jb, kc);
      }
    }
  }

  jac.Initialize(3, 3);
  jac(0, 0) = dx._x; jac(0, 1) = dy._x; jac(0, 2) = dz._x;
  jac(1, 0) = dx._y; jac(1, 1) = dy._y; jac(1, 2) = dz._y;
  jac(2, 0) = dx._z; jac(2, 1) = dy._z; jac(2, 2) = dz._z;
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformation3D
::EvaluateJacobian(irtkMatrix &jac, int i, int j, int k) const
{
  if (_FFD.IsInside(i, j, k)) ::EvaluateJacobian(&_CPImage, jac, i, j, k);
  else                        ::EvaluateJacobian( _CPValue, jac, i, j, k);
}

// -----------------------------------------------------------------------------
template <class CPImage>
void EvaluateJacobian(const CPImage *coeff, irtkMatrix &jac, double x, double y)
{
  typedef irtkBSplineFreeFormTransformation3D::Kernel Kernel;

  int i = static_cast<int>(floor(x));
  int j = static_cast<int>(floor(y));

  const int A = Kernel::VariableToIndex(x - i);
  const int B = Kernel::VariableToIndex(y - j);

  typename CPImage::VoxelType dx, dy;
  double                      wx[2], wy[2];
  int                         ia, jb;

  --i, --j;
  for (int b = 0; b < 4; ++b) {
    jb = j + b;
    wy[0] = Kernel::LookupTable  [B][b];
    wy[1] = Kernel::LookupTable_I[B][b];
    for (int a = 0; a < 4; ++a) {
      ia = i + a;
      wx[0] = Kernel::LookupTable  [A][a];
      wx[1] = Kernel::LookupTable_I[A][a];
      dx += (wx[1] * wy[0]) * coeff->Get(ia, jb);
      dy += (wx[0] * wy[1]) * coeff->Get(ia, jb);
    }
  }

  jac.Initialize(3, 3);
  jac(0, 0) = dx._x; jac(0, 1) = dy._x;
  jac(1, 0) = dx._y; jac(1, 1) = dy._y;
  jac(2, 0) = dx._z; jac(2, 1) = dy._z;
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformation3D
::EvaluateJacobian(irtkMatrix &jac, double x, double y) const
{
  if (_FFD.IsInside(x, y)) ::EvaluateJacobian(&_CPImage, jac, x, y);
  else                     ::EvaluateJacobian( _CPValue, jac, x, y);
}

// -----------------------------------------------------------------------------
template <class CPImage>
void EvaluateJacobian(const CPImage *coeff, irtkMatrix &jac, double x, double y, double z)
{
  typedef irtkBSplineFreeFormTransformation3D::Kernel Kernel;

  int i = static_cast<int>(floor(x));
  int j = static_cast<int>(floor(y));
  int k = static_cast<int>(floor(z));

  const int A = Kernel::VariableToIndex(x - i);
  const int B = Kernel::VariableToIndex(y - j);
  const int C = Kernel::VariableToIndex(z - k);

  typename CPImage::VoxelType dx, dy, dz;
  double                      wx[2], wy[2], wz[2];
  int                         ia, jb, kc;

  --i, --j, --k;
  for (int c = 0; c < 4; ++c) {
    kc = k + c;
    wz[0] = Kernel::LookupTable  [C][c];
    wz[1] = Kernel::LookupTable_I[C][c];
    for (int b = 0; b < 4; ++b) {
      jb = j + b;
      wy[0] = Kernel::LookupTable  [B][b];
      wy[1] = Kernel::LookupTable_I[B][b];
      for (int a = 0; a < 4; ++a) {
        ia = i + a;
        wx[0] = Kernel::LookupTable  [A][a];
        wx[1] = Kernel::LookupTable_I[A][a];
        dx += (wx[1] * wy[0] * wz[0]) * coeff->Get(ia, jb, kc);
        dy += (wx[0] * wy[1] * wz[0]) * coeff->Get(ia, jb, kc);
        dz += (wx[0] * wy[0] * wz[1]) * coeff->Get(ia, jb, kc);
      }
    }
  }

  jac.Initialize(3, 3);
  jac(0, 0) = dx._x; jac(0, 1) = dy._x; jac(0, 2) = dz._x;
  jac(1, 0) = dx._y; jac(1, 1) = dy._y; jac(1, 2) = dz._y;
  jac(2, 0) = dx._z; jac(2, 1) = dy._z; jac(2, 2) = dz._z;
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformation3D
::EvaluateJacobian(irtkMatrix &jac, double x, double y, double z) const
{
  if (_FFD.IsInside(x, y, z)) ::EvaluateJacobian(&_CPImage, jac, x, y, z);
  else                        ::EvaluateJacobian( _CPValue, jac, x, y, z);
}

// -----------------------------------------------------------------------------
template <class CPImage>
void EvaluateHessian(const CPImage *coeff, irtkMatrix hessian[3], int i, int j)
{
  typedef irtkBSplineFreeFormTransformation3D::Kernel Kernel;

  const double *w[3] = {
    Kernel::LatticeWeights,
    Kernel::LatticeWeights_I,
    Kernel::LatticeWeights_II
  };

  typename CPImage::VoxelType dxx, dxy, dyy;
  int                         ia, jb;

  --i, --j;
  for (int b = 0; b < 3; ++b) {
    jb = j + b;
    for (int a = 0; a < 3; ++a) {
      ia = i + a;
      dxx += (w[2][a] * w[0][b]) * coeff->Get(ia, jb);
      dxy += (w[1][a] * w[1][b]) * coeff->Get(ia, jb);
      dyy += (w[0][a] * w[2][b]) * coeff->Get(ia, jb);
    }
  }

  irtkMatrix &hx = hessian[0];
  hx.Initialize(3, 3);
  hx(0, 0) = dxx._x; hx(0, 1) = dxy._x;
  hx(1, 0) = dxy._x; hx(1, 1) = dyy._x;

  irtkMatrix &hy = hessian[1];
  hy.Initialize(3, 3);
  hy(0, 0) = dxx._y; hy(0, 1) = dxy._y;
  hy(1, 0) = dxy._y; hy(1, 1) = dyy._y;

  irtkMatrix &hz = hessian[2];
  hz.Initialize(3, 3);
  hz(0, 0) = dxx._z; hz(0, 1) = dxy._z;
  hz(1, 0) = dxy._z; hz(1, 1) = dyy._z;
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformation3D
::EvaluateHessian(irtkMatrix hessian[3], int i, int j) const
{
  if (_FFD.IsInside(i, j)) ::EvaluateHessian(&_CPImage, hessian, i, j);
  else                     ::EvaluateHessian( _CPValue, hessian, i, j);
}

// -----------------------------------------------------------------------------
template <class CPImage>
void EvaluateHessian(const CPImage *coeff, irtkMatrix hessian[3], int i, int j, int k)
{
  typedef irtkBSplineFreeFormTransformation3D::Kernel Kernel;

  const double *w[3] = {
    Kernel::LatticeWeights,
    Kernel::LatticeWeights_I,
    Kernel::LatticeWeights_II
  };

  typename CPImage::VoxelType dxx, dxy, dxz, dyy, dyz, dzz;
  int                         ia, jb, kc;

  --i, --j, --k;
  for (int c = 0; c < 3; ++c) {
    kc = k + c;
    for (int b = 0; b < 3; ++b) {
      jb = j + b;
      for (int a = 0; a < 3; ++a) {
        ia = i + a;
        dxx += (w[2][a] * w[0][b] * w[0][c]) * coeff->Get(ia, jb, kc);
        dxy += (w[1][a] * w[1][b] * w[0][c]) * coeff->Get(ia, jb, kc);
        dxz += (w[1][a] * w[0][b] * w[1][c]) * coeff->Get(ia, jb, kc);
        dyy += (w[0][a] * w[2][b] * w[0][c]) * coeff->Get(ia, jb, kc);
        dyz += (w[0][a] * w[1][b] * w[1][c]) * coeff->Get(ia, jb, kc);
        dzz += (w[0][a] * w[0][b] * w[2][c]) * coeff->Get(ia, jb, kc);
      }
    }
  }

  irtkMatrix &hx = hessian[0];
  hx.Initialize(3, 3);
  hx(0, 0) = dxx._x; hx(0, 1) = dxy._x; hx(0, 2) = dxz._x;
  hx(1, 0) = dxy._x; hx(1, 1) = dyy._x; hx(1, 2) = dyz._x;
  hx(2, 0) = dxz._x; hx(2, 1) = dyz._x; hx(2, 2) = dzz._x;

  irtkMatrix &hy = hessian[1];
  hy.Initialize(3, 3);
  hy(0, 0) = dxx._y; hy(0, 1) = dxy._y; hy(0, 2) = dxz._y;
  hy(1, 0) = dxy._y; hy(1, 1) = dyy._y; hy(1, 2) = dyz._y;
  hy(2, 0) = dxz._y; hy(2, 1) = dyz._y; hy(2, 2) = dzz._y;

  irtkMatrix &hz = hessian[2];
  hz.Initialize(3, 3);
  hz(0, 0) = dxx._z; hz(0, 1) = dxy._z; hz(0, 2) = dxz._z;
  hz(1, 0) = dxy._z; hz(1, 1) = dyy._z; hz(1, 2) = dyz._z;
  hz(2, 0) = dxz._z; hz(2, 1) = dyz._z; hz(2, 2) = dzz._z;
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformation3D
::EvaluateHessian(irtkMatrix hessian[3], int i, int j, int k) const
{
  if (_FFD.IsInside(i, j, k)) ::EvaluateHessian(&_CPImage, hessian, i, j, k);
  else                        ::EvaluateHessian( _CPValue, hessian, i, j, k);
}

// -----------------------------------------------------------------------------
template <class CPImage>
void EvaluateHessian(const CPImage *coeff, irtkMatrix hessian[3], double x, double y)
{
  typedef irtkBSplineFreeFormTransformation3D::Kernel Kernel;

  int i = static_cast<int>(floor(x));
  int j = static_cast<int>(floor(y));

  const int A = Kernel::VariableToIndex(x - i);
  const int B = Kernel::VariableToIndex(y - j);

  typename CPImage::VoxelType dxx, dxy, dyy;
  double                      wx[3], wy[3];
  int                         ia, jb;

  --i, --j;
  for (int b = 0; b < 4; ++b) {
    jb = j + b;
    wy[0] = Kernel::LookupTable   [B][b];
    wy[1] = Kernel::LookupTable_I [B][b];
    wy[2] = Kernel::LookupTable_II[B][b];
    for (int a = 0; a < 4; ++a) {
      ia = i + a;
      wx[0] = Kernel::LookupTable   [A][a];
      wx[1] = Kernel::LookupTable_I [A][a];
      wx[2] = Kernel::LookupTable_II[A][a];
      dxx += (wx[2] * wy[0]) * coeff->Get(ia, jb);
      dxy += (wx[1] * wy[1]) * coeff->Get(ia, jb);
      dyy += (wx[0] * wy[2]) * coeff->Get(ia, jb);
    }
  }

  irtkMatrix &hx = hessian[0];
  hx.Initialize(3, 3);
  hx(0, 0) = dxx._x; hx(0, 1) = dxy._x;
  hx(1, 0) = dxy._x; hx(1, 1) = dyy._x;

  irtkMatrix &hy = hessian[1];
  hy.Initialize(3, 3);
  hy(0, 0) = dxx._y; hy(0, 1) = dxy._y;
  hy(1, 0) = dxy._y; hy(1, 1) = dyy._y;

  irtkMatrix &hz = hessian[2];
  hz.Initialize(3, 3);
  hz(0, 0) = dxx._z; hz(0, 1) = dxy._z;
  hz(1, 0) = dxy._z; hz(1, 1) = dyy._z;
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformation3D
::EvaluateHessian(irtkMatrix hessian[3], double x, double y) const
{
  if (_FFD.IsInside(x, y)) ::EvaluateHessian(&_CPImage, hessian, x, y);
  else                     ::EvaluateHessian( _CPValue, hessian, x, y);
}

// -----------------------------------------------------------------------------
template <class CPImage>
void EvaluateHessian(const CPImage *coeff, irtkMatrix hessian[3], double x, double y, double z)
{
  typedef irtkBSplineFreeFormTransformation3D::Kernel Kernel;

  int i = static_cast<int>(floor(x));
  int j = static_cast<int>(floor(y));
  int k = static_cast<int>(floor(z));

  const int A = Kernel::VariableToIndex(x - i);
  const int B = Kernel::VariableToIndex(y - j);
  const int C = Kernel::VariableToIndex(z - k);

  typename CPImage::VoxelType dxx, dxy, dxz, dyy, dyz, dzz;
  double                      wx[3], wy[3], wz[3];
  int                         ia, jb, kc;

  --i, --j, --k;
  for (int c = 0; c < 4; ++c) {
    kc = k + c;
    wz[0] = Kernel::LookupTable   [C][c];
    wz[1] = Kernel::LookupTable_I [C][c];
    wz[2] = Kernel::LookupTable_II[C][c];
    for (int b = 0; b < 4; ++b) {
      jb = j + b;
      wy[0] = Kernel::LookupTable   [B][b];
      wy[1] = Kernel::LookupTable_I [B][b];
      wy[2] = Kernel::LookupTable_II[B][b];
      for (int a = 0; a < 4; ++a) {
        ia = i + a;
        wx[0] = Kernel::LookupTable   [A][a];
        wx[1] = Kernel::LookupTable_I [A][a];
        wx[2] = Kernel::LookupTable_II[A][a];
        dxx += (wx[2] * wy[0] * wz[0]) * coeff->Get(ia, jb, kc);
        dxy += (wx[1] * wy[1] * wz[0]) * coeff->Get(ia, jb, kc);
        dxz += (wx[1] * wy[0] * wz[1]) * coeff->Get(ia, jb, kc);
        dyy += (wx[0] * wy[2] * wz[0]) * coeff->Get(ia, jb, kc);
        dyz += (wx[0] * wy[1] * wz[1]) * coeff->Get(ia, jb, kc);
        dzz += (wx[0] * wy[0] * wz[2]) * coeff->Get(ia, jb, kc);
      }
    }
  }

  irtkMatrix &hx = hessian[0];
  hx.Initialize(3, 3);
  hx(0, 0) = dxx._x; hx(0, 1) = dxy._x; hx(0, 2) = dxz._x;
  hx(1, 0) = dxy._x; hx(1, 1) = dyy._x; hx(1, 2) = dyz._x;
  hx(2, 0) = dxz._x; hx(2, 1) = dyz._x; hx(2, 2) = dzz._x;

  irtkMatrix &hy = hessian[1];
  hy.Initialize(3, 3);
  hy(0, 0) = dxx._y; hy(0, 1) = dxy._y; hy(0, 2) = dxz._y;
  hy(1, 0) = dxy._y; hy(1, 1) = dyy._y; hy(1, 2) = dyz._y;
  hy(2, 0) = dxz._y; hy(2, 1) = dyz._y; hy(2, 2) = dzz._y;

  irtkMatrix &hz = hessian[2];
  hz.Initialize(3, 3);
  hz(0, 0) = dxx._z; hz(0, 1) = dxy._z; hz(0, 2) = dxz._z;
  hz(1, 0) = dxy._z; hz(1, 1) = dyy._z; hz(1, 2) = dyz._z;
  hz(2, 0) = dxz._z; hz(2, 1) = dyz._z; hz(2, 2) = dzz._z;
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformation3D
::EvaluateHessian(irtkMatrix hessian[3], double x, double y, double z) const
{
  if (_FFD.IsInside(x, y, z)) ::EvaluateHessian(&_CPImage, hessian, x, y, z);
  else                        ::EvaluateHessian( _CPValue, hessian, x, y, z);
}

// -----------------------------------------------------------------------------
template <class CPImage>
void EvaluateLaplacian(const CPImage *coeff, double laplacian[3], int i, int j, int k)
{
  typedef irtkBSplineFreeFormTransformation3D::Kernel Kernel;

  const double *w[3] = {
    Kernel::LatticeWeights,
    Kernel::LatticeWeights_I,
    Kernel::LatticeWeights_II
  };

  typename CPImage::VoxelType v = .0;
  int                         ia, jb, kc;

  --i, --j, --k;
  for (int c = 0; c < 3; ++c) {
    kc = k + c;
    for (int b = 0; b < 3; ++b) {
      jb = j + b;
      for (int a = 0; a < 3; ++a) {
        ia = i + a;
        v += (w[2][a] * w[0][b] * w[0][c]) * coeff->Get(ia, jb, kc);
        v += (w[0][a] * w[2][b] * w[0][c]) * coeff->Get(ia, jb, kc);
        v += (w[0][a] * w[0][b] * w[2][c]) * coeff->Get(ia, jb, kc);
      }
    }
  }

  laplacian[0] = v._x;
  laplacian[1] = v._y;
  laplacian[2] = v._z;
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformation3D
::EvaluateLaplacian(double laplacian[3], int i, int j, int k) const
{
  if (_FFD.IsInside(i, j, k)) ::EvaluateLaplacian(&_CPImage, laplacian, i, j, k);
  else                        ::EvaluateLaplacian( _CPValue, laplacian, i, j, k);
}

// -----------------------------------------------------------------------------
template <class CPImage>
void EvaluateLaplacian(const CPImage *coeff, double &x, double &y, double &z)
{
  typedef irtkBSplineFreeFormTransformation3D::Kernel Kernel;

  int i = static_cast<int>(floor(x));
  int j = static_cast<int>(floor(y));
  int k = static_cast<int>(floor(z));

  const int A = Kernel::VariableToIndex(x - i);
  const int B = Kernel::VariableToIndex(y - j);
  const int C = Kernel::VariableToIndex(z - k);

  typename CPImage::VoxelType v;
  double                      wx[2], wy[2], wz[2];
  int                         ia, jb, kc;

  --i, --j, --k;
  for (int c = 0; c < 4; ++c) {
    kc = k + c;
    wz[0] = Kernel::LookupTable   [C][c];
    wz[1] = Kernel::LookupTable_II[C][c];
    for (int b = 0; b < 4; ++b) {
      jb = j + b;
      wy[0] = Kernel::LookupTable   [B][b];
      wy[1] = Kernel::LookupTable_II[B][b];
      for (int a = 0; a < 4; ++a) {
        ia = i + a;
        wx[0] = Kernel::LookupTable   [A][a];
        wx[1] = Kernel::LookupTable_II[A][a];
        v += (wx[1] * wy[0] * wz[0]) * coeff->Get(ia, jb, kc);
        v += (wx[0] * wy[1] * wz[0]) * coeff->Get(ia, jb, kc);
        v += (wx[0] * wy[0] * wz[1]) * coeff->Get(ia, jb, kc);
      }
    }
  }

  x = v._x, y = v._y, z = v._z;
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformation3D
::EvaluateLaplacian(double laplacian[3], double x, double y, double z) const
{
  if (_FFD.IsInside(x, y, z)) ::EvaluateLaplacian(&_CPImage, x, y, z);
  else                        ::EvaluateLaplacian( _CPValue, x, y, z);
  laplacian[0] = x, laplacian[1] = y, laplacian[2] = z;
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformation3D
::EvaluateLaplacian(double &x, double &y, double &z) const
{
  if (_FFD.IsInside(x, y, z)) ::EvaluateLaplacian(&_CPImage, x, y, z);
  else                        ::EvaluateLaplacian( _CPValue, x, y, z);
}

// =============================================================================
// Point transformation
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkBSplineFreeFormTransformation3D::CanModifyDisplacement(int) const
{
  return true;
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformation3D
::DisplacementAfterDOFChange(int dof, double dv, irtkGenericImage<double> &disp,
           /* unused --> */ double, double, const irtkWorldCoordsImage *) const
{
  // Corresponding control point index and dimension of displacement field
  const int cp  = this->DOFToIndex    (dof);
  const int dim = this->DOFToDimension(dof);

  // Bounding box of control point in image coordinates
  double                x1, y1, z1, x2, y2, z2;
  this->BoundingBox(cp, x1, y1, z1, x2, y2, z2);
  disp.WorldToImage(x1, y1, z1);
  disp.WorldToImage(x2, y2, z2);

  if (x1 > x2) swap(x1, x2);
  if (y1 > y2) swap(y1, y2);
  if (z1 > z2) swap(z1, z2);

  // Calculate increment for voxel offset to kernel function weight lookup table
  const double dx = (Kernel::LookupTableSize - 1) / (x2 - x1);
  const double dy = (Kernel::LookupTableSize - 1) / (y2 - y1);
  const double dz = (Kernel::LookupTableSize - 1) / (z2 - z1);

  // Bounding box of control point in voxel indices
  int                          i1, j1, k1, i2, j2, k2;
  this->BoundingBox(&disp, cp, i1, j1, k1, i2, j2, k2, 1.0 / _SpeedupFactor);

  // Get pointer to displacement of lower-left voxel in bounding box and
  // increments for update of pointer while looping over the voxels
  const int s1 = 1;
  const int s2 = (disp.X() - (i2 - i1 + 1));
  const int s3 = (disp.Y() - (j2 - j1 + 1)) * disp.X();
  double    *d = disp.Data(i1, j1, k1, dim);

  // Loop over voxels in bounding box of control point and add additional
  // displacement along dim induced by change of control point parameter
  double di, dj, dk; // displacement change = (B_i * (B_j * (B_k * dv)))
  for (int k = k1; k <= k2; ++k, d += s3) {
    dk = Kernel::WeightLookupTable[iround((k - z1) * dz)] * dv;
    for (int j = j1; j <= j2; ++j, d += s2) {
      dj = Kernel::WeightLookupTable[iround((j - y1) * dy)] * dk;
      for (int i = i1; i <= i2; ++i, d += s1) {
        di = Kernel::WeightLookupTable[iround((i - x1) * dx)] * dj;
        (*d) += di;
      }
    }
  }
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformation3D
::JacobianDetDerivative(irtkMatrix *detdev, int x, int y, int z) const
{
  // Values of the B-spline basis functions and its 1st derivative
  const double *w[2] = {
    Kernel::LatticeWeights,
    Kernel::LatticeWeights_I
  };

  // 1D B-splines
  double x_0 = 0, y_0 = 0, z_0 = 0;
  double x_1 = 0, y_1 = 0, z_1 = 0;

  if (-1 <= x && x <= 1) {
    x_0 = w[0][x + 1];
    x_1 = w[1][x + 1];
  }
  if (-1 <= y && y <= 1) {
    y_0 = w[0][y + 1];
    y_1 = w[1][y + 1];
  }
  if (-1 <= z && z <= 1) {
    z_0 = w[0][z + 1];
    z_1 = w[1][z + 1];
  }

  // B-spline tensor product
  double b_i = x_1 * y_0 * z_0;
  double b_j = x_0 * y_1 * z_0;
  double b_k = x_0 * y_0 * z_1;

  // Return
  for (int i = 0; i < 3; ++i) {
    detdev[i].Initialize(3, 3);
    // w.r.t lattice coordinates
    detdev[i](i, 0) = b_i;
    detdev[i](i, 1) = b_j;
    detdev[i](i, 2) = b_k;
    // w.r.t world coordinates
    JacobianToWorld(detdev[i]);
  }
}

// =============================================================================
// Properties
// =============================================================================

// -----------------------------------------------------------------------------
int irtkBSplineFreeFormTransformation3D::KernelSize() const
{
  return 4;
}

// -----------------------------------------------------------------------------
double irtkBSplineFreeFormTransformation3D
::BendingEnergy(double x, double y, double z, double, double, bool wrt_world) const
{
  // Convert to lattice coordinates
  this->WorldToLattice(x, y, z);
  // Calculate 2nd order derivatives
  irtkMatrix hessian[3];
  if (_z == 1) EvaluateHessian(hessian, x, y);
  else         EvaluateHessian(hessian, x, y, z);
  // Convert derivatives to world coordinates
  if (wrt_world) HessianToWorld(hessian);
  // Calculate bending energy
  return Bending3D(hessian);
}

// -----------------------------------------------------------------------------
double irtkBSplineFreeFormTransformation3D::BendingEnergy(bool incl_passive, bool wrt_world) const
{
  irtkMatrix hessian[3];
  double     bending = .0;
  int        nactive = 0;

  for (int k = 0; k < _z; ++k) {
    for (int j = 0; j < _y; ++j) {
      for (int i = 0; i < _x; ++i) {
        if (incl_passive || IsActive(i, j, k)) {
          if (_z == 1) EvaluateHessian(hessian, i, j);
          else         EvaluateHessian(hessian, i, j, k);
          if (wrt_world) HessianToWorld(hessian);
          bending += Bending3D(hessian);
          ++nactive;
        }
      }
    }
  }

  if (nactive) bending /= nactive;
  return bending;
}

// -----------------------------------------------------------------------------
double irtkBSplineFreeFormTransformation3D
::BendingEnergy(const irtkImageAttributes &attr, double, bool wrt_world) const
{
  const int nvox = attr.NumberOfSpatialPoints();
  if (nvox == 0) return .0;

  double     bending = .0;
  double     x, y, z;
  irtkMatrix hessian[3];

  for (int k = 0; k < attr._z; ++k) {
    for (int j = 0; j < attr._y; ++j) {
      for (int i = 0; i < attr._x; ++i) {
        x = i, y = j, z = k;
        attr .LatticeToWorld(x, y, z);
        this->WorldToLattice(x, y, z);
        if (_z == 1) EvaluateHessian(hessian, x, y);
        else         EvaluateHessian(hessian, x, y, z);
        if (wrt_world) HessianToWorld(hessian);
        bending += Bending3D(hessian);
      }
    }
  }

  return bending / nvox;
}


namespace irtkBSplineFreeFormTransformation3DUtils {

// -----------------------------------------------------------------------------
/// Voxel function which evaluates the 2nd order derivatives of one component
/// of a B-spline FFD at control point lattice coordinates
struct Evaluate2ndOrderBSplineFFDDerivatives : public irtkVoxelFunction
{
  typedef irtkBSplineFreeFormTransformation3D::CPExtrapolator Extrapolator;
  typedef irtkBSplineFreeFormTransformation3D::Vector         Vector;
  typedef irtkBSplineFreeFormTransformation3D::Kernel         Kernel;

  const irtkBSplineFreeFormTransformation3D *_Transformation; ///< B-spline free-form deformation
  const Extrapolator                        *_CPValue;        ///< Coefficients of B-spline FFD
  bool                                       _WrtWorld;       ///< Whether to compute derivatives
                                                              ///< w.r.t world or lattice coordinates

  /// Derivatives of 2nd order derivatives w.r.t control point parameters
  static double LookupTable_2D[9][3];

  /// Derivatives of 2nd order derivatives w.r.t control point parameters
  static double LookupTable_3D[27][6];

  /// Initialize static lookup tables of 2nd order derivatives w.r.t control
  /// point parameters evaluated for each voxel in the 2D kernel support region
  static void InitializeLookupTable2D();

  /// Initialize static lookup tables of 2nd order derivatives w.r.t control
  /// point parameters evaluated for each voxel in the 3D kernel support region
  static void InitializeLookupTable3D();

  /// Constructor
  Evaluate2ndOrderBSplineFFDDerivatives(const irtkBSplineFreeFormTransformation3D *ffd,
                                        bool wrt_world = false)
  :
    _Transformation(ffd), _CPValue(ffd->Extrapolator()), _WrtWorld(wrt_world)
  {}

  /// Evaluate 2nd order derivatives of 2D FFD at given lattice coordinate
  void operator()(int i, int j, int k, int, Vector *dxx, Vector *dxy, Vector *dyy)
  {
    // Note: Derivatives are evaluated on a lattice that has an
    //       additional boundary margin of one voxel. Therefore,
    //       CP indices I and J are shifted by an offset of -1.
    int n = 0;
    for (int J = j-2; J <= j; ++J) {
      for (int I = i-2; I <= i; ++I, ++n) {
        *dxx += _CPValue->Get(I, J) * LookupTable_2D[n][0];
        *dxy += _CPValue->Get(I, J) * LookupTable_2D[n][1];
        *dyy += _CPValue->Get(I, J) * LookupTable_2D[n][2];
      }
    }
    // Apply product and chain rule to convert derivatives to ones w.r.t the world
    if (_WrtWorld) {
      _Transformation->HessianToWorld(dxx->_x, dxy->_x, dyy->_x);
      _Transformation->HessianToWorld(dxx->_y, dxy->_y, dyy->_y);
      _Transformation->HessianToWorld(dxx->_z, dxy->_z, dyy->_z);
    }
  }

  /// Evaluate 2nd order derivatives of 3D FFD at given lattice coordinate
  void operator()(int i, int j, int k, int,
                  Vector *dxx, Vector *dxy, Vector *dxz, Vector *dyy, Vector *dyz, Vector *dzz)
  {
    // Note: Derivatives are evaluated on a lattice that has an
    //       additional boundary margin of one voxel. Therefore,
    //       CP indices I, J and K are shifted by an offset of -1.
    int n = 0;
    for (int K = k-2; K <= k; ++K) {
      for (int J = j-2; J <= j; ++J) {
        for (int I = i-2; I <= i; ++I, ++n) {
          *dxx += _CPValue->Get(I, J, K) * LookupTable_3D[n][0];
          *dxy += _CPValue->Get(I, J, K) * LookupTable_3D[n][1];
          *dxz += _CPValue->Get(I, J, K) * LookupTable_3D[n][2];
          *dyy += _CPValue->Get(I, J, K) * LookupTable_3D[n][3];
          *dyz += _CPValue->Get(I, J, K) * LookupTable_3D[n][4];
          *dzz += _CPValue->Get(I, J, K) * LookupTable_3D[n][5];
        }
      }
    }
    // Apply product and chain rule to convert derivatives to ones w.r.t the world
    if (_WrtWorld) {
      _Transformation->HessianToWorld(dxx->_x, dxy->_x, dxz->_x, dyy->_x, dyz->_x, dzz->_x);
      _Transformation->HessianToWorld(dxx->_y, dxy->_y, dxz->_y, dyy->_y, dyz->_y, dzz->_y);
      _Transformation->HessianToWorld(dxx->_z, dxy->_z, dxz->_z, dyy->_z, dyz->_z, dzz->_z);
    }
  }
};

// -----------------------------------------------------------------------------
double Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_2D[9][3] = {{.0}};
void Evaluate2ndOrderBSplineFFDDerivatives::InitializeLookupTable2D()
{
  static bool initialized = false;
  if (initialized) return;

  const double *w[3] = {
    Kernel::LatticeWeights,
    Kernel::LatticeWeights_I,
    Kernel::LatticeWeights_II
  };

  int n = 0;
  for (int b = 0; b < 3; ++b)
  for (int a = 0; a < 3; ++a, ++n) {
    LookupTable_2D[n][0] = w[2][a] * w[0][b];
    LookupTable_2D[n][1] = w[1][a] * w[1][b];
    LookupTable_2D[n][2] = w[0][a] * w[2][b];
  }

  initialized = true;
}

// -----------------------------------------------------------------------------
double Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[27][6] = {{.0}};
void Evaluate2ndOrderBSplineFFDDerivatives::InitializeLookupTable3D()
{
  static bool initialized = false;
  if (initialized) return;

  const double *w[3] = {
    Kernel::LatticeWeights,
    Kernel::LatticeWeights_I,
    Kernel::LatticeWeights_II
  };

  int n = 0;
  for (int c = 0; c < 3; ++c)
  for (int b = 0; b < 3; ++b)
  for (int a = 0; a < 3; ++a, ++n) {
    LookupTable_3D[n][0] = w[2][a] * w[0][b] * w[0][c];
    LookupTable_3D[n][1] = w[1][a] * w[1][b] * w[0][c];
    LookupTable_3D[n][2] = w[1][a] * w[0][b] * w[1][c];
    LookupTable_3D[n][3] = w[0][a] * w[2][b] * w[0][c];
    LookupTable_3D[n][4] = w[0][a] * w[1][b] * w[1][c];
    LookupTable_3D[n][5] = w[0][a] * w[0][b] * w[2][c];
  }

  initialized = true;
}

} // namespace irtkBSplineFreeFormTransformation3DUtils
using namespace irtkBSplineFreeFormTransformation3DUtils;


// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformation3D::BendingEnergyGradient(double *gradient, double weight, bool incl_passive, bool wrt_world) const
{
  Vector sum;
  int    m, n;

  IRTK_START_TIMING();

  // Pre-multiply weight by derivative of square function (2) and normalization factor
  const int ncps = this->NumberOfActiveCPs();
  if (ncps == 0) return;
  weight *= 2.0 / ncps;

  // ---------------------------------------------------------------------------
  // Bending energy of 2D FFD
  if (_z == 1) {

    // Initialize static lookup table
    Evaluate2ndOrderBSplineFFDDerivatives::InitializeLookupTable2D();

    // Add a layer of boundary voxels to avoid additional boundary conditions
    irtkImageAttributes attr = this->Attributes();
    attr._x += 2, attr._y += 2;

    // Compute 2nd order derivatives w.r.t control point lattice coordinates,
    // evaluated at each control point of the lattice
    irtkGenericImage<Vector> dxx(attr);
    irtkGenericImage<Vector> dxy(attr);
    irtkGenericImage<Vector> dyy(attr);

    Evaluate2ndOrderBSplineFFDDerivatives eval(this, wrt_world);
    ParallelForEachVoxel(attr, dxx, dxy, dyy, eval);

    // Compute 3rd order derivatives, twice w.r.t. lattice or world coordinate
    // and once w.r.t. transformation parameters of control point
    double w[9][4];
    if (wrt_world) {
      // Loop over support region (3x3) of a control point
      //
      // Note that the following terms are independent of the transformation
      // parameters and therefore can be pre-computed here as they are identical for
      // all control point positions at which the bending gradient is evaluated.
      for (n = 0; n < 9; ++n) {
        m = 0; // index of 2nd order derivative ( d^2 T(x[n]) / dx_i dx_j )
               // which is multiplied by weight w[n][m] according to the chain rule
        for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j, ++m) {
            w[n][m]  = Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_2D[n][0] * _matW2L(0, i) * _matW2L(0, j); // ( d^3 T(x[n]) / dudu dPhi ) * du/dx_i * du/dx_j
            w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_2D[n][1] * _matW2L(0, i) * _matW2L(1, j); // ( d^3 T(x[n]) / dudv dPhi ) * du/dx_i * dv/dx_j
            w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_2D[n][1] * _matW2L(1, i) * _matW2L(0, j); // ( d^3 T(x[n]) / dvdu dPhi ) * dv/dx_i * du/dx_j
            w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_2D[n][2] * _matW2L(1, i) * _matW2L(1, j); // ( d^3 T(x[n]) / dvdv dPhi ) * dv/dx_i * dv/dx_j
        }
        // Add weights of ( d^2 T(x[n]) / dx_i dx_j ) and ( d^2 T(x[n]) / dx_j dx_i )
        // for i != j as these derivatives are identical and re-order remaining weights.
        w[n][1] = w[n][1] + w[n][2];
        w[n][2] = w[n][3];
      }
    } else {
      for (n = 0; n < 9; ++n) {
        w[n][0] =       Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_2D[n][0];
        w[n][1] = 2.0 * Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_2D[n][1];
        w[n][2] =       Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_2D[n][2];
      }
    }

    // Compute derivative of bending energy w.r.t each control point
    int xdof, ydof, zdof;
    for (int cj = 0; cj < _y; ++cj)
    for (int ci = 0; ci < _x; ++ci) {
      if (incl_passive || IsActive(ci, cj)) {
        sum = .0;
        // Loop over support region (3x3) of control point
        //
        // Note: Derivatives were evaluated on a lattice that has an
        //       additional boundary margin of one voxel. Therefore,
        //       indices i and j are shifted by an offset of +1.
        n = 0;
        for (int j = cj; j <= cj+2; ++j)
        for (int i = ci; i <= ci+2; ++i, ++n) {
            sum += dxx(i, j) * w[n][0];
            sum += dxy(i, j) * w[n][1];
            sum += dyy(i, j) * w[n][2];
        }
        this->IndexToDOFs(this->LatticeToIndex(ci, cj), xdof, ydof, zdof);
        gradient[xdof] += weight * sum._x;
        gradient[ydof] += weight * sum._y;
        gradient[zdof] += weight * sum._z;
      }
    }

  // ---------------------------------------------------------------------------
  // Bending energy of 3D FFD
  } else {

    // Initialize static lookup table
    Evaluate2ndOrderBSplineFFDDerivatives::InitializeLookupTable3D();

    // Add a layer of boundary voxels to avoid additional boundary conditions
    irtkImageAttributes attr = this->Attributes();
    attr._x += 2, attr._y += 2, attr._z += 2;

    // Compute 2nd order derivatives w.r.t. lattice or world coordinates,
    // respectively, evaluated at each control point of the lattice
    irtkGenericImage<Vector> dxx(attr);
    irtkGenericImage<Vector> dxy(attr);
    irtkGenericImage<Vector> dxz(attr);
    irtkGenericImage<Vector> dyy(attr);
    irtkGenericImage<Vector> dyz(attr);
    irtkGenericImage<Vector> dzz(attr);

    Evaluate2ndOrderBSplineFFDDerivatives eval(this, wrt_world);
    ParallelForEachVoxel(attr, dxx, dxy, dxz, dyy, dyz, dzz, eval);

    // Compute 3rd order derivatives, twice w.r.t. lattice or world coordinate
    // and once w.r.t. transformation parameters of control point
    double w[27][9];
    if (wrt_world) {
      // Loop over support region (3x3x3) of a control point
      //
      // Note that the following terms are independent of the transformation
      // parameters and therefore can be pre-computed here as they are identical for
      // all control point positions at which the bending gradient is evaluated.
      for (n = 0; n < 27; ++n) {
        m = 0; // index of 2nd order derivative ( d^2 T(x[n]) / dx_i dx_j )
               // which is multiplied by weight w[n][m] according to the chain rule
        for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j, ++m) {
          w[n][m]  = Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][0] * _matW2L(0, i) * _matW2L(0, j); // ( d^3 T(x[n]) / dudu dPhi ) * du/dx_i * du/dx_j
          w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][1] * _matW2L(0, i) * _matW2L(1, j); // ( d^3 T(x[n]) / dudv dPhi ) * du/dx_i * dv/dx_j
          w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][2] * _matW2L(0, i) * _matW2L(2, j); // ( d^3 T(x[n]) / dudw dPhi ) * du/dx_i * dw/dx_j
          w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][1] * _matW2L(1, i) * _matW2L(0, j); // ( d^3 T(x[n]) / dvdu dPhi ) * dv/dx_i * du/dx_j
          w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][3] * _matW2L(1, i) * _matW2L(1, j); // ( d^3 T(x[n]) / dvdv dPhi ) * dv/dx_i * dv/dx_j
          w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][4] * _matW2L(1, i) * _matW2L(2, j); // ( d^3 T(x[n]) / dvdw dPhi ) * dv/dx_i * dw/dx_j
          w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][2] * _matW2L(2, i) * _matW2L(0, j); // ( d^3 T(x[n]) / dwdu dPhi ) * dw/dx_i * du/dx_j
          w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][4] * _matW2L(2, i) * _matW2L(1, j); // ( d^3 T(x[n]) / dwdv dPhi ) * dw/dx_i * dv/dx_j
          w[n][m] += Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][5] * _matW2L(2, i) * _matW2L(2, j); // ( d^3 T(x[n]) / dwdw dPhi ) * dw/dx_i * dw/dx_j
        }
        // Add weights of ( d^2 T(x[n]) / dx_i dx_j ) and ( d^2 T(x[n]) / dx_j dx_i )
        // for i != j as these derivatives are identical and re-order remaining weights.
        w[n][1] = w[n][1] + w[n][3];
        w[n][2] = w[n][2] + w[n][6];
        w[n][3] = w[n][4];
        w[n][4] = w[n][5] + w[n][7];
        w[n][5] = w[n][8];
      }
    } else {
      for (n = 0; n < 27; ++n) {
        w[n][0] =       Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][0];
        w[n][1] = 2.0 * Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][1];
        w[n][2] = 2.0 * Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][2];
        w[n][3] =       Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][3];
        w[n][4] = 2.0 * Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][4];
        w[n][5] =       Evaluate2ndOrderBSplineFFDDerivatives::LookupTable_3D[n][5];
      }
    }

    // Compute derivative of bending energy w.r.t each control point
    int xdof, ydof, zdof;
    for (int ck = 0; ck < _z; ++ck)
    for (int cj = 0; cj < _y; ++cj)
    for (int ci = 0; ci < _x; ++ci) {
      if (incl_passive || IsActive(ci, cj, ck)) {
        sum = .0;
        // Loop over support region (3x3x3) of control point
        //
        // Note: Derivatives were evaluated on a lattice that has an
        //       additional boundary margin of one voxel. Therefore,
        //       indices i, j and k are shifted by an offset of +1.
        n = 0;
        for (int k = ck; k <= ck+2; ++k)
        for (int j = cj; j <= cj+2; ++j)
        for (int i = ci; i <= ci+2; ++i, ++n) {
          sum += dxx(i, j, k) * w[n][0];
          sum += dxy(i, j, k) * w[n][1];
          sum += dxz(i, j, k) * w[n][2];
          sum += dyy(i, j, k) * w[n][3];
          sum += dyz(i, j, k) * w[n][4];
          sum += dzz(i, j, k) * w[n][5];
        }
        this->IndexToDOFs(this->LatticeToIndex(ci, cj, ck), xdof, ydof, zdof);
        gradient[xdof] += weight * sum._x;
        gradient[ydof] += weight * sum._y;
        gradient[zdof] += weight * sum._z;
      }
    }

  }
  IRTK_DEBUG_TIMING(2, "bending gradient computation");
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformation3D::Print(irtkIndent indent) const
{
  cout << indent << "3D B-spline FFD:" << endl;
  irtkFreeFormTransformation3D::Print(indent + 1);
}

// -----------------------------------------------------------------------------
bool irtkBSplineFreeFormTransformation3D::CanRead(irtkTransformationType format) const
{
  switch (format) {
    case IRTKTRANSFORMATION_BSPLINE_FFD_2D_v1:
    case IRTKTRANSFORMATION_BSPLINE_FFD_3D_v1:
    case IRTKTRANSFORMATION_BSPLINE_FFD_3D_v2:
    case IRTKTRANSFORMATION_BSPLINE_FFD_3D_v3:
    case IRTKTRANSFORMATION_BSPLINE_FFD_3D_v4:
      return true;
    default:
      return false;
  }
}

// =============================================================================
// Deprecated
// =============================================================================

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformation3D::FFD2D(double &x, double &y) const
{
  if (x < -2 || y < -2 || x >= _x || y >= _y) {
    x = y = .0;
    return;
  }

  int i = static_cast<int>(floor(x));
  int j = static_cast<int>(floor(y));

  const int A = Kernel::VariableToIndex(x - i);
  const int B = Kernel::VariableToIndex(y - j);

  --i, --j;

  CPValue v;
  int ia, jb;

  for (int b = 0; b <= 3; ++b) {
    jb = j + b;
    if (0 <= jb && jb < _y) {
      for (int a = 0; a <= 3; ++a) {
        ia = i + a;
        if (0 <= ia && ia < _x) {
          v += Kernel::LookupTable[A][a] *
               Kernel::LookupTable[B][b] * _CPImage(ia, jb);
        }
      }
    }
  }

  x = v._x, y = v._y;
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformation3D::FFD3D(double &x, double &y, double &z) const
{
  if (x < -2 || y < -2 || z < -2 || x >= _x || y >= _y || z >= _z) {
    x = y = z = .0;
    return;
  }

  int i = static_cast<int>(floor(x));
  int j = static_cast<int>(floor(y));
  int k = static_cast<int>(floor(z));

  const int A = Kernel::VariableToIndex(x - i);
  const int B = Kernel::VariableToIndex(y - j);
  const int C = Kernel::VariableToIndex(z - k);

  --i, --j, --k;

  CPValue v;
  double  w;
  int ia, jb, kc;

  for (int c = 0; c <= 3; ++c) {
    kc = k + c;
    if (0 <= kc && kc < _z) {
      for (int b = 0; b <= 3; ++b) {
        jb = j + b;
        if (0 <= jb && jb < _y) {
          w = Kernel::LookupTable[B][b] * Kernel::LookupTable[C][c];
          for (int a = 0; a <= 3; ++a) {
            ia = i + a;
            if (0 <= ia && ia < _x) {
              v += Kernel::LookupTable[A][a] * w * _CPImage(ia, jb, kc);
            }
          }
        }
      }
    }
  }

  x = v._x, y = v._y, z = v._z;
}

// -----------------------------------------------------------------------------
inline double irtkBSplineFreeFormTransformation3D::_Bending2D(int i, int j) const
{
  int I, J;
  double B_J, B_I, B_J_I, B_I_I, B_J_II, B_I_II, v;

  double y_jj=0, x_jj=0, y_ii=0, x_ii=0, y_ij=0, x_ij=0;

  double b[3]    = {1.0/6.0, 2.0/3.0, 1.0/6.0};
  double b_i[3]  = {-0.5, 0, 0.5};
  double b_ii[3] = {1.0, -2.0, 1.0};

  for (J = j-1; J < j+2; J++) {
    if (0 <= J && J < _y) {
      B_J    = b   [J-(j-1)];
      B_J_I  = b_i [J-(j-1)];
      B_J_II = b_ii[J-(j-1)];
      for (I = i-1; I < i+2; I++) {
        if (0 <= I && I < _x) {
          B_I    = b   [I-(i-1)];
          B_I_I  = b_i [I-(i-1)];
          B_I_II = b_ii[I-(i-1)];

          v = B_I_II * B_J;
          y_ii += _CPImage(I, J)._y * v;
          x_ii += _CPImage(I, J)._x * v;

          v = B_I * B_J_II;
          y_jj += _CPImage(I, J)._y * v;
          x_jj += _CPImage(I, J)._x * v;

          v = B_I_I * B_J_I;
          y_ij += _CPImage(I, J)._y * v;
          x_ij += _CPImage(I, J)._x * v;
        }
      }
    }
  }

  return (x_ii*x_ii + x_jj*x_jj + y_ii*y_ii + y_jj*y_jj + 2*(x_ij*x_ij + y_ij*y_ij));
}

// -----------------------------------------------------------------------------
inline double irtkBSplineFreeFormTransformation3D::_Bending2D(double x, double y) const
{
  int    I, J;
  double B_J, B_I, B_J_I, B_I_I, B_J_II, B_I_II, v;
  double y_jj=0, x_jj=0, y_ii=0, x_ii=0, y_ij=0, x_ij=0;

  double z = 0;
  this->WorldToLattice(x, y, z);

  const int l = static_cast<int>(floor(x));
  const int m = static_cast<int>(floor(y));

  const int S = Kernel::VariableToIndex(x - l);
  const int T = Kernel::VariableToIndex(y - m);

  for (int j = 0; j < 4; ++j) {
    J = j + m - 1;
    if (0 <= J && J < _y) {
      B_J    = Kernel::LookupTable   [T][j];
      B_J_I  = Kernel::LookupTable_I [T][j];
      B_J_II = Kernel::LookupTable_II[T][j];
      for (int i = 0; i < 4; ++i) {
        I = i + l - 1;
        if (0 <= I && I < _x) {
          B_I    = Kernel::LookupTable   [S][i];
          B_I_I  = Kernel::LookupTable_I [S][i];
          B_I_II = Kernel::LookupTable_II[S][i];

          v = B_I * B_J_II;
          y_jj += _CPImage(I, J)._y * v;
          x_jj += _CPImage(I, J)._x * v;

          v = B_I_II * B_J;
          y_ii += _CPImage(I, J)._y * v;
          x_ii += _CPImage(I, J)._x * v;

          v = B_I_I * B_J_I;
          y_ij += _CPImage(I, J)._y * v;
          x_ij += _CPImage(I, J)._x * v;
        }
      }
    }
  }

  return (x_ii*x_ii + x_jj*x_jj + y_ii*y_ii + y_jj*y_jj + 2*(x_ij*x_ij + y_ij*y_ij));
}

// -----------------------------------------------------------------------------
inline double irtkBSplineFreeFormTransformation3D::_Bending3D(int i, int j, int k) const
{
  int I, J, K;
  double B_K, B_J, B_I, B_K_I, B_J_I, B_I_I, B_K_II, B_J_II, B_I_II, v;

  double z_kk=0, y_kk=0, x_kk=0, z_jj=0, y_jj=0, x_jj=0, z_ii=0, y_ii=0, x_ii=0;
  double z_ij=0, y_ij=0, x_ij=0, z_ik=0, y_ik=0, x_ik=0, z_jk=0, y_jk=0, x_jk=0;

  double b[3]    = {1.0/6.0, 2.0/3.0, 1.0/6.0};
  double b_i[3]  = {-0.5, 0, 0.5};
  double b_ii[3] = {1.0, -2.0, 1.0};

  for (K = k-1; K < k+2; K++) {
    if (0 <= K && K < _z) {
      B_K    = b   [K-(k-1)];
      B_K_I  = b_i [K-(k-1)];
      B_K_II = b_ii[K-(k-1)];
      for (J = j-1; J < j+2; J++) {
        if (0 <= J && J < _y) {
          B_J    = b   [J-(j-1)];
          B_J_I  = b_i [J-(j-1)];
          B_J_II = b_ii[J-(j-1)];
          for (I = i-1; I < i+2; I++) {
            if (0 <= I && I < _x) {
              B_I    = b   [I-(i-1)];
              B_I_I  = b_i [I-(i-1)];
              B_I_II = b_ii[I-(i-1)];

              v = B_I * B_J * B_K_II;
              z_kk += _CPImage(I, J, K)._z * v;
              y_kk += _CPImage(I, J, K)._y * v;
              x_kk += _CPImage(I, J, K)._x * v;

              v = B_I * B_J_II * B_K;
              z_jj += _CPImage(I, J, K)._z * v;
              y_jj += _CPImage(I, J, K)._y * v;
              x_jj += _CPImage(I, J, K)._x * v;

              v = B_I_II * B_J * B_K;
              z_ii += _CPImage(I, J, K)._z * v;
              y_ii += _CPImage(I, J, K)._y * v;
              x_ii += _CPImage(I, J, K)._x * v;

              v = B_I_I * B_J_I * B_K;
              z_ij += _CPImage(I, J, K)._z * v;
              y_ij += _CPImage(I, J, K)._y * v;
              x_ij += _CPImage(I, J, K)._x * v;

              v = B_I_I * B_J * B_K_I;
              z_ik += _CPImage(I, J, K)._z * v;
              y_ik += _CPImage(I, J, K)._y * v;
              x_ik += _CPImage(I, J, K)._x * v;

              v = B_I * B_J_I * B_K_I;
              z_jk += _CPImage(I, J, K)._z * v;
              y_jk += _CPImage(I, J, K)._y * v;
              x_jk += _CPImage(I, J, K)._x * v;
            }
          }
        }
      }
    }
  }

  return (x_ii*x_ii + x_jj*x_jj + x_kk*x_kk +
          y_ii*y_ii + y_jj*y_jj + y_kk*y_kk +
          z_ii*z_ii + z_jj*z_jj + z_kk*z_kk +
          2*(x_ij*x_ij + y_ij*y_ij + z_ij*z_ij +
             x_ik*x_ik + y_ik*y_ik + z_ik*z_ik +
             x_jk*x_jk + y_jk*y_jk + z_jk*z_jk));
}

// -----------------------------------------------------------------------------
inline double irtkBSplineFreeFormTransformation3D::_Bending3D(double x, double y, double z) const
{
  int    I, J, K;
  double B_K, B_J, B_I, B_K_I, B_J_I, B_I_I, B_K_II, B_J_II, B_I_II, v;
  double z_kk=0, y_kk=0, x_kk=0, z_jj=0, y_jj=0, x_jj=0, z_ii=0, y_ii=0, x_ii=0;
  double z_ij=0, y_ij=0, x_ij=0, z_ik=0, y_ik=0, x_ik=0, z_jk=0, y_jk=0, x_jk=0;

  this->WorldToLattice(x, y, z);

  const int l = static_cast<int>(floor(x));
  const int m = static_cast<int>(floor(y));
  const int n = static_cast<int>(floor(z));

  const int S = Kernel::VariableToIndex(x - l);
  const int T = Kernel::VariableToIndex(y - m);
  const int U = Kernel::VariableToIndex(z - n);

  for (int k = 0; k < 4; ++k) {
    K = k + n - 1;
    if (0 <= K && K < _z) {
      B_K    = Kernel::LookupTable   [U][k];
      B_K_I  = Kernel::LookupTable_I [U][k];
      B_K_II = Kernel::LookupTable_II[U][k];
      for (int j = 0; j < 4; ++j) {
        J = j + m - 1;
        if (0 <= J && J < _y) {
          B_J    = Kernel::LookupTable   [T][j];
          B_J_I  = Kernel::LookupTable_I [T][j];
          B_J_II = Kernel::LookupTable_II[T][j];
          for (int i = 0; i < 4; ++i) {
            I = i + l - 1;
            if (0 <= I && I < _x) {
              B_I    = Kernel::LookupTable   [S][i];
              B_I_I  = Kernel::LookupTable_I [S][i];
              B_I_II = Kernel::LookupTable_II[S][i];

              v = B_I * B_J * B_K_II;
              z_kk += _CPImage(I, J, K)._z * v;
              y_kk += _CPImage(I, J, K)._y * v;
              x_kk += _CPImage(I, J, K)._x * v;

              v = B_I * B_J_II * B_K;
              z_jj += _CPImage(I, J, K)._z * v;
              y_jj += _CPImage(I, J, K)._y * v;
              x_jj += _CPImage(I, J, K)._x * v;

              v = B_I_II * B_J * B_K;
              z_ii += _CPImage(I, J, K)._z * v;
              y_ii += _CPImage(I, J, K)._y * v;
              x_ii += _CPImage(I, J, K)._x * v;

              v = B_I_I * B_J_I * B_K;
              z_ij += _CPImage(I, J, K)._z * v;
              y_ij += _CPImage(I, J, K)._y * v;
              x_ij += _CPImage(I, J, K)._x * v;

              v = B_I_I * B_J * B_K_I;
              z_ik += _CPImage(I, J, K)._z * v;
              y_ik += _CPImage(I, J, K)._y * v;
              x_ik += _CPImage(I, J, K)._x * v;

              v = B_I * B_J_I * B_K_I;
              z_jk += _CPImage(I, J, K)._z * v;
              y_jk += _CPImage(I, J, K)._y * v;
              x_jk += _CPImage(I, J, K)._x * v;
            }
          }
        }
      }
    }
  }

  return (x_ii*x_ii + x_jj*x_jj + x_kk*x_kk +
          y_ii*y_ii + y_jj*y_jj + y_kk*y_kk +
          z_ii*z_ii + z_jj*z_jj + z_kk*z_kk +
          2*(x_ij*x_ij + y_ij*y_ij + z_ij*z_ij +
             x_ik*x_ik + y_ik*y_ik + z_ik*z_ik +
             x_jk*x_jk + y_jk*y_jk + z_jk*z_jk));
}

// -----------------------------------------------------------------------------
double irtkBSplineFreeFormTransformation3D::Bending(double x, double y, double z) const
{
  return (_z == 1) ? _Bending2D(x, y) : _Bending3D(x, y, z);
}

// -----------------------------------------------------------------------------
double irtkBSplineFreeFormTransformation3D::Bending() const
{
  double bending = .0;

  if (_z == 1) {
    for (int j = 0; j < _y; ++j) {
      for (int i = 0; i < _x; ++i) {
        bending += _Bending2D(i, j);
      }
    }
  } else {
    for (int k = 0; k < _z; ++k) {
      for (int j = 0; j < _y; ++j) {
        for (int i = 0; i < _x; ++i) {
          bending += _Bending3D(i, j, k);
        }
      }
    }
  }

  return bending;
}

// -----------------------------------------------------------------------------
inline void irtkBSplineFreeFormTransformation3D::_BendingGradient2D(double *gradient) const
{
  double B_J, B_I, B_J_I, B_I_I, B_J_II, B_I_II, v;

  Vector **d_ii = Allocate<Vector>(_x, _y);
  Vector **d_ij = Allocate<Vector>(_x, _y);
  Vector **d_jj = Allocate<Vector>(_x, _y);

  for (int j = 0; j < _y; ++j) {
    for (int i = 0; i < _x; ++i) {
      for (int J = j-1; J < j+2; ++J) {
        if (0 <= J && J < _y) {
          B_J    = Kernel::LatticeWeights   [J-(j-1)];
          B_J_I  = Kernel::LatticeWeights_I [J-(j-1)];
          B_J_II = Kernel::LatticeWeights_II[J-(j-1)];
          for (int I = i-1; I < i+2; ++I) {
            if (0 <= I && I < _x) {
              B_I    = Kernel::LatticeWeights   [I-(i-1)];
              B_I_I  = Kernel::LatticeWeights_I [I-(i-1)];
              B_I_II = Kernel::LatticeWeights_II[I-(i-1)];

              v = B_I * B_J_II;
              d_jj[j][i]._y += 2 * _CPImage(I, J)._y * v;
              d_jj[j][i]._x += 2 * _CPImage(I, J)._x * v;

              v = B_I_II * B_J;
              d_ii[j][i]._y += 2 * _CPImage(I, J)._y * v;
              d_ii[j][i]._x += 2 * _CPImage(I, J)._x * v;

              v = B_I_I * B_J_I;
              d_ij[j][i]._y += 4 * _CPImage(I, J)._y * v;
              d_ij[j][i]._x += 4 * _CPImage(I, J)._x * v;
            }
          }
        }
      }
    }
  }

  double gx, gy;
  int    xdof, ydof;

  for (int j = 0; j < _y; ++j) {
    for (int i = 0; i < _x; ++i) {
      gx = gy = .0;
      for (int J = j-1; J < j+2; ++J) {
        if (0 <= J && J < _y) {
          B_J    = Kernel::LatticeWeights   [J-(j-1)];
          B_J_I  = Kernel::LatticeWeights_I [J-(j-1)];
          B_J_II = Kernel::LatticeWeights_II[J-(j-1)];
          for (int I = i-1; I < i+2; ++I) {
            if (0 <= I && I < _x) {
              B_I    = Kernel::LatticeWeights   [I-(i-1)];
              B_I_I  = Kernel::LatticeWeights_I [I-(i-1)];
              B_I_II = Kernel::LatticeWeights_II[I-(i-1)];

              v = B_I * B_J_II;
              gx += d_jj[J][I]._x * v;
              gy += d_jj[J][I]._y * v;

              v = B_I_II * B_J;
              gx += d_ii[J][I]._x * v;
              gy += d_ii[J][I]._y * v;

              v = B_I_I * B_J_I;
              gx += d_ij[J][I]._x * v;
              gy += d_ij[J][I]._y * v;
            }
          }
        }
      }
      this->IndexToDOFs(this->LatticeToIndex(i, j), xdof, ydof);
      gradient[xdof] += - gx;
      gradient[ydof] += - gy;
    }
  }

  Deallocate(d_ii);
  Deallocate(d_ij);
  Deallocate(d_jj);
}

// -----------------------------------------------------------------------------
inline void irtkBSplineFreeFormTransformation3D::_BendingGradient3D(double *gradient) const
{
  double B_K, B_J, B_I, B_K_I, B_J_I, B_I_I, B_K_II, B_J_II, B_I_II, v;

  Vector ***d_ii = Allocate<Vector>(_x, _y, _z);
  Vector ***d_ij = Allocate<Vector>(_x, _y, _z);
  Vector ***d_ik = Allocate<Vector>(_x, _y, _z);
  Vector ***d_jj = Allocate<Vector>(_x, _y, _z);
  Vector ***d_jk = Allocate<Vector>(_x, _y, _z);
  Vector ***d_kk = Allocate<Vector>(_x, _y, _z);

  for (int k = 0; k < _z; ++k) {
    for (int j = 0; j < _y; ++j) {
      for (int i = 0; i < _x; ++i) {
        for (int K = k-1; K < k+2; ++K) {
          if (0 <= K && K < _z) {
            B_K    = Kernel::LatticeWeights   [K-(k-1)];
            B_K_I  = Kernel::LatticeWeights_I [K-(k-1)];
            B_K_II = Kernel::LatticeWeights_II[K-(k-1)];
            for (int J = j-1; J < j+2; ++J) {
              if (0 <= J && J < _y) {
                B_J    = Kernel::LatticeWeights   [J-(j-1)];
                B_J_I  = Kernel::LatticeWeights_I [J-(j-1)];
                B_J_II = Kernel::LatticeWeights_II[J-(j-1)];
                for (int I = i-1; I < i+2; ++I) {
                  if (0 <= I && I < _x) {
                    B_I    = Kernel::LatticeWeights   [I-(i-1)];
                    B_I_I  = Kernel::LatticeWeights_I [I-(i-1)];
                    B_I_II = Kernel::LatticeWeights_II[I-(i-1)];

                    v = B_I * B_J * B_K_II;
                    d_kk[k][j][i]._z += 2 * _CPImage(I, J, K)._z * v;
                    d_kk[k][j][i]._y += 2 * _CPImage(I, J, K)._y * v;
                    d_kk[k][j][i]._x += 2 * _CPImage(I, J, K)._x * v;

                    v = B_I * B_J_II * B_K;
                    d_jj[k][j][i]._z += 2 * _CPImage(I, J, K)._z * v;
                    d_jj[k][j][i]._y += 2 * _CPImage(I, J, K)._y * v;
                    d_jj[k][j][i]._x += 2 * _CPImage(I, J, K)._x * v;

                    v = B_I_II * B_J * B_K;
                    d_ii[k][j][i]._z += 2 * _CPImage(I, J, K)._z * v;
                    d_ii[k][j][i]._y += 2 * _CPImage(I, J, K)._y * v;
                    d_ii[k][j][i]._x += 2 * _CPImage(I, J, K)._x * v;

                    v = B_I_I * B_J_I * B_K;
                    d_ij[k][j][i]._z += 4 * _CPImage(I, J, K)._z * v;
                    d_ij[k][j][i]._y += 4 * _CPImage(I, J, K)._y * v;
                    d_ij[k][j][i]._x += 4 * _CPImage(I, J, K)._x * v;

                    v = B_I_I * B_J * B_K_I;
                    d_ik[k][j][i]._z += 4 * _CPImage(I, J, K)._z * v;
                    d_ik[k][j][i]._y += 4 * _CPImage(I, J, K)._y * v;
                    d_ik[k][j][i]._x += 4 * _CPImage(I, J, K)._x * v;

                    v = B_I * B_J_I * B_K_I;
                    d_jk[k][j][i]._z += 4 * _CPImage(I, J, K)._z * v;
                    d_jk[k][j][i]._y += 4 * _CPImage(I, J, K)._y * v;
                    d_jk[k][j][i]._x += 4 * _CPImage(I, J, K)._x * v;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  double gx, gy, gz;
  int    xdof, ydof, zdof;

  for (int k = 0; k < _z; ++k) {
    for (int j = 0; j < _y; ++j) {
      for (int i = 0; i < _x; ++i) {
        gx = gy = gz = .0;
        for (int K = k-1; K < k+2; ++K) {
          if (0 <= K && K < _z) {
            B_K    = Kernel::LatticeWeights   [K-(k-1)];
            B_K_I  = Kernel::LatticeWeights_I [K-(k-1)];
            B_K_II = Kernel::LatticeWeights_II[K-(k-1)];
            for (int J = j-1; J < j+2; ++J) {
              if (0 <= J && J < _y) {
                B_J    = Kernel::LatticeWeights   [J-(j-1)];
                B_J_I  = Kernel::LatticeWeights_I [J-(j-1)];
                B_J_II = Kernel::LatticeWeights_II[J-(j-1)];
                for (int I = i-1; I < i+2; ++I) {
                  if (0 <= I && I < _x) {
                    B_I    = Kernel::LatticeWeights   [I-(i-1)];
                    B_I_I  = Kernel::LatticeWeights_I [I-(i-1)];
                    B_I_II = Kernel::LatticeWeights_II[I-(i-1)];

                    v = B_I * B_J * B_K_II;
                    gx += d_kk[K][J][I]._x * v;
                    gy += d_kk[K][J][I]._y * v;
                    gz += d_kk[K][J][I]._z * v;

                    v = B_I * B_J_II * B_K;
                    gx += d_jj[K][J][I]._x * v;
                    gy += d_jj[K][J][I]._y * v;
                    gz += d_jj[K][J][I]._z * v;

                    v = B_I_II * B_J * B_K;
                    gx += d_ii[K][J][I]._x * v;
                    gy += d_ii[K][J][I]._y * v;
                    gz += d_ii[K][J][I]._z * v;

                    v = B_I_I * B_J_I * B_K;
                    gx += d_ij[K][J][I]._x * v;
                    gy += d_ij[K][J][I]._y * v;
                    gz += d_ij[K][J][I]._z * v;

                    v = B_I_I * B_J * B_K_I;
                    gx += d_ik[K][J][I]._x * v;
                    gy += d_ik[K][J][I]._y * v;
                    gz += d_ik[K][J][I]._z * v;

                    v = B_I * B_J_I * B_K_I;
                    gx += d_jk[K][J][I]._x * v;
                    gy += d_jk[K][J][I]._y * v;
                    gz += d_jk[K][J][I]._z * v;
                  }
                }
              }
            }
          }
        }
        this->IndexToDOFs(this->LatticeToIndex(i, j, k), xdof, ydof, zdof);
        gradient[xdof] += - gx;
        gradient[ydof] += - gy;
        gradient[zdof] += - gz;
      }
    }
  }

  Deallocate(d_ii);
  Deallocate(d_ij);
  Deallocate(d_ik);
  Deallocate(d_jj);
  Deallocate(d_jk);
  Deallocate(d_kk);
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformation3D::BendingGradient(double *gradient) const
{
  if (_z == 1) _BendingGradient2D(gradient);
  else         _BendingGradient3D(gradient);
}
