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

#ifdef HAS_VTK

extern bool interactiveVTK;
extern bool displayVTK;
extern bool firstVTK;

extern void update   (irtkPointSet& pset, int nx, int ny, int nz);
extern void visualize(irtkPointSet& pset, int nx, int ny, int nz);

#endif

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
void irtkRigidTransformation::UpdateRotationSineCosine()
{
  const double radians = M_PI / 180.0;
  _sinrx = sin(_Param[RX] * radians);
  _sinry = sin(_Param[RY] * radians);
  _sinrz = sin(_Param[RZ] * radians);
  _cosrx = cos(_Param[RX] * radians);
  _cosry = cos(_Param[RY] * radians);
  _cosrz = cos(_Param[RZ] * radians);
}

// -----------------------------------------------------------------------------
irtkRigidTransformation::irtkRigidTransformation(int ndofs)
:
  irtkHomogeneousTransformation(ndofs)
{
  UpdateRotationSineCosine();
}

// -----------------------------------------------------------------------------
irtkRigidTransformation::irtkRigidTransformation(const irtkRigidTransformation &t, int ndofs)
:
  irtkHomogeneousTransformation(t, ndofs),
  _sinrx(t._sinrx),
  _sinry(t._sinry),
  _sinrz(t._sinrz),
  _cosrx(t._cosrx),
  _cosry(t._cosry),
  _cosrz(t._cosrz)
{
}

// -----------------------------------------------------------------------------
irtkRigidTransformation::irtkRigidTransformation()
:
  irtkHomogeneousTransformation(6)
{
  UpdateRotationSineCosine();
}

// -----------------------------------------------------------------------------
irtkRigidTransformation::irtkRigidTransformation(const irtkRigidTransformation &t)
:
  irtkHomogeneousTransformation(t),
  _sinrx(t._sinrx),
  _sinry(t._sinry),
  _sinrz(t._sinrz),
  _cosrx(t._cosrx),
  _cosry(t._cosry),
  _cosrz(t._cosrz)
{
}

// -----------------------------------------------------------------------------
irtkRigidTransformation::~irtkRigidTransformation()
{
}

// =============================================================================
// Approximation
// =============================================================================

// -----------------------------------------------------------------------------
void irtkRigidTransformation
::ApproximateDOFs(const double *x,  const double *y,  const double *z, const double *,
                  const double *dx, const double *dy, const double *dz, int no)
{
  if (no < 4) {
    cerr << "irtkRigidTransformation::ApproximateDOFs: Must have at least four points" << endl;
    exit(1);
  }

  // Calculate centroids
  irtkPoint c1, c2;
  for (int n = 0; n < no; ++n) {
    c1._x += x[n] / no;
    c1._y += y[n] / no;
    c1._z += z[n] / no;
    c2._x += (x[n] + dx[n]) / no;
    c2._y += (y[n] + dy[n]) / no;
    c2._z += (z[n] + dz[n]) / no;
  }

  // Subtract centroids
  irtkPointSet target(no), source(no);
  for (int n = 0; n < no; ++n) {
    irtkPoint &p1 = target(n);
    irtkPoint &p2 = source(n);
    p1._x = x[n] - c1._x;
    p1._y = y[n] - c1._y;
    p1._z = z[n] - c1._z;
    p2._x = x[n] + dx[n] - c2._x;
    p2._y = y[n] + dy[n] - c2._y;
    p2._z = z[n] + dz[n] - c2._z;
  }

  // Calculate covariance matrix
  irtkMatrix h(3, 3);
  irtkMatrix a(3, 1);
  irtkMatrix b(1, 3);
  for (int n = 0; n < no; ++n) {
    a(0, 0) = target(n)._x;
    a(1, 0) = target(n)._y;
    a(2, 0) = target(n)._z;
    b(0, 0) = source(n)._x;
    b(0, 1) = source(n)._y;
    b(0, 2) = source(n)._z;
    h += a * b;
  }

  // Calculate SVD
  irtkMatrix u, v;
  irtkVector w;
  h.SVD(u, w, v);

  // Calculate rotation matrix
  u.Transpose();
  irtkMatrix r = v * u;

  // Reflect axis with most singular value if necessary
  if (r.Det () < 0.999) {

    // Search for most singular value
    int i = 0;
    if (w(0) < w(1) && w(0) < w(2)) i = 0;
    if (w(1) < w(0) && w(1) < w(2)) i = 1;
    if (w(2) < w(1) && w(2) < w(0)) i = 2;

    // Multiply column with most singular value by -1
    v(0, i) *= -1.0;
    v(1, i) *= -1.0;
    v(2, i) *= -1.0;

    // Recalculate rotation matrix
    r = v * u;
  }

  // Calculate rotated centroid
  const double cx = c1._x, cy = c1._y, cz = c1._z;
  c1._x = r(0, 0) * cx + r(0, 1) * cy + r(0, 2) * cz;
  c1._y = r(1, 0) * cx + r(1, 1) * cy + r(1, 2) * cz;
  c1._z = r(2, 0) * cx + r(2, 1) * cy + r(2, 2) * cz;

  // Set homogeneous transformation matrix
  irtkMatrix m(4, 4);
  m(0, 0) = r(0, 0);
  m(1, 0) = r(1, 0);
  m(2, 0) = r(2, 0);
  m(3, 0) = 0.0;
  m(0, 1) = r(0, 1);
  m(1, 1) = r(1, 1);
  m(2, 1) = r(2, 1);
  m(3, 1) = 0.0;
  m(0, 2) = r(0, 2);
  m(1, 2) = r(1, 2);
  m(2, 2) = r(2, 2);
  m(3, 2) = 0.0;
  m(0, 3) = c2._x - c1._x;
  m(1, 3) = c2._y - c1._y;
  m(2, 3) = c2._z - c1._z;
  m(3, 3) = 1.0;
  this->PutMatrix(m);
}

// =============================================================================
// Transformation parameters
// =============================================================================

// -----------------------------------------------------------------------------
irtkMatrix irtkRigidTransformation::DOFs2Matrix(const double *param)
{
  irtkMatrix matrix(4, 4);
  const double radians = M_PI / 180.0;
  RigidParametersToMatrix(param[TX], param[TY], param[TZ],
                          param[RX] * radians,
                          param[RY] * radians,
                          param[RZ] * radians,
                          matrix);
  return matrix;
}

// -----------------------------------------------------------------------------
void irtkRigidTransformation::Matrix2DOFs(const irtkMatrix &matrix, double *param)
{
  // Get rigid parameters
  MatrixToRigidParameters(matrix, param[TX], param[TY], param[TZ],
                                  param[RX], param[RY], param[RZ]);

  // Convert angles to degrees
  const double degrees = 180.0 / M_PI;
  param[RX] *= degrees;
  param[RY] *= degrees;
  param[RZ] *= degrees;
}

// -----------------------------------------------------------------------------
void irtkRigidTransformation::UpdateMatrix()
{
  // Recompute sines and cosines
  UpdateRotationSineCosine();

  // Update transformation matrix
  RigidParametersToMatrix(_Param[TX], _Param[TY], _Param[TZ],
                          _cosrx, _cosry, _cosrz,
                          _sinrx, _sinry, _sinrz, _matrix);
  _inverse = _matrix.Inverse();
}

// -----------------------------------------------------------------------------
void irtkRigidTransformation::UpdateDOFs()
{
  // Get rigid parameters
  MatrixToRigidParameters(_matrix, _Param[TX], _Param[TY], _Param[TZ],
                                   _Param[RX], _Param[RY], _Param[RZ]);

  // Update sines and cosines
  _cosrx = cos(_Param[RX]);
  _cosry = cos(_Param[RY]);
  _cosrz = cos(_Param[RZ]);
  _sinrx = sin(_Param[RX]);
  _sinry = sin(_Param[RY]);
  _sinrz = sin(_Param[RZ]);

  // Convert angles to degrees
  const double degrees = 180.0 / M_PI;
  _Param[RX] *= degrees;
  _Param[RY] *= degrees;
  _Param[RZ] *= degrees;
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
void irtkRigidTransformation::JacobianDOFs(double jac[3], int dof, double x, double y, double z, double, double) const
{
  switch (dof) {
    case TX: {
      jac[0] = 1.0;
      jac[1] = .0;
      jac[2] = .0;
    } break;
    case TY: {
      jac[0] = .0;
      jac[1] = 1.0;
      jac[2] = .0;
    } break;
    case TZ: {
      jac[0] = .0;
      jac[1] = .0;
      jac[2] = 1.0;
    } break;
    case RX: {
      const double deg2rad = M_PI / 180.0;

      jac[0] = .0;
      jac[1] = (+ (_cosrx * _sinry * _cosrz + _sinrx * _sinrz) * x
                + (_cosrx * _sinry * _sinrz - _sinrx * _cosrz) * y
                + (_cosrx * _cosry                           ) * z) * deg2rad;
      jac[2] = (+ (_cosrx * _sinrz - _sinrx * _sinry * _cosrz) * x
                - (_cosrx * _cosrz + _sinrx * _sinry * _sinrz) * y
                - (_sinrx * _cosry                           ) * z) * deg2rad;
    } break;
    case RY: {
      const double deg2rad = M_PI / 180.0;
      jac[0] = (- (_sinry * _cosrz         ) * x
                - (_sinry * _sinrz         ) * y
                - (_cosry                  ) * z) * deg2rad;
      jac[1] = (+ (_sinrx * _cosry * _cosrz) * x
                + (_sinrx * _cosry * _sinrz) * y
                - (_sinry * _sinrx         ) * z) * deg2rad;
      jac[2] = (+ (_cosrx * _cosry * _cosrz) * x
                + (_cosrx * _cosry * _sinrz) * y
                - (_cosrx * _sinry         ) * z) * deg2rad;
    } break;
    case RZ: {
      const double deg2rad = M_PI / 180.0;
      jac[0] = (- (_sinrz * _cosry                             ) * x
                + (_cosry * _cosrz                             ) * y) * deg2rad;
      jac[1] = (+ (- _sinrx * _sinry * _sinrz - _cosrx * _cosrz) * x
                + (  _sinrx * _sinry * _cosrz - _cosrx * _sinrz) * y) * deg2rad;
      jac[2] = (+ (- _cosrx * _sinry * _sinrz + _sinrx * _cosrz) * x
                + (  _cosrx * _sinry * _cosrz + _sinrx * _sinrz) * y) * deg2rad;
    } break;
    default:
      cerr << this->NameOfClass() << "::JacobianDOFs(): No such parameter: " << dof << endl;
      exit(1);
  }
}

// -----------------------------------------------------------------------------
void irtkRigidTransformation::DeriveJacobianWrtDOF(irtkMatrix &dJdp, int dof, double, double, double, double, double) const
{
  dJdp.Initialize(3, 3);

  switch (dof) {
    case TX: case TY: case TZ: {
      //Do nothing
      return;
    } break;
    case RX: {
      const double deg2rad = M_PI / 180.0;

      dJdp(0, 0) = .0;
      dJdp(1, 0) = + (_cosrx * _sinry * _cosrz + _sinrx * _sinrz);
      dJdp(2, 0) = + (_cosrx * _sinrz - _sinrx * _sinry * _cosrz);
      dJdp(0, 1) = .0;
      dJdp(1, 1) = + (_cosrx * _sinry * _sinrz - _sinrx * _cosrz);
      dJdp(2, 1) = - (_cosrx * _cosrz + _sinrx * _sinry * _sinrz);
      dJdp(0, 2) = .0;
      dJdp(1, 2) = + (_cosrx * _cosry);
      dJdp(2, 2) = - (_sinrx * _cosry);

      dJdp *= deg2rad;
    } break;
    case RY: {
      const double deg2rad = M_PI / 180.0;

      dJdp(0, 0) = - (_sinry * _cosrz);
      dJdp(1, 0) = + (_sinrx * _cosry * _cosrz);
      dJdp(2, 0) = + (_cosrx * _cosry * _cosrz);
      dJdp(0, 1) = - (_sinry * _sinrz);
      dJdp(1, 1) = + (_sinrx * _cosry * _sinrz);
      dJdp(2, 1) = + (_cosrx * _cosry * _sinrz);
      dJdp(0, 2) = - (_cosry);
      dJdp(1, 2) = - (_sinry * _sinrx);
      dJdp(2, 2) = - (_cosrx * _sinry);

      dJdp *= deg2rad;
    } break;
    case RZ: {
      const double deg2rad = M_PI / 180.0;

      dJdp(0, 0) = - (_sinrz * _cosry);
      dJdp(1, 0) = + (-_sinrx * _sinry * _sinrz - _cosrx * _cosrz);
      dJdp(2, 0) = + (-_cosrx * _sinry * _sinrz + _sinrx * _cosrz);
      dJdp(0, 1) = + (_cosry * _cosrz);
      dJdp(1, 1) = + (_sinrx * _sinry * _cosrz - _cosrx * _sinrz);
      dJdp(2, 1) = + (_cosrx * _sinry * _cosrz + _sinrx * _sinrz);
      dJdp(0, 2) = .0;
      dJdp(1, 2) = .0;
      dJdp(2, 2) = .0;

      dJdp *= deg2rad;
    } break;
    default:
      cerr << this->NameOfClass() << "::DeriveJacobianWrtDOF(): No such parameter: " << dof << endl;
      exit(1);
  }
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void irtkRigidTransformation::Print(irtkIndent indent) const
{
  cout.setf(ios::right);
  cout.setf(ios::fixed);
  streamsize previous_precision = cout.precision(4);

  if (_Status[TX] == _Active || !fequal(_Param[TX], .0, 1e-4) ||
      _Status[TY] == _Active || !fequal(_Param[TY], .0, 1e-4) ||
      _Status[TZ] == _Active || !fequal(_Param[TZ], .0, 1e-4)) {
    cout << indent;
    if (_Status[TX] == _Active) cout << "tx  = " << setw(8) << _Param[TX] << "  ";
    if (_Status[TY] == _Active) cout << "ty  = " << setw(8) << _Param[TY] << "  ";
    if (_Status[TZ] == _Active) cout << "tz  = " << setw(8) << _Param[TZ];
    cout << endl;
  }
  if (_Status[RX] == _Active || !fequal(_Param[RX], .0, 1e-4) ||
      _Status[RY] == _Active || !fequal(_Param[RY], .0, 1e-4) ||
      _Status[RZ] == _Active || !fequal(_Param[RZ], .0, 1e-4)) {
    cout << indent;
    if (_Status[RX] == _Active) cout << "rx  = " << setw(8) << _Param[RX] << "  ";
    if (_Status[RY] == _Active) cout << "ry  = " << setw(8) << _Param[RY] << "  ";
    if (_Status[RZ] == _Active) cout << "rz  = " << setw(8) << _Param[RZ];
    cout << endl;
  }

  cout.precision(previous_precision);
  cout.unsetf(ios::right);
  cout.unsetf(ios::fixed);
}
