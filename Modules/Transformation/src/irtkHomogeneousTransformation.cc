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

#include <irtkPointSamples.h>
#include <irtkAdaptiveLineSearch.h>
#include <irtkConjugateGradientDescent.h>
#include <irtkTransformationApproximationError.h>


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
// Note: DoFs are parameters (e.g., rotation) from which 4x4 matrix is constructed
irtkHomogeneousTransformation::irtkHomogeneousTransformation(int ndofs)
:
  irtkTransformation(ndofs),
  _matrix (4, 4),
  _inverse(4, 4)
{
  _matrix .Ident();
  _inverse.Ident();
}

// -----------------------------------------------------------------------------
// Note: DoFs are parameters (e.g., rotation) from which 4x4 matrix is constructed
irtkHomogeneousTransformation::irtkHomogeneousTransformation(const irtkHomogeneousTransformation &t, int ndofs)
:
  irtkTransformation(t, ndofs),
  _matrix (t._matrix),
  _inverse(t._inverse)
{
}

// -----------------------------------------------------------------------------
// Note: DoFs are the matrix elements
irtkHomogeneousTransformation::irtkHomogeneousTransformation()
:
  irtkTransformation(16),
  _inverse(4, 4)
{
  // Use memory allocated by irtkTransformation when possible
  if (sizeof(DOFValue) == sizeof(irtkMatrix::ElementType)) {
    _matrix.Initialize(4, 4, _Param);
  } else {
    _matrix.Initialize(4, 4);
  }
  _matrix .Ident();
  _inverse.Ident();
  // Last row should always be "0 0 0 1"
  _Status[12] = Passive;
  _Status[13] = Passive;
  _Status[14] = Passive;
  _Status[15] = Passive;
}

// -----------------------------------------------------------------------------
// Note: DoFs are the matrix elements
irtkHomogeneousTransformation::irtkHomogeneousTransformation(const irtkMatrix &matrix)
:
  irtkTransformation(16),
  _matrix (4, 4, _Param), // uses memory allocated by irtkTransformation
  _inverse(4, 4)
{
  _matrix  = matrix;
  _inverse = matrix.Inverse();
  // Last row should always be "0 0 0 1"
  _Status[12] = Passive;
  _Status[13] = Passive;
  _Status[14] = Passive;
  _Status[15] = Passive;
}

// -----------------------------------------------------------------------------
// Note: DoFs are either matrix elements if other transformation is an instance
//       of irtkHomogeneousTransformation itself and not a subclass or
//       parameters (e.g., rotation) from which 4x4 matrix is constructed
irtkHomogeneousTransformation::irtkHomogeneousTransformation(const irtkHomogeneousTransformation &t)
:
  irtkTransformation(t)
{
  if (t.GetMatrix().GetPointerToElements() == t._Param) {
    _matrix.Initialize(4, 4, _Param);
  } else {
    _matrix = t._matrix;
  }
  _inverse = t._inverse;
}

// -----------------------------------------------------------------------------
irtkHomogeneousTransformation::~irtkHomogeneousTransformation()
{
}

// =============================================================================
// Approximation
// =============================================================================

// -----------------------------------------------------------------------------
double irtkHomogeneousTransformation
::Approximate(const irtkImageAttributes &domain,
              double *dx, double *dy, double *dz,
              int niter, double max_error)
{
  // Check input arguments
  const int no = domain.NumberOfPoints();
  if (no <= 0) return .0;
  double error = numeric_limits<double>::infinity();
  if (niter < 1) return error;

  // Compute world coordinates of lattice points
  double *x = Allocate<double>(no);
  double *y = Allocate<double>(no);
  double *z = Allocate<double>(no);
  double *t = Allocate<double>(no);
  domain.LatticeToWorld(x, y, z, t);

  // Allocate memory for transformed world coordinates
  double *wx = Allocate<double>(no);
  double *wy = Allocate<double>(no);
  double *wz = Allocate<double>(no);

  // Copy original displacements
  double *tx = Allocate<double>(no);
  double *ty = Allocate<double>(no);
  double *tz = Allocate<double>(no);
  memcpy(tx, dx, no * sizeof(double));
  memcpy(ty, dy, no * sizeof(double));
  memcpy(tz, dz, no * sizeof(double));

  // Evaluate approximation error and residual displacements
  if (this->RequiresCachingOfDisplacements()) {
    error = EvaluateRMSError(domain, dx, dy, dz);
  } else {
    error = EvaluateRMSError(x, y, z, t, dx, dy, dz, no);
  }

  // Repeat approximation n times or until error drops below a threshold
  for (int iter = 0; iter < niter && error > max_error; ++iter) {

    // Transform world coordinates
    memcpy(wx, x, no * sizeof(double));
    memcpy(wy, y, no * sizeof(double));
    memcpy(wz, z, no * sizeof(double));
    this->Transform(no, wx, wy, wz, t);

    // Approximate residual displacements by new parameters
    const irtkMatrix matrix = _matrix;
    this->ApproximateDOFs(wx, wy, wz, t, dx, dy, dz, no);

    // Compose with previous transformation
    this->PutMatrix(_matrix * matrix);

    // Evaluate approximation error and residual displacements
    memcpy(dx, tx, no * sizeof(double));
    memcpy(dy, ty, no * sizeof(double));
    memcpy(dz, tz, no * sizeof(double));
    if (this->RequiresCachingOfDisplacements()) {
      error = EvaluateRMSError(domain, dx, dy, dz);
    } else {
      error = EvaluateRMSError(x, y, z, t, dx, dy, dz, no);
    }
  }

  // Free memory
  Deallocate(tx);
  Deallocate(ty);
  Deallocate(tz);
  Deallocate(wx);
  Deallocate(wy);
  Deallocate(wz);
  Deallocate(x);
  Deallocate(y);
  Deallocate(z);
  Deallocate(t);

  return error;
}

// -----------------------------------------------------------------------------
double irtkHomogeneousTransformation
::Approximate(const double *x,  const double *y,  const double *z,
              double       *dx, double       *dy, double       *dz, int no,
              int niter, double max_error)
{
  // Check input arguments
  if (no <= 0) return .0;
  double error = numeric_limits<double>::infinity();
  if (niter < 1) return error;

  // Compute (fixed) time coordinates
  const double t0 = .0;
  double *t = CAllocate<double>(no, &t0);

  // Allocate memory for transformed world coordinates
  double *wx = Allocate<double>(no);
  double *wy = Allocate<double>(no);
  double *wz = Allocate<double>(no);

  // Copy original displacements
  double *tx = Allocate<double>(no);
  double *ty = Allocate<double>(no);
  double *tz = Allocate<double>(no);
  memcpy(tx, dx, no * sizeof(double));
  memcpy(ty, dy, no * sizeof(double));
  memcpy(tz, dz, no * sizeof(double));

  // Evaluate approximation error and residual displacements
  error = EvaluateRMSError(x, y, z, t0, dx, dy, dz, no);

  // Repeat approximation n times or until error drops below a threshold
  for (int iter = 0; iter < niter && error > max_error; ++iter) {

    // Transform world coordinates
    memcpy(wx, x, no * sizeof(double));
    memcpy(wy, y, no * sizeof(double));
    memcpy(wz, z, no * sizeof(double));
    this->Transform(no, wx, wy, wz, t0);

    // Approximate residual displacements by new parameters
    const irtkMatrix matrix = _matrix;
    this->ApproximateDOFs(wx, wy, wz, t, dx, dy, dz, no);

    // Compose with previous transformation
    this->PutMatrix(_matrix * matrix);

    // Evaluate error of approximation and residual displacements
    memcpy(dx, tx, no * sizeof(double));
    memcpy(dy, ty, no * sizeof(double));
    memcpy(dz, tz, no * sizeof(double));
    error = EvaluateRMSError(x, y, z, t0, dx, dy, dz, no);
  }

  // Free memory
  Deallocate(tx);
  Deallocate(ty);
  Deallocate(tz);
  Deallocate(wx);
  Deallocate(wy);
  Deallocate(wz);
  Deallocate(t);

  return error;
}

// -----------------------------------------------------------------------------
double irtkHomogeneousTransformation
::Approximate(const double *x,  const double *y,  const double *z,  const double *t,
              double       *dx, double       *dy, double       *dz, int no,
              int niter, double max_error)
{
  // Check input arguments
  if (no <= 0) return .0;
  double error = numeric_limits<double>::infinity();
  if (niter < 1) return error;

  // Allocate memory for transformed world coordinates
  double *wx = Allocate<double>(no);
  double *wy = Allocate<double>(no);
  double *wz = Allocate<double>(no);

  // Copy original displacements
  double *tx = Allocate<double>(no);
  double *ty = Allocate<double>(no);
  double *tz = Allocate<double>(no);
  memcpy(tx, dx, no * sizeof(double));
  memcpy(ty, dy, no * sizeof(double));
  memcpy(tz, dz, no * sizeof(double));

  // Evaluate error of approximation and residual displacements
  error = EvaluateRMSError(x, y, z, t, dx, dy, dz, no);

  // Repeat approximation n times or until error drops below a threshold
  for (int iter = 0; iter < niter; ++iter) {

    // Transform world coordinates
    memcpy(wx, x, no * sizeof(double));
    memcpy(wy, y, no * sizeof(double));
    memcpy(wz, z, no * sizeof(double));
    this->Transform(no, wx, wy, wz, t);

    // Approximate residual displacements by new parameters
    const irtkMatrix matrix = _matrix;
    this->ApproximateDOFs(wx, wy, wz, t, dx, dy, dz, no);

    // Compose with previous transformation
    this->PutMatrix(_matrix * matrix);

    // Evaluate error of approximation and residual displacements
    memcpy(dx, tx, no * sizeof(double));
    memcpy(dy, ty, no * sizeof(double));
    memcpy(dz, tz, no * sizeof(double));
    error = EvaluateRMSError(x, y, z, t, dx, dy, dz, no);
  }

  // Free memory
  Deallocate(tx);
  Deallocate(ty);
  Deallocate(tz);
  Deallocate(wx);
  Deallocate(wy);
  Deallocate(wz);

  return error;
}

// -----------------------------------------------------------------------------
double irtkHomogeneousTransformation::Approximate(const irtkMatrix &matrix)
{
  return this->ApproximateAsNew(matrix);
}

// -----------------------------------------------------------------------------
double irtkHomogeneousTransformation::ApproximateAsNew(const irtkMatrix &matrix)
{
  if (this->NumberOfActiveDOFs() == 12) {
    this->PutMatrix(matrix);
    return .0;
  }

  const int    nsamples = 100;
  const double t0       =  .0;

  irtkPointSamples samples(nsamples, /*fixed seed=*/ 0);
  samples.SampleSphere(.0, 100.0);

  double *x  =  Allocate<double>(nsamples);
  double *y  =  Allocate<double>(nsamples);
  double *z  =  Allocate<double>(nsamples);
  double *t  = CAllocate<double>(nsamples, &t0);
  double *dx =  Allocate<double>(nsamples);
  double *dy =  Allocate<double>(nsamples);
  double *dz =  Allocate<double>(nsamples);

  for (int i = 0; i < nsamples; ++i) {
    x [i] = samples(i)._x;
    y [i] = samples(i)._y;
    z [i] = samples(i)._z;
    dx[i] = matrix(0, 0) * x[i] + matrix(0, 1) * y[i] + matrix(0, 2) * z[i] + matrix(0, 3) - x[i];
    dy[i] = matrix(1, 0) * x[i] + matrix(1, 1) * y[i] + matrix(1, 2) * z[i] + matrix(1, 3) - y[i];
    dz[i] = matrix(2, 0) * x[i] + matrix(2, 1) * y[i] + matrix(2, 2) * z[i] + matrix(2, 3) - z[i];
  }

  double rms = this->ApproximateAsNew(x, y, z, t, dx, dy, dz, nsamples);

  Deallocate(x);
  Deallocate(y);
  Deallocate(z);
  Deallocate(t);
  Deallocate(dx);
  Deallocate(dy);
  Deallocate(dz);

  return rms;
}

// -----------------------------------------------------------------------------
void irtkHomogeneousTransformation
::ApproximateDOFs(const double *x,  const double *y,  const double *z, const double *t,
                  const double *dx, const double *dy, const double *dz, int no)
{
  IRTK_START_TIMING();

  // Note: Do not reset transformation parameters as this function is usually
  //       called by subclass implementations after these were initialized!
  irtkMatrix pre(4, 4), post(4, 4);

  // Mean squared error function
  irtkTransformationApproximationError error(this, x, y, z, t, dx, dy, dz, no);

  // Center point sets
  if (_Status[TX] == Active && _Status[TY] == Active && _Status[TZ] == Active) {
    error.CenterPoints();
  }

  // Adjust transformation
  pre.Ident();
  pre(0, 3)  = + error.TargetCenter()._x;
  pre(1, 3)  = + error.TargetCenter()._y;
  pre(2, 3)  = + error.TargetCenter()._z;

  post.Ident();
  post(0, 3) = - error.SourceCenter()._x;
  post(1, 3) = - error.SourceCenter()._y;
  post(2, 3) = - error.SourceCenter()._z;

  this->PutMatrix(post * this->GetMatrix() * pre);

  // Optimization method
  irtkAdaptiveLineSearch linesearch;
  linesearch.MaxRejectedStreak(0);
  linesearch.StrictStepLengthRange(false);
  linesearch.ReusePreviousStepLength(true);
  linesearch.NumberOfIterations(20);
  linesearch.MinStepLength(1e-6);
  linesearch.MaxStepLength(10.0);

  irtkConjugateGradientDescent optimizer;
  optimizer.Function(&error);
  optimizer.LineSearch(&linesearch);
  optimizer.NumberOfSteps(100);
  optimizer.Epsilon(1e-6);
  optimizer.Delta(1e-12);

  // Find transformation parameters which minimize the approximation error
  optimizer.Run();
  optimizer.ConjugateGradientOff();
  optimizer.Run();

  // Include centering transformations in final transformation
  pre.Ident();
  pre(0, 3)  = - error.TargetCenter()._x;
  pre(1, 3)  = - error.TargetCenter()._y;
  pre(2, 3)  = - error.TargetCenter()._z;

  post.Ident();
  post(0, 3) = + error.SourceCenter()._x;
  post(1, 3) = + error.SourceCenter()._y;
  post(2, 3) = + error.SourceCenter()._z;

  this->PutMatrix(post * this->GetMatrix() * pre);
  IRTK_DEBUG_TIMING(5, "irtkHomogeneousTransformation::ApproximateDOFs");
}

// =============================================================================
// Transformation parameters
// =============================================================================

// -----------------------------------------------------------------------------
void irtkHomogeneousTransformation::UpdateMatrix()
{
  if (_matrix.GetPointerToElements() != _Param) {
    if (NumberOfDOFs() != 16) {
      cerr << "irtkHomogeneousTransformation::UpdateMatrix: Override in subclass!" << endl;
      exit(1);
    }
    memcpy(_matrix.GetPointerToElements(), _Param, 16 * sizeof(double));
  }
  _inverse = _matrix.Inverse();
}

// -----------------------------------------------------------------------------
void irtkHomogeneousTransformation::UpdateDOFs()
{
  if (_matrix.GetPointerToElements() != _Param) {
    if (NumberOfDOFs() != 16) {
      cerr << "irtkHomogeneousTransformation::UpdateDOFs: Override in subclass!" << endl;
      exit(1);
    }
    memcpy(_Param, _matrix.GetPointerToElements(), 16 * sizeof(double));
  }
}

// ---------------------------------------------------------------------------
bool irtkHomogeneousTransformation::CopyFrom(const irtkTransformation *other)
{
  const irtkHomogeneousTransformation *linear;
  if ((linear = dynamic_cast<const irtkHomogeneousTransformation *>(other))) {
    this->Reset();
    const int ndofs = min(_NumberOfDOFs, other->NumberOfDOFs());
    for (int dof = 0; dof < ndofs; ++dof) {
      if (_Status[dof] == Active) _Param[dof] = other->Get(dof);
    }
    this->Update(MATRIX);
    return true;
  } else {
    return false;
  }
}

// -----------------------------------------------------------------------------
void irtkHomogeneousTransformation::Reset()
{
  _matrix.Ident();
  this->Update(DOFS);
}

// -----------------------------------------------------------------------------
void irtkHomogeneousTransformation::Invert()
{
  irtkMatrix matrix = _matrix;
  _matrix  = _inverse;
  _inverse = matrix;
  this->Update(DOFS);
}

// =============================================================================
// Properties
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkHomogeneousTransformation::IsIdentity() const
{
  for (int i = 0; i < _matrix.Rows(); ++i) {
    for (int j = 0; j < _matrix.Cols(); ++j) {
      if (_matrix(i, j) != static_cast<double>(i == j)) return false;
    }
  }
  return true;
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void irtkHomogeneousTransformation::Print(irtkIndent indent) const
{
  _matrix.Print(indent);
}

// ---------------------------------------------------------------------------
irtkCifstream &irtkHomogeneousTransformation::ReadDOFs(irtkCifstream &from, irtkTransformationType format)
{
  irtkTransformation::ReadDOFs(from, format);
  this->Update(MATRIX);
  return from;
}
