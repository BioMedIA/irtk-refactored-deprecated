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
#include <irtkDisplacementToVelocityField.h>
#include <irtkImageToInterpolationCoefficients.h>

#include "irtkFreeFormTransformationIntegration.h"


// =============================================================================
// Integration methods
// =============================================================================

IRTK_FFDIM2(RKE1,   irtkBSplineFreeFormTransformationTD);
IRTK_FFDIM2(RKEH12, irtkBSplineFreeFormTransformationTD);
IRTK_FFDIM2(RKE2,   irtkBSplineFreeFormTransformationTD);
IRTK_FFDIM2(RKH2,   irtkBSplineFreeFormTransformationTD);
IRTK_FFDIM2(RKBS23, irtkBSplineFreeFormTransformationTD);
IRTK_FFDIM2(RK4,    irtkBSplineFreeFormTransformationTD);
IRTK_FFDIM2(RKF45,  irtkBSplineFreeFormTransformationTD);
IRTK_FFDIM2(RKCK45, irtkBSplineFreeFormTransformationTD);
IRTK_FFDIM2(RKDP45, irtkBSplineFreeFormTransformationTD);

// =============================================================================
// Forward declaration of CUDA accelerated functions
// =============================================================================
#ifdef USE_CUDA

// -----------------------------------------------------------------------------
namespace irtkCUBSplineFreeFormTransformationTD {
  template <class VoxelType>
  void Displacement(const irtkBSplineFreeFormTransformationTD &h_ffd,
                    irtkGenericImage<VoxelType>               &h_disp,
                    FFDIntegrationMethod                       im,
                    double t1, double t2, double mindt, double maxdt, double tol);
}

#endif
// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkBSplineFreeFormTransformationTD::irtkBSplineFreeFormTransformationTD()
:
  _IntegrationMethod(FFDIM_RKE2),
  _MinTimeStep      (0.01),
  _MaxTimeStep      (0.1),
  _Tolerance        (1.e-3)
{
  _ExtrapolationMode = Extrapolation_NN;
}

// -----------------------------------------------------------------------------
irtkBSplineFreeFormTransformationTD
::irtkBSplineFreeFormTransformationTD(const irtkImageAttributes &attr,
                                      double dx, double dy, double dz, double dt)
:
  _IntegrationMethod(FFDIM_RKE2),
  _MinTimeStep      (0.01),
  _MaxTimeStep      (0.1),
  _Tolerance        (1.e-3)
{
  _ExtrapolationMode = Extrapolation_NN;
  Initialize(attr, dx, dy, dz, dt);
}

// -----------------------------------------------------------------------------
irtkBSplineFreeFormTransformationTD
::irtkBSplineFreeFormTransformationTD(const irtkBaseImage &target,
                                      double dx, double dy, double dz, double dt)
:
  _IntegrationMethod(FFDIM_RKE2),
  _MinTimeStep      (0.01),
  _MaxTimeStep      (0.1),
  _Tolerance        (1.e-3)
{
  _ExtrapolationMode = Extrapolation_NN;
  Initialize(target.Attributes(), dx, dy, dz, dt);
}

// -----------------------------------------------------------------------------
irtkBSplineFreeFormTransformationTD
::irtkBSplineFreeFormTransformationTD(const irtkBSplineFreeFormTransformationTD &ffd)
:
  irtkBSplineFreeFormTransformation4D(ffd),
  _IntegrationMethod(ffd._IntegrationMethod),
  _MinTimeStep      (ffd._MinTimeStep),
  _MaxTimeStep      (ffd._MaxTimeStep),
  _Tolerance        (ffd._Tolerance)
{
}

// -----------------------------------------------------------------------------
irtkBSplineFreeFormTransformationTD
::~irtkBSplineFreeFormTransformationTD()
{
}

// =============================================================================
// Approximation/Interpolation
// =============================================================================

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationTD
::ApproximateDOFs(const irtkGenericImage<double> * const *disp,
                  const double *t1, const double *t2, int no,
                  bool smooth, int nterms, int niter)
{
  // TODO: Refactor this method taking advantage of recent changes related
  //       to velocity field based transformation parameterization.
  //       - Use Vector as voxel type for 4D vector field image.
  //       - Use fast inter-/extrapolators which also handle vector images.
  //       - Possibly use cubic B-spline interpolator for BCH filter,
  //         similar to the irtkBSplineFreeFormTransformationSV::Add method.
  irtkInterpolateImageFunction *f;

  irtkDisplacementToVelocityFieldBCH<double> dtov;
  dtov.SetSmoothVelocities  (smooth);
  dtov.SetNumberOfTerms     (nterms);
  dtov.SetNumberOfIterations(niter);

  irtkGenericImage<double> vx(_attr);
  irtkGenericImage<double> vy(_attr);
  irtkGenericImage<double> vz(_attr);

  double *vw = new double[_t];
  memset(vw, 0, _t * sizeof(double));

  irtkGenericImage<double> d(_attr, 3);
  irtkGenericImage<double> v(_attr, 3);
 
  for (int n = 0; n < no; n++) {
    // Get control point frames covered by n-th displacement field
    int l1 = static_cast<int>(floor(this->TimeToLattice(t1[n])));
    int l2 = static_cast<int>(ceil (this->TimeToLattice(t2[n])));

    if (l1 > l2) {
      int l = l1;
      l1 = l2;
      l2 = l;
    }

    // Check if within temporal control point domain
    if ((l1 < 0 || l1 >= _t) && (l2 < 0 || l2 >= _t)) continue;

    // Sample displacement field at control points using linear interpolation
    if (disp[n]->GetImageAttributes() == _attr) {
      for (int k = 0; k < _z; ++k) {
        for (int j = 0; j < _y; ++j) {
          for (int i = 0; i < _x; ++i) {
            d.Put(i, j, k, 0, disp[n]->Get(i, j, k, 0));
            d.Put(i, j, k, 1, disp[n]->Get(i, j, k, 1));
            d.Put(i, j, k, 2, disp[n]->Get(i, j, k, 2));
          }
        }
      }
    } else {
      double x, y, z, vec[3] = {.0, .0, .0};

      irtkGenericImage<double> *input = const_cast<irtkGenericImage<double> *>(disp[n]);
      f = irtkInterpolateImageFunction::New(Interpolation_Linear,
                                            Extrapolation_NN,
                                            input);
      f->SetInput(input);
      f->Initialize();

      for (int k = 0; k < _z; k++) {
        for (int j = 0; j < _y; j++) {
          for (int i = 0; i < _x; i++) {
            x = i, y = j, z = k;
            d.ImageToWorld(x, y, z);
            disp[n]->WorldToImage(x, y, z);
            f->Evaluate(vec, x, y, z);
            d.Put(i, j, k, 0, vec[0]);
            d.Put(i, j, k, 1, vec[1]);
            d.Put(i, j, k, 2, vec[2]);
          }
        }
      }
    }

    // Compute stationary velocity field
    const double T = t2[n] - t1[n];

    dtov.SetNumberOfSteps(static_cast<int>(round(T / _MinTimeStep)));
    dtov.SetInput (&d);
    dtov.SetOutput(&v);
    dtov.Run();

    // Add stationary velocities to frames in [l1, l2]
    for (int l = l1; l <= l2; l++) {
      if (0 <= l && l < _t) {
        for (int k = 0; k < _z; k++) {
          for (int j = 0; j < _y; j++) {
            for (int i = 0; i < _x; i++) {
              vx(i, j, k, l) += v(i, j, k, 0);
              vy(i, j, k, l) += v(i, j, k, 1);
              vz(i, j, k, l) += v(i, j, k, 2);
            }
          }
        }
        vw[l] += T;
      }
    }
  }
 
  // Compute average velocities at control points
  for (int l = 0; l < _t; l++) {
    if (vw[l] == .0) vw[l] = 1.0;
    for (int k = 0; k < _z; k++) {
      for (int j = 0; j < _y; j++) {
        for (int i = 0; i < _x; i++) {
          _CPImage(i, j, k, l) = vx(i, j, k, l) / vw[l];
          _CPImage(i, j, k, l) = vy(i, j, k, l) / vw[l];
          _CPImage(i, j, k, l) = vz(i, j, k, l) / vw[l];
        }
      }
    }
  }

  // Convert to B-spline coefficients
  ConvertToSplineCoefficients(3, _CPImage);

  // Clean up
  delete[] vw;
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationTD
::ApproximateDOFs(const double *, const double *, const double *, const double *,
                  const double *, const double *, const double *, int)
{
  cerr << this->NameOfClass() << "::ApproximateDOFs: Not implemented" << endl;
  cerr << "  --> Try the other overloaded approximation function" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationTD
::ApproximateDOFsGradient(const double *, const double *, const double *, const double *,
                          const double *, const double *, const double *, int,
                          double *, double) const
{
  cerr << this->NameOfClass() << "::ApproximateDOFsGradient: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
double irtkBSplineFreeFormTransformationTD
::ApproximateAsNew(irtkGenericImage<double> **disp,
                   const double *t1, const double *t2, int no,
                   bool smooth, int nterms, int niter)
{
  // Approximate 3D+t velocity field
  this->ApproximateDOFs(disp, t1, t2, no, smooth, nterms, niter);

  // Evaluate RMS of approximation error
  double error  = .0;
  int    ntotal = 0;
  double dx, dy, dz;

  for (int n = 0; n < no; ++n) {
    irtkGenericImage<double> &d = *disp[n];
    for (int k = 0; k < d.Z(); ++k) {
      for (int j = 0; j < d.Y(); ++j) {
        for (int i = 0; i < d.X(); ++i) {
          dx = i, dy = j, dz = k;
          d.ImageToWorld(dx, dy, dz);
          this->Displacement(dx, dy, dz, t1[n], t2[n]);
          d(i, j, k, 0) -= dx;
          d(i, j, k, 1) -= dy;
          d(i, j, k, 2) -= dz;
          error += sqrt(d(i, j, k, 0) * d(i, j, k, 0) +
                        d(i, j, k, 1) * d(i, j, k, 1) +
                        d(i, j, k, 2) * d(i, j, k, 2));
          ++ntotal;
        }
      }
    }
  }
  if (ntotal > 0) error /= ntotal;
  return error;
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationTD::Interpolate(const double *, const double *, const double *)
{
  cerr << this->NameOfClass() << "::Interpolate: Not implemented" << endl;
  exit(1);
}

// =============================================================================
// Parameters (non-DoFs)
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkBSplineFreeFormTransformationTD::Set(const char *name, const char *value)
{

  if (strcmp(name, "Integration method")          == 0 ||
      strcmp(name, "Velocity integration method") == 0) {
    return FromString(value, _IntegrationMethod);
  } else if (strcmp(name, "Length of integration steps") == 0) {
    double dt;
    if (!FromString(value, dt) || dt <= .0) return false;
    _MinTimeStep = _MaxTimeStep = dt;
    return true;
  } else if (strcmp(name, "Minimum length of integration steps") == 0) {
    return FromString(value, _MinTimeStep) && _MinTimeStep > .0;
  } else if (strcmp(name, "Maximum length of integration steps") == 0) {
    return FromString(value, _MaxTimeStep) && _MaxTimeStep > .0;
  } else if (strcmp(name, "Integration tolerance")          == 0 ||
             strcmp(name, "Velocity integration tolerance") == 0) {
    return FromString(value, _Tolerance) && _Tolerance >= .0;
  }
  return irtkBSplineFreeFormTransformation4D::Set(name, value);
}

// -----------------------------------------------------------------------------
irtkParameterList irtkBSplineFreeFormTransformationTD::Parameter() const
{
  irtkParameterList params = irtkBSplineFreeFormTransformation4D::Parameter();
  Insert(params, "Integration method",                  ToString(_IntegrationMethod));
  Insert(params, "Minimum length of integration steps", ToString(_MinTimeStep));
  Insert(params, "Maximum length of integration steps", ToString(_MaxTimeStep));
  Insert(params, "Integration tolerance",               ToString(_Tolerance));
  return params;
}

// =============================================================================
// Point transformation
// =============================================================================

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationTD::LocalTransform(double &x, double &y, double &z, double t, double t0) const
{
  if      (_IntegrationMethod == FFDIM_RKE1)   RKE1  ::Transform(this, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RKE2)   RKE2  ::Transform(this, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RKH2)   RKH2  ::Transform(this, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RK4)    RK4   ::Transform(this, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RKEH12) RKEH12::Transform(this, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKBS23) RKBS23::Transform(this, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKF45)  RKF45 ::Transform(this, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKCK45) RKCK45::Transform(this, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKDP45) RKDP45::Transform(this, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else {
    cerr << "irtkBSplineFreeFormTransformationTD::LocalTransform: Unknown integration method: " << _IntegrationMethod << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationTD::TransformAndJacobian(irtkMatrix &jac, double &x, double &y, double &z, double t, double t0) const
{
  if      (_IntegrationMethod == FFDIM_RKE1)   RKE1  ::Jacobian(this, jac, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RKE2)   RKE2  ::Jacobian(this, jac, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RKH2)   RKH2  ::Jacobian(this, jac, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RK4)    RK4   ::Jacobian(this, jac, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RKEH12) RKEH12::Jacobian(this, jac, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKBS23) RKBS23::Jacobian(this, jac, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKF45)  RKF45 ::Jacobian(this, jac, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKCK45) RKCK45::Jacobian(this, jac, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKDP45) RKDP45::Jacobian(this, jac, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else {
    cerr << "irtkBSplineFreeFormTransformationTD::TransformAndJacobian: Unknown integration method: " << _IntegrationMethod << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationTD::TransformAndJacobianDOFs(irtkMatrix &jac, int ci, int cj, int ck, int cl, double &x, double &y, double &z, double t, double t0) const
{
  if      (_IntegrationMethod == FFDIM_RKE1)   RKE1  ::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RKE2)   RKE2  ::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RKH2)   RKH2  ::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RK4)    RK4   ::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RKEH12) RKEH12::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKBS23) RKBS23::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKF45)  RKF45 ::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKCK45) RKCK45::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKDP45) RKDP45::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else {
    cerr << "irtkBSplineFreeFormTransformationTD::TransformAndJacobianDOFs: Unknown integration method: " << _IntegrationMethod << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationTD::TransformAndJacobianDOFs(irtkMatrix &jac, int cp, double &x, double &y, double &z, double t, double t0) const
{
  int ci, cj, ck, cl;
  this->IndexToLattice(cp, ci, cj, ck, cl);
  this->TransformAndJacobianDOFs(jac, ci, cj, ck, cl, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationTD::TransformAndJacobianDOFs(irtkTransformationJacobian &jac, double &x, double &y, double &z, double t, double t0) const
{
  if      (_IntegrationMethod == FFDIM_RKE1)   RKE1  ::JacobianDOFs(this, jac, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RKE2)   RKE2  ::JacobianDOFs(this, jac, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RKH2)   RKH2  ::JacobianDOFs(this, jac, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RK4)    RK4   ::JacobianDOFs(this, jac, x, y, z, t0, t, _MinTimeStep);
  else if (_IntegrationMethod == FFDIM_RKEH12) RKEH12::JacobianDOFs(this, jac, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKBS23) RKBS23::JacobianDOFs(this, jac, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKF45)  RKF45 ::JacobianDOFs(this, jac, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKCK45) RKCK45::JacobianDOFs(this, jac, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else if (_IntegrationMethod == FFDIM_RKDP45) RKDP45::JacobianDOFs(this, jac, x, y, z, t0, t, _MinTimeStep, _MaxTimeStep, _Tolerance);
  else {
    cerr << "irtkBSplineFreeFormTransformationTD::TransformAndJacobianDOFs: Unknown integration method: " << _IntegrationMethod << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationTD::Displacement(irtkGenericImage<double> &disp, double t, double t0, const irtkWorldCoordsImage *i2w) const
{
  if (disp.GetT() < 2 || disp.GetT() > 3) {
    cerr << this->NameOfClass() << "::Displacement: Input/output image must have either 2 or 3 vector components (_t)" << endl;
    exit(1);
  }

#ifdef USE_CUDA
  if (use_gpu) {
    irtkCUBSplineFreeFormTransformationTD
        ::Displacement(*this, disp, _IntegrationMethod, t0, t,
                                    _MinTimeStep, _MaxTimeStep, _Tolerance);
  } else
#endif
  {
    // Use (multi-threaded) superclass implementation instead
    irtkBSplineFreeFormTransformation4D::Displacement(disp, t, t0, i2w);
  }
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationTD::Displacement(irtkGenericImage<float> &disp, double t, double t0, const irtkWorldCoordsImage *i2w) const
{
  if (disp.GetT() < 2 || disp.GetT() > 3) {
    cerr << this->NameOfClass() << "::Displacement: Input/output image must have either 2 or 3 vector components (_t)" << endl;
    exit(1);
  }

#ifdef USE_CUDA
  if (use_gpu) {
    irtkCUBSplineFreeFormTransformationTD
        ::Displacement(*this, disp, _IntegrationMethod, t0, t,
                                    _MinTimeStep, _MaxTimeStep, _Tolerance);
  } else
#endif
  {
    // Use (multi-threaded) superclass implementation instead
    irtkBSplineFreeFormTransformation4D::Displacement(disp, t, t0, i2w);
  }
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationTD
::ParametricGradient(const irtkPointSet &pos, const irtkVector3D<double> *in,
                     double *out, double t, double t0, double w) const
{
  // Cannot use irtk(BSpline)FreeFormTransformation3D::ParametricGradient
  // computation which assumes a local support region of the interpolation kernel
  // because point trajectories may enter support regions of more than one kernel
  irtkFreeFormTransformation::ParametricGradient(pos, in, out, t, t0, w);
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationTD::Print(irtkIndent indent) const
{
  cout << indent << "B-spline TD FFD:" << endl;
  ++indent;
  irtkFreeFormTransformation4D::Print(indent);
  cout << indent << "Numerical integration:" << endl;
  ++indent;
  cout << indent << "Method:      " << ToString(_IntegrationMethod) << endl;
  if (_MinTimeStep < _MaxTimeStep && _IntegrationMethod >= FFDIM_RKEH12) {
    cout << indent << "Step size:   [" << _MinTimeStep << ", " << _MaxTimeStep << "]" << endl;
    cout << indent << "Local error: " << _Tolerance << endl;
  } else {
    cout << indent << "Step size:   " << _MinTimeStep << endl;  
  }
}

// -----------------------------------------------------------------------------
bool irtkBSplineFreeFormTransformationTD::CanRead(irtkTransformationType format) const
{
  switch (format) {
    case IRTKTRANSFORMATION_BSPLINE_FFD_TD_v1:
    case IRTKTRANSFORMATION_BSPLINE_FFD_TD_v2:
    case IRTKTRANSFORMATION_BSPLINE_FFD_TD_v3:
    case IRTKTRANSFORMATION_BSPLINE_FFD_TD_v4:
      return true;
    default:
      return false;
  }
}

// -----------------------------------------------------------------------------
irtkCifstream &irtkBSplineFreeFormTransformationTD::ReadDOFs(irtkCifstream &from, irtkTransformationType format)
{
  // Read FFD data
  switch (format) {
    case IRTKTRANSFORMATION_BSPLINE_FFD_TD_v1:
    case IRTKTRANSFORMATION_BSPLINE_FFD_TD_v2:
      irtkBSplineFreeFormTransformation4D::ReadDOFs(from, IRTKTRANSFORMATION_BSPLINE_FFD_4D_v1);
      break;
    case IRTKTRANSFORMATION_BSPLINE_FFD_TD_v3:
      irtkBSplineFreeFormTransformation4D::ReadDOFs(from, IRTKTRANSFORMATION_BSPLINE_FFD_4D_v2);
      break;
    default:
      irtkBSplineFreeFormTransformation4D::ReadDOFs(from, IRTKTRANSFORMATION_BSPLINE_FFD_4D);
  }

  // Read integration method
  if (format == IRTKTRANSFORMATION_BSPLINE_FFD_4D_v1) {
    _IntegrationMethod = FFDIM_RKE1;
  } else {
    unsigned int tmp = FFDIM_UNKNOWN;
    if (!from.ReadAsUInt(&tmp, 1)) {
      cerr << "irtkBSplineFreeFormTransformationTD::ReadDOFs: Failed to read integration method" << endl;
      exit(1);
    }
    _IntegrationMethod = static_cast<FFDIntegrationMethod>(tmp);
  }

  // Read minimum/maximum time step length
  if (!from.ReadAsDouble(&_MinTimeStep, 1) || !from.ReadAsDouble(&_MaxTimeStep, 1)) {
    cerr << "irtkBSplineFreeFormTransformationTD::ReadDOFs: Failed to read min/max time step" << endl;
    exit(1);
  }

  // Read local error tolerance
  if (format == IRTKTRANSFORMATION_BSPLINE_FFD_TD_v1) {
    _Tolerance = 1.e-3;
  } else {
    if (!from.ReadAsDouble(&_Tolerance, 1)) {
      cerr << "irtkBSplineFreeFormTransformationTD::ReadDOFs: Failed to read integration tolerance" << endl;
      exit(1);
    }
  }

  return from;
}

// -----------------------------------------------------------------------------
irtkCofstream &irtkBSplineFreeFormTransformationTD::WriteDOFs(irtkCofstream &to) const
{
  // Write FFD data
  irtkBSplineFreeFormTransformation4D::WriteDOFs(to);

  // Write integration method
  unsigned int tmp = _IntegrationMethod;
  to.WriteAsUInt(&tmp, 1);

  // Write minimum/maximum time step length
  to.WriteAsDouble(&_MinTimeStep, 1);
  to.WriteAsDouble(&_MaxTimeStep, 1);

  // Write local error tolerance
  to.WriteAsDouble(&_Tolerance, 1);

  return to;
}

// =============================================================================
// Others
// =============================================================================

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationTD::Invert()
{
  Vector tmp;
  for (int l = 0; l < _t / 2; ++l) {
    int L = _t - l - 1;
    for (int k = 0; k < _z; ++k) {
      for (int j = 0; j < _y; ++j) {
        for (int i = 0; i < _x; ++i) {
          tmp = _CPImage(i, j, k, l);
          _CPImage(i, j, k, l) = - _CPImage(i, j, k, L);
          _CPImage(i, j, k, L) = - tmp;
        }
      }
    }
  }
}
