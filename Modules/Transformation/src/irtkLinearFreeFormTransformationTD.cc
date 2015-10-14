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


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkLinearFreeFormTransformationTD::irtkLinearFreeFormTransformationTD()
:
  _MinTimeStep(0.01),
  _MaxTimeStep(0.5)
{
  _ExtrapolationMode = Extrapolation_NN;
}

// -----------------------------------------------------------------------------
irtkLinearFreeFormTransformationTD
::irtkLinearFreeFormTransformationTD(const irtkImageAttributes &attr,
                                     double dx, double dy, double dz, double dt)
:
  _MinTimeStep(0.01),
  _MaxTimeStep(0.5)
{
  _ExtrapolationMode = Extrapolation_NN;
  Initialize(attr, dx, dy, dz, dt);
}

// -----------------------------------------------------------------------------
irtkLinearFreeFormTransformationTD
::irtkLinearFreeFormTransformationTD(const irtkBaseImage &target,
                                     double dx, double dy, double dz, double dt)
:
_MinTimeStep(0.01),
_MaxTimeStep(0.5)
{
  _ExtrapolationMode = Extrapolation_NN;
  Initialize(target.Attributes(), dx, dy, dz, dt);
}

// -----------------------------------------------------------------------------
irtkLinearFreeFormTransformationTD
::irtkLinearFreeFormTransformationTD(const irtkBSplineFreeFormTransformationTD &t)
:
  irtkLinearFreeFormTransformation4D(t),
  _MinTimeStep(t.MinTimeStep()),
  _MaxTimeStep(t.MaxTimeStep())
{
}

// -----------------------------------------------------------------------------
irtkLinearFreeFormTransformationTD
::irtkLinearFreeFormTransformationTD(const irtkLinearFreeFormTransformationTD &t)
:
  irtkLinearFreeFormTransformation4D(t),
  _MinTimeStep(t._MinTimeStep),
  _MaxTimeStep(t._MaxTimeStep)
{
}

// -----------------------------------------------------------------------------
irtkLinearFreeFormTransformationTD
::~irtkLinearFreeFormTransformationTD()
{
}

// =============================================================================
// Approximation/Interpolation
// =============================================================================

// -----------------------------------------------------------------------------
void irtkLinearFreeFormTransformationTD
::ApproximateDOFs(const irtkGenericImage<double> * const *disp,
                  const double *t1, const double *t2, int no,
                  bool smooth, int nterms, int niter)
{
  if (no == 0) return;
 
  // Get last displacement field
  int last = 0;
  for (int n = 1; n < no; n++) {
    if (t2[n] > t2[last]) last = n;
  }

  // Determine number of 3D+t samples
  int npoints = 0;
  for (int n = 0; n < no; n++) {
    npoints += disp[n]->GetX() * disp[n]->GetY() * disp[n]->GetZ();
  }

  // Allocate memory for computed stationary velocity fields
  double *px = Allocate<double>(npoints);
  double *py = Allocate<double>(npoints);
  double *pz = Allocate<double>(npoints);
  double *pt = Allocate<double>(npoints);
  double *vx = Allocate<double>(npoints);
  double *vy = Allocate<double>(npoints);
  double *vz = Allocate<double>(npoints);

  // Compute stationary velocity fields from displacement fields
  // and store them in 3D+t vector field with associated time t
  irtkDisplacementToVelocityFieldBCH<double> dtov;
  dtov.SetSmoothVelocities  (smooth);
  dtov.SetNumberOfTerms     (nterms);
  dtov.SetNumberOfIterations(niter);

  int p = 0;
  for (int n = 0; n < no; n++) {
    irtkGenericImage<double> v;
    dtov.SetT(t2[n] - t1[n]);
    dtov.SetNumberOfSteps(static_cast<int>(round(dtov.GetT() / _MinTimeStep)));
    dtov.SetInput(const_cast<irtkGenericImage<double> *>(disp[n]));
    dtov.SetOutput(&v);
    dtov.Run();
    for (int k = 0; k < v.GetZ(); k++) {
      for (int j = 0; j < v.GetY(); j++) {
        for (int i = 0; i < v.GetX(); i++) {
          px[p] = i;
          py[p] = j;
          pz[p] = k;
          v.ImageToWorld(px[p], py[p], pz[p]);
          pt[p] = t1[n];
          vx[p] = v.Get(i, j, k, 0);
          vy[p] = v.Get(i, j, k, 1);
          vz[p] = v.Get(i, j, k, 2);
          p++;
        }
        if (n == last) {
          px[p] = px[p - 1];
          py[p] = py[p - 1];
          pz[p] = pz[p - 1];
          pt[p] = t2[last ];
          vx[p] = vx[p - 1];
          vy[p] = vy[p - 1];
          vz[p] = vz[p - 1];
          p++;
        }
      }
    }
  }

  // Approximate 3D+t velocity field by linear FFD
  irtkLinearFreeFormTransformation4D::ApproximateDOFs(px, py, pz, pt, vx, vy, vz, npoints);

  // Free memory
  Deallocate(px);
  Deallocate(py);
  Deallocate(pz);
  Deallocate(pt);
  Deallocate(vx);
  Deallocate(vy);
  Deallocate(vz);
}

// -----------------------------------------------------------------------------
void irtkLinearFreeFormTransformationTD
::ApproximateDOFs(const double *, const double *, const double *, const double *,
                  const double *, const double *, const double *, int)
{
  cerr << this->NameOfClass() << "::ApproximateDOFs: Not implemented" << endl;
  cerr << "  --> Try the other overloaded approximation function" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
void irtkLinearFreeFormTransformationTD
::ApproximateDOFsGradient(const double *, const double *, const double *, const double *,
                          const double *, const double *, const double *, int,
                          double *, double) const
{
  cerr << this->NameOfClass() << "::ApproximateDOFsGradient: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
double irtkLinearFreeFormTransformationTD
::ApproximateAsNew(irtkGenericImage<double> **disp,
                   const double *t1, const double *t2, int no,
                   bool smooth, int nterms, int niter)
{
  // Approximate 3D+t velocity field
  this->ApproximateDOFs(disp, t1, t2, no, smooth, nterms, niter);

  // Evaluate RMS of approximation error
  double error  = 0;
  int    ntotal = 0;
  double dx, dy, dz;

  for (int n = 0; n < no; n++) {
    irtkGenericImage<double> &d = *disp[n];
    const int X = d.GetX();
    const int Y = d.GetY();
    const int Z = d.GetZ();
    for (int k = 0; k < Z; k++) {
      for (int j = 0; j < Y; j++) {
        for (int i = 0; i < X; i++) {
          dx = i;
          dy = j;
          dz = k;
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
void irtkLinearFreeFormTransformationTD::Interpolate(const double *, const double *, const double *)
{
  cerr << this->NameOfClass() << "::Interpolate: Not implemented" << endl;
  exit(1);
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void irtkLinearFreeFormTransformationTD::Print(irtkIndent indent) const
{
  cout << indent << "Linear TD FFD:" << endl;
  ++indent;
  irtkFreeFormTransformation4D::Print(indent);
  cout << indent << "Minimum length of steps: " << _MinTimeStep << endl;
  cout << indent << "Maximum length of steps: " << _MaxTimeStep << endl;
}

// -----------------------------------------------------------------------------
bool irtkLinearFreeFormTransformationTD::CanRead(irtkTransformationType format) const
{
  switch (format) {
    case IRTKTRANSFORMATION_LINEAR_FFD_TD_v1:
    case IRTKTRANSFORMATION_LINEAR_FFD_TD_v2:
    case IRTKTRANSFORMATION_LINEAR_FFD_TD_v3:
      return true;
    default:
      return false;
  }
}

// -----------------------------------------------------------------------------
irtkCifstream &irtkLinearFreeFormTransformationTD::ReadDOFs(irtkCifstream &from, irtkTransformationType format)
{
  // Read FFD data
  if (format < IRTKTRANSFORMATION_LINEAR_FFD_TD_v2) {
    irtkLinearFreeFormTransformation4D::ReadDOFs(from, IRTKTRANSFORMATION_LINEAR_FFD_4D_v1);
  } else {
    irtkLinearFreeFormTransformation4D::ReadDOFs(from, IRTKTRANSFORMATION_LINEAR_FFD_4D_v2);
  }

  // Read minimum/maximum time step length
  from.ReadAsDouble(&_MinTimeStep, 1);
  from.ReadAsDouble(&_MaxTimeStep, 1);

  return from;
}

// -----------------------------------------------------------------------------
irtkCofstream &irtkLinearFreeFormTransformationTD::WriteDOFs(irtkCofstream &to) const
{
  // Write FFD data
  irtkLinearFreeFormTransformation4D::WriteDOFs(to);

  // Write minimum/maximum time step length
  to.WriteAsDouble(&_MinTimeStep, 1);
  to.WriteAsDouble(&_MaxTimeStep, 1);

  return to;
}
