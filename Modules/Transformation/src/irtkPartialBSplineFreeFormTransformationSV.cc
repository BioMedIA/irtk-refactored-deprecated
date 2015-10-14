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


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkPartialBSplineFreeFormTransformationSV
::irtkPartialBSplineFreeFormTransformationSV(irtkBSplineFreeFormTransformationSV *t, double f)
:
  _Transformation(t),
  _Fraction      (f)
{
}

// -----------------------------------------------------------------------------
irtkPartialBSplineFreeFormTransformationSV
::~irtkPartialBSplineFreeFormTransformationSV()
{
}

// =============================================================================
// Transformation parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkPartialBSplineFreeFormTransformationSV::CopyFrom(const irtkTransformation *t)
{
  const irtkPartialBSplineFreeFormTransformationSV *part = NULL;
  if ((part = dynamic_cast<const irtkPartialBSplineFreeFormTransformationSV *>(t))) {
    return _Transformation->CopyFrom(part->Transformation());
  } else {
    return _Transformation->CopyFrom(t);
  }
}

// -----------------------------------------------------------------------------
int irtkPartialBSplineFreeFormTransformationSV::NumberOfDOFs() const
{
  return _Transformation->NumberOfDOFs();
}

// -----------------------------------------------------------------------------
void irtkPartialBSplineFreeFormTransformationSV::Put(int dof, double x)
{
  _Transformation->Put(dof, x);
}

// -----------------------------------------------------------------------------
double irtkPartialBSplineFreeFormTransformationSV::Get(int dof) const
{
  return _Transformation->Get(dof);
}

// -----------------------------------------------------------------------------
void irtkPartialBSplineFreeFormTransformationSV::Put(const DOFValue *x)
{
  _Transformation->Put(x);
}

// -----------------------------------------------------------------------------
void irtkPartialBSplineFreeFormTransformationSV::Add(const DOFValue *dx)
{
  _Transformation->Add(dx);
}

// -----------------------------------------------------------------------------
double irtkPartialBSplineFreeFormTransformationSV::Update(const DOFValue *dx)
{
  return _Transformation->Update(dx);
}

// -----------------------------------------------------------------------------
void irtkPartialBSplineFreeFormTransformationSV::Get(DOFValue *x) const
{
  _Transformation->Get(x);
}

// -----------------------------------------------------------------------------
void irtkPartialBSplineFreeFormTransformationSV::PutStatus(int dof, DOFStatus s)
{
  _Transformation->PutStatus(dof, s);
}

// -----------------------------------------------------------------------------
DOFStatus irtkPartialBSplineFreeFormTransformationSV::GetStatus(int dof) const
{
  return _Transformation->GetStatus(dof);
}

// -----------------------------------------------------------------------------
bool irtkPartialBSplineFreeFormTransformationSV::HasSameDOFsAs(const irtkTransformation *t) const
{
  return HaveSameDOFs(_Transformation, t);
}

// -----------------------------------------------------------------------------
bool irtkPartialBSplineFreeFormTransformationSV::IsIdentity() const
{
  return (_Fraction == .0) || _Transformation->IsIdentity();
}

// =============================================================================
// Parameters (non-DoFs)
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkPartialBSplineFreeFormTransformationSV::Set(const char *param, const char *value)
{
  return _Transformation->Set(param, value);
}

// -----------------------------------------------------------------------------
irtkParameterList irtkPartialBSplineFreeFormTransformationSV::Parameter() const
{
  return _Transformation->Parameter();
}

// =============================================================================
// Point transformation
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkPartialBSplineFreeFormTransformationSV::RequiresCachingOfDisplacements() const
{
  return _Transformation->RequiresCachingOfDisplacements();
}

// -----------------------------------------------------------------------------
inline double irtkPartialBSplineFreeFormTransformationSV::UpperIntegrationLimit(double t, double t0) const
{
  return _Fraction * _Transformation->UpperIntegrationLimit(t, t0);
}

// -----------------------------------------------------------------------------
void irtkPartialBSplineFreeFormTransformationSV::GlobalTransform(double &x, double &y, double &z, double t, double t0) const
{
  _Transformation->GlobalTransform(x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
void irtkPartialBSplineFreeFormTransformationSV::LocalTransform(double &x, double &y, double &z, double t, double t0) const
{
  const double T = UpperIntegrationLimit(t, t0);
  if (T) _Transformation->LocalTransform(x, y, z, T, .0);
}

// -----------------------------------------------------------------------------
void irtkPartialBSplineFreeFormTransformationSV::Transform(double &x, double &y, double &z, double t, double t0) const
{
  const double T = UpperIntegrationLimit(t, t0);
  if (T) _Transformation->Transform(x, y, z, T, .0);
}

// -----------------------------------------------------------------------------
void irtkPartialBSplineFreeFormTransformationSV::GlobalInverse(double &x, double &y, double &z, double t, double t0) const
{
  _Transformation->GlobalInverse(x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
bool irtkPartialBSplineFreeFormTransformationSV::LocalInverse(double &x, double &y, double &z, double t, double t0) const
{
  const double T = UpperIntegrationLimit(t, t0);
  if (T) _Transformation->LocalInverse(x, y, z, T, .0);
  return true;
}

// -----------------------------------------------------------------------------
bool irtkPartialBSplineFreeFormTransformationSV::Inverse(double &x, double &y, double &z, double t, double t0) const
{
  const double T = UpperIntegrationLimit(t, t0);
  if (T) _Transformation->Inverse(x, y, z, T, .0);
  return true;
}

// ---------------------------------------------------------------------------
void irtkPartialBSplineFreeFormTransformationSV::Displacement(irtkGenericImage<double> &disp, double t, double t0, const irtkWorldCoordsImage *wc) const
{
  const double T = UpperIntegrationLimit(t, t0);
  if (T) _Transformation->Displacement(disp, T, .0, wc);
}

// ---------------------------------------------------------------------------
void irtkPartialBSplineFreeFormTransformationSV::Displacement(irtkGenericImage<float> &disp, double t, double t0, const irtkWorldCoordsImage *wc) const
{
  const double T = UpperIntegrationLimit(t, t0);
  if (T) _Transformation->Displacement(disp, T, .0, wc);
}

// ---------------------------------------------------------------------------
int irtkPartialBSplineFreeFormTransformationSV::InverseDisplacement(irtkGenericImage<double> &disp, double t, double t0, const irtkWorldCoordsImage *wc) const
{
  const double T = UpperIntegrationLimit(t, t0);
  if (T) _Transformation->InverseDisplacement(disp, T, .0, wc);
  return 0;
}

// ---------------------------------------------------------------------------
int irtkPartialBSplineFreeFormTransformationSV::InverseDisplacement(irtkGenericImage<float> &disp, double t, double t0, const irtkWorldCoordsImage *wc) const
{
  const double T = UpperIntegrationLimit(t, t0);
  if (T) _Transformation->InverseDisplacement(disp, T, .0, wc);
  return 0;
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
void irtkPartialBSplineFreeFormTransformationSV::GlobalJacobian(irtkMatrix &jac, double x, double y, double z, double t, double t0) const
{
  _Transformation->GlobalJacobian(jac, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
void irtkPartialBSplineFreeFormTransformationSV::LocalJacobian(irtkMatrix &jac, double x, double y, double z, double t, double t0) const
{
  const double T = UpperIntegrationLimit(t, t0);
  if (T) _Transformation->LocalJacobian(jac, x, y, z, T, .0);
}

// -----------------------------------------------------------------------------
void irtkPartialBSplineFreeFormTransformationSV::Jacobian(irtkMatrix &jac, double x, double y, double z, double t, double t0) const
{
  const double T = UpperIntegrationLimit(t, t0);
  if (T) _Transformation->Jacobian(jac, x, y, z, T, .0);
}

// -----------------------------------------------------------------------------
void irtkPartialBSplineFreeFormTransformationSV::GlobalHessian(irtkMatrix h[3], double x, double y, double z, double t, double t0) const
{
  _Transformation->GlobalHessian(h, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
void irtkPartialBSplineFreeFormTransformationSV::LocalHessian(irtkMatrix h[3], double x, double y, double z, double t, double t0) const
{
  const double T = UpperIntegrationLimit(t, t0);
  if (T) _Transformation->LocalHessian(h, x, y, z, T, .0);
}

// -----------------------------------------------------------------------------
void irtkPartialBSplineFreeFormTransformationSV::Hessian(irtkMatrix h[3], double x, double y, double z, double t, double t0) const
{
  const double T = UpperIntegrationLimit(t, t0);
  if (T) _Transformation->Hessian(h, x, y, z, T, .0);
}

// -----------------------------------------------------------------------------
void irtkPartialBSplineFreeFormTransformationSV::JacobianDOFs(double jac[3], int dof, double x, double y, double z, double t, double t0) const
{
  const double T = UpperIntegrationLimit(t, t0);
  if (T) _Transformation->JacobianDOFs(jac, dof, x, y, z, T, .0);
}

// -----------------------------------------------------------------------------
void irtkPartialBSplineFreeFormTransformationSV
::ParametricGradient(const irtkGenericImage<double> *in, double *out,
                     const irtkWorldCoordsImage *i2w, const irtkWorldCoordsImage *wc,
                     double t0, double w) const
{
  const double T = UpperIntegrationLimit(in->GetTOrigin(), t0);
  if (T) _Transformation->ParametricGradient(in, out, i2w, wc, T, .0, w);
}

// -----------------------------------------------------------------------------
void irtkPartialBSplineFreeFormTransformationSV
::ParametricGradient(const irtkGenericImage<double> **in, int n, double *out,
                     const irtkWorldCoordsImage *i2w, const irtkWorldCoordsImage *wc,
                     const double *t0, double w) const
{
  double t, T;
  for (int i = 0; i < n; ++i) {
    t = in[i]->GetTOrigin();
    T = UpperIntegrationLimit(t, t0 ? t0[i] : t);
    if (T) _Transformation->ParametricGradient(in[i], out, i2w, wc, T, .0, w);
  }
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void irtkPartialBSplineFreeFormTransformationSV::Print(irtkIndent indent) const
{
  cout << indent << "Partial B-spline SV FFD" << endl;
  _Transformation->Print(indent+1);
  cout << indent << "Fraction of transformation applied: " << _Fraction << endl;
}

// -----------------------------------------------------------------------------
irtkCifstream &irtkPartialBSplineFreeFormTransformationSV::Read(irtkCifstream &from)
{
  return _Transformation->Read(from);
}

// -----------------------------------------------------------------------------
irtkCofstream &irtkPartialBSplineFreeFormTransformationSV::Write(irtkCofstream &to) const
{
  return _Transformation->Write(to);
}
