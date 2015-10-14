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
irtkPartialMultiLevelStationaryVelocityTransformation
::irtkPartialMultiLevelStationaryVelocityTransformation(irtkMultiLevelStationaryVelocityTransformation *t, double f)
:
  _Transformation(t),
  _Fraction      (f)
{
}

// -----------------------------------------------------------------------------
irtkPartialMultiLevelStationaryVelocityTransformation
::~irtkPartialMultiLevelStationaryVelocityTransformation()
{
}

// =============================================================================
// Levels
// =============================================================================

// -----------------------------------------------------------------------------
int irtkPartialMultiLevelStationaryVelocityTransformation::NumberOfLevels() const
{
  return _Transformation->NumberOfLevels();
}

// -----------------------------------------------------------------------------
irtkAffineTransformation *irtkPartialMultiLevelStationaryVelocityTransformation::GetGlobalTransformation()
{
  return _Transformation->GetGlobalTransformation();
}

// -----------------------------------------------------------------------------
const irtkAffineTransformation *irtkPartialMultiLevelStationaryVelocityTransformation::GetGlobalTransformation() const
{
  return _Transformation->GetGlobalTransformation();
}

// -----------------------------------------------------------------------------
irtkFreeFormTransformation *irtkPartialMultiLevelStationaryVelocityTransformation::GetLocalTransformation(int pos)
{
  return _Transformation->GetLocalTransformation(pos);
}

// -----------------------------------------------------------------------------
const irtkFreeFormTransformation *irtkPartialMultiLevelStationaryVelocityTransformation::GetLocalTransformation(int pos) const
{
  return _Transformation->GetLocalTransformation(pos);
}

// ---------------------------------------------------------------------------
irtkFreeFormTransformation *irtkPartialMultiLevelStationaryVelocityTransformation::PutLocalTransformation(irtkFreeFormTransformation *transformation, int pos)
{
  return _Transformation->PutLocalTransformation(transformation, pos);
}

// ---------------------------------------------------------------------------
void irtkPartialMultiLevelStationaryVelocityTransformation::PushLocalTransformation(irtkFreeFormTransformation *transformation)
{
  _Transformation->PushLocalTransformation(transformation);
}

// ---------------------------------------------------------------------------
void irtkPartialMultiLevelStationaryVelocityTransformation::InsertLocalTransformation(irtkFreeFormTransformation *transformation, int pos)
{
  _Transformation->InsertLocalTransformation(transformation);
}

// ---------------------------------------------------------------------------
irtkFreeFormTransformation *irtkPartialMultiLevelStationaryVelocityTransformation::PopLocalTransformation()
{
  return _Transformation->PopLocalTransformation();
}

// ---------------------------------------------------------------------------
irtkFreeFormTransformation *irtkPartialMultiLevelStationaryVelocityTransformation::RemoveLocalTransformation(int pos)
{
  return _Transformation->RemoveLocalTransformation(pos);
}

// ---------------------------------------------------------------------------
void irtkPartialMultiLevelStationaryVelocityTransformation::CombineLocalTransformation()
{
  _Transformation->CombineLocalTransformation();
}

// ---------------------------------------------------------------------------
void irtkPartialMultiLevelStationaryVelocityTransformation::MergeGlobalIntoLocalDisplacement()
{
  _Transformation->MergeGlobalIntoLocalDisplacement();
}

// =============================================================================
// Transformation parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkPartialMultiLevelStationaryVelocityTransformation::CopyFrom(const irtkTransformation *t)
{
  const irtkPartialMultiLevelStationaryVelocityTransformation *part = NULL;
  if ((part = dynamic_cast<const irtkPartialMultiLevelStationaryVelocityTransformation *>(t))) {
    return _Transformation->CopyFrom(part->Transformation());
  } else {
    return _Transformation->CopyFrom(t);
  }
}

// -----------------------------------------------------------------------------
int irtkPartialMultiLevelStationaryVelocityTransformation::NumberOfDOFs() const
{
  return _Transformation->NumberOfDOFs();
}

// -----------------------------------------------------------------------------
void irtkPartialMultiLevelStationaryVelocityTransformation::Put(int dof, double x)
{
  _Transformation->Put(dof, x);
}

// -----------------------------------------------------------------------------
double irtkPartialMultiLevelStationaryVelocityTransformation::Get(int dof) const
{
  return _Transformation->Get(dof);
}

// -----------------------------------------------------------------------------
void irtkPartialMultiLevelStationaryVelocityTransformation::Put(const DOFValue *x)
{
  _Transformation->Put(x);
}

// -----------------------------------------------------------------------------
void irtkPartialMultiLevelStationaryVelocityTransformation::Add(const DOFValue *dx)
{
  _Transformation->Add(dx);
}

// -----------------------------------------------------------------------------
double irtkPartialMultiLevelStationaryVelocityTransformation::Update(const DOFValue *dx)
{
  return _Transformation->Update(dx);
}

// -----------------------------------------------------------------------------
void irtkPartialMultiLevelStationaryVelocityTransformation::Get(DOFValue *x) const
{
  _Transformation->Get(x);
}

// -----------------------------------------------------------------------------
void irtkPartialMultiLevelStationaryVelocityTransformation::PutStatus(int dof, DOFStatus s)
{
  _Transformation->PutStatus(dof, s);
}

// -----------------------------------------------------------------------------
DOFStatus irtkPartialMultiLevelStationaryVelocityTransformation::GetStatus(int dof) const
{
  return _Transformation->GetStatus(dof);
}

// -----------------------------------------------------------------------------
bool irtkPartialMultiLevelStationaryVelocityTransformation::HasSameDOFsAs(const irtkTransformation *t) const
{
  return HaveSameDOFs(_Transformation, t);
}

// -----------------------------------------------------------------------------
bool irtkPartialMultiLevelStationaryVelocityTransformation::IsIdentity() const
{
  return (_Fraction == .0) || _Transformation->IsIdentity();
}

// =============================================================================
// Parameters (non-DoFs)
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkPartialMultiLevelStationaryVelocityTransformation::Set(const char *param, const char *value)
{
  return _Transformation->Set(param, value);
}

// -----------------------------------------------------------------------------
irtkParameterList irtkPartialMultiLevelStationaryVelocityTransformation::Parameter() const
{
  return _Transformation->Parameter();
}

// =============================================================================
// Point transformation
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkPartialMultiLevelStationaryVelocityTransformation::RequiresCachingOfDisplacements() const
{
  return _Transformation->RequiresCachingOfDisplacements();
}

// -----------------------------------------------------------------------------
void irtkPartialMultiLevelStationaryVelocityTransformation::GlobalTransform(double &x, double &y, double &z, double t, double t0) const
{
  _Transformation->GlobalTransform(x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
void irtkPartialMultiLevelStationaryVelocityTransformation::LocalTransform(int m, int n, double &x, double &y, double &z, double t, double t0) const
{
  _Transformation->RKE1((m < 0 ? 0 : m), n, x, y, z, + _Fraction * _Transformation->UpperIntegrationLimit(t, t0));
}

// -----------------------------------------------------------------------------
void irtkPartialMultiLevelStationaryVelocityTransformation::Transform(int m, int n, double &x, double &y, double &z, double t, double t0) const
{
  _Transformation->RKE1(m, n, x, y, z, + _Fraction * _Transformation->UpperIntegrationLimit(t, t0));
}

// -----------------------------------------------------------------------------
void irtkPartialMultiLevelStationaryVelocityTransformation::GlobalInverse(double &x, double &y, double &z, double t, double t0) const
{
  _Transformation->GlobalInverse(x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
bool irtkPartialMultiLevelStationaryVelocityTransformation::LocalInverse(int m, int n, double &x, double &y, double &z, double t, double t0) const
{
  _Transformation->RKE1((m < 0 ? 0 : m), n, x, y, z, - _Fraction * _Transformation->UpperIntegrationLimit(t, t0));
  return true;
}

// -----------------------------------------------------------------------------
bool irtkPartialMultiLevelStationaryVelocityTransformation::Inverse(int m, int n, double &x, double &y, double &z, double t, double t0) const
{
  _Transformation->RKE1(m, n, x, y, z, - _Fraction * _Transformation->UpperIntegrationLimit(t, t0));
  return true;
}

// ---------------------------------------------------------------------------
void irtkPartialMultiLevelStationaryVelocityTransformation::Displacement(int m, int n, irtkGenericImage<double> &disp, double t, double t0, const irtkWorldCoordsImage *wc) const
{
  double T = _Fraction * _Transformation->UpperIntegrationLimit(t, t0);
  if (T) _Transformation->Displacement(m, n, disp, T, .0, wc);
}

// ---------------------------------------------------------------------------
void irtkPartialMultiLevelStationaryVelocityTransformation::Displacement(int m, int n, irtkGenericImage<float> &disp, double t, double t0, const irtkWorldCoordsImage *wc) const
{
  double T = _Fraction * _Transformation->UpperIntegrationLimit(t, t0);
  if (T) _Transformation->Displacement(m, n, disp, T, .0, wc);
}

// ---------------------------------------------------------------------------
int irtkPartialMultiLevelStationaryVelocityTransformation::InverseDisplacement(int m, int n, irtkGenericImage<double> &disp, double t, double t0, const irtkWorldCoordsImage *wc) const
{
  double T = _Fraction * _Transformation->UpperIntegrationLimit(t, t0);
  if (T) _Transformation->InverseDisplacement(m, n, disp, T, .0, wc);
  return 0;
}

// ---------------------------------------------------------------------------
int irtkPartialMultiLevelStationaryVelocityTransformation::InverseDisplacement(int m, int n, irtkGenericImage<float> &disp, double t, double t0, const irtkWorldCoordsImage *wc) const
{
  double T = _Fraction * _Transformation->UpperIntegrationLimit(t, t0);
  if (T) _Transformation->InverseDisplacement(m, n, disp, T, .0, wc);
  return 0;
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
void irtkPartialMultiLevelStationaryVelocityTransformation
::ParametricGradient(const irtkGenericImage<double> *in, double *out,
                     const irtkWorldCoordsImage *i2w, const irtkWorldCoordsImage *wc,
                     double t0, double w) const
{
  double t = in->GetTOrigin();
  double T = _Fraction * _Transformation->UpperIntegrationLimit(t, t0);
  if (T) _Transformation->ParametricGradient(in, out, i2w, wc, T, .0, w);
}

// -----------------------------------------------------------------------------
void irtkPartialMultiLevelStationaryVelocityTransformation
::ParametricGradient(const irtkGenericImage<double> **in, int n, double *out,
                     const irtkWorldCoordsImage *i2w, const irtkWorldCoordsImage *wc,
                     const double *t0, double w) const
{
  double t, T;
  for (int i = 0; i < n; ++i) {
    t = in[i]->GetTOrigin();
    T = _Fraction * _Transformation->UpperIntegrationLimit(t, t0 ? t0[i] : t);
    if (T) _Transformation->ParametricGradient(in[i], out, i2w, wc, T, .0, w);
  }
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void irtkPartialMultiLevelStationaryVelocityTransformation::Print(irtkIndent indent) const
{
  cout << indent << "Partial Multi-level B-spline SV FFD" << endl;
  _Transformation->Print(indent+1);
  cout << indent << "Fraction of transformation applied: " << _Fraction << endl;
}

// -----------------------------------------------------------------------------
irtkCifstream &irtkPartialMultiLevelStationaryVelocityTransformation::Read(irtkCifstream &from)
{
  return _Transformation->Read(from);
}

// -----------------------------------------------------------------------------
irtkCofstream &irtkPartialMultiLevelStationaryVelocityTransformation::Write(irtkCofstream &to) const
{
  return _Transformation->Write(to);
}
