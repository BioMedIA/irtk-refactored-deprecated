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
// Declaration of inverse transformation (cf. irtkTransformationInverse.cc)
// =============================================================================

// -----------------------------------------------------------------------------
bool EvaluateInverse(const irtkMultiLevelTransformation *, int, int,
                     double &, double &, double &, double, double);

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkMultiLevelTransformation::irtkMultiLevelTransformation()
:
  _NumberOfLevels(0)
{
  for (int l = 0; l < MAX_TRANS; ++l) {
    _LocalTransformation      [l] = NULL;
    _LocalTransformationStatus[l] = Passive;
  }
}

// -----------------------------------------------------------------------------
irtkMultiLevelTransformation::irtkMultiLevelTransformation(const irtkRigidTransformation &t)
:
  _GlobalTransformation(t),
  _NumberOfLevels(0)
{
  for (int l = 0; l < MAX_TRANS; ++l) {
    _LocalTransformation      [l] = NULL;
    _LocalTransformationStatus[l] = Passive;
  }
}

// -----------------------------------------------------------------------------
irtkMultiLevelTransformation::irtkMultiLevelTransformation(const irtkAffineTransformation &t)
:
  _GlobalTransformation(t),
  _NumberOfLevels(0)
{
  for (int l = 0; l < MAX_TRANS; ++l) {
    _LocalTransformation      [l] = NULL;
    _LocalTransformationStatus[l] = Passive;
  }
}

// -----------------------------------------------------------------------------
irtkMultiLevelTransformation::irtkMultiLevelTransformation(const irtkMultiLevelTransformation &t)
:
  irtkTransformation(t),
  _NumberOfLevels(t._NumberOfLevels)
{
  for (int l = 0; l < _NumberOfLevels; ++l) {
    _LocalTransformation[l] = dynamic_cast<irtkFreeFormTransformation *>(irtkTransformation::New(t._LocalTransformation[l]));
    if (_LocalTransformation[l] == NULL) {
      cerr << "irtkMultiLevelTransformation::irtkMultiLevelTransformation: Failed to copy local transformation at level ";
      cerr << l << " and of type " << t._LocalTransformation[l]->NameOfClass() << endl;
      exit(1);
    }
    _LocalTransformationStatus[l] = t._LocalTransformationStatus[l];
  }
}

// -----------------------------------------------------------------------------
irtkMultiLevelTransformation::~irtkMultiLevelTransformation()
{
  for (int l = 0; l < MAX_TRANS; ++l) Delete(_LocalTransformation[l]);
}

// =============================================================================
// Approximation
// =============================================================================

// -----------------------------------------------------------------------------
double irtkMultiLevelTransformation
::Approximate(const irtkImageAttributes &, double *, double *, double *,
              int, double)
{
  cerr << this->NameOfClass() << "::Approximate: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
double irtkMultiLevelTransformation
::Approximate(const double *, const double *, const double *,
              double       *, double       *, double       *, int,
              int, double)
{
  cerr << this->NameOfClass() << "::Approximate: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
double irtkMultiLevelTransformation
::Approximate(const double *, const double *, const double *, const double *,
              double       *, double       *, double       *, int,
              int, double)
{
  cerr << this->NameOfClass() << "::Approximate: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
double irtkMultiLevelTransformation
::ApproximateAsNew(const irtkImageAttributes &, double *, double *, double *,
                   int, double)
{
  cerr << this->NameOfClass() << "::ApproximateAsNew: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
double irtkMultiLevelTransformation
::ApproximateAsNew(const double *, const double *, const double *,
                   double       *, double       *, double       *, int,
                   int, double)
{
  cerr << this->NameOfClass() << "::ApproximateAsNew: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
double irtkMultiLevelTransformation
::ApproximateAsNew(const double *, const double *, const double *, const double *,
                   double       *, double       *, double       *, int,
                   int, double)
{
  cerr << this->NameOfClass() << "::ApproximateAsNew: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
void irtkMultiLevelTransformation
::ApproximateDOFs(const double *, const double *, const double *, const double *,
                  const double *, const double *, const double *, int)
{
  cerr << this->NameOfClass() << "::ApproximateDOFs: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
void irtkMultiLevelTransformation
::ApproximateDOFsGradient(const double *, const double *, const double *, const double *,
                          const double *, const double *, const double *, int,
                          double *, double) const
{
  cerr << this->NameOfClass() << "::ApproximateDOFsGradient: Not implemented" << endl;
  exit(1);
}

// =============================================================================
// Transformation parameters (DOFs)
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkMultiLevelTransformation::CopyFrom(const irtkTransformation *other)
{
  const irtkHomogeneousTransformation *lin  = NULL;
  const irtkMultiLevelTransformation  *mffd = NULL;
  const irtkFreeFormTransformation    *affd = NULL;

  (lin  = dynamic_cast<const irtkHomogeneousTransformation *>(other)) ||
  (mffd = dynamic_cast<const irtkMultiLevelTransformation  *>(other)) ||
  (affd = dynamic_cast<const irtkFreeFormTransformation    *>(other));

  if (mffd && mffd->NumberOfLevels() == 0) {
    lin  = mffd->GetGlobalTransformation();
    mffd = NULL;
  }

  if (lin) {
    this->Reset();
    _GlobalTransformation.CopyFrom(lin);
    return true;
  }
  if (affd) {
    if (_NumberOfLevels == 0) return false;
    for (int i = 0; i < _NumberOfLevels; ++i) {
      if (_LocalTransformation[i]->CopyFrom(affd)) {
        _GlobalTransformation.Reset();
        for (int j = 0; j < _NumberOfLevels; ++j) {
          if (i != j) _LocalTransformation[j]->Reset();
        }
        return true;
      }
    }
    return false;
  }
  if (mffd && mffd->NumberOfLevels() == _NumberOfLevels
           && strcmp(this->NameOfClass(), mffd->NameOfClass()) == 0) {
    if (!_GlobalTransformation.CopyFrom(mffd->GetGlobalTransformation())) {
      return false;
    }
    for (int i = 0; i < _NumberOfLevels; ++i) {
      if (!_LocalTransformation[i]->CopyFrom(mffd->GetLocalTransformation(i))) {
        this->Reset();
        return false;
      }
    }
    return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
double irtkMultiLevelTransformation::DOFGradientNorm(const double *gradient) const
{
  double norm, max = .0;
  for (int l = 0; l < this->NumberOfLevels(); ++l) {
    if (!this->LocalTransformationIsActive(l)) continue;
    const irtkFreeFormTransformation *ffd = this->GetLocalTransformation(l);
    norm = ffd->DOFGradientNorm(gradient);
    if (norm > max) max = norm;
  }
  return max;
}

// -----------------------------------------------------------------------------
int irtkMultiLevelTransformation::NumberOfDOFs() const
{
  int ndofs = 0;
  for (int l = 0; l < this->NumberOfLevels(); ++l) {
    if (!this->LocalTransformationIsActive(l)) continue;
    ndofs += this->GetLocalTransformation(l)->NumberOfDOFs();
  }
  return ndofs;
}

// -----------------------------------------------------------------------------
void irtkMultiLevelTransformation::Put(int idx, double value)
{
  irtkFreeFormTransformation *ffd;
  int                         dof;
  DOFIndexToLocalTransformation(this, idx, ffd, dof);
  ffd->Put(dof, value);
  this->Changed(true);
}

// -----------------------------------------------------------------------------
double irtkMultiLevelTransformation::Get(int idx) const
{
  const irtkFreeFormTransformation *ffd;
  int                               dof;
  DOFIndexToLocalTransformation(this, idx, ffd, dof);
  return ffd->Get(dof);
}

// -----------------------------------------------------------------------------
void irtkMultiLevelTransformation::Put(const DOFValue *x)
{
  irtkFreeFormTransformation *ffd;
  for (int l = 0; l < this->NumberOfLevels(); ++l) {
    if (!this->LocalTransformationIsActive(l)) continue;
    ffd = this->GetLocalTransformation(l);
    ffd->Put(x);
    x += ffd->NumberOfDOFs();
  }
  this->Changed(true);
}

// -----------------------------------------------------------------------------
void irtkMultiLevelTransformation::Add(const DOFValue *dx)
{
  irtkFreeFormTransformation *ffd;
  for (int l = 0; l < this->NumberOfLevels(); ++l) {
    if (!this->LocalTransformationIsActive(l)) continue;
    ffd = this->GetLocalTransformation(l);
    ffd->Add(dx);
    dx += ffd->NumberOfDOFs();
  }
  this->Changed(true);
}

// -----------------------------------------------------------------------------
double irtkMultiLevelTransformation::Update(const DOFValue *dx)
{
  double max_delta = .0;
  irtkFreeFormTransformation *ffd;
  for (int l = 0; l < this->NumberOfLevels(); ++l) {
    if (!this->LocalTransformationIsActive(l)) continue;
    ffd = this->GetLocalTransformation(l);
    max_delta = max(max_delta, ffd->Update(dx));
    dx += ffd->NumberOfDOFs();
  }
  if (max_delta > .0) this->Changed(true);
  return max_delta;
}

// -----------------------------------------------------------------------------
void irtkMultiLevelTransformation::Get(DOFValue *x) const
{
  const irtkFreeFormTransformation *ffd;
  for (int l = 0; l < this->NumberOfLevels(); ++l) {
    if (!this->LocalTransformationIsActive(l)) continue;
    ffd = this->GetLocalTransformation(l);
    ffd->Get(x);
    x += ffd->NumberOfDOFs();
  }
}

// -----------------------------------------------------------------------------
void irtkMultiLevelTransformation::PutStatus(int idx, DOFStatus status)
{
  irtkFreeFormTransformation *ffd;
  int                         dof;
  DOFIndexToLocalTransformation(this, idx, ffd, dof);
  ffd->PutStatus(dof, status);
  this->Changed(true);
}

// -----------------------------------------------------------------------------
DOFStatus irtkMultiLevelTransformation::GetStatus(int idx) const
{
  const irtkFreeFormTransformation *ffd;
  int                               dof;
  DOFIndexToLocalTransformation(this, idx, ffd, dof);
  return ffd->GetStatus(dof);
}

// -----------------------------------------------------------------------------
bool irtkMultiLevelTransformation::HasSameDOFsAs(const irtkTransformation *t) const
{
  const irtkMultiLevelTransformation *mffd;
  mffd = dynamic_cast<const irtkMultiLevelTransformation *>(t);

  if (mffd) {
    if (this->NumberOfLevels() != mffd->NumberOfLevels()) return false;
    for (int l = 0; l < this->NumberOfLevels(); ++l) {
      if (this->LocalTransformationStatus(l) != mffd->LocalTransformationStatus(l)) return false;
      const irtkFreeFormTransformation *ffd1 = this->GetLocalTransformation(l);
      const irtkFreeFormTransformation *ffd2 = mffd->GetLocalTransformation(l);
      if (!ffd1->HasSameDOFsAs(ffd2)) return false;
    }
    return true;
  } else {
    return (this->NumberOfLevels() == 1) &&
            this->LocalTransformationIsActive(0) &&
            this->GetLocalTransformation(0)->HasSameDOFsAs(t);
  }
}

// -----------------------------------------------------------------------------
bool irtkMultiLevelTransformation::IsIdentity() const
{
  if (!_GlobalTransformation.IsIdentity()) return false;
  for (int l = 0; l < this->NumberOfLevels(); ++l) {
    if (!this->GetLocalTransformation(l)->IsIdentity()) return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
void irtkMultiLevelTransformation::Reset()
{
  _GlobalTransformation.Reset();
  for (int l = 0; l < _NumberOfLevels; ++l) {
    _LocalTransformation[l]->Reset();
  }
}

// -----------------------------------------------------------------------------
void irtkMultiLevelTransformation::Clear()
{
  _GlobalTransformation.Reset();
  for (int l = 0; l < _NumberOfLevels; ++l) {
    Delete(_LocalTransformation[l]);
  }
  _NumberOfLevels = 0;
}

// =============================================================================
// Parameters (non-DoFs)
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkMultiLevelTransformation::Set(const char *name, const char *value)
{
  bool known = this->GetGlobalTransformation()->Set(name, value);
  for (int l = 0; l < this->NumberOfLevels(); ++l) {
    known = this->GetLocalTransformation(l)->Set(name, value) || known;
  }
  return known;
}

// -----------------------------------------------------------------------------
irtkParameterList irtkMultiLevelTransformation::Parameter() const
{
  irtkParameterList params = this->GetGlobalTransformation()->Parameter();
  for (int l = 0; l < this->NumberOfLevels(); ++l) {
    Insert(params, this->GetLocalTransformation(l)->Parameter());
  }
  return params;
}

// =============================================================================
// Levels
// =============================================================================

// -----------------------------------------------------------------------------
int irtkMultiLevelTransformation::NumberOfCPs(bool excl_passive) const
{
  int n = 0;
  if (excl_passive) {
    for (int l = 0; l < _NumberOfLevels; ++l) {
      if (_LocalTransformationStatus[l] == Active) {
        n += _LocalTransformation[l]->NumberOfCPs();
      }
    }
  } else {
    for (int l = 0; l < _NumberOfLevels; ++l) {
      n += _LocalTransformation[l]->NumberOfCPs();
    }
  }
  return n;
}

// -----------------------------------------------------------------------------
int irtkMultiLevelTransformation::NumberOfActiveCPs() const
{
  int n = 0;
  for (int l = 0; l < _NumberOfLevels; ++l) {
    if (_LocalTransformationStatus[l] == Active) {
      n += _LocalTransformation[l]->NumberOfActiveCPs();
    }
  }
  return n;
}

// -----------------------------------------------------------------------------
void irtkMultiLevelTransformation::CheckTransformation(irtkFreeFormTransformation *transformation) const
{
  if (transformation == NULL) {
    cerr << this->NameOfClass() << "::CheckTransformation: NULL pointer given" << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
irtkFreeFormTransformation *irtkMultiLevelTransformation::PutLocalTransformation(irtkFreeFormTransformation *transformation, int i)
{
  this->CheckTransformation(transformation);
  if (i < 0 || i >= _NumberOfLevels) {
    cerr << this->NameOfClass() << "::PutLocalTransformation: No such transformation: " << i << endl;
    exit(1);
  }
  irtkFreeFormTransformation * const old = _LocalTransformation[i];
  // Copy status of transformation if it is part of stack at different
  // location. This is the case when this method is used to rearrange the
  // local transformations within the stack.
  for (int l = 0; l < _NumberOfLevels; ++l) {
    if (transformation == _LocalTransformation[l]) {
      _LocalTransformationStatus[i] = _LocalTransformationStatus[l];
    }
  }
  // Set new transformation
  _LocalTransformation[i] = transformation;
  return old;
}

// -----------------------------------------------------------------------------
void irtkMultiLevelTransformation::PushLocalTransformation(irtkFreeFormTransformation *transformation)
{
  this->CheckTransformation(transformation);
  if (_NumberOfLevels == MAX_TRANS) {
    cerr << "irtkMultiLevelTransformation::PushLocalTransformation: Stack overflow" << endl;
    exit(1);
  }
  // Change status of transformation that is currently on top of stack
  // to passive if no more than one transformation was active before
  int nactive = 0;
  for (int l = 0; l < _NumberOfLevels; ++l) {
    if (_LocalTransformationStatus[l] == Active) ++nactive;
  }
  if (nactive == 1) {
    _LocalTransformationStatus[_NumberOfLevels-1] = Passive;
  }
  // Add transformation at top of stack
  _LocalTransformation      [_NumberOfLevels] = transformation;
  _LocalTransformationStatus[_NumberOfLevels] = Active;
  _NumberOfLevels++;
}

// -----------------------------------------------------------------------------
void irtkMultiLevelTransformation::InsertLocalTransformation(irtkFreeFormTransformation *transformation, int pos)
{
  this->CheckTransformation(transformation);
  if (_NumberOfLevels == MAX_TRANS) {
    cerr << "irtkMultiLevelTransformation::InsertLocalTransformation: Stack overflow" << endl;
    exit(1);
  }
  _NumberOfLevels++;
  for (int l = pos; l < _NumberOfLevels; ++l) {
    _LocalTransformation      [l+1] = _LocalTransformation      [l];
    _LocalTransformationStatus[l+1] = _LocalTransformationStatus[l];
  }
  _LocalTransformation      [pos] = transformation;
  _LocalTransformationStatus[pos] = Passive;
}

// -----------------------------------------------------------------------------
irtkFreeFormTransformation *irtkMultiLevelTransformation::PopLocalTransformation()
{
  irtkFreeFormTransformation *localTransformation = NULL;
  if (_NumberOfLevels > 0) {
    // Change status of transformation that is now on top of stack
    // to active if no more than one transformation was active before
    int nactive = 0;
    for (int l = 0; l < _NumberOfLevels; ++l) {
      if (_LocalTransformationStatus[l] == Active) ++nactive;
    }
    if (_NumberOfLevels > 1 && nactive == 1) {
      _LocalTransformationStatus[_NumberOfLevels-2] = Active;
    }
    // Remove transformation at top of stack
    _NumberOfLevels--;
    localTransformation = _LocalTransformation[_NumberOfLevels];
    _LocalTransformation      [_NumberOfLevels] = NULL;
    _LocalTransformationStatus[_NumberOfLevels] = Passive;
  }
  return localTransformation;
}

// -----------------------------------------------------------------------------
irtkFreeFormTransformation *irtkMultiLevelTransformation::RemoveLocalTransformation(int pos)
{
  irtkFreeFormTransformation *localTransformation;
  if (0 <= pos && pos < _NumberOfLevels) {
    localTransformation = _LocalTransformation[pos];
    for (int l = pos; l < _NumberOfLevels; ++l) {
      _LocalTransformation      [l] = _LocalTransformation      [l+1];
      _LocalTransformationStatus[l] = _LocalTransformationStatus[l+1];
    }
    _NumberOfLevels--;
    _LocalTransformation      [_NumberOfLevels] = NULL;
    _LocalTransformationStatus[_NumberOfLevels] = Passive;
  } else {
    cerr << "irtkMultiLevelTransformation::RemoveLocalTransformation: No such transformation: " << pos << endl;
    exit(1);
  }
  return localTransformation;
}

// -----------------------------------------------------------------------------
void irtkMultiLevelTransformation::CombineLocalTransformation()
{
  cerr << this->NameOfClass() << "::CombineLocalTransformation: Not implemented for this transformation" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
void irtkMultiLevelTransformation::MergeGlobalIntoLocalDisplacement()
{
  cerr << this->NameOfClass() << "::MergeGlobalIntoLocalDisplacement: Not implemented for this transformation" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
void irtkMultiLevelTransformation::InterpolateGlobalDisplacement(irtkFreeFormTransformation *f)
{
  double x, y, z;

  const int no = f->NumberOfDOFs() / 3;

  double *dx = new double[no];
  double *dy = new double[no];
  double *dz = new double[no];

  int idx = 0;
  for (int l = 0; l < f->GetT(); l++) {
    for (int k = 0; k < f->GetZ(); k++){
      for (int j = 0; j < f->GetY(); j++){
        for (int i = 0; i < f->GetX(); i++){
          x = i;
          y = j;
          z = k;
          f->LatticeToWorld(x, y, z);
          _GlobalTransformation.Displacement(x, y, z);
          dx[idx] = x;
          dy[idx] = y;
          dz[idx] = z;
          idx++;
        }
      }
    }
  }

  f->Interpolate(dx, dy, dz);

  delete[] dx;
  delete[] dy;
  delete[] dz;

  if (f->GetX() < 4 || f->GetY() < 4 || f->GetZ() < 4) {
    cerr << "irtkMultiLevelTransformation::InterpolateGlobalDisplacement: ";
    cerr << "Very small lattice for interpolation. Result likely to be inaccurate." << endl;
    return;
  }

  double totVol = f->GetX() * f->GetXSpacing() + f->GetY() * f->GetYSpacing() + f->GetZ() * f->GetZSpacing();
  double effVol = (f->GetX()-4) * f->GetXSpacing() + (f->GetY()-4) * f->GetYSpacing() + (f->GetZ()-4) * f->GetZSpacing();
/*
  Not sure if the temporal lattice should be considered here as the global transformation is only defined in 3D.
  Probably this method is best only applied to 3D MFFD's, but not 3D+t MFFD's.
  -as12321

  if (f->GetT() > 1) {
    totVol +=  f->GetT()    * f->GetTSpacing();
    effVol += (f->GetT()-4) * f->GetTSpacing();
  }
*/
  cout << "irtkMultiLevelTransformation::InterpolateGlobalDisplacement: ";
  cout << "Accurate interpolation of affine transformation over ";
  printf("% .1f %% of lattice volume\n", 100.0 * effVol / totVol);
}

// =============================================================================
// Point transformation
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkMultiLevelTransformation
::Inverse(int m, int n, double &x, double &y, double &z, double t, double t0) const
{
  return EvaluateInverse(this, m, n, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
class irtkMultiLevelTransformationToDisplacementField
{
public:

  const irtkMultiLevelTransformation *_Transformation;
  irtkGenericImage<VoxelType>        *_Output;
  double                              _TargetTime;
  double                              _SourceTime;
  int                                 _M, _N;

  // ---------------------------------------------------------------------------
  /// Evaluate transformation at displacement field voxels in 2D
  void operator() (const blocked_range2d<int>& r) const
  {
    double x,  y, z;
    double dx, dy;

    for (int j = r.rows().begin(); j != r.rows().end(); ++j)
    for (int i = r.cols().begin(); i != r.cols().end(); ++i) {
      // Transform point into world coordinates
      x = i, y = j, z = 0;
      _Output->ImageToWorld(x, y, z);
      // Apply current displacement
      dx = _Output->Get(i, j, 0, 0);
      dy = _Output->Get(i, j, 0, 1);
      x += dx, y += dy;
      // Calculate displacement
      _Transformation->Displacement(_M, _N, x, y, z, _SourceTime, _TargetTime);
      // Update displacement
      _Output->Put(i, j, 0, 0, x + dx);
      _Output->Put(i, j, 0, 1, y + dy);
    }
  }

  // ---------------------------------------------------------------------------
  /// Evaluate transformation at displacement field voxels in 3D
  void operator() (const blocked_range3d<int>& r) const
  {
    double x,  y,  z;
    double dx, dy, dz;

    for (int k = r.pages().begin(); k != r.pages().end(); ++k)
    for (int j = r.rows ().begin(); j != r.rows ().end(); ++j)
    for (int i = r.cols ().begin(); i != r.cols ().end(); ++i) {
      // Transform point into world coordinates
      x = i, y = j, z = k;
      _Output->ImageToWorld(x, y, z);
      // Apply current displacement
      dx = _Output->Get(i, j, k, 0);
      dy = _Output->Get(i, j, k, 1);
      dz = _Output->Get(i, j, k, 2);
      x += dx, y += dy, z += dz;
      // Calculate displacement
      _Transformation->Displacement(_M, _N, x, y, z, _SourceTime, _TargetTime);
      // Update displacement
      _Output->Put(i, j, k, 0, x + dx);
      _Output->Put(i, j, k, 1, y + dy);
      _Output->Put(i, j, k, 2, z + dz);
    }
  }

};

// -----------------------------------------------------------------------------
template <class VoxelType>
class irtkMultiLevelTransformationToInverseDisplacementField
{
public:

  const irtkMultiLevelTransformation *_Transformation;
  irtkGenericImage<VoxelType>        *_Output;
  double                              _TargetTime;
  double                              _SourceTime;
  int                                 _M, _N;
  int                                 _NumberOfSingularPoints;

  // ---------------------------------------------------------------------------
  /// Default constructor
  irtkMultiLevelTransformationToInverseDisplacementField()
  :
    _Transformation(NULL),
    _Output(NULL),
    _TargetTime(.0),
    _SourceTime(.0),
    _M(0), _N(0),
    _NumberOfSingularPoints(0)
  {}

  // ---------------------------------------------------------------------------
  /// Split constructor
  irtkMultiLevelTransformationToInverseDisplacementField(
    const irtkMultiLevelTransformationToInverseDisplacementField &other, split
  ) :
    _Transformation(other._Transformation),
    _Output(other._Output),
    _TargetTime(other._TargetTime),
    _SourceTime(other._SourceTime),
    _M(other._M), _N(other._N),
    _NumberOfSingularPoints(0)
  {}

  // ---------------------------------------------------------------------------
  /// Join results
  void join(const irtkMultiLevelTransformationToInverseDisplacementField &other)
  {
    _NumberOfSingularPoints += other._NumberOfSingularPoints;
  }

  // ---------------------------------------------------------------------------
  /// Evaluate transformation at displacement field voxels in 2D
  void operator() (const blocked_range2d<int>& r)
  {
    double x,  y, z;
    double dx, dy;

    for (int j = r.rows().begin(); j != r.rows().end(); ++j)
    for (int i = r.cols().begin(); i != r.cols().end(); ++i) {
      // Transform point into world coordinates
      x = i, y = j, z = 0;
      _Output->ImageToWorld(x, y, z);
      // Apply current displacement
      dx = _Output->Get(i, j, 0, 0);
      dy = _Output->Get(i, j, 0, 1);
      x += dx, y += dy;
      // Calculate inverse displacement
      if (!_Transformation->InverseDisplacement(_M, _N, x, y, z, _SourceTime, _TargetTime)) {
        ++_NumberOfSingularPoints;
      }
      // Update displacement
      _Output->Put(i, j, 0, 0, x + dx);
      _Output->Put(i, j, 0, 1, y + dy);
    }
  }

  // ---------------------------------------------------------------------------
  /// Evaluate transformation at displacement field voxels in 3D
  void operator() (const blocked_range3d<int>& r)
  {
    double x,  y,  z;
    double dx, dy, dz;

    for (int k = r.pages().begin(); k != r.pages().end(); ++k)
    for (int j = r.rows ().begin(); j != r.rows ().end(); ++j)
    for (int i = r.cols ().begin(); i != r.cols ().end(); ++i) {
      // Transform point into world coordinates
      x = i, y = j, z = k;
      _Output->ImageToWorld(x, y, z);
      // Apply current displacement
      dx = _Output->Get(i, j, k, 0);
      dy = _Output->Get(i, j, k, 1);
      dz = _Output->Get(i, j, k, 2);
      x += dx, y += dy, z += dz;
      // Calculate inverse displacement
      if (!_Transformation->InverseDisplacement(_M, _N, x, y, z, _SourceTime, _TargetTime)) {
        ++_NumberOfSingularPoints;
      }
      // Update displacement
      _Output->Put(i, j, k, 0, x + dx);
      _Output->Put(i, j, k, 1, y + dy);
      _Output->Put(i, j, k, 2, z + dz);
    }
  }

};

// -----------------------------------------------------------------------------
void irtkMultiLevelTransformation::Displacement(int m, int n, irtkGenericImage<double> &disp, double t, double t0, const irtkWorldCoordsImage *) const
{
  if (disp.GetT() < 2 || disp.GetT() > 3) {
    cerr << "irtkMultiLevelTransformation::Displacement: Input/output image must have either 2 or 3 vector components (_t)" << endl;
    exit(1);
  }

  irtkMultiLevelTransformationToDisplacementField<double> eval;
  eval._Transformation = this;
  eval._Output         = &disp;
  eval._TargetTime     = t0;
  eval._SourceTime     = t;
  eval._M              = m;
  eval._N              = n;

  if (disp.GetZ() <= 1) {
    blocked_range2d<int> voxels(0, disp.GetY(), 0, disp.GetX());
    parallel_for(voxels, eval);
  } else {
    blocked_range3d<int> voxels(0, disp.GetZ(), 0, disp.GetY(), 0, disp.GetX());
    parallel_for(voxels, eval);
  }
}

// -----------------------------------------------------------------------------
void irtkMultiLevelTransformation::Displacement(int m, int n, irtkGenericImage<float> &disp, double t, double t0, const irtkWorldCoordsImage *) const
{
  if (disp.GetT() < 2 || disp.GetT() > 3) {
    cerr << "irtkMultiLevelTransformation::Displacement: Input/output image must have either 2 or 3 vector components (_t)" << endl;
    exit(1);
  }

  irtkMultiLevelTransformationToDisplacementField<float> eval;
  eval._Transformation = this;
  eval._Output         = &disp;
  eval._TargetTime     = t0;
  eval._SourceTime     = t;
  eval._M              = m;
  eval._N              = n;

  if (disp.GetZ() <= 1) {
    blocked_range2d<int> voxels(0, disp.GetY(), 0, disp.GetX());
    parallel_for(voxels, eval);
  } else {
    blocked_range3d<int> voxels(0, disp.GetZ(), 0, disp.GetY(), 0, disp.GetX());
    parallel_for(voxels, eval);
  }
}

// -----------------------------------------------------------------------------
int irtkMultiLevelTransformation::InverseDisplacement(int m, int n, irtkGenericImage<double> &disp, double t, double t0, const irtkWorldCoordsImage *) const
{
  if (disp.GetT() < 2 || disp.GetT() > 3) {
    cerr << "irtkMultiLevelTransformation::InverseDisplacement: Input/output image must have either 2 or 3 vector components (_t)" << endl;
    exit(1);
  }

  irtkMultiLevelTransformationToInverseDisplacementField<double> eval;
  eval._Transformation = this;
  eval._Output         = &disp;
  eval._TargetTime     = t0;
  eval._SourceTime     = t;
  eval._M              = m;
  eval._N              = n;

  if (disp.GetZ() <= 1) {
    blocked_range2d<int> voxels(0, disp.GetY(), 0, disp.GetX());
    parallel_reduce(voxels, eval);
  } else {
    blocked_range3d<int> voxels(0, disp.GetZ(), 0, disp.GetY(), 0, disp.GetX());
    parallel_reduce(voxels, eval);
  }

  return eval._NumberOfSingularPoints;
}

// -----------------------------------------------------------------------------
int irtkMultiLevelTransformation::InverseDisplacement(int m, int n, irtkGenericImage<float> &disp, double t, double t0, const irtkWorldCoordsImage *) const
{
  if (disp.GetT() < 2 || disp.GetT() > 3) {
    cerr << "irtkMultiLevelTransformation::InverseDisplacement: Input/output image must have either 2 or 3 vector components (_t)" << endl;
    exit(1);
  }

  irtkMultiLevelTransformationToInverseDisplacementField<float> eval;
  eval._Transformation = this;
  eval._Output         = &disp;
  eval._TargetTime     = t0;
  eval._SourceTime     = t;
  eval._M              = m;
  eval._N              = n;

  if (disp.GetZ() <= 1) {
    blocked_range2d<int> voxels(0, disp.GetY(), 0, disp.GetX());
    parallel_reduce(voxels, eval);
  } else {
    blocked_range3d<int> voxels(0, disp.GetZ(), 0, disp.GetY(), 0, disp.GetX());
    parallel_reduce(voxels, eval);
  }

  return eval._NumberOfSingularPoints;
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
void irtkMultiLevelTransformation
::JacobianDOFs(double jac[3], int idx, double x, double y, double z, double t, double t0) const
{
  const irtkFreeFormTransformation *ffd;
  int                               dof;
  DOFIndexToLocalTransformation(this, idx, ffd, dof);
  ffd->JacobianDOFs(jac, dof, x, y, z, t, t0);
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void irtkMultiLevelTransformation::Print(irtkIndent indent) const
{
  irtkIndent subindent(indent + 1);

  // Print global transformation
  cout << indent << "Global transformation:" << endl;
  this->GetGlobalTransformation()->Print(subindent);

  // Print local transformations
  for (int l = 0; l < _NumberOfLevels; ++l) {
    cout << indent << "Local transformation";
    cout << (this->LocalTransformationIsActive(l) ? " (active)" : " (passive)");
    cout << ":" << endl;
    this->GetLocalTransformation(l)->Print(subindent);
  }
}

// -----------------------------------------------------------------------------
irtkCifstream &irtkMultiLevelTransformation::ReadDOFs(irtkCifstream &from, irtkTransformationType)
{
  unsigned int magic_no, trans_type;

  // Delete old local transformations
  for (int l = 0; l < _NumberOfLevels; ++l) {
    Delete(_LocalTransformation[l]);
  }

  // Read number of local transformations
  from.ReadAsInt(&_NumberOfLevels, 1);

  // Read global transformation
  this->GetGlobalTransformation()->Read(from);

  // Read local transformations
  for (int l = 0; l < _NumberOfLevels; ++l) {

    // Remember current file location
    int offset = from.Tell();

    // Read magic no. for transformations
    from.ReadAsUInt(&magic_no, 1);
    if (magic_no != IRTKTRANSFORMATION_MAGIC) {
      cerr << this->NameOfClass() << "::Read: Not a valid transformation found at file offset " << offset << endl;
      exit(1);
    }

    // Read transformation type
    from.ReadAsUInt(&trans_type, 1);

    // Instantiate new local transformation
    irtkTransformation *ffd = irtkTransformation::New(static_cast<irtkTransformationType>(trans_type));
    _LocalTransformation[l] = dynamic_cast<irtkFreeFormTransformation *>(ffd);
    if (_LocalTransformation[l] == NULL) {
      delete ffd;
      cerr << this->NameOfClass() << "::Read: Not a valid FFD (ID " << trans_type << ") found at file offset " << offset << endl;
      exit(1);
    }

    // Read local transformation
    from.Seek(offset);
    _LocalTransformation[l]->Read(from);
  }

  return from;
}

// -----------------------------------------------------------------------------
irtkCofstream &irtkMultiLevelTransformation::WriteDOFs(irtkCofstream &to) const
{
  // Write number of levels
  to.WriteAsInt(&_NumberOfLevels, 1);

  // Write global transformation
  _GlobalTransformation.Write(to);

  // Write local transformations
  for (int l = 0; l < _NumberOfLevels; ++l) _LocalTransformation[l]->Write(to);

  return to;
}
