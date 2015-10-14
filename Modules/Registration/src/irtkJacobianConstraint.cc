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

#include <irtkJacobianConstraint.h>


// =============================================================================
// Auxiliaries
// =============================================================================

namespace irtkJacobianConstraintUtils {


// -----------------------------------------------------------------------------
/// Body of Update function for parallel execution
struct irtkJacobianConstraintUpdate
{
private:

  const irtkJacobianConstraint     *_This;
  const irtkFreeFormTransformation *_FFD;
  irtkMatrix                       *_AdjJacobian;
  double                           *_DetJacobian;

public:

  void operator()(const blocked_range<int> &re) const
  {
    int    i, j, k, l;
    double x, y, z, t;
    for (int cp = re.begin(); cp != re.end(); ++cp) {
      _FFD->IndexToLattice(cp, i, j, k, l);
      x = i, y = j, z = k, t = _FFD->LatticeToTime(l);
      _FFD->LatticeToWorld(x, y, z);
      _DetJacobian[cp] = _This->Jacobian(_FFD, x, y, z, t, _AdjJacobian[cp]);
    }
  }

  static void Run(const irtkJacobianConstraint     *obj,
                  const irtkFreeFormTransformation *ffd,
                  double                           *det,
                  irtkMatrix                       *adj)
  {
    irtkJacobianConstraintUpdate body;
    body._This        = obj;
    body._FFD         = ffd;
    body._AdjJacobian = adj;
    body._DetJacobian = det;
    blocked_range<int> cps(0, ffd->NumberOfCPs());
    parallel_for(cps, body);
  }
};


} // namespace irtkJacobianConstraintUtils
using namespace irtkJacobianConstraintUtils;

// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkJacobianConstraint::irtkJacobianConstraint(const char *name)
:
  irtkTransformationConstraint(name),
  _DetJacobian(NULL),
  _AdjJacobian(NULL),
  _NumberOfCPs(0)
{
}

// -----------------------------------------------------------------------------
irtkJacobianConstraint::~irtkJacobianConstraint()
{
  Deallocate(_DetJacobian);
  Deallocate(_AdjJacobian);
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void irtkJacobianConstraint::Update(bool)
{
  typedef irtkJacobianConstraintUpdate Body;

  const irtkMultiLevelTransformation *mffd = NULL;
  const irtkFreeFormTransformation   *ffd  = NULL;

  (mffd = MFFD()) || (ffd = FFD());

  int ncps = 0;
  if      (ffd)  ncps =  ffd->NumberOfCPs();
  else if (mffd) ncps = mffd->NumberOfCPs(true);
  if (ncps == 0) return;

  if (ncps != _NumberOfCPs) {
    Deallocate(_DetJacobian);
    Deallocate(_AdjJacobian);
    _NumberOfCPs = ncps;
    _DetJacobian = Allocate<double>    (ncps);
    _AdjJacobian = Allocate<irtkMatrix>(ncps);
  }

  if (mffd) {
    double     *det = _DetJacobian;
    irtkMatrix *adj = _AdjJacobian;
    for (int n = 0; n < mffd->NumberOfLevels(); ++n) {
      if (mffd->LocalTransformationIsActive(n)) {
        ffd = mffd->GetLocalTransformation(n);
        Body::Run(this, ffd, det, adj);
        det += ffd->NumberOfCPs();
        adj += ffd->NumberOfCPs();
      }
    }
  } else if (ffd) {
    Body::Run(this, ffd, _DetJacobian, _AdjJacobian);
  }
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
void irtkJacobianConstraint
::WriteJacobian(const char *fname, const irtkFreeFormTransformation *ffd, const double *det) const
{
  irtkGenericImage<double> jacobian(ffd->Attributes(), const_cast<double *>(det));
  jacobian.Write(fname);
}

// -----------------------------------------------------------------------------
void irtkJacobianConstraint::WriteDataSets(const char *p, const char *suffix, bool) const
{
  if (debug < 3) return;

  const int   sz = 1024;
  char        fname[sz];
  string      _prefix = Prefix(p);
  const char *prefix  = _prefix.c_str();

  const irtkMultiLevelTransformation *mffd = NULL;
  const irtkFreeFormTransformation   *ffd  = NULL;

  (mffd = MFFD()) || (ffd = FFD());

  if (mffd) {
    double     *det = _DetJacobian;
    for (int n = 0; n < mffd->NumberOfLevels(); ++n) {
      if (mffd->LocalTransformationIsActive(n)) {
        ffd = mffd->GetLocalTransformation(n);
        if (mffd->NumberOfActiveLevels() == 1) {
          snprintf(fname, sz, "%sjacobian_determinant%s.nii.gz", prefix, suffix);
        } else {
          snprintf(fname, sz, "%sjacobian_determinant_of_ffd_at_level_%d%s.nii.gz", prefix, n+1, suffix);
        }
        WriteJacobian(fname, ffd, det);
        det += ffd->NumberOfCPs();
      }
    }
  } else if (ffd) {
    snprintf(fname, sz, "%sjacobian_determinant%s.nii.gz", prefix, suffix);
    WriteJacobian(fname, ffd, _DetJacobian);
  }
}
