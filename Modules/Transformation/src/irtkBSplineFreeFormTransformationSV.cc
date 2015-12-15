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

#include <irtkScalingAndSquaring.h>
#include <irtkDisplacementToVelocityField.h>
#include <irtkLieBracketImageFilter.h>
#include <irtkVoxelFunction.h>
#include <irtkImageToInterpolationCoefficients.h>

#include "irtkFreeFormTransformationIntegration.h"
#include <memory>

// =============================================================================
// Integration methods
// =============================================================================

// Tolerance of embedded Runge-Kutta methods with automatic step length control
static const double SVFFD_RKTOL = 1.e-3;

IRTK_FFDIM2(RKE1,   irtkBSplineFreeFormTransformationSV);
IRTK_FFDIM2(RKEH12, irtkBSplineFreeFormTransformationSV);
IRTK_FFDIM2(RKE2,   irtkBSplineFreeFormTransformationSV);
IRTK_FFDIM2(RKH2,   irtkBSplineFreeFormTransformationSV);
IRTK_FFDIM2(RKBS23, irtkBSplineFreeFormTransformationSV);
IRTK_FFDIM2(RK4,    irtkBSplineFreeFormTransformationSV);
IRTK_FFDIM2(RKF45,  irtkBSplineFreeFormTransformationSV);
IRTK_FFDIM2(RKCK45, irtkBSplineFreeFormTransformationSV);
IRTK_FFDIM2(RKDP45, irtkBSplineFreeFormTransformationSV);

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
inline double DefaultMaximumScaledVelocity(double dx, double dy, double dz)
{
  double ds = .0;
  if (dx > .0 && (ds == .0 || dx < ds)) ds = dx;
  if (dy > .0 && (ds == .0 || dy < ds)) ds = dy;
  if (dz > .0 && (ds == .0 || dz < ds)) ds = dz;
  return 0.4 * ds;
}

// -----------------------------------------------------------------------------
irtkBSplineFreeFormTransformationSV
::irtkBSplineFreeFormTransformationSV()
:
  _T                (1.0),
  _TimeUnit         (1.0),
  _NumberOfSteps    (64),
  _MaxScaledVelocity(-1.0),
  _IntegrationMethod(FFDIM_FastSS),
  _LieDerivative    (false),
  _NumberOfBCHTerms (4),
  _JacobianDOFs     (NULL),
  _JacobianDOFsIntervalLength(.0)
{
  _ExtrapolationMode = Extrapolation_NN;
}

// -----------------------------------------------------------------------------
irtkBSplineFreeFormTransformationSV
::irtkBSplineFreeFormTransformationSV(const irtkImageAttributes &attr,
                                      double dx, double dy, double dz)
:
  _T                (1.0),
  _TimeUnit         (1.0),
  _NumberOfSteps    (64),
  _MaxScaledVelocity(-1.0),
  _IntegrationMethod(FFDIM_FastSS),
  _LieDerivative    (false),
  _NumberOfBCHTerms (4),
  _JacobianDOFs     (NULL),
  _JacobianDOFsIntervalLength(.0)
{
  _ExtrapolationMode = Extrapolation_NN;
  Initialize(attr, dx, dy, dz);
}

// -----------------------------------------------------------------------------
irtkBSplineFreeFormTransformationSV
::irtkBSplineFreeFormTransformationSV(const irtkBaseImage &target,
                                      double dx, double dy, double dz)
:
  _T                (1.0),
  _TimeUnit         (1.0),
  _NumberOfSteps    (64),
  _MaxScaledVelocity(-1.0),
  _IntegrationMethod(FFDIM_FastSS),
  _LieDerivative    (false),
  _NumberOfBCHTerms (4),
  _JacobianDOFs     (NULL),
  _JacobianDOFsIntervalLength(.0)
{
  _ExtrapolationMode = Extrapolation_NN;
  Initialize(target.Attributes(), dx, dy, dz);
}

// -----------------------------------------------------------------------------
irtkBSplineFreeFormTransformationSV
::irtkBSplineFreeFormTransformationSV(const irtkGenericImage<double> &image, bool disp)
:
  _T                (1.0),
  _TimeUnit         (1.0),
  _NumberOfSteps    (64),
  _MaxScaledVelocity(-1.0),
  _IntegrationMethod(FFDIM_FastSS),
  _LieDerivative    (false),
  _NumberOfBCHTerms (4),
  _JacobianDOFs     (NULL),
  _JacobianDOFsIntervalLength(.0)
{
  Initialize(image, disp);
}

// -----------------------------------------------------------------------------
irtkBSplineFreeFormTransformationSV
::irtkBSplineFreeFormTransformationSV(const irtkBSplineFreeFormTransformationSV &ffd)
:
  irtkBSplineFreeFormTransformation3D(ffd),
  _T                (ffd._T),
  _TimeUnit         (ffd._TimeUnit),
  _NumberOfSteps    (ffd._NumberOfSteps),
  _MaxScaledVelocity(ffd._MaxScaledVelocity),
  _IntegrationMethod(ffd._IntegrationMethod),
  _LieDerivative    (ffd._LieDerivative),
  _NumberOfBCHTerms (ffd._NumberOfBCHTerms),
  _JacobianDOFs     (NULL),
  _JacobianDOFsIntervalLength(.0)
{
}

// -----------------------------------------------------------------------------
irtkBSplineFreeFormTransformationSV
::~irtkBSplineFreeFormTransformationSV()
{
  delete _JacobianDOFs;
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationSV::Initialize(const irtkImageAttributes &attr)
{
  irtkBSplineFreeFormTransformation3D::Initialize(attr);
  if (_MaxScaledVelocity < .0) _MaxScaledVelocity = DefaultMaximumScaledVelocity(_dx, _dy, _dz);
  Delete(_JacobianDOFs);
  _JacobianDOFsIntervalLength = .0;
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationSV::Initialize(const irtkImageAttributes &attr,
                                                     double dx, double dy, double dz,
                                                     const irtkTransformation *dof)
{
  // Initialize free-form deformation (for extended image grid)
  //
  // Ensure that for all target voxels the displacement can be recovered
  // without requiring any extrapolation of the velocity field during
  // computation of the trajectory (integration, i.e., exponentiation)
  this->Initialize(ApproximationDomain(attr, dof), dx, dy, dz);

  // Approximate given transformation
  this->ApproximateAsNew(dof);
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationSV::Subdivide(bool subdivide_x, bool subdivide_y, bool subdivide_z, bool subdivide_t)
{
  irtkBSplineFreeFormTransformation3D::Subdivide(subdivide_x, subdivide_y, subdivide_z, subdivide_t);
  if (_MaxScaledVelocity > .0) _MaxScaledVelocity /= 2;
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationSV::Changed(bool changed)
{
  irtkBSplineFreeFormTransformation3D::Changed(changed);
  if (changed) _JacobianDOFsIntervalLength = .0;
}

// =============================================================================
// Auxiliary functions
// =============================================================================

// -----------------------------------------------------------------------------
/// Voxel function used to evaluate Lie bracket at each lattice point using
/// the Lie derivative definition which is based on the Jacobian matrices
/// of the vector fields. Uses nearest neighbor extrapolation of the velocity field.
class SVFFDEvaluateLieBracket : public irtkVoxelFunction
{
private:

  typedef irtkBSplineFreeFormTransformationSV::Vector Vector;
  typedef irtkBSplineFreeFormTransformationSV::Kernel Kernel;

  const irtkBSplineFreeFormTransformationSV *_FFD;

  double                          _tau; ///< Scaling of left/first vector field
  const irtkGenericImage<Vector> &_v;   ///< Left/First vector field
  const irtkGenericImage<Vector> &_w;   ///< Right/Second vector field

public:

  // ---------------------------------------------------------------------------
  /// Constructor
  SVFFDEvaluateLieBracket(const irtkBSplineFreeFormTransformationSV *ffd,
                          const irtkGenericImage<Vector>            &v,
                          const irtkGenericImage<Vector>            &w)
  :
    _FFD(ffd), _tau(1.0), _v(v), _w(w)
  {}

  // ---------------------------------------------------------------------------
  /// Constructor
  SVFFDEvaluateLieBracket(const irtkBSplineFreeFormTransformationSV *ffd,
                          double                                     tau,
                          const irtkGenericImage<Vector>            &v,
                          const irtkGenericImage<Vector>            &w)
  :
    _FFD(ffd), _tau(tau), _v(v), _w(w)
  {}

  // ---------------------------------------------------------------------------
  /// Evaluate velocity at lattice point
  void Evaluate(Vector &x, double s, const irtkGenericImage<Vector> &v, int i, int j, int k) const
  {
    double B_I, B_J, B_K;
    int    II,  JJ,  KK;
  
    x = .0;
    for (int K = k-1; K <= k+1; K++) {
      B_K = Kernel::LatticeWeights[K - (k-1)];
      if      (K <  0)        KK = 0;
      else if (K >= v.GetZ()) KK = v.GetZ()-1;
      else                    KK = K;
      for (int J = j-1; J <= j+1; J++) {
        B_J = Kernel::LatticeWeights[J - (j-1)];
        if      (J <  0)        JJ = 0;
        else if (J >= v.GetY()) JJ = v.GetY()-1;
        else                    JJ = J;
        for (int I = i-1; I <= i+1; I++) {
          B_I = Kernel::LatticeWeights[I - (i-1)];
          if      (I <  0)        II = 0;
          else if (I >= v.GetX()) II = v.GetX()-1;
          else                    II = I;
          x += B_I * B_J * B_K * s * v(II, JJ, KK);
        }
      }
    }
  }

  // ---------------------------------------------------------------------------
  /// Evaluate Jacobian of velocity field at lattice point
  void Jacobian(irtkMatrix &jac, double s, const irtkGenericImage<Vector> &v, int i, int j, int k) const
  {
    int    I, J, K, II, JJ, KK;
    double B_I, B_J, B_K, B_I_I, B_J_I, B_K_I;
    Vector dx, dy, dz;

    for (K = k-1; K <= k+1; K++) {
      B_K   = Kernel::LatticeWeights  [K - (k-1)];
      B_K_I = Kernel::LatticeWeights_I[K - (k-1)];
      if      (K <  0)        KK = 0;
      else if (K >= v.GetZ()) KK = v.GetZ()-1;
      else                    KK = K;
      for (J = j-1; J <= j+1; J++) {
        B_J   = Kernel::LatticeWeights  [J - (j-1)];
        B_J_I = Kernel::LatticeWeights_I[J - (j-1)];
        if      (J <  0)        JJ = 0;
        else if (J >= v.GetY()) JJ = v.GetY()-1;
        else                    JJ = J;
        for (I = i-1; I <= i+1; I++) {
          B_I   = Kernel::LatticeWeights  [I - (i-1)];
          B_I_I = Kernel::LatticeWeights_I[I - (i-1)];
          if      (I <  0)        II = 0;
          else if (I >= v.GetX()) II = v.GetX()-1;
          else                    II = I;
          const Vector &coeff = v(II, JJ, KK);
          dx += B_I_I * B_J   * B_K   * s * coeff;
          dy += B_I   * B_J_I * B_K   * s * coeff;
          dz += B_I   * B_J   * B_K_I * s * coeff;
        }
      }
    }

    jac.Initialize(3, 3);
    jac(0, 0) = dx._x;
    jac(0, 1) = dy._x;
    jac(0, 2) = dz._x;
    jac(1, 0) = dx._y;
    jac(1, 1) = dy._y;
    jac(1, 2) = dz._y;
    jac(2, 0) = dx._z;
    jac(2, 1) = dy._z;
    jac(2, 2) = dz._z;

    _FFD->JacobianToWorld(jac);
  }

  // ---------------------------------------------------------------------------
  /// Compute product of 3x3 matrix and 3D column vector
  static Vector MatrixProduct(const irtkMatrix &jac, const Vector &vel)
  {
    Vector v;
    v._x  = jac(0, 0) * vel._x + jac(0, 1) * vel._y + jac(0, 2) * vel._z;
    v._y  = jac(1, 0) * vel._x + jac(1, 1) * vel._y + jac(1, 2) * vel._z;
    v._z  = jac(2, 0) * vel._x + jac(2, 1) * vel._y + jac(2, 2) * vel._z;
    return v;
  }

  // ---------------------------------------------------------------------------
  /// Evaluate Lie bracket at given lattice point, u = [v, w]
  void operator ()(int i, int j, int k, int, Vector *u) const
  {
    irtkMatrix jac(3, 3);
    Vector     vel;
    // u = J_w * v
    Jacobian(jac,  1.0, _w, i, j, k);
    Evaluate(vel, _tau, _v, i, j, k);
    *u  = MatrixProduct(jac, vel);
    // u = J_w * v - J_v * w
    Jacobian(jac, _tau, _v, i, j, k);
    Evaluate(vel,  1.0, _w, i, j, k);
    *u -= MatrixProduct(jac, vel);
  }
};

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationSV
::EvaluateBCHFormula(int nterms, CPImage &u, double tau, const CPImage &v, const CPImage &w, bool minus_v) const
{
  IRTK_START_TIMING();
  irtkGenericImage<Vector> l1, l2, l3, l4;
  const irtkImageAttributes &lattice = u.Attributes();

  // Calculate required Lie brackets...
  if (_LieDerivative) {
    // ...using Lie derivative
    if (nterms >= 3) {
      // - [v, w]
      l1.Initialize(lattice, 3);
      ParallelForEachVoxel(SVFFDEvaluateLieBracket(this, tau, v, w), lattice, l1);
      ConvertToCubicBSplineCoefficients(l1);
      if (nterms >= 4) {
        // - [v, [v, w]]
        l2.Initialize(lattice, 3);
        ParallelForEachVoxel(SVFFDEvaluateLieBracket(this, tau, v, l1), lattice, l2);
        ConvertToCubicBSplineCoefficients(l2);
        if (nterms >= 5) {
          // - [[v, w], w]
          l3.Initialize(lattice, 3);
          ParallelForEachVoxel(SVFFDEvaluateLieBracket(this, l1, w), lattice, l3);
          ConvertToCubicBSplineCoefficients(l3);
          if (nterms >= 6) {
            // - [[v, [v, w]], w]
            l4.Initialize(lattice, 3);
            ParallelForEachVoxel(SVFFDEvaluateLieBracket(this, l2, w), lattice, l4);
            ConvertToCubicBSplineCoefficients(l4);
            // - [[w, [v, w]], v] == [[v, [v, w]], w]
          }
        }
      }
    }
  } else {
    // ...using composition of vector fields
    if (nterms >= 3) {
      irtkDifferenceOfCompositionLieBracketImageFilter3D<Vector> lb;
      lb.Interpolation(Interpolation_CubicBSpline);
      lb.Extrapolation(Extrapolation_NN);
      lb.ComputeInterpolationCoefficients(false);
      // - [v, w]
      lb.SetInput  (0, const_cast<irtkGenericImage<Vector> *>(&v));
      lb.SetInput  (1, const_cast<irtkGenericImage<Vector> *>(&w));
      lb.SetOutput (&l1);
      lb.SetScaling(0, tau);
      lb.Run();
      lb.SetScaling(0, 1.0);
      ConvertToCubicBSplineCoefficients(l1);
      if (nterms >= 4) {
        // - [v, [v, w]]
        lb.SetInput  (0, const_cast<irtkGenericImage<Vector> *>(&v));
        lb.SetInput  (1, &l1);
        lb.SetOutput (&l2);
        lb.SetScaling(0, tau);
        lb.Run();
        lb.SetScaling(0, 1.0);
        ConvertToCubicBSplineCoefficients(l2);
        if (nterms >= 5) {
          // - [[v, w], w]
          lb.SetInput (0, &l1);
          lb.SetInput (1, const_cast<irtkGenericImage<Vector> *>(&w));
          lb.SetOutput(&l3);
          lb.Run();
          ConvertToCubicBSplineCoefficients(l3);
          if (nterms >= 6) {
            // - [[v, [v, w]], w]
            lb.SetInput (0, &l2);
            lb.SetInput (1, const_cast<irtkGenericImage<Vector> *>(&w));
            lb.SetOutput(&l4);
            lb.Run();
            ConvertToCubicBSplineCoefficients(l4);
            // - [[w, [v, w]], v] == [[v, [v, w]], w]
          }
        }
      }
    }
  }

  // Evaluate BCH formula given all pre-computed terms and their respective weights
  const double weight1[] = {0.0, 1.0, 1.0/2.0, 1.0/12.0, 1.0/12.0, 1.0/48.0, 1.0/48.0};
  const double weight2[] = {1.0, 1.0, 1.0/2.0, 1.0/12.0, 1.0/12.0, 1.0/48.0, 1.0/48.0};

  irtkNaryVoxelFunction::VoxelWiseWeightedSum bch;
  if (minus_v) bch._Weight = weight1;
  else         bch._Weight = weight2;

  switch (nterms) {
    case 1: { if (minus_v) { u = .0; } else { u = v; }                 break; }
    case 2: { ParallelForEachScalar(v, w,                     u, bch); break; }
    case 3: { ParallelForEachScalar(v, w, l1,                 u, bch); break; }
    case 4: { ParallelForEachScalar(v, w, l1, l2,             u, bch); break; }
    case 5: { ParallelForEachScalar(v, w, l1, l2, l3,         u, bch); break; }
    case 6: { ParallelForEachScalar(v, w, l1, l2, l3, l4,     u, bch); break; }
    case 7: { ParallelForEachScalar(v, w, l1, l2, l3, l4, l4, u, bch); break; }
    default:
      cerr << "irtkBSplineFreeFormTransformationSV::EvaluateBCHFormula: Invalid number of terms " << nterms << endl;
      exit(1);
  };
  IRTK_DEBUG_TIMING(3, "evaluation of BCH formula");
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationSV
::EvaluateBCHFormula(int nterms, CPImage &u, const CPImage &v, const CPImage &w, bool minus_v) const
{
  EvaluateBCHFormula(nterms, u, 1.0, v, w, minus_v);
}

// =============================================================================
// Approximation/Interpolation
// =============================================================================

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationSV
::ApproximateDOFs(const double *x,  const double *y,  const double *z, const double *t,
                  const double *dx, const double *dy, const double *dz, int no)
{
  // FIXME The initial approximation of the displacements on the control point grid
  //       smoothes the displacement field too much and hence introduces quite some
  //       error. Use the overloaded ApproximateAsNew(disp) method when possible.
  //       -as12312

  // Approximate displacements at control points
  double *rx = Allocate<double>(no);
  double *ry = Allocate<double>(no);
  double *rz = Allocate<double>(no);

  memcpy(rx, dx, no * sizeof(double));
  memcpy(ry, dy, no * sizeof(double));
  memcpy(rz, dz, no * sizeof(double));

  irtkBSplineFreeFormTransformation3D::ApproximateDOFs(x, y, z, t, rx, ry, rz, no);

  Deallocate(rx);
  Deallocate(ry);
  Deallocate(rz);

  // Find stationary velocity field which approximates the displacements
  irtkGenericImage<double> disp(this->Attributes(), 3);
  for (int k = 0; k < _z; ++k) {
    for (int j = 0; j < _y; ++j) {
      for (int i = 0; i < _x; ++i) {
        disp(i, j, k, 0) = _CPImage(i, j, k)._x;
        disp(i, j, k, 1) = _CPImage(i, j, k)._y;
        disp(i, j, k, 2) = _CPImage(i, j, k)._z;
      }
    }
  }
  this->ApproximateAsNew(disp);
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationSV
::ApproximateDOFsGradient(const double *, const double *, const double *, const double *,
                          const double *, const double *, const double *, int,
                          double *, double) const
{
  cerr << this->NameOfClass() << "::ApproximateDOFsGradient: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
irtkImageAttributes irtkBSplineFreeFormTransformationSV
::ApproximationDomain(const irtkImageAttributes &attr, const irtkTransformation *dof)
{
  if (!dof) return attr;

  irtkImageAttributes grid(attr);

  // Ensure that for all target voxels the displacement can be recovered
  // without requiring any extrapolation of the velocity field during
  // computation of the trajectory (integration, i.e., exponentiation)
  const irtkMatrix i2w = grid.GetImageToWorldMatrix();
  const irtkMatrix w2i = grid.GetWorldToImageMatrix();

  double margin_top    = .0;
  double margin_bottom = .0;
  double margin_left   = .0;
  double margin_right  = .0;
  double margin_front  = .0;
  double margin_back   = .0;

  double x, y, z, wx, wy, wz;
  for (int k = 0; k < attr._z; ++k) {
    for (int j = 0; j < attr._y; ++j) {
      for (int i = 0; i < attr._x; ++i) {
        wx = i2w(0, 0) * i + i2w(0, 1) * j + i2w(0, 2) * k + i2w(0, 3);
        wy = i2w(1, 0) * i + i2w(1, 1) * j + i2w(1, 2) * k + i2w(1, 3);
        wz = i2w(2, 0) * i + i2w(2, 1) * j + i2w(2, 2) * k + i2w(2, 3);
        dof->Transform(wx, wy, wz);
        x = w2i(0, 0) * wx + w2i(0, 1) * wy + w2i(0, 2) * wz + w2i(0, 3);
        y = w2i(1, 0) * wx + w2i(1, 1) * wy + w2i(1, 2) * wz + w2i(1, 3);
        z = w2i(2, 0) * wx + w2i(2, 1) * wy + w2i(2, 2) * wz + w2i(2, 3);
        if (x <  0           && -x              > margin_left)   margin_left   = -x;
        if (y <  0           && -y              > margin_bottom) margin_bottom = -y;
        if (z <  0           && -z              > margin_front)  margin_front  = -z;
        if (x >= grid._x - 1 && x - grid._x - 1 > margin_right)  margin_right  = x - grid._x - 1;
        if (y >= grid._y - 1 && y - grid._y - 1 > margin_top)    margin_top    = y - grid._y - 1;
        if (z >= grid._z - 1 && z - grid._z - 1 > margin_back)   margin_back   = z - grid._z - 1;
      }
    }
  }

  // Account for inter-/extrapolation error on boundary of FFD lattice and
  // therefore make lattice a bit bigger than otherwise needed
  const double margin_safety = 1.5;
  margin_left   = ceil(margin_left   * margin_safety);
  margin_right  = ceil(margin_right  * margin_safety);
  margin_bottom = ceil(margin_bottom * margin_safety);
  margin_top    = ceil(margin_top    * margin_safety);
  margin_front  = ceil(margin_front  * margin_safety);
  margin_back   = ceil(margin_back   * margin_safety);

  // Compute offsets by which lattice origin must be moved such that
  // the lattice origin is the center of the extended lattice again
  const double ox = (margin_right - margin_left)   * grid._dx / 2.0;
  const double oy = (margin_top   - margin_bottom) * grid._dy / 2.0;
  const double oz = (margin_back  - margin_front)  * grid._dz / 2.0;

  // Initialize free-form deformation (for extended image grid)
  grid._x       += margin_left + margin_right;
  grid._y       += margin_bottom + margin_top;
  grid._z       += margin_front + margin_back;
  grid._xorigin += grid._xaxis[0] * ox + grid._yaxis[0] * oy + grid._zaxis[0] * oz;
  grid._yorigin += grid._xaxis[1] * ox + grid._yaxis[1] * oy + grid._zaxis[1] * oz;
  grid._zorigin += grid._xaxis[2] * ox + grid._yaxis[2] * oy + grid._zaxis[2] * oz;

  return grid;
}

// -----------------------------------------------------------------------------
double irtkBSplineFreeFormTransformationSV
::ApproximateAsNew(const irtkImageAttributes &domain, const irtkTransformation *dof,
                   int niter, double max_error)
{
  const irtkHomogeneousTransformation       *lin   = NULL;
  const irtkBSplineFreeFormTransformationSV *svffd = NULL;

  (lin   = dynamic_cast<const irtkHomogeneousTransformation       *>(dof)) ||
  (svffd = dynamic_cast<const irtkBSplineFreeFormTransformationSV *>(dof));

  // Approximate any other transformation using the base class implementation
  // which simply evaluates the displacement of the transformation at each
  // control point and then calls Interpolate in order to interpolate these
  // control point displacements
  if (!lin && !svffd) {
    return irtkBSplineFreeFormTransformation3D::ApproximateAsNew(domain, dof, niter, max_error);
  }

  double error = numeric_limits<double>::infinity();
  if (niter < 1) return error;

  for (int iter = 0; iter < niter && error > max_error; ++iter) {

    // Compute velocities at control points using log map of affine matrix
    if (lin) {
      EvaluateGlobalSVFFD logA(logm(lin->GetMatrix()), &_CPImage);
      ParallelForEachVoxel(_CPImage.Attributes(), &_CPImage, logA);
    // Evaluate velocities of other SV FFD at control points of this SV FFD
    } else {
      Vector *v = _CPImage.GetPointerToVoxels();
      for (int k = 0; k < _z; ++k) {
        for (int j = 0; j < _y; ++j) {
          for (int i = 0; i < _x; ++i, ++v) {
            v->_x = i, v->_y = j, v->_z = k;
            svffd->Evaluate(v->_x, v->_y, v->_z);
          }
        }
      }
    }

    // Convert velocities to B-spline coefficients
    ConvertToSplineCoefficients(3, _CPImage);

    // Evaluate approximation error
    error = EvaluateRMSError(domain, dof);
  }

  return error;
}

// -----------------------------------------------------------------------------
//double irtkBSplineFreeFormTransformationSV
//::ApproximateFFDAsNew(irtkGenericImage<double> &disp, double T, int niter)
//{
//  irtkGenericImage<double>  dv;
//  irtkWorldCoordinatesImage wc;
//  disp.ImageToWorld(wc);
//  for (int iter = 0; iter < niter; ++iter) {
//    // Interpolate velocity field
//    const irtkImageAttributes &grid = d.GetImageAttributes();
//    v->Initialize(grid, d.GetT());
//    switch (d.GetT()) {
//      case 2: { ParallelForEachVoxel(EvaluateBSplineSVFFD2D(this, v), grid, v); break; }
//      case 3: { ParallelForEachVoxel(EvaluateBSplineSVFFD3D(this, v), grid, v); break; }
//      default: {
//        cerr << "irtkBSplineFreeFormTransformationSV::ScalingAndSquaring: Vector field must have 2 or 3 components (_t)" << endl;
//        exit(1);
//      }
//    }
//    // Exponentiate velocity field
//    irtkVelocityToDisplacementFieldSS<double> exp;
//    exp.T(T);
//    exp.NumberOfSteps    (NumberOfStepsForIntervalLength(T));
//    exp.MaxScaledVelocity(_MaxScaledVelocity);
//    exp.Interpolation    (Interpolation_BSpline);
//    exp.Upsample         (false); // better, but too expensive
//    exp.ExternalCache    (cache); // avoid frequent allocation/deallocation
//    exp.SetInput         (0, v);  // velocity field to be exponentiated
//    exp.SetInput         (1, &d); // input displacement field (may be zero)
//    exp.SetOutput        (v);     // result is exp(v) o d
//    exp.Run();
//    d.CopyFrom(*v);
//    // dv = exp(-v) o disp
//    dv = disp;
//    ScalingAndSquaring(dv, T, &wc);
//    // Approximate dv by B-spline coefficients at control points
//    irtkBSplineFreeFormTransformation3D ffd(this->GetFFDAttributes());
//    ffd.Approximate(wc.GetPointerToVoxels(0, 0, 0, 0),
//                    wc.GetPointerToVoxels(0, 0, 0, 1),
//                    wc.GetPointerToVoxels(0, 0, 0, 2),
//                    dv.GetPointerToVoxels(0, 0, 0, 0),
//                    dv.GetPointerToVoxels(0, 0, 0, 1),
//                    dv.GetPointerToVoxels(0, 0, 0, 2));
//    dv.Initialize(_x, _y, _z, 1, 3);
//    for (int k = 0; k < _z; ++k) {
//      for (int j = 0; j < _y; ++j) {
//        for (int i = 0; i < _x; ++i) {
//          dv(i, j, k, 0) = ffd._data[k][j][i]._x;
//          dv(i, j, k, 1) = ffd._data[k][j][i]._x;
//          dv(i, j, k, 2) = ffd._data[k][j][i]._x;
//        }
//      }
//    }
//    // Calculate required Lie brackets
//    if (_NumberOfBCHTerms >= 3) {
//      irtkDifferenceOfCompositionLieBracketImageFilter3D<double> lb;
//      lb.Interpolation(Interpolation_CubicBSpline);
//      lb.Extrapolation(Extrapolation_NN);
//      lb.ComputeInterpolationCoefficients(false);
//      // - [v, dv]
//      lb.SetInput (0, v);
//      lb.SetInput (1, d);
//      lb.SetOutput(l1);
//      lb.Run();
//      ConvertToCubicBSplineCoefficients(*l1);
//      if (_NumberOfBCHTerms >= 4) {
//        // - [v, [v, dv]]
//        lb.SetInput (0, v);
//        lb.SetInput (1, l1);
//        lb.SetOutput(l2);
//        lb.Run();
//        ConvertToCubicBSplineCoefficients(*l2);
//        if (_NumberOfBCHTerms >= 5) {
//          // - [dv, [v, dv]]
//          lb.SetInput (0, d);
//          lb.SetInput (1, l1);
//          lb.SetOutput(l3);
//          lb.Run();
//          ConvertToCubicBSplineCoefficients(*l3);
//          if (_NumberOfBCHTerms >= 6) {
//            // - [dv, [v, [v, dv]]]
//            lb.SetInput (0, d);
//            lb.SetInput (1, l2);
//            lb.SetOutput(l4);
//            lb.Run();
//            ConvertToCubicBSplineCoefficients(*l4);
//          }
//        }
//      }
//    }
//    // Update velocity coefficients using BCH formula
//    irtkNaryVoxelFunction::EvaluateBCHFormula bch;
//    switch (_NumberOfBCHTerms) {
//      case 2: { /* nothing else to do */                              break; }
//      case 3: { ParallelForEachScalar(v, dv, l1,             v, bch); break; }
//      case 4: { ParallelForEachScalar(v, dv, l1, l2,         v, bch); break; }
//      case 5: { ParallelForEachScalar(v, dv, l1, l2, l3,     v, bch); break; }
//      case 6: { ParallelForEachScalar(v, dv, l1, l2, l3, l4, v, bch); break; }
//    };
//  }
//}

// -----------------------------------------------------------------------------
double irtkBSplineFreeFormTransformationSV
::ApproximateAsNew(irtkGenericImage<double> &disp, int niter, double max_error)
{
  return this->ApproximateAsNew(disp, false, 3, niter * 8);
}

// -----------------------------------------------------------------------------
double irtkBSplineFreeFormTransformationSV
::ApproximateAsNew(irtkGenericImage<double> &disp, bool smooth, int nterms, int niter)
{
  // TODO: Refactor/review implementation again after update of velocities
  //       from gradient is now implemented and working.

  irtkImageAttributes grid = this->Attributes();
  grid._t  = 3;
  grid._dt = disp.GetTSize(); // ignore difference in _dt

  // Sample displacement field at control points using linear interpolation
  irtkGenericImage<double> *d = NULL;

  if (disp.GetImageAttributes() == grid) {
    d = &disp;
  } else {    
    double x, y, z, vec[3] = {.0, .0, .0};

    std::unique_ptr<irtkInterpolateImageFunction> f(irtkInterpolateImageFunction::New(Interpolation_Linear, Extrapolation_NN, &disp));
    f->SetInput(&disp);
    f->Initialize();

    d = new irtkGenericImage<double>(grid);

    for (int k = 0; k < grid._z; k++) {
      for (int j = 0; j < grid._y; j++) {
        for (int i = 0; i < grid._x; i++) {
          x = i;
          y = j;
          z = k;
          d->ImageToWorld(x, y, z);
          disp.WorldToImage(x, y, z);
          f->Evaluate(vec, x, y, z);
          d->Put(i, j, k, 0, vec[0]);
          d->Put(i, j, k, 1, vec[1]);
          d->Put(i, j, k, 2, vec[2]);
        }
      }
    }
  }

  // Compute stationary velocity field at control points
  irtkGenericImage<double>                   v;
  irtkDisplacementToVelocityFieldBCH<double> dtov;

  dtov.SetInput (d);
  dtov.SetOutput(&v);

  dtov.SetT                 (UpperIntegrationLimit(0, 1));
  dtov.SetNumberOfIterations(niter);
  dtov.SetNumberOfTerms     (nterms);
  dtov.SetNumberOfSteps     (NumberOfStepsForIntervalLength(dtov.GetT()));
  dtov.SetSmoothVelocities  (smooth);

  dtov.Run();

  // Free temporary displacement field
  if (d != &disp) {
    delete d;
    d = NULL;
  }

  // Interpolate velocities by B-spline function
  irtkBSplineFreeFormTransformation3D::Interpolate(v.GetPointerToVoxels(0, 0, 0, 0),
                                                   v.GetPointerToVoxels(0, 0, 0, 1),
                                                   v.GetPointerToVoxels(0, 0, 0, 2));

  // Evaluate RMS of approximation error
  double error = .0;

  v = .0;
  this->Displacement(v, 0, 1);

  for (int k = 0; k < disp.GetZ(); k++) {
    for (int j = 0; j < disp.GetY(); j++) {
      for (int i = 0; i < disp.GetX(); i++) {
        disp(i, j, k, 0) -= v(i, j, k, 0);
        disp(i, j, k, 1) -= v(i, j, k, 1);
        disp(i, j, k, 2) -= v(i, j, k, 2);
        error += sqrt(disp(i, j ,k, 0) * disp(i, j, k, 0) +
                      disp(i, j ,k, 1) * disp(i, j, k, 1) +
                      disp(i, j ,k, 2) * disp(i, j, k, 2));
      }
    }
  }
  error /= static_cast<double>(disp.GetX() * disp.GetY() * disp.GetZ());

  return error;
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationSV
::Interpolate(const double *, const double *, const double *)
{
  cerr << this->NameOfClass() << "::Interpolate: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationSV::CombineWith(const irtkTransformation *dof)
{
  // Convert transformation into SV FFD
  const irtkBSplineFreeFormTransformationSV *svffd = NULL;
  svffd = dynamic_cast<const irtkBSplineFreeFormTransformationSV *>(dof);
  if (!svffd) {
    irtkBSplineFreeFormTransformationSV *tmp;
    tmp = new irtkBSplineFreeFormTransformationSV();
    tmp->Initialize(this->Attributes());
    tmp->ApproximateAsNew(dof);
    svffd = tmp;
  }
  // Compute coefficients of composite SV FFD using BCH formula
  EvaluateBCHFormula(4, _CPImage, _CPImage, svffd->_CPImage);
  // Delete temporary SV FFD
  if (svffd != dof) delete svffd;
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationSV::Invert()
{
  _CPImage *= -1.0;
}

// =============================================================================
// Parameters (non-DoFs)
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkBSplineFreeFormTransformationSV::Set(const char *name, const char *value)
{
  if (strcmp(name, "Cross-sectional time interval") == 0 ||
      strcmp(name, "Cross sectional time interval") == 0) {
    return FromString(value, _T);
  } else if (strcmp(name, "Time unit of integration interval") == 0) {
    return FromString(value, _TimeUnit);
  } else if (strcmp(name,    "No. of integration steps") == 0 ||
             strcmp(name, "Number of integration steps") == 0) {
    return FromString(value, _NumberOfSteps) && _NumberOfSteps > 0;
  } else if (strcmp(name,    "No. of squaring steps")    == 0 ||
             strcmp(name, "Number of squaring steps") == 0) {
    if (!FromString(value, _NumberOfSteps) || _NumberOfSteps <= 0) return false;
    _NumberOfSteps = pow(2.0, _NumberOfSteps);
    if (_IntegrationMethod != FFDIM_SS && _IntegrationMethod != FFDIM_FastSS) {
      _IntegrationMethod = FFDIM_FastSS;
    }
  } else if (strcmp(name, "Maximum scaled velocity") == 0) {
    return FromString(value, _MaxScaledVelocity);

  } else if (strcmp(name, "Use Lie derivative") == 0) {
    return FromString(value, _LieDerivative);
  } else if (strcmp(name,    "No. of BCH terms")    == 0 ||
             strcmp(name, "Number of BCH terms") == 0) {
    return FromString(value, _NumberOfBCHTerms) && _NumberOfBCHTerms <= 6;
  } else if (strcmp(name, "Integration method") == 0) {
    return FromString(value, _IntegrationMethod) && _IntegrationMethod != FFDIM_UNKNOWN;
  // deprecated parameters
  } else if (strcmp(name, "Use scaling and squaring") == 0) {
    bool useSS = false;
    if (!FromString(value, useSS)) return false;
    if (useSS) {
      if (_IntegrationMethod != FFDIM_SS && _IntegrationMethod != FFDIM_FastSS) {
        _IntegrationMethod = FFDIM_FastSS;
      }
    } else {
      if (_IntegrationMethod == FFDIM_SS || _IntegrationMethod == FFDIM_FastSS) {
        _IntegrationMethod = FFDIM_RKE1;
      }
    }
    return true;
  } else if (strcmp(name, "Fast scaling and squaring") == 0) {
    bool fastSS;
    if (!FromString(value, fastSS)) return false;
    if (_IntegrationMethod == FFDIM_SS && fastSS) _IntegrationMethod = FFDIM_FastSS;
    return true;
  }
  return irtkBSplineFreeFormTransformation3D::Set(name, value);
}

// -----------------------------------------------------------------------------
irtkParameterList irtkBSplineFreeFormTransformationSV::Parameter() const
{
  irtkParameterList params = irtkBSplineFreeFormTransformation3D::Parameter();
  Insert(params, "Integration method",                ToString(_IntegrationMethod));
  Insert(params, "Cross-sectional time interval",     ToString(_T));
  Insert(params, "Time unit of integration interval", ToString(_TimeUnit));
  Insert(params, "No. of integration steps",          ToString(_NumberOfSteps));
  Insert(params, "Maximum scaled velocity",           ToString(_MaxScaledVelocity));
  Insert(params, "Use Lie derivative",                ToString(_LieDerivative));
  Insert(params, "No. of BCH terms",                  ToString(_NumberOfBCHTerms));
  return params;
}

// =============================================================================
// Point transformation
// =============================================================================

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationSV::IntegrateVelocities(double &x, double &y, double &z, double T) const
{
  const double dt = StepLengthForIntervalLength(T);
  if (dt) {
    if      (_IntegrationMethod == FFDIM_FastSS ||
             _IntegrationMethod == FFDIM_SS     ||
             _IntegrationMethod == FFDIM_RKE1)   RKE1  ::Transform(this, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RKE2)   RKE2  ::Transform(this, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RKH2)   RKH2  ::Transform(this, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RK4)    RK4   ::Transform(this, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RKEH12) RKEH12::Transform(this, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKBS23) RKBS23::Transform(this, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKF45)  RKF45 ::Transform(this, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKCK45) RKCK45::Transform(this, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKDP45) RKDP45::Transform(this, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else {
      cerr << "irtkBSplineFreeFormTransformationSV::IntegrateVelocities: Unknown integration method: " << _IntegrationMethod << endl;
      exit(1);
    }
  }
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationSV
::LocalTransform(double &x, double &y, double &z, double t, double t0) const
{
  IntegrateVelocities(x, y, z, + UpperIntegrationLimit(t, t0));
}

// -----------------------------------------------------------------------------
bool irtkBSplineFreeFormTransformationSV
::LocalInverse(double &x, double &y, double &z, double t, double t0) const
{
  IntegrateVelocities(x, y, z, - UpperIntegrationLimit(t, t0));
  return true;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkBSplineFreeFormTransformationSV
::ScalingAndSquaring(irtkGenericImage<VoxelType> *d,
                     double T, const irtkWorldCoordsImage *wc) const
{
  ScalingAndSquaring<VoxelType>(d->Attributes(), d, NULL, NULL, NULL, NULL, T, wc);
}

template void irtkBSplineFreeFormTransformationSV
::ScalingAndSquaring(irtkGenericImage<float> *d,
                     double T, const irtkWorldCoordsImage *) const;

template void irtkBSplineFreeFormTransformationSV
::ScalingAndSquaring(irtkGenericImage<double> *d,
                     double T, const irtkWorldCoordsImage *) const;

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkBSplineFreeFormTransformationSV
::ScalingAndSquaring(const irtkImageAttributes   &a,
                     irtkGenericImage<VoxelType> *d,
                     irtkGenericImage<VoxelType> *dx,
                     irtkGenericImage<VoxelType> *dj,
                     irtkGenericImage<VoxelType> *lj,
                     irtkGenericImage<VoxelType> *dv,
                     double T, const irtkWorldCoordsImage *) const
{
  // Whether to perform a fast scaling-and-squaring on the control point lattice
  const bool fast = (_IntegrationMethod == FFDIM_FastSS);
  // Attributes of output images
  irtkImageAttributes attr(a);
  if (!attr) {
    if      (d ) attr = d ->Attributes();
    else if (dx) attr = dx->Attributes();
    else if (dj) attr = dj->Attributes();
    else if (lj) attr = lj->Attributes();
    else if (dv) attr = dv->Attributes();
  }
  attr._t = 1, attr._dt = .0;
  if (!attr) return;
  // Copy input displacement field
  irtkGenericImage<VoxelType> *din = (d ? new irtkGenericImage<VoxelType>(*d) : NULL);
  // TODO: Improve running time of irtkScalingAndSquaring filter. The previously
  //       used irtkVelocityToDisplacementFieldSS image filter has a considerably
  //       shorter run time (e.g., single-threaded 12s vs 23s, and multi-threaded
  //       3.3s vs 6.7s on MacBook Pro Retina Early 2013, 2.7 GHz Intel Core i7).
  //       The most time consuming step (for the fast scaling and squaring) is
  //       irtkScalingAndSquaring::Resample. Once the running time of the
  //       irtkScalingAndSquaring filter is acceptable, remove the if block
  //       and use the else block only.
  if (d && !dx && !dj && !lj && !dv) {

    irtkGenericImage<VoxelType> v;
    // Use only vector fields defined at control points with B-spline interpolation.
    // This results in an approximate solution due to the error at each squaring.
    if (fast) {
      v.Initialize(this->Attributes(), 3);
      VoxelType *vx = v.Data(0, 0, 0, 0);
      VoxelType *vy = v.Data(0, 0, 0, 1);
      VoxelType *vz = v.Data(0, 0, 0, 2);
      const Vector *vp = _CPImage.Data();
      for (int idx = 0; idx < _CPImage.NumberOfVoxels(); ++idx, ++vx, ++vy, ++vz, ++vp) {
        *vx = static_cast<VoxelType>(vp->_x);
        *vy = static_cast<VoxelType>(vp->_y);
        *vz = static_cast<VoxelType>(vp->_z);
      }
    // Evaluate velocities at output voxels beforehand and use
    // linear interpolation of dense vector fields during squaring
    } else {
      v.Initialize(attr, 3);
      ParallelForEachVoxel(EvaluateBSplineSVFFD3D(this, &v), attr, v);
    }
    // Exponentiate velocity field
    irtkVelocityToDisplacementFieldSS<VoxelType> exp;
    exp.T                (T);
    exp.NumberOfSteps    (NumberOfStepsForIntervalLength(T));
    exp.MaxScaledVelocity(_MaxScaledVelocity);
    exp.Interpolation    (fast ? Interpolation_FastCubicBSpline : Interpolation_Linear);
    exp.Upsample         (false);  // better, but too expensive
    exp.SetInput         (0, &v);  // velocity field to be exponentiated
    exp.SetInput         (1, din); // input displacement field (may be zero)
    exp.SetOutput        (d);      // result is exp(v) o d
    exp.ComputeInterpolationCoefficients(!fast);
    exp.Run();

  } else {

    // Copy B-spline coefficients of velocity field
    irtkGenericImage<VoxelType> v;
    v.Initialize(this->Attributes(), 3);
    VoxelType *vx = v.Data(0, 0, 0, 0);
    VoxelType *vy = v.Data(0, 0, 0, 1);
    VoxelType *vz = v.Data(0, 0, 0, 2);
    const Vector *vp = _CPImage.Data();
    for (int idx = 0; idx < _CPImage.NumberOfVoxels(); ++idx, ++vx, ++vy, ++vz, ++vp) {
      *vx = static_cast<VoxelType>(vp->_x);
      *vy = static_cast<VoxelType>(vp->_y);
      *vz = static_cast<VoxelType>(vp->_z);
    }
    // Exponentiate velocity field
    irtkScalingAndSquaring<VoxelType> exp;
    exp.IntegrationLimit  (T);
    exp.NumberOfSteps     (NumberOfStepsForIntervalLength(T));
    exp.MaxScaledVelocity (_MaxScaledVelocity);
    exp.Interpolation     (fast ? Interpolation_FastCubicBSpline : Interpolation_Linear);
    exp.InterimAttributes (fast ? this->Attributes()             : attr);
    exp.OutputAttributes  (attr);
    exp.Upsample          (false); // better, but too computationally expensive
    exp.InputVelocity     (&v);    // velocity field to be exponentiated
    exp.InputDisplacement (din);   // input displacement field (may be zero)
    exp.OutputDisplacement(d);     // i.e., d = exp(v) o din
    exp.OutputJacobian    (dx);    // i.e., Jacobian
    exp.OutputDetJacobian (dj);    // i.e., det(Jacobian)
    exp.OutputLogJacobian (lj);    // i.e., log(det(Jacobian)
    exp.OutputJacobianDOFs(dv);    // i.e., Jacobian w.r.t. v
    exp.ComputeInterpolationCoefficients(false); // v contains B-spline coefficients
    exp.Run();

  }
  // Free copy of input displacement field
  Delete(din);
}

template void irtkBSplineFreeFormTransformationSV
::ScalingAndSquaring(const irtkImageAttributes &,
                     irtkGenericImage<float> *,
                     irtkGenericImage<float> *,
                     irtkGenericImage<float> *,
                     irtkGenericImage<float> *,
                     irtkGenericImage<float> *,
                     double, const irtkWorldCoordsImage *) const;

template void irtkBSplineFreeFormTransformationSV
::ScalingAndSquaring(const irtkImageAttributes &,
                     irtkGenericImage<double> *,
                     irtkGenericImage<double> *,
                     irtkGenericImage<double> *,
                     irtkGenericImage<double> *,
                     irtkGenericImage<double> *,
                     double, const irtkWorldCoordsImage *) const;

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationSV
::Displacement(irtkGenericImage<float> &d, double t, double t0, const irtkWorldCoordsImage *wc) const
{
  double T;
  if ((T = + UpperIntegrationLimit(t, t0))) {
    IRTK_START_TIMING();
    // Use scaling and squaring method to efficiently compute displacements when possible
    if ((_IntegrationMethod == FFDIM_SS || _IntegrationMethod == FFDIM_FastSS) &&
        ((_z <= 1 && d.GetZ() <= 1) || (_z > 1 && d.GetZ() > 1))) {
      ScalingAndSquaring(&d, T, wc);
    // Evaluate transformation at each voxel separately using numerical integration
    } else {
      irtkTransformation::Displacement(d, t, t0, wc);
    }
    IRTK_DEBUG_TIMING(3, "computation of exp(" << T << "*v)");
  }
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationSV
::Displacement(irtkGenericImage<double> &d, double t, double t0, const irtkWorldCoordsImage *wc) const
{
  double T;
  if ((T = + UpperIntegrationLimit(t, t0))) {
    IRTK_START_TIMING();
    // Use scaling and squaring method to efficiently compute displacements when possible
    if ((_IntegrationMethod == FFDIM_SS || _IntegrationMethod == FFDIM_FastSS) &&
        ((_z <= 1 && d.GetZ() <= 1) || (_z > 1 && d.GetZ() > 1))) {
      ScalingAndSquaring(&d, T, wc);
    // Evaluate transformation at each voxel separately using numerical integration
    } else {
      irtkTransformation::Displacement(d, t, t0, wc);
    }
    IRTK_DEBUG_TIMING(3, "computation of exp(" << T << "*v)");
  }
}

// -----------------------------------------------------------------------------
int irtkBSplineFreeFormTransformationSV
::InverseDisplacement(irtkGenericImage<float> &d, double t, double t0, const irtkWorldCoordsImage *wc) const
{
  double T;
  if ((T = - UpperIntegrationLimit(t, t0))) {
    IRTK_START_TIMING();
    // Use scaling and squaring method to efficiently compute displacements when possible
    if ((_IntegrationMethod == FFDIM_SS || _IntegrationMethod == FFDIM_FastSS) &&
        ((_z <= 1 && d.GetZ() <= 1) || (_z > 1 && d.GetZ() > 1))) {
      ScalingAndSquaring(&d, T, wc);
    // Evaluate transformation at each voxel separately using numerical integration
    } else {
      irtkTransformation::InverseDisplacement(d, t, t0, wc);
    }
    IRTK_DEBUG_TIMING(3, "computation of exp(" << T << "*v)");
  }
  return 0;
}

// -----------------------------------------------------------------------------
int irtkBSplineFreeFormTransformationSV
::InverseDisplacement(irtkGenericImage<double> &d, double t, double t0, const irtkWorldCoordsImage *wc) const
{
  double T;
  if ((T = - UpperIntegrationLimit(t, t0))) {
    IRTK_START_TIMING();
    // Use scaling and squaring method to efficiently compute displacements when possible
    if ((_IntegrationMethod == FFDIM_SS || _IntegrationMethod == FFDIM_FastSS) &&
        ((_z <= 1 && d.GetZ() <= 1) || (_z > 1 && d.GetZ() > 1))) {
      ScalingAndSquaring(&d, T, wc);
    // Evaluate transformation at each voxel separately using numerical integration
    } else {
      irtkTransformation::InverseDisplacement(d, t, t0, wc);
    }
    IRTK_DEBUG_TIMING(3, "computation of exp(" << T << "*v)");
  }
  return 0;
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationSV::LocalJacobian(irtkMatrix &jac, double x, double y, double z, double t, double t0) const
{
  jac.Initialize(3, 3);
  jac.Ident();
  double dt, T;
  if ((dt = StepLengthForIntervalLength(T = UpperIntegrationLimit(t, t0)))) {
    if      (_IntegrationMethod == FFDIM_SS     ||
             _IntegrationMethod == FFDIM_FastSS ||
             _IntegrationMethod == FFDIM_RKE1)   RKE1  ::Jacobian(this, jac, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RKE2)   RKE2  ::Jacobian(this, jac, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RKH2)   RKH2  ::Jacobian(this, jac, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RK4)    RK4   ::Jacobian(this, jac, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RKEH12) RKEH12::Jacobian(this, jac, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKBS23) RKBS23::Jacobian(this, jac, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKF45)  RKF45 ::Jacobian(this, jac, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKCK45) RKCK45::Jacobian(this, jac, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKDP45) RKDP45::Jacobian(this, jac, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else {
      cerr << "irtkBSplineFreeFormTransformationSV::Jacobian: Unknown integration method: " << _IntegrationMethod << endl;
      exit(1);
    }
  }
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationSV::LocalHessian(irtkMatrix [3], double, double, double, double, double) const
{
  cerr << this->NameOfClass() << "::LocalHessian: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationSV::JacobianDOFs(irtkMatrix &jac, int cp, double x, double y, double z, double t, double t0) const
{
  jac.Initialize(3, 3);
  jac.Ident();
  double dt, T;
  if ((dt = StepLengthForIntervalLength(T = UpperIntegrationLimit(t, t0)))) {
    int ci, cj, ck, cl;
    this->IndexToLattice(cp, ci, cj, ck, cl);
    if      (_IntegrationMethod == FFDIM_SS     ||
             _IntegrationMethod == FFDIM_FastSS ||
             _IntegrationMethod == FFDIM_RKE1)   RKE1  ::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RKE2)   RKE2  ::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RKH2)   RKH2  ::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RK4)    RK4   ::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RKEH12) RKEH12::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKBS23) RKBS23::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKF45)  RKF45 ::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKCK45) RKCK45::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKDP45) RKDP45::JacobianDOFs(this, jac, ci, cj, ck, cl, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else {
      cerr << "irtkBSplineFreeFormTransformationSV::JacobianDOFs: Unknown integration method: " << _IntegrationMethod << endl;
      exit(1);
    }
  }
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationSV::JacobianDOFs(double jac[3], int dof, double x, double y, double z, double t, double t0) const
{
  irtkMatrix dTdp(3, 3);
  this->JacobianDOFs(dTdp, dof / 3, x, y, z, t, t0);
  const int c = dof % 3;
  jac[0] = dTdp(0, c);
  jac[1] = dTdp(1, c);
  jac[2] = dTdp(2, c);
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationSV
::EvaluateJacobianDOFs(irtkTransformationJacobian &jac, double x, double y) const
{
  int i = static_cast<int>(floor(x));
  int j = static_cast<int>(floor(y));

  const int A = Kernel::VariableToIndex(x - i);
  const int B = Kernel::VariableToIndex(y - j);

  --i, --j;

  double wxy, wy;
  int    ci, cj, xdof, ydof;

  for (int b = 0; b <= 3; ++b) {
    cj = j + b;
    if (cj < 0 || cj >= _y) continue;
    wy = Kernel::LookupTable[B][b];
    for (int a = 0; a <= 3; ++a) {
      ci = i + a;
      if (ci < 0 || ci >= _x) continue;
      wxy = Kernel::LookupTable[A][a] * wy;
      IndexToDOFs(LatticeToIndex(ci, cj), xdof, ydof);
      jac(xdof)._x = jac(ydof)._y = wxy;
    }
  }
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationSV
::EvaluateJacobianDOFs(irtkTransformationJacobian &jac, double x, double y, double z) const
{
  int i = static_cast<int>(floor(x));
  int j = static_cast<int>(floor(y));
  int k = static_cast<int>(floor(z));

  const int A = Kernel::VariableToIndex(x - i);
  const int B = Kernel::VariableToIndex(y - j);
  const int C = Kernel::VariableToIndex(z - k);

  --i, --j, --k;

  double wxyz, wyz, wz;
  int    ci, cj, ck, xdof, ydof, zdof;

  for (int c = 0; c <= 3; ++c) {
    ck = k + c;
    if (ck < 0 || ck >= _z) continue;
    wz = Kernel::LookupTable[C][c];
    for (int b = 0; b <= 3; ++b) {
      cj = j + b;
      if (cj < 0 || cj >= _y) continue;
      wyz = Kernel::LookupTable[B][b] * wz;
      for (int a = 0; a <= 3; ++a) {
        ci = i + a;
        if (ci < 0 || ci >= _x) continue;
        wxyz = Kernel::LookupTable[A][a] * wyz;
        IndexToDOFs(LatticeToIndex(ci, cj, ck), xdof, ydof, zdof);
        jac(xdof)._x = jac(ydof)._y = jac(zdof)._z = wxyz;
      }
    }
  }
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationSV::JacobianDOFs(irtkTransformationJacobian &jac, double x, double y, double z, double t, double t0) const
{
  jac.Clear();
  double dt, T;
  if ((dt = StepLengthForIntervalLength(T = UpperIntegrationLimit(t, t0)))) {
    if      (_IntegrationMethod == FFDIM_SS     ||
             _IntegrationMethod == FFDIM_FastSS ||
             _IntegrationMethod == FFDIM_RKE1)   RKE1  ::JacobianDOFs(this, jac, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RKE2)   RKE2  ::JacobianDOFs(this, jac, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RKH2)   RKH2  ::JacobianDOFs(this, jac, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RK4)    RK4   ::JacobianDOFs(this, jac, x, y, z, .0, T, dt);
    else if (_IntegrationMethod == FFDIM_RKEH12) RKEH12::JacobianDOFs(this, jac, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKBS23) RKBS23::JacobianDOFs(this, jac, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKF45)  RKF45 ::JacobianDOFs(this, jac, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKCK45) RKCK45::JacobianDOFs(this, jac, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else if (_IntegrationMethod == FFDIM_RKDP45) RKDP45::JacobianDOFs(this, jac, x, y, z, .0, T, 0.5 * dt, 2.0 * dt, SVFFD_RKTOL);
    else {
      cerr << "irtkBSplineFreeFormTransformationSV::JacobianDOFs: Unknown integration method: " << _IntegrationMethod << endl;
      exit(1);
    }
  }
}

namespace irtkBSplineFreeFormTransformationSVUtils {

// -----------------------------------------------------------------------------
struct MultiplyDerivatives : public irtkVoxelFunction
{
  const int x, y, z, xx, xy, xz, yx, yy, yz, zx, zy, zz; // offsets

  MultiplyDerivatives(int n)
  :
    x(0), y(x+n), z(y+n),
    xx(   0), xy(xx+n), xz(xy+n),
    yx(xz+n), yy(yx+n), yz(yy+n),
    zx(yz+n), zy(zx+n), zz(zy+n)
  {}

  void operator ()(const irtkGenericImage<double> &, int, const double *in, const double *d, double *out) const
  {
    // Attention: out can be equal d, therefore use temporary variables
    double gx = in[x] * d[xx] + in[x] * d[xy] + in[x] * d[xz];
    double gy = in[y] * d[yx] + in[y] * d[yy] + in[y] * d[yz];
    double gz = in[z] * d[zx] + in[z] * d[zy] + in[z] * d[zz];
    out[x] = gx, out[y] = gy, out[z] = gz;
  }
};

// -----------------------------------------------------------------------------
struct MultiplyApproximateDerivatives : public irtkVoxelFunction
{
  typedef irtkGenericFastCubicBSplineInterpolateImageFunction<irtkGenericImage<double> > JacobianDOFs;

  enum { xx, xy, xz, yx, yy, yz, zx, zy, zz }; // offsets
  const int x, y, z;                           // offsets

  irtkGenericImage<double> *_Output;
  const JacobianDOFs       *_JacobianDOFs;

  MultiplyApproximateDerivatives(const JacobianDOFs       *dv,
                                 irtkGenericImage<double> *out)
  :
    x(0), y(x+out->NumberOfSpatialVoxels()), z(y+out->NumberOfSpatialVoxels()),
    _Output(out), _JacobianDOFs(dv)
  {}

  void operator ()(int i, int j, int k, int, const double *in, double *out) const
  {
    irtkPoint p;
    double d[9];
    p._x = i, p._y = j, p._z = k;
    _Output      ->ImageToWorld(p);
    _JacobianDOFs->WorldToImage(p);
    _JacobianDOFs->Evaluate(d,  p._x, p._y, p._z);
    out[x] = in[x] * d[xx] + in[x] * d[xy] + in[x] * d[xz];
    out[y] = in[y] * d[yx] + in[y] * d[yy] + in[y] * d[yz];
    out[z] = in[z] * d[zx] + in[z] * d[zy] + in[z] * d[zz];
  }

  template <class TCoord>
  void operator ()(int i, int j, int k, int, const TCoord *wc, const double *in, double *out) const
  {
    irtkPoint p;
    double d[9];
    p._x = wc[x], p._y = wc[y], p._z = wc[z];
    _JacobianDOFs->WorldToImage(p);
    _JacobianDOFs->Evaluate(d,  p._x, p._y, p._z);
    out[x] = in[x] * d[xx] + in[x] * d[xy] + in[x] * d[xz];
    out[y] = in[y] * d[yx] + in[y] * d[yy] + in[y] * d[yz];
    out[z] = in[z] * d[zx] + in[z] * d[zy] + in[z] * d[zz];
  }
};

// -----------------------------------------------------------------------------
struct MultiplyPointWiseDerivatives
{
  enum { xx, xy, xz, yx, yy, yz, zx, zy, zz };

  const irtkPointSet         *_PointSet;
  const irtkVector3D<double> *_Input;
  irtkVector3D<double>       *_Output;

  const irtkGenericFastCubicBSplineInterpolateImageFunction<irtkGenericImage<double> > *_JacobianDOFs;

  void operator ()(const blocked_range<int> &ids) const
  {
    irtkPoint p;
    double d[9];
    for (int id = ids.begin(); id != ids.end(); ++id) {
      _PointSet->GetPoint(id, p);
      _JacobianDOFs->WorldToImage(p);
      _JacobianDOFs->Evaluate(d, p._x, p._y, p._z);
      _Output[id]._x = _Input[id]._x * d[xx] + _Input[id]._x * d[xy] + _Input[id]._x * d[xz];
      _Output[id]._y = _Input[id]._y * d[yx] + _Input[id]._y * d[yy] + _Input[id]._y * d[yz];
      _Output[id]._z = _Input[id]._z * d[zx] + _Input[id]._z * d[zy] + _Input[id]._z * d[zz];
    }
  }
};

} // namespace irtkBSplineFreeFormTransformationSVUtils
using namespace irtkBSplineFreeFormTransformationSVUtils;

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationSV
::ParametricGradient(const irtkGenericImage<double> *in, double *out,
                     const irtkWorldCoordsImage *i2w, const irtkWorldCoordsImage *wc,
                     double t, double t0, double w) const
{
  // Upper integration limit for given interval
  const double T = UpperIntegrationLimit(t, t0);
  if (T == .0) return;

  // ---------------------------------------------------------------------------
  // BCH based velocity update computation
  if (_NumberOfBCHTerms > 1) {

    IRTK_START_TIMING();
    // Apply chain rule to compute gradient w.r.t control point displacements
    CPImage d(this->Attributes());
    DOFValue * const grd = reinterpret_cast<DOFValue *>(d.Data());
    irtkBSplineFreeFormTransformation3D::ParametricGradient(in, grd, i2w, wc, t0, 1.0);
    // Convert displacement gradient to velocity update
    EvaluateBCHFormula(_NumberOfBCHTerms, d, T, _CPImage, d, true);
    // Adjust weight as update field is computed for tau * v, i.e.,
    //   exp(tau * v_{i+1}) = exp(tau v_i) o exp(\delta u)
    //   ==> v_{i+1} = log(exp(tau * v_{i+1})) / tau
    w /= T;
    // Add weighted gradient to total energy gradient
    for (int dof = 0; dof < this->NumberOfDOFs(); ++dof) out[dof] += w * grd[dof];
    IRTK_DEBUG_TIMING(2, "parametric gradient computation (BCH)");

  // ---------------------------------------------------------------------------
  // Scaling and squaring based gradient computation
  } else if (_IntegrationMethod == FFDIM_FastSS) {

    IRTK_START_TIMING();

    // Compute derivative of transformation T = exp(v) w.r.t. v
    if (_JacobianDOFsIntervalLength != T || !_JacobianDOFs) {
      if (!_JacobianDOFs) _JacobianDOFs = new irtkGenericImage<double>();
      ScalingAndSquaring<double>(_attr, NULL, NULL, NULL, NULL, _JacobianDOFs, T);
      _JacobianDOFsIntervalLength = T;
    }

    // Initialize interpolator for evaluation of derivatives at non-CP locations
    irtkGenericFastCubicBSplineInterpolateImageFunction<irtkGenericImage<double> > dv;
    dv.Input(_JacobianDOFs);
    dv.Initialize();

    // Multiply input derivatives w.r.t. T by the derivative of T w.r.t. v
    irtkGenericImage<double> *tmp = new irtkGenericImage<double>(in->Attributes());
    MultiplyApproximateDerivatives mul(&dv, tmp);
    if (wc) {
      ParallelForEachVoxel(in->Attributes(), wc, in, tmp, mul);
    } else {
      ParallelForEachVoxel(in->Attributes(),     in, tmp, mul);
    }

    // Multiply resulting vectors by derivative of v w.r.t. the DoFs
    irtkBSplineFreeFormTransformation3D::ParametricGradient(tmp, out, i2w, wc, t0, w);
    delete tmp;
    IRTK_DEBUG_TIMING(2, "parametric gradient computation (FastSS)");

  } else if (_IntegrationMethod == FFDIM_SS) {

    IRTK_START_TIMING();

    // Compute derivative of transformation T = exp(v) w.r.t. v
    irtkGenericImage<double> dv;
    ScalingAndSquaring<double>(in->Attributes(), NULL, NULL, NULL, NULL, &dv, T, wc);

    // Multiply input derivatives w.r.t. T by the derivative of T w.r.t. v
    MultiplyDerivatives mul(in->NumberOfSpatialVoxels());
    ParallelForEachVoxel(in, &dv, &dv, mul);

    // Multiply resulting vectors by derivative of v w.r.t. the DoFs
    irtkBSplineFreeFormTransformation3D::ParametricGradient(&dv, out, wc, t0, w);
    IRTK_DEBUG_TIMING(2, "parametric gradient computation (SS)");

  // ---------------------------------------------------------------------------
  // Runge-Kutta integration based gradient computation similar to TD FFD
  // transformation parameterized by non-stationary velocity field
  } else {

    // Note: T = in->GetTOrigin() - t0
    irtkFreeFormTransformation::ParametricGradient(in, out, i2w, wc, in->GetTOrigin() - T, w);

  }
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationSV
::ParametricGradient(const irtkPointSet &pos, const irtkVector3D<double> *in,
                     double *out, double t, double t0, double w) const
{
  // ---------------------------------------------------------------------------
  // Scaling and squaring based gradient computation for dense point clouds
  if (_IntegrationMethod == FFDIM_FastSS || _IntegrationMethod == FFDIM_SS) {

    IRTK_START_TIMING();
 
    // Upper integration limit for given interval
    const double T = UpperIntegrationLimit(t, t0);
    if (T == .0) return;

    // Compute derivative of transformation T = exp(v) w.r.t. v
    if (_JacobianDOFsIntervalLength != T || !_JacobianDOFs) {
      irtkImageAttributes attr = _attr;
      if (_IntegrationMethod == FFDIM_SS) {
        // FIXME: Should be more adaptive and not specific to typical image
        //        resolution and size encountered in MR brain imaging
        attr._dx = (_dx > 1.0 ? 1.0 : .0);
        attr._dy = (_dy > 1.0 ? 1.0 : .0);
        attr._dz = (_dz > 1.0 ? 1.0 : .0);
        attr._x  = (attr._dx > .0 ? static_cast<int>(ceil(_x * _dx / attr._dx)) : 1);
        attr._y  = (attr._dy > .0 ? static_cast<int>(ceil(_y * _dy / attr._dy)) : 1);
        attr._z  = (attr._dz > .0 ? static_cast<int>(ceil(_z * _dz / attr._dz)) : 1);
        if (attr._x > 256) attr._x = 256;
        if (attr._y > 256) attr._y = 256;
        if (attr._z > 256) attr._z = 256;
        attr._dx = (_x * _dx / attr._x);
        attr._dy = (_y * _dy / attr._y);
        attr._dz = (_z * _dz / attr._z);
      }
      if (!_JacobianDOFs) _JacobianDOFs = new irtkGenericImage<double>();
      _JacobianDOFsIntervalLength = T;
      ScalingAndSquaring<double>(attr, NULL, NULL, NULL, NULL, _JacobianDOFs, T);
    }

    // Initialize interpolator for evaluation of derivatives at non-CP locations
    irtkGenericFastCubicBSplineInterpolateImageFunction<irtkGenericImage<double> > dv;
    dv.Input(_JacobianDOFs);
    dv.Initialize();
 
    // Multiply input derivatives w.r.t. T by the derivative of T w.r.t. v
    MultiplyPointWiseDerivatives mul;
    mul._PointSet     = &pos;
    mul._Input        = in;
    mul._JacobianDOFs = &dv;
    mul._Output       = new irtkVector3D<double>[pos.Size()];
    parallel_for(blocked_range<int>(0, pos.Size()), mul);

    // Multiply resulting vectors by derivative of v w.r.t. the DoFs
    irtkBSplineFreeFormTransformation3D::ParametricGradient(pos, mul._Output, out, t, t0, w);
    delete[] mul._Output;

    IRTK_DEBUG_TIMING(2, "point-wise parametric gradient computation (" << ToString(_IntegrationMethod) << ")");

  // ---------------------------------------------------------------------------
  // Runge-Kutta integration based gradient computation similar to TD FFD
  // transformation parameterized by non-stationary velocity field
  } else {

    irtkFreeFormTransformation::ParametricGradient(pos, in, out, t, t0, w);

  }
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationSV::Print(irtkIndent indent) const
{
  cout << indent << "B-spline SV FFD:" << endl;
  indent++;
  // Print FFD attributes
  irtkFreeFormTransformation3D::Print(indent);
  // Change output stream settings
  const streamsize         w = cout.width    (0);
  const streamsize         p = cout.precision(2);
  const ios_base::fmtflags f = cout.flags();
  cout.unsetf(ios_base::floatfield);
  // Print SV FFD parameters
  cout << indent << "Integration method:                " << setw(6) << ToString(_IntegrationMethod) << endl;
  cout << indent << "Cross-sectional time interval:     " << setw(6) << _T << endl;
  cout << indent << "Time unit of integration interval: " << setw(6) << _TimeUnit << endl;
  cout << indent << "Maximum scaled velocity:           " << setw(6) << _MaxScaledVelocity << endl;
  cout << indent << "No. of integration steps per unit: " << setw(6) << _NumberOfSteps << endl;
  cout << indent << "No. of cross-sectional steps:      " << setw(6) << NumberOfStepsForIntervalLength(_T) << endl;
  cout << indent << "No. of BCH terms:                  " << setw(6) << _NumberOfBCHTerms << endl;
  cout << indent << "Use Lie derivative:                " << setw(6) << ToString(_LieDerivative) << endl;
  // Restore output stream settings
  cout.width    (w);
  cout.precision(p);
  cout.flags    (f);
}

// -----------------------------------------------------------------------------
bool irtkBSplineFreeFormTransformationSV::CanRead(irtkTransformationType format) const
{
  switch (format) {
    case IRTKTRANSFORMATION_BSPLINE_FFD_SV_v1:
    case IRTKTRANSFORMATION_BSPLINE_FFD_SV_v2:
    case IRTKTRANSFORMATION_BSPLINE_FFD_SV_v3:
    case IRTKTRANSFORMATION_BSPLINE_FFD_SV_v4:
    case IRTKTRANSFORMATION_BSPLINE_FFD_SV_v5:
    case IRTKTRANSFORMATION_BSPLINE_FFD_SV_v6:
    case IRTKTRANSFORMATION_BSPLINE_FFD_SV_v7:
    case IRTKTRANSFORMATION_BSPLINE_FFD_SV_v8:
      return true;
    default:
      return false;
  }
}

// -----------------------------------------------------------------------------
irtkCifstream &irtkBSplineFreeFormTransformationSV::ReadDOFs(irtkCifstream &from, irtkTransformationType format)
{
  // Read FFD data
  switch (format) {
    case IRTKTRANSFORMATION_BSPLINE_FFD_SV_v1:
    case IRTKTRANSFORMATION_BSPLINE_FFD_SV_v2:
    case IRTKTRANSFORMATION_BSPLINE_FFD_SV_v3:
    case IRTKTRANSFORMATION_BSPLINE_FFD_SV_v4:
    case IRTKTRANSFORMATION_BSPLINE_FFD_SV_v5:
    case IRTKTRANSFORMATION_BSPLINE_FFD_SV_v6:
      irtkBSplineFreeFormTransformation3D::ReadDOFs(from, IRTKTRANSFORMATION_BSPLINE_FFD_3D_v2);
    case IRTKTRANSFORMATION_BSPLINE_FFD_SV_v7:
      irtkBSplineFreeFormTransformation3D::ReadDOFs(from, IRTKTRANSFORMATION_BSPLINE_FFD_3D_v3);
    default:
      irtkBSplineFreeFormTransformation3D::ReadDOFs(from, IRTKTRANSFORMATION_BSPLINE_FFD_3D);
  }

  // Read number of integration steps
  from.ReadAsInt(&_NumberOfSteps, 1);

  if (format == IRTKTRANSFORMATION_BSPLINE_FFD_SV_v1) return from;

  // Read upper integration limit
  from.ReadAsDouble(&_T, 1);
  // Read number of BCH terms to use for update
  from.ReadAsInt(&_NumberOfBCHTerms, 1);

  if (format == IRTKTRANSFORMATION_BSPLINE_FFD_SV_v2) return from;

  // Read time unit of integration interval
  from.ReadAsDouble(&_TimeUnit, 1);

  if (format == IRTKTRANSFORMATION_BSPLINE_FFD_SV_v3) return from;

  if (format <= IRTKTRANSFORMATION_BSPLINE_FFD_SV_v6) {

    // Whether to use scaling and squaring
    char useSS;
    from.ReadAsChar(&useSS, 1);

    // Maximum scaled velocity
    from.ReadAsDouble(&_MaxScaledVelocity, 1);

    if (format == IRTKTRANSFORMATION_BSPLINE_FFD_SV_v4) return from;

    // Whether to use fast scaling and squaring
    char fastSS;
    from.ReadAsChar(&fastSS, 1);

    // Set integration method
    _IntegrationMethod = (useSS ? (fastSS ? FFDIM_FastSS : FFDIM_SS) : FFDIM_RKE1);

  } else {

    // Integration method
    unsigned int integration_method;
    from.ReadAsUInt(&integration_method, 1);
    _IntegrationMethod = static_cast<FFDIntegrationMethod>(integration_method);

    // Maximum scaled velocity
    from.ReadAsDouble(&_MaxScaledVelocity, 1);

  }

  return from;
}

// -----------------------------------------------------------------------------
irtkCofstream &irtkBSplineFreeFormTransformationSV::WriteDOFs(irtkCofstream &to) const
{
  // Write FFD data
  irtkBSplineFreeFormTransformation3D::WriteDOFs(to);

  // Write number of integration steps
  to.WriteAsInt(&_NumberOfSteps, 1);
  // Write upper integration limit
  to.WriteAsDouble(&_T, 1);
  // Write number of BCH terms to use for update
  to.WriteAsInt(&_NumberOfBCHTerms, 1);
  // Write time unit of integration interval
  to.WriteAsDouble(&_TimeUnit, 1);
  // Integration method
  const unsigned int integration_method = _IntegrationMethod;
  to.WriteAsUInt(&integration_method, 1);
  // Maximum scaled velocity
  to.WriteAsDouble(&_MaxScaledVelocity, 1);

  return to;
}
