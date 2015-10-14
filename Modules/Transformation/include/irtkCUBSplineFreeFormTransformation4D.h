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

#ifndef IRTKCUBSPLINEFREEFORMTRANSFORMATION4D_H_
#define IRTKCUBSPLINEFREEFORMTRANSFORMATION4D_H_


#include <irtkBSplineFunction.h>
#include <irtkCUJacobianDOFs.h>


class irtkBSplineFreeFormTransformation4D;


/**
 * Spatio-temporal B-spline FFD representation for use by device code.
 *
 * \note No polymorphism, i.e., virtual methods can be used by device code!
 *
 * \note Keep number of data members at a minimum to reduce the memory
 *       footprint of the structure so it can be passed by value to a
 *       CUDA kernel. Current size limit for all kernel arguments is 256 bytes.
 */
class irtkCUBSplineFreeFormTransformation4D : public irtkCUFreeFormTransformation4D
{
public:

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  IRTKCU_HOST_API irtkCUBSplineFreeFormTransformation4D();

  /// Constructor
  IRTKCU_HOST_API irtkCUBSplineFreeFormTransformation4D(const irtkBSplineFreeFormTransformation4D &);

  /// Copy constructor, makes a shallow copy only
  IRTKCU_HOST_API irtkCUBSplineFreeFormTransformation4D(const irtkCUBSplineFreeFormTransformation4D &);

  /// Destructor
  virtual IRTKCU_HOST_API ~irtkCUBSplineFreeFormTransformation4D();

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Evaluate FFD at given lattice coordinates without boundary checks
  IRTKCU_API void EvaluateFFD(real3 &, realt) const;

  /// Evaluate FFD at given lattice coordinates using NN extrapolation
  IRTKCU_API void FFD(real3 &, realt) const;

  /// Compute Jacobian of FFD at given lattice coordinates
  IRTKCU_API void FFDJacobian(real3x3 &, real4 p) const;

  /// Compute Jacobian of FFD at given lattice coordinates
  IRTKCU_API void FFDJacobian(real3x3 &, real3 p, realt t) const;

  /// Compute Jacobian of FFD in word coordinates at given lattice coordinates
  IRTKCU_API void FFDJacobianWorld(real3x3 &, real4 p) const;

  /// Compute Jacobian of FFD in word coordinates at given lattice coordinates
  IRTKCU_API void FFDJacobianWorld(real3x3 &, real3 p, realt t) const;

  /// Compute Jacobian of FFD w.r.t control point at given lattice coordinates
  IRTKCU_API void FFDJacobianDOFs(real3 &, uint4 cp, real4 p) const;

  /// Compute Jacobian of FFD w.r.t control point at given lattice coordinates
  IRTKCU_API void FFDJacobianDOFs(real3 &, uint4 cp, real3 p, realt t) const;

  /// Compute Jacobian of FFD w.r.t parameters at given lattice coordinates
  IRTKCU_API void FFDJacobianDOFs(irtkCUJacobianDOFs &, real4 p) const;

  /// Compute Jacobian of FFD w.r.t parameters at given lattice coordinates
  IRTKCU_API void FFDJacobianDOFs(irtkCUJacobianDOFs &, real3 p, realt t) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
inline irtkCUBSplineFreeFormTransformation4D::irtkCUBSplineFreeFormTransformation4D()
:
  irtkCUFreeFormTransformation4D()
{
}

// -----------------------------------------------------------------------------
inline irtkCUBSplineFreeFormTransformation4D
::irtkCUBSplineFreeFormTransformation4D(const irtkBSplineFreeFormTransformation4D &ffd)
:
  irtkCUFreeFormTransformation4D(ffd)
{
}

// -----------------------------------------------------------------------------
inline irtkCUBSplineFreeFormTransformation4D
::irtkCUBSplineFreeFormTransformation4D(const irtkCUBSplineFreeFormTransformation4D &other)
:
  irtkCUFreeFormTransformation4D(other)
{
}

// -----------------------------------------------------------------------------
inline irtkCUBSplineFreeFormTransformation4D::~irtkCUBSplineFreeFormTransformation4D()
{
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
inline void irtkCUBSplineFreeFormTransformation4D::EvaluateFFD(real3 &p, realt t) const
{
  realt B_L, B_K, B_J;

  const int3  ip = make_int3(p);
  const int   it = int      (t);
  const real3 fp = p - make_float3(ip);
  const realt ft = t - realt      (it);

  const real3 *iter   = GetPointerToControlPoints(ip-1, it-1);
  const uint4  stride = GetControlPointsStride();

  p = make_real3(.0f);
  for (int l = 0; l < 4; l++) {
    B_L = irtkBSpline<realt>::B(l, ft);
    for (int k = 0; k < 4; k++) {
      B_K = irtkBSpline<realt>::B(k, fp.z) * B_L;
      for (int j = 0; j < 4; j++) {
        B_J = irtkBSpline<realt>::B(j, fp.y) * B_K;
        // begin: innermost loop unrolled
        p += (*iter) * irtkBSpline<realt>::B0(fp.x) * B_J; iter += stride.x;
        p += (*iter) * irtkBSpline<realt>::B1(fp.x) * B_J; iter += stride.x;
        p += (*iter) * irtkBSpline<realt>::B2(fp.x) * B_J; iter += stride.x;
        p += (*iter) * irtkBSpline<realt>::B3(fp.x) * B_J; iter += stride.x;
        // end: innermost loop unrolled
        iter += stride.y;
      }
      iter += stride.z;
    }
    iter += stride.w;
  }
}

// -----------------------------------------------------------------------------
inline void irtkCUBSplineFreeFormTransformation4D::FFD(real3 &p, realt t) const
{
  if      (p.x <    0.0f) p.x = 0.0f;
  else if (p.x >= _dim.x) p.x = _dim.x - 1.0f;
  if      (p.y <    0.0f) p.y = 0.0f;
  else if (p.y >= _dim.y) p.y = _dim.y - 1.0f;
  if      (p.z <    0.0f) p.z = 0.0f;
  else if (p.z >= _dim.z) p.z = _dim.z - 1.0f;
  if      (t   <    0.0f) t   = 0.0f;
  else if (t   >= _dim.w) t   = _dim.w - 1.0f;
  EvaluateFFD(p, t);
}

// -----------------------------------------------------------------------------
inline void irtkCUBSplineFreeFormTransformation4D::FFDJacobian(real3x3 &jac, real4 p) const
{
  jac = .0f;

  if (p.x < -2 || p.x > _dim.x+1 ||
      p.y < -2 || p.y > _dim.y+1 ||
      p.z < -2 || p.z > _dim.z+1 ||
      p.w < -2 || p.w > _dim.w+1) return;

  // Partial derivatives with respect to time are not provided by this function.
  // I.e. none of dTx/dt , dTy/dt or dTz/dt is provided (they should be zero, right?).

  // We are only returning the first three columns of the Jacobian
  // (the full Jacobian is a 3x4 matrix and we return a 3x3 one)

  const int4  c = make_int4(floor(p));
  const real4 u = p - make_real4(c);

  const real3 *iter   = GetPointerToControlPoints(c - 1);
  const uint4  stride = GetControlPointsStride();

  realt b, B_I, B_J, B_K, B_L, B_I_I, B_J_I, B_K_I;
  for (int l = 0; l < 4; l++) {
    B_L = irtkBSpline<realt>::B(l, u.w);
    for (int k = 0; k < 4; k++) {
      B_K   = irtkBSpline<realt>::B  (k, u.z);
      B_K_I = irtkBSpline<realt>::B_I(k, u.z);
      for (int j = 0; j < 4; j++) {
        B_J   = irtkBSpline<realt>::B  (j, u.y);
        B_J_I = irtkBSpline<realt>::B_I(j, u.y);
        for (int i = 0; i < 4; i++) {
          B_I   = irtkBSpline<realt>::B  (i, u.x);
          B_I_I = irtkBSpline<realt>::B_I(i, u.x);
          b = B_I_I * B_J * B_K * B_L;
          jac.a.x += b * iter->x;
          jac.b.x += b * iter->y;
          jac.c.x += b * iter->z;
          b = B_I * B_J_I * B_K * B_L;
          jac.a.y += b * iter->x;
          jac.b.y += b * iter->y;
          jac.c.y += b * iter->z;
          b = B_I * B_J * B_K_I * B_L;
          jac.a.z += b * iter->x;
          jac.b.z += b * iter->y;
          jac.c.z += b * iter->z;
          iter += stride.x;
        }
        iter += stride.y;
      }
      iter += stride.z;
    }
    iter += stride.w;
  }
}

// -----------------------------------------------------------------------------
inline void irtkCUBSplineFreeFormTransformation4D::FFDJacobian(real3x3 &jac, real3 p, realt t) const
{
  FFDJacobian(jac, make_real4(p.x, p.y, p.z, t));
}

// -----------------------------------------------------------------------------
inline void irtkCUBSplineFreeFormTransformation4D::FFDJacobianWorld(real3x3 &jac, real4 p) const
{
  FFDJacobian(jac, p);
  jac = jac * make_real3x3(_matW2L);
}

// -----------------------------------------------------------------------------
inline void irtkCUBSplineFreeFormTransformation4D::FFDJacobianWorld(real3x3 &jac, real3 p, realt t) const
{
  FFDJacobianWorld(jac, make_real4(p.x, p.y, p.z, t));
}

// -----------------------------------------------------------------------------
inline void irtkCUBSplineFreeFormTransformation4D::FFDJacobianDOFs(real3 &jac, uint4 cp, real4 p) const
{
  jac.x = jac.y = jac.z =   irtkBSpline<realt>::B(p.x - cp.x)
                          * irtkBSpline<realt>::B(p.y - cp.y)
                          * irtkBSpline<realt>::B(p.z - cp.z)
                          * irtkBSpline<realt>::B(p.w - cp.w);
}

// -----------------------------------------------------------------------------
inline void irtkCUBSplineFreeFormTransformation4D::FFDJacobianDOFs(real3 &jac, uint4 cp, real3 p, realt t) const
{
  FFDJacobianDOFs(jac, cp, make_real4(p.x, p.y, p.z, t));
}

// -----------------------------------------------------------------------------
inline void irtkCUBSplineFreeFormTransformation4D::FFDJacobianDOFs(irtkCUJacobianDOFs &jac, real4 p) const
{
  const uint numcp = _dim.x * _dim.y * _dim.z * _dim.w;
  const int4 o     = make_int4(floor(p)) - 1;

  int4  c;
  real4 b;

  for (int dw = 0; dw < 4; dw++) {
    c.w = o.w + dw;
    if (0 <= c.w && c.w < _dim.w) {
      b.w = irtkBSpline<realt>::B(p.w - c.w);
      for (int dz = 0; dz < 4; dz++) {
        c.z = o.z + dz;
        if (0 <= c.z && c.z < _dim.z) {
          b.z = irtkBSpline<realt>::B(p.z - c.z) * b.w;
          for (int dy = 0; dy < 4; dy++) {
            c.y = o.y + dy;
            if (0 <= c.y && c.y < _dim.y) {
              b.y = irtkBSpline<realt>::B(p.y - c.y) * b.z;
              for (int dx = 0; dx < 4; dx++) {
                c.x = o.x + dx;
                if (0 <= c.x && c.x < _dim.x) {
                  const int cp = this->LatticeToIndex(c);
                  b.x = irtkBSpline<realt>::B(p.x - c.x) * b.y;
                  jac(cp).x = jac(cp + numcp).y = jac(cp + (numcp<<1)).z = b.x;
                }
              }
            }
          }
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
inline void irtkCUBSplineFreeFormTransformation4D::FFDJacobianDOFs(irtkCUJacobianDOFs &jac, real3 p, realt t) const
{
  FFDJacobianDOFs(jac, make_real4(p.x, p.y, p.z, t));
}


#endif
