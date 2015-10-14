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

#ifndef IRTKCUFREEFORMTRANSFORMATIONRUNGEKUTTA_H_
#define IRTKCUFREEFORMTRANSFORMATIONRUNGEKUTTA_H_


// TODO Change velocities to be given in lattice units rather than world units.
//      This reduces the number of conversions from world to lattice coordinates
//      during the numerical integration.
//      -as12312

#include <irtkCUCommon.h>
#include <irtkCUMath.h>
#include <irtkCUJacobianDOFs.h>


// =============================================================================
// Base class of Runge-Kutta integration methods
// =============================================================================

/**
 * Base class of integration methods used by FFDs which are parameterized by a
 * velocity field to compute the displacements.
 */
class irtkCUFreeFormTransformationRungeKutta
{
protected:

  // ---------------------------------------------------------------------------
  /// Helper for computation of Jacobian w.r.t control point
  static IRTKCU_API void dkdp(real3x3 &dk, const real3x3 &Dv, const real3x3 &dx, const real3 &dv, realt h)
  {
    dk.a.x = (Dv.a.x * dx.a.x + Dv.a.y * dx.b.x + Dv.a.z * dx.c.x + dv.x) * h;
    dk.a.y = (Dv.a.x * dx.a.y + Dv.a.y * dx.b.y + Dv.a.z * dx.c.y       ) * h;
    dk.a.z = (Dv.a.x * dx.a.z + Dv.a.y * dx.b.z + Dv.a.z * dx.c.z       ) * h;
    dk.b.x = (Dv.b.x * dx.a.x + Dv.b.y * dx.b.x + Dv.b.z * dx.c.x       ) * h;
    dk.b.y = (Dv.b.x * dx.a.y + Dv.b.y * dx.b.y + Dv.b.z * dx.c.y + dv.y) * h;
    dk.b.z = (Dv.b.x * dx.a.z + Dv.b.y * dx.b.z + Dv.b.z * dx.c.z       ) * h;
    dk.c.x = (Dv.c.x * dx.a.x + Dv.c.y * dx.b.x + Dv.c.z * dx.c.x       ) * h;
    dk.c.y = (Dv.c.x * dx.a.y + Dv.c.y * dx.b.y + Dv.c.z * dx.c.y       ) * h;
    dk.c.z = (Dv.c.x * dx.a.z + Dv.c.y * dx.b.z + Dv.c.z * dx.c.z + dv.z) * h;
  }

  // -------------------------------------------------------------------------
  /// Helper for computation of Jacobian w.r.t control point
  static IRTKCU_API void dkdp(irtkCUJacobianDOFs &dk, const real3x3 &Dv, const irtkCUJacobianDOFs &dx, const irtkCUJacobianDOFs &dv, realt h)
  {
    dk  = dx;
    dk *= Dv; //  Dv * dx
    dk += dv; //  Dv * dx + dv
    dk *= h;  // (Dv * dx + dv) * h
  }

}; // irtkCUFreeFormTransformationRungeKutta


// =============================================================================
// Generic implementation of explicit Runge-Kutta
// =============================================================================

/**
 * Explicit Runge-Kutta integration method for FFD parameterized by velocity field.
 * 
 * The first template argument is the type of the FFD which has to implement
 * the following methods:
 *
 * - FFDEvaluateWorld
 * - FFDJacobianWorld
 * - FFDJacobianDOFs
 *
 * The second template argument is a structure of the respective Butcher tableau
 * of the Runge-Kutta method which has been defined using the macro
 * IRTK_DEFINE_FFDRK_EXPLICIT.
 */
template <class TFreeFormTransformation, class TButcherTableau>
class irtkCUFreeFormTransformationExplicitRungeKutta
:
  public irtkCUFreeFormTransformationRungeKutta
{
public:

  typedef TFreeFormTransformation FFD; ///< Short-hand for FFD type
  typedef TButcherTableau         BT;  ///< Short-hand for Butcher tableau

  // -----------------------------------------------------------------------------
  static IRTKCU_API void Transform(const FFD &v, real3 &p, realt t1, realt t2, realt dt)
  {
    if (t1 == t2) return;

    real3        k[BT::s];                    // Intermediate evaluations
    const realt  d = copysign(1.0f, t2 - t1); // Direction of integration
    realt        h = d * fabs(dt);            // Initial step size
    unsigned int i, j;                        // Butcher tableau indices
    realt        l;                           // Temporal lattice coordinate

    // Integrate from t=t1 to t=t2
    realt t = t1;
    while (d * t < d * t2) {
      // Ensure that last step ends at t2
      if (d * (t + h) > d * t2) h = t2 - t;
      // Evaluate velocities at intermediate steps
      if (BT::fsal && t != t1) k[0] = k[BT::s - 1], i = 1;
      else                                          i = 0;
      for (/*i = 0|1*/; i < BT::s; i++) {
        k[i] = p;
        for (j = 0; j < i; j++) k[i] += k[j] * BT::a(i,j);
        v.WorldToLattice(k[i]);
        l = v.TimeToLattice(t + BT::c(i) * h);
        v.FFD(k[i], l);
        k[i] *= h;
      }
      // Perform step
      for (i = 0; i < BT::s; i++) {
        p += k[i] * BT::b(i);
      }
      t += h;
    }
  }

  // -------------------------------------------------------------------------
  static IRTKCU_API void Jacobian(const FFD &, real3x3 &, real3 &, realt, realt, realt)
  {
    cerr << "irtkCUFreeFormTransformationExplicitRungeKutta::Jacobian: Not implemented" << endl;
    exit(1);
  }

  // -------------------------------------------------------------------------
  static IRTKCU_API void JacobianDOFs(const FFD &v, real3x3 &jac, uint4 cp, real3 &p, realt t1, realt t2, realt dt)
  {
    if (t1 == t2) return;

    real3        k [BT::s];                   // Intermediate evaluations
    real3x3      dk[BT::s];                   // Derivative of k_i w.r.t control point
    real3x3      dx;                          // Derivative of intermediate location w.r.t control point
    real3x3      Dv[BT::s];                   // Partial derivative of velocity field w.r.t spatial coordinates
    real3        dv[BT::s];                   // Partial derivative of velocity field w.r.t control point
    const realt  d = copysign(1.0f, t2 - t1); // Direction of integration
    realt        h = d * fabs(dt);            // Initial step size
    unsigned int i, j;                        // Butcher tableau indices
    realt        l;                           // Temporal lattice coordinate

    // Integrate from t=t1 to t=t2
    realt t = t1;
    while (d * t < d * t2) {
      // Ensure that last step ends at t2
      if (d * (t + h) > d * t2) h = t2 - t;
      // Evaluate at intermediate steps
      if (BT::fsal && t != t1) {
        i     = BT::s - 1;
        k [0] = k [i];
        Dv[0] = Dv[i];
        dv[0] = dv[i];
        i     = 1;
      } else {
        i     = 0;
      }
      for (/*i = 0|1*/; i < BT::s; i++) {
        // World coordinates of current step
        k[i] = p;
        for (j = 0; j < i; j++) k[i] += k[j] * BT::a(i, j);
        // Convert to lattice coordinates
        v.WorldToLattice(k[i]);
        l = v.TimeToLattice(t + BT::c(i) * h);
        // Evaluate partial derivatives of velocity
        v.FFDJacobianWorld(Dv[i],     k[i], l);
        v.FFDJacobianDOFs (dv[i], cp, k[i], l);
        // Evaluate velocity
        v.FFD(k[i], l);
        k[i] *= h;
      }
      // Calculate derivatives of k_i
      for (i = 0; i < BT::s; i++) {
        dx = jac;
        for (j = 0; j < i; j++) dx += dk[j] * BT::a(i, j);
        dkdp(dk[i], Dv[i], dx, dv[i], h); // dk = (Dv * dx + dv) * h
      }
      // Perform step with local extrapolation
      for (i = 0; i < BT::s; i++) {
        p   +=  k[i] * BT::b(i);
        jac += dk[i] * BT::b(i);
      }
      t += h;
    }
  }

  // -------------------------------------------------------------------------
  static IRTKCU_API void JacobianDOFs(const FFD &v, irtkCUJacobianDOFs &jac, real3 &p, realt t1, realt t2, realt dt)
  {
    if (t1 == t2) return;

    real3              k [BT::s];                   // Intermediate evaluations
    irtkCUJacobianDOFs dk[BT::s];                   // Derivative of k_i w.r.t control point and temporarily used
                                                    // for partial derivative of velocity field w.r.t control point
    irtkCUJacobianDOFs dx;                          // Derivative of intermediate location w.r.t control point
    real3x3            Dv;                          // Partial derivative of velocity field w.r.t spatial coordinates
    const realt        d = copysign(1.0f, t2 - t1); // Direction of integration
    realt              h = d * fabs(dt);            // Initial step size
    unsigned int       i, j;                        // Butcher tableau indices
    realt              l;                           // Temporal lattice coordinate
 
    // Integrate from t=t1 to t=t2
    realt t = t1;
    while (d * t < d * t2) {
      // Ensure that last step ends at t2
      if (d * (t + h) > d * t2) h = t2 - t;
      // Evaluate at intermediate steps
      for (i = 0; i < BT::s; i++) {
        // World coordinates of current step
        k[i] = p;
        for (j = 0; j < i; j++) k[i] += k[j] * BT::a(i, j);
        // Convert to lattice coordinates
        v.WorldToLattice(k[i]);
        l = v.TimeToLattice(t + BT::c(i) * h);
        // Evaluate partial derivatives of velocity
        v.FFDJacobianWorld(Dv,           k[i], l);
        v.FFDJacobianDOFs (dk[i]/*=dv*/, k[i], l);
        // Evaluate velocity
        v.FFD(k[i], l);
        k[i] *= h;
        // Jacobian of current step
        dx = jac;
        for (j = 0; j < i; j++) dx.add(dk[j], BT::a(i, j));
        // Calculate derivatives of k_i
        dx *= Dv * h;     //  Dv * dx * h
        dx.add(dk[i], h); // (Dv * dx + dv) * h
        dk[i] = dx;
      }
      // Perform step with local extrapolation
      for (i = 0; i < BT::s; i++) {
        p += k[i] * BT::b(i);
        jac.add(dk[i], BT::b(i));
      }
      t += h;
    }
  }

};

// =============================================================================
// Generic implementation of embedded Runge-Kutta
// =============================================================================

/**
 * Embedded Runge-Kutta integration method for FFD parameterized by velocity field.
 * 
 * The first template argument is the type of the FFD which has to implement
 * the following methods:
 *
 * - FFDEvaluateWorld
 * - FFDJacobianWorld
 * - FFDJacobianDOFs
 *
 * The second template argument is a structure of the respective Butcher tableau
 * of the Runge-Kutta method which has been defined using the macro
 * IRTK_DEFINE_FFDRK_EMBEDDED.
 */
template <class TFreeFormTransformation, class TButcherTableau>
class irtkCUFreeFormTransformationEmbeddedRungeKutta
:
  public irtkCUFreeFormTransformationRungeKutta
{
public:

  typedef TFreeFormTransformation FFD; ///< Short-hand for FFD type
  typedef TButcherTableau         BT;  ///< Short-hand for Butcher tableau

  // -------------------------------------------------------------------------
  static IRTKCU_API void Transform(const FFD &v, real3 &p, realt t1, realt t2, realt mindt, realt maxdt, realt tol)
  {
    if (t1 == t2) return;

    real3        k[BT::s];                             // k_i = h * v(t + c_i * h, x + sum_{j=0}^{i-1} a_j * k_j)
    real3        tmp;                                  // Solution of order p-1
    real3        next;                                 // Solution of order p
    realt        error;                                // Local error estimate
    realt        h = t2 - t1;                          // Initial step size
    realt        hnext;                                // Next step size
    const realt  d = copysign(1.0f, h);                // Direction of integration
    const realt  e = 1.0f / static_cast<realt>(BT::p); // Exponent for step size scaling factor
    unsigned int i, j;                                 // Butcher tableau indices
    realt        l;                                    // Temporal lattice coordinate

    // Decrease initial step size if necessary
    if (fabs(h) > maxdt) h = copysign(maxdt, d);

    // Integrate from t=t1 to t=t2
    realt t = t1;
    while (d * t < d * t2) {
      // Ensure that last step ends at t2
      if (d * (t + h) > d * t2) h = t2 - t;
      // Evaluate velocities at intermediate steps
      if (BT::fsal && t != t1) k[0] = k[BT::s - 1], i = 1;
      else                                          i = 0;
      for (/*i = 0|1*/; i < BT::s; i++) {
        k[i] = p;
        for (j = 0; j < i; j++) k[i] += k[j] * BT::a(i, j);
        v.WorldToLattice(k[i]);
        l = v.TimeToLattice(t + BT::c(i) * h);
        v.FFD(k[i], l);
        k[i] *= h;
      }
      // Calculate solution of order p
      next = p;
      for (i = 0; i < BT::s; i++) next += k[i] * BT::b(i);
      // Adapt step size
      if (mindt < maxdt) {
        // Calculate solution of order p-1
        tmp = p;
        for (i = 0; i < BT::s; i++) tmp += k[i] * BT::b0(i);
        // Estimate local error
        error = max(fabs(next - tmp));
        // If local error exceeds tolerance...
        if (fabs(h) > mindt && error > tol) {
          // ...decrease step size
          h *= 0.8f * pow(tol / error, e);
          if (fabs(h) < mindt) h = copysign(mindt, d);
          // ...and redo current step
          continue;
        // Otherwise, increase step size
        } else {
          hnext = 0.8f * pow(tol / error, e) * h;
          if      (fabs(hnext) < mindt) hnext = copysign(mindt, d);
          else if (fabs(hnext) > maxdt) hnext = copysign(maxdt, d);
        }
      } else {
        hnext = h;
      }
      // Perform step with local extrapolation
      p  = next;
      t += h;
      // Update step size
      h = hnext;
    }
  }

  // -------------------------------------------------------------------------
  static IRTKCU_API void Jacobian(const FFD &, real3x3 &, real3 &, realt, realt, realt, realt, realt)
  {
    cerr << "irtkCUFreeFormTransformationEmbeddedRungeKutta::Jacobian: Not implemented" << endl;
    exit(1);
  }

  // -------------------------------------------------------------------------
  static IRTKCU_API void JacobianDOFs(const FFD &v, real3x3 &jac, uint4 cp, real3 &p, realt t1, realt t2, realt mindt, realt maxdt, realt tol)
  {
    if (t1 == t2) return;

    real3        k [BT::s];                            // Intermediate evaluations
    real3x3      dk[BT::s];                            // Derivative of k_i w.r.t control point
    real3x3      dx;                                   // Derivative of intermediate location w.r.t control point
    real3x3      Dv[BT::s];                            // Partial derivative of velocity field w.r.t spatial coordinates
    real3        dv[BT::s];                            // Partial derivative of velocity field w.r.t control point
    realt        h = t2 - t1;                          // Initial step size
    realt        hnext;                                // Next step size
    const realt  d = copysign(1.0f, h);                // Direction of integration
    real3        tmp;                                  // Solution of order p-1
    real3        next;                                 // Solution of order p
    realt        error;                                // Local error estimate
    const realt  e = 1.0f / static_cast<realt>(BT::p); // Exponent for step size scaling factor
    unsigned int i, j;                                 // Butcher tableau indices
    realt        l;                                    // Temporal lattice coordinate

    // Decrease initial step size if necessary
    if (fabs(h) > maxdt) h = copysign(maxdt, d);

    // Integrate from t=t1 to t=t2
    realt t = t1;
    while (d * t < d * t2) {
      // Ensure that last step ends at t2
      if (d * (t + h) > d * t2) h = t2 - t;
      // Evaluate at intermediate steps
      if (BT::fsal && t != t1) {
        i     = BT::s - 1;
        k [0] = k [i];
        Dv[0] = Dv[i];
        dv[0] = dv[i];
        i     = 1;
      } else {
        i     = 0;
      }
      for (/*i = 0|1*/; i < BT::s; i++) {
        // Lattice coordinates
        k[i] = p;
        for (j = 0; j < i; j++) k[i] += k[j] * BT::a(i, j);
        v.WorldToLattice(k[i]);
        l = v.TimeToLattice(t + BT::c(i) * h);
        // Evaluate partial derivatives of velocity
        v.FFDJacobianWorld(Dv[i],     k[i], l);
        v.FFDJacobianDOFs (dv[i], cp, k[i], l);
        // Evaluate velocity
        v.FFD(k[i], l);
        k[i] *= h;
      }
      // Calculate solution of order p
      next = p;
      for (i = 0; i < BT::s; i++) next += k[i] * BT::b(i);
      // Adapt step size
      if (mindt < maxdt) {
        // Calculate solution of order p-1
        tmp = p;
        for (i = 0; i < BT::s; i++) tmp += k[i] * BT::b0(i);
        // Estimate local error
        error = max(fabs(next - tmp));
        // If local error exceeds tolerance...
        if (fabs(h) > mindt && error > tol) {
          // ...decrease step size
          h *= 0.8f * pow(tol / error, e);
          if (fabs(h) < mindt) h = copysign(mindt, d);
          // ...and redo current step
          continue;
        // Otherwise, increase step size
        } else {
          hnext = 0.8f * pow(tol / error, e) * h;
          if      (fabs(hnext) < mindt) hnext = copysign(mindt, d);
          else if (fabs(hnext) > maxdt) hnext = copysign(maxdt, d);
        }
      } else {
        hnext = h;
      }
      // Calculate derivatives of k_i
      for (i = 0; i < BT::s; i++) {
        dx = jac;
        for (j = 0; j < i; j++) dx += dk[j] * BT::a(i, j);
        dkdp(dk[i], Dv[i], dx, dv[i], h); // dk = (Dv * dx + dv) * h
      }
      // Update Jacobian
      for (i = 0; i < BT::s; i++) jac += dk[i] * BT::b(i);
      // Perform step with local extrapolation
      p  = next;
      t += h;
      // Update step size
      h = hnext;
    }
  }

  // -------------------------------------------------------------------------
  static IRTKCU_API void JacobianDOFs(const FFD &v, irtkCUJacobianDOFs &jac, real3 &p, realt t1, realt t2, realt mindt, realt maxdt, realt tol)
  {
    if (t1 == t2) return;

    real3              k [BT::s];                            // Intermediate evaluations
    irtkCUJacobianDOFs dk[BT::s];                            // Derivative of k_i w.r.t control point
    irtkCUJacobianDOFs dx;                                   // Derivative of intermediate location w.r.t control point
    real3x3            Dv[BT::s];                            // Partial derivative of velocity field w.r.t spatial coordinates
    irtkCUJacobianDOFs dv[BT::s];                            // Partial derivative of velocity field w.r.t control point
    realt              h = t2 - t1;                          // Initial step size
    realt              hnext;                                // Next step size
    const realt        d = copysign(1.0f, h);                // Direction of integration
    real3              tmp;                                  // Solution of order p-1
    real3              next;                                 // Solution of order p
    realt              error;                                // Local error estimate
    const realt        e = 1.0f / static_cast<realt>(BT::p); // Exponent for step size scaling factor
    unsigned int       i, j;                                 // Butcher tableau indices
    realt              l;                                    // Temporal lattice coordinate

    // Decrease initial step size if necessary
    if (fabs(h) > maxdt) h = copysign(maxdt, d);

    // Integrate from t=t1 to t=t2
    realt t = t1;
    while (d * t < d * t2) {
      // Ensure that last step ends at t2
      if (d * (t + h) > d * t2) h = t2 - t;
      // Evaluate at intermediate steps
      if (BT::fsal && t != t1) {
        i     = BT::s - 1;
        k [0] = k [i];
        Dv[0] = Dv[i];
        dv[0] = dv[i];
        i     = 1;
      } else {
        i     = 0;
      }
      for (/*i = 0|1*/; i < BT::s; i++) {
        // Lattice coordinates
        k[i] = p;
        for (j = 0; j < i; j++) k[i] += k[j] * BT::a(i, j);
        v.WorldToLattice(k[i]);
        l = v.TimeToLattice(t + BT::c(i) * h);
        // Evaluate partial derivatives of velocity
        v.FFDJacobianWorld(Dv[i], k[i], l);
        v.FFDJacobianDOFs (dv[i], k[i], l);
        // Evaluate velocity
        v.FFD(k[i], l);
        k[i] *= h;
      }
      // Calculate solution of order p
      next = p;
      for (i = 0; i < BT::s; i++) next += k[i] * BT::b(i);
      // Adapt step size
      if (mindt < maxdt) {
        // Calculate solution of order p-1
        tmp = p;
        for (i = 0; i < BT::s; i++) tmp += k[i] * BT::b0(i);
        // Estimate local error
        error = max(fabs(next - tmp));
        // If local error exceeds tolerance...
        if (fabs(h) > mindt && error > tol) {
          // ...decrease step size
          h *= 0.8f * pow(tol / error, e);
          if (fabs(h) < mindt) h = copysign(mindt, d);
          // ...and redo current step
          continue;
        // Otherwise, increase step size
        } else {
          hnext = 0.8f * pow(tol / error, e) * h;
          if      (fabs(hnext) < mindt) hnext = copysign(mindt, d);
          else if (fabs(hnext) > maxdt) hnext = copysign(maxdt, d);
        }
      } else {
        hnext = h;
      }
      // Calculate derivatives of k_i
      for (i = 0; i < BT::s; i++) {
        dx = jac;
        for (j = 0; j < i; j++) dx += dk[j] * BT::a(i, j);
        dkdp(dk[i], Dv[i], dx, dv[i], h); // dk = (Dv * dx + dv) * h
      }
      // Update Jacobian
      for (i = 0; i < BT::s; i++) jac += dk[i] * BT::b(i);
      // Perform step with local extrapolation
      p  = next;
      t += h;
      // Update step size
      h = hnext;
    }
  }
};

// =============================================================================
// Auxiliary macros
// =============================================================================

// -----------------------------------------------------------------------------
/// Define Butcher tableau of explicit Runge-Kutta method
#define IRTKCU_DEFINE_FFDRK_EXPLICIT(NAME, S, P, FSAL, C, A, B)                \
  struct irtkCUFreeFormTransformationButcherTableau##NAME                      \
  {                                                                            \
    static const int   s    = S;                                               \
    static const int   p    = P;                                               \
    static const bool  fsal = FSAL;                                            \
                                                                               \
    IRTKCU_API static realt a(int i, int j)                                    \
    {                                                                          \
      const realt a[S][S] = A;                                                 \
      return a[i][j];                                                          \
    }                                                                          \
                                                                               \
    /* Coefficients for higher order approximation of integral */              \
    IRTKCU_API static realt b(int i)                                           \
    {                                                                          \
      const realt b[S] = B;                                                    \
      return b[i];                                                             \
    }                                                                          \
                                                                               \
    /* Coefficients for lower order approximation of integral */               \
    IRTKCU_API static realt b0(int i)                                          \
    {                                                                          \
      return b(i);                                                             \
    }                                                                          \
                                                                               \
    IRTKCU_API static realt c(int i)                                           \
    {                                                                          \
      const realt c[S] = C;                                                    \
      return c[i];                                                             \
    }                                                                          \
  }

// -----------------------------------------------------------------------------
/// Define Butcher tableau of embedded Runge-Kutta method
#define IRTKCU_DEFINE_FFDRK_EMBEDDED(NAME, S, P, FSAL, C, A, B0, B)            \
  struct irtkCUFreeFormTransformationButcherTableau##NAME                      \
  {                                                                            \
    static const int   s    = S;                                               \
    static const int   p    = P;                                               \
    static const bool  fsal = FSAL;                                            \
                                                                               \
    IRTKCU_API static realt a(int i, int j)                                    \
    {                                                                          \
      const realt a[S][S] = A;                                                 \
      return a[i][j];                                                          \
    }                                                                          \
                                                                               \
    /* Coefficients for higher order approximation of integral */              \
    IRTKCU_API static realt b(int i)                                           \
    {                                                                          \
      const realt b[S] = B;                                                    \
      return b[i];                                                             \
    }                                                                          \
                                                                               \
    /* Coefficients for lower order approximation of integral */               \
    IRTKCU_API static realt b0(int i)                                          \
    {                                                                          \
      const realt b[S] = B0;                                                   \
      return b[i];                                                             \
    }                                                                          \
                                                                               \
    IRTKCU_API static realt c(int i)                                           \
    {                                                                          \
      const realt c[S] = C;                                                    \
      return c[i];                                                             \
    }                                                                          \
  }

// ---------------------------------------------------------------------------
/// Declare explicit Runge-Kutta method
#define IRTKCU_DECLARE_FFDRK_EXPLICIT(NAME)                                    \
  template <class TFFD>                                                        \
  class irtkCUFreeFormTransformationIntegration##NAME                          \
    : public irtkCUFreeFormTransformationExplicitRungeKutta                    \
             <TFFD, irtkCUFreeFormTransformationButcherTableau##NAME>          \
  { }

// ---------------------------------------------------------------------------
/// Declare embedded Runge-Kutta method
#define IRTKCU_DECLARE_FFDRK_EMBEDDED(NAME)                                    \
  template <class TFFD>                                                        \
  class irtkCUFreeFormTransformationIntegration##NAME                          \
    : public irtkCUFreeFormTransformationEmbeddedRungeKutta                    \
             <TFFD, irtkCUFreeFormTransformationButcherTableau##NAME>          \
  { }

// ============================================================================
// Runge-Kutta integration methods
// ============================================================================

#define BT_R(...) {__VA_ARGS__}
#define BT_A(...) {__VA_ARGS__}
#define BT_B(b)    b
#define BT_C(c)    c


// -----------------------------------------------------------------------------
// Forward Euler method
IRTKCU_DEFINE_FFDRK_EXPLICIT(RKE1, 1, 1, false,
  BT_C(BT_R(0.0f)),
  /* ------------- */
  BT_A(BT_R(0.0f)),
  /* ------------- */
  BT_B(BT_R(1.0f))
);

// -----------------------------------------------------------------------------
// Modified Euler method (Heun's method, explicit midpoint rule)
IRTKCU_DEFINE_FFDRK_EXPLICIT(RKE2, 2, 2, false,
  BT_C(BT_R(0.0f, 0.5f)),
  /* ----------------- */
  BT_A(BT_R(0.0f, 0.0f),
       BT_R(0.5f, 0.0f)),
  /* ----------------- */
  BT_B(BT_R(0.0f, 1.0f))
);

// -----------------------------------------------------------------------------
// Improved Euler method (Heun's method, explicit trapezoidal rule)
IRTKCU_DEFINE_FFDRK_EXPLICIT(RKH2, 2, 2, false,
  BT_C(BT_R(0.0f, 0.5f)),
  /* ----------------- */
  BT_A(BT_R(0.0f, 0.0f),
       BT_R(1.0f, 0.0f)),
  /* ----------------- */
  BT_B(BT_R(0.5f, 0.5f))
);

// -----------------------------------------------------------------------------
// Classical Runge-Kutta method
IRTKCU_DEFINE_FFDRK_EXPLICIT(RK4, 4, 4, false,
  BT_C(BT_R(0.0f,      1.0f/2.0f, 1.0f/2.0f, 1.0f)),
  /* ------------------------------------------------ */
  BT_A(BT_R(0.0f,      0.0f,      0.0f,      0.0f),
       BT_R(1.0f/2.0f, 0.0f,      0.0f,      0.0f),
       BT_R(0.0f,      1.0f/2.0f, 0.0f,      0.0f),
       BT_R(0.0f,      0.0f,      1.0f,      0.0f)),
  /* ------------------------------------------------ */
  BT_B(BT_R(1.0f/6.0f, 1.0f/3.0f, 1.0f/3.0f, 1.0f/6.0f))
);

// ---------------------------------------------------------------------------
// Euler-Heun method of order 2(1)
IRTKCU_DEFINE_FFDRK_EMBEDDED(RKEH12, 2, 2, false,
  BT_C(BT_R(0.0, 1.0)),
  /* ----------------- */
  BT_A(BT_R(0.0, 0.0),
       BT_R(1.0, 0.0)),
  /* ----------------- */
  BT_B(BT_R(0.5, 0.5)),
  BT_B(BT_R(1.0, 0.0))
);

// ---------------------------------------------------------------------------
// Bogacki-Shampine method of order 3(2)
IRTKCU_DEFINE_FFDRK_EMBEDDED(RKBS23, 4, 3, true,
  BT_C(BT_R(0.0,      1.0/2.0,  3.0/4.0, 1.0    )),
  /* ---------------------------------------------- */
  BT_A(BT_R(0.0,      0.0,      0.0,     0.0    ),
       BT_R(1.0/2.0,  0.0,      0.0,     0.0    ),
       BT_R(0.0,      3.0/4.0,  0.0,     0.0    ),
       BT_R(2.0/9.0,  1.0/3.0,  4.0/9.0, 0.0    )),
  /* ---------------------------------------------- */
  BT_B(BT_R(2.0/9.0,  2.0/3.0,  4.0/9.0, 0.0    )),
  BT_B(BT_R(7.0/24.0, 1.0/4.0,  1.0/3.0, 1.0/8.0))
);

// ---------------------------------------------------------------------------
// Fehlberg method of order 5(4)
IRTKCU_DEFINE_FFDRK_EMBEDDED(RKF45, 6, 5, false,
  BT_C(BT_R(   0.0,            1.0/4.0,        3.0/8.0,        12.0/13.0,     1.0,       1.0/2.0 )),
  /* ----------------------------------------------------------------------------------------------- */
  BT_A(BT_R(   0.0,            0.0,            0.0,             0.0,          0.0,       0.0     ),
       BT_R(   1.0/4.0,        0.0,            0.0,             0.0,          0.0,       0.0     ),
       BT_R(   3.0/32.0,       9.0/32.0,       0.0,             0.0,          0.0,       0.0     ),
       BT_R(1932.0/2197.0, -7200.0/2197.0,  7296.0/2197.0,      0.0,          0.0,       0.0     ),
       BT_R( 439.0/216.0,     -8.0,         3680.0/513.0,    -845.0/4104.0,   0.0,       0.0     ),
       BT_R(  -8.0/27.0,       2.0,        -3544.0/2565.0,   1859.0/4104.0, -11.0/40.0,  0.0     )),
  /* ----------------------------------------------------------------------------------------------- */
  BT_B(BT_R(  25.0/216.0,      0.0,         1408.0/2565.0,   2197.0/4104.0,  -1.0/5.0,   0.0     )),
  BT_B(BT_R(  16.0/135.0,      0.0,         6656.0/12825.0, 28561.0/56430.0, -9.0/50.0,  2.0/55.0))
);

// ---------------------------------------------------------------------------
// Dormand–Prince method of order 5(4)
IRTKCU_DEFINE_FFDRK_EMBEDDED(RKDP45, 7, 5, true,
  BT_C(BT_R(     0.0,             1.0/5.0,        3.0/10.0,        4.0/5.0,        8.0/9.0,        1.0,        1.0     )),
  /* --------------------------------------------------------------------------------------------------------------------- */
  BT_A(BT_R(     0.0,             0.0,            0.0,             0.0,            0.0,            0.0,        0.0     ),
       BT_R(     1.0/5.0,         0.0,            0.0,             0.0,            0.0,            0.0,        0.0     ),
       BT_R(     3.0/40.0,        9.0/40.0,       0.0,             0.0,            0.0,            0.0,        0.0     ),
       BT_R(    44.0/45.0,      -56.0/15.0,      32.0/9.0,         0.0,            0.0,            0.0,        0.0     ),
       BT_R( 19372.0/6561.0, -25360.0/2187.0, 64448.0/6561.0,   -212.0/729.0,      0.0,            0.0,        0.0     ),
       BT_R(  9017.0/3168.0,   -355.0/33.0,   46732.0/5247.0,     49.0/176.0,  -5103.0/18656.0,    0.0,        0.0     ),
       BT_R(    35.0/384.0,       0.0,          500.0/1113.0,    125.0/192.0,  -2187.0/6784.0,    11.0/84.0,   0.0     )),
  /* --------------------------------------------------------------------------------------------------------------------- */
  BT_B(BT_R(    35.0/384.0,       0.0,          500.0/1113.0,    125.0/192.0,  -2187.0/6784.0,    11.0/84.0,   0.0     )),
  BT_B(BT_R(  5179.0/57600.0,     0.0,         7571.0/16695.0,   393.0/640.0, -92097.0/339200.0, 187.0/2100.0, 1.0/40.0))
);

// ---------------------------------------------------------------------------
// Cash-Karp method of order 5(4)
IRTKCU_DEFINE_FFDRK_EMBEDDED(RKCK45, 6, 5, false,
  BT_C(BT_R(   0.0,           1.0/5.0,     3.0/10.0,        3.0/5.0,        1.0,           7.0/8.0   )),
  /* --------------------------------------------------------------------------------------------------- */
  BT_A(BT_R(   0.0,           0.0,         0.0,             0.0,            0.0,           0.0       ),
       BT_R(   1.0/5.0,       0.0,         0.0,             0.0,            0.0,           0.0       ),
       BT_R(   3.0/40.0,      9.0/40.0,    0.0,             0.0,            0.0,           0.0       ),
       BT_R(   3.0/10.0,     -9.0/10.0,    6.0/5.0,         0.0,            0.0,           0.0       ),
       BT_R( -11.0/54.0,      5.0/2.0,   -70.0/27.0,       35.0/27.0,       0.0,           0.0       ),
       BT_R(1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0,    0.0       )),
  /* --------------------------------------------------------------------------------------------------- */
  BT_B(BT_R(  37.0/378.0,     0.0,       250.0/621.0,     125.0/594.0,      0.0,         512.0/1771.0)),
  BT_B(BT_R(2825.0/27648.0,   0.0,     18575.0/48384.0, 13525.0/55296.0,  277.0/14336.0,   1.0/4.0   ))
);


#undef BT_A
#undef BT_B
#undef BT_C
#undef BT_R

// ============================================================================
// Runge-Kutta integration methods
// ============================================================================

IRTKCU_DECLARE_FFDRK_EXPLICIT(RKE1);   ///< Forward Euler method
IRTKCU_DECLARE_FFDRK_EMBEDDED(RKEH12); ///< Euler-Heun method of order 2(1)
IRTKCU_DECLARE_FFDRK_EXPLICIT(RKE2);   ///< Modified Euler method (Heun's method, explicit midpoint rule)
IRTKCU_DECLARE_FFDRK_EXPLICIT(RKH2);   ///< Improved Euler method (Heun's method, explicit trapezoidal rule)
IRTKCU_DECLARE_FFDRK_EMBEDDED(RKBS23); ///< Bogacki-Shampine method of order 3(2)
IRTKCU_DECLARE_FFDRK_EXPLICIT(RK4);    ///< Classical Runge-Kutta method
IRTKCU_DECLARE_FFDRK_EMBEDDED(RKF45);  ///< Fehlberg method of order 5(4)
IRTKCU_DECLARE_FFDRK_EMBEDDED(RKCK45); ///< Cash-Karp method of order 5(4)
IRTKCU_DECLARE_FFDRK_EMBEDDED(RKDP45); ///< Dormand–Prince method of order 5(4)

// ============================================================================
// Auxiliary function
// ============================================================================

// -----------------------------------------------------------------------------
template <unsigned int FFDIM, class FFD>
IRTKCU_API real3 RungeKutta(const FFD &ffd, real3 p, realt t1, realt t2, realt mindt, realt maxdt = .0f, realt tol = .0f)
{
  // Note: if statements are evaluated at compile time, not runtime!
  if      (FFDIM == FFDIM_RKE1)   irtkCUFreeFormTransformationIntegrationRKE1  <FFD>::Transform(ffd, p, t1, t2, mindt);
  else if (FFDIM == FFDIM_RKE2)   irtkCUFreeFormTransformationIntegrationRKE2  <FFD>::Transform(ffd, p, t1, t2, mindt);
  else if (FFDIM == FFDIM_RKH2)   irtkCUFreeFormTransformationIntegrationRKH2  <FFD>::Transform(ffd, p, t1, t2, mindt);
  else if (FFDIM == FFDIM_RK4)    irtkCUFreeFormTransformationIntegrationRK4   <FFD>::Transform(ffd, p, t1, t2, mindt);
  else if (FFDIM == FFDIM_RKEH12) irtkCUFreeFormTransformationIntegrationRKEH12<FFD>::Transform(ffd, p, t1, t2, mindt, maxdt, tol);
  else if (FFDIM == FFDIM_RKBS23) irtkCUFreeFormTransformationIntegrationRKBS23<FFD>::Transform(ffd, p, t1, t2, mindt, maxdt, tol);
  else if (FFDIM == FFDIM_RKF45)  irtkCUFreeFormTransformationIntegrationRKF45 <FFD>::Transform(ffd, p, t1, t2, mindt, maxdt, tol);
  else if (FFDIM == FFDIM_RKCK45) irtkCUFreeFormTransformationIntegrationRKCK45<FFD>::Transform(ffd, p, t1, t2, mindt, maxdt, tol);
  else if (FFDIM == FFDIM_RKDP45) irtkCUFreeFormTransformationIntegrationRKDP45<FFD>::Transform(ffd, p, t1, t2, mindt, maxdt, tol);
  return p;
}

// -----------------------------------------------------------------------------
template <unsigned int FFDIM, class FFD>
IRTKCU_API real3 RungeKutta(const FFD &ffd, real3x3 &J, uint4 cp, real3 p, realt t1, realt t2, realt mindt, realt maxdt = .0f, realt tol = .0f)
{
  // Note: if statements are evaluated at compile time, not runtime!
  if      (FFDIM == FFDIM_RKE1)   irtkCUFreeFormTransformationIntegrationRKE1  <FFD>::JacobianDOFs(ffd, J, cp, p, t1, t2, mindt);
  else if (FFDIM == FFDIM_RKE2)   irtkCUFreeFormTransformationIntegrationRKE2  <FFD>::JacobianDOFs(ffd, J, cp, p, t1, t2, mindt);
  else if (FFDIM == FFDIM_RKH2)   irtkCUFreeFormTransformationIntegrationRKH2  <FFD>::JacobianDOFs(ffd, J, cp, p, t1, t2, mindt);
  else if (FFDIM == FFDIM_RK4)    irtkCUFreeFormTransformationIntegrationRK4   <FFD>::JacobianDOFs(ffd, J, cp, p, t1, t2, mindt);
  else if (FFDIM == FFDIM_RKEH12) irtkCUFreeFormTransformationIntegrationRKEH12<FFD>::JacobianDOFs(ffd, J, cp, p, t1, t2, mindt, maxdt, tol);
  else if (FFDIM == FFDIM_RKBS23) irtkCUFreeFormTransformationIntegrationRKBS23<FFD>::JacobianDOFs(ffd, J, cp, p, t1, t2, mindt, maxdt, tol);
  else if (FFDIM == FFDIM_RKF45)  irtkCUFreeFormTransformationIntegrationRKF45 <FFD>::JacobianDOFs(ffd, J, cp, p, t1, t2, mindt, maxdt, tol);
  else if (FFDIM == FFDIM_RKCK45) irtkCUFreeFormTransformationIntegrationRKCK45<FFD>::JacobianDOFs(ffd, J, cp, p, t1, t2, mindt, maxdt, tol);
  else if (FFDIM == FFDIM_RKDP45) irtkCUFreeFormTransformationIntegrationRKDP45<FFD>::JacobianDOFs(ffd, J, cp, p, t1, t2, mindt, maxdt, tol);
  return p;
}

// -----------------------------------------------------------------------------
template <unsigned int FFDIM, class FFD>
IRTKCU_API real3 RungeKutta(const FFD &ffd, irtkCUJacobianDOFs &J, real3 p, realt t1, realt t2, realt mindt, realt maxdt = .0f, realt tol = .0f)
{
  // Note: if statements are evaluated at compile time, not runtime!
  if      (FFDIM == FFDIM_RKE1)   irtkCUFreeFormTransformationIntegrationRKE1  <FFD>::JacobianDOFs(ffd, J, p, t1, t2, mindt);
  else if (FFDIM == FFDIM_RKE2)   irtkCUFreeFormTransformationIntegrationRKE2  <FFD>::JacobianDOFs(ffd, J, p, t1, t2, mindt);
  else if (FFDIM == FFDIM_RKH2)   irtkCUFreeFormTransformationIntegrationRKH2  <FFD>::JacobianDOFs(ffd, J, p, t1, t2, mindt);
  else if (FFDIM == FFDIM_RK4)    irtkCUFreeFormTransformationIntegrationRK4   <FFD>::JacobianDOFs(ffd, J, p, t1, t2, mindt);
  else if (FFDIM == FFDIM_RKEH12) irtkCUFreeFormTransformationIntegrationRKEH12<FFD>::JacobianDOFs(ffd, J, p, t1, t2, mindt, maxdt, tol);
  else if (FFDIM == FFDIM_RKBS23) irtkCUFreeFormTransformationIntegrationRKBS23<FFD>::JacobianDOFs(ffd, J, p, t1, t2, mindt, maxdt, tol);
  else if (FFDIM == FFDIM_RKF45)  irtkCUFreeFormTransformationIntegrationRKF45 <FFD>::JacobianDOFs(ffd, J, p, t1, t2, mindt, maxdt, tol);
  else if (FFDIM == FFDIM_RKCK45) irtkCUFreeFormTransformationIntegrationRKCK45<FFD>::JacobianDOFs(ffd, J, p, t1, t2, mindt, maxdt, tol);
  else if (FFDIM == FFDIM_RKDP45) irtkCUFreeFormTransformationIntegrationRKDP45<FFD>::JacobianDOFs(ffd, J, p, t1, t2, mindt, maxdt, tol);
  return p;
}


#endif
