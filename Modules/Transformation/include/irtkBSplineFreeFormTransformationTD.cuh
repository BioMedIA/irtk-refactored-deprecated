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

#ifndef _IRTKBSPLINEFREEFORMTRANSFORMATION_CUH
#define _IRTKBSPLINEFREEFORMTRANSFORMATION_CUH

#include <irtkTransformation.h>

#include <irtkCUCommon.h>


/**
 * B-spline TD FFD for use by CUDA accelerated device code.
 */
class irtkCUBSplineFreeFormTransformationTD
{
  // ---------------------------------------------------------------------------
  // Data members
protected:

  real4   _origin;  ///< Origin of control point lattice.
  real4   _spacing; ///< Spacing between control points in each dimension.
  uint4   _size;    ///< Number of control points in each dimension.
  real3x4 _matL2W;  ///< Transform spatial lattice to world coordinates.
  real3x4 _matW2L;  ///< Transform spatial world to lattice coordinates.
  real3  *_cp;      ///< Control point data in global device memory.

public:

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Construct FFD usable by device code from host transformation
  IRTKCU_HOST_API irtkCUBSplineFreeFormTransformationTD(const irtkFreeFormTransformation4D &);

  /// Copy constructor
  IRTKCU_HOST_API irtkCUBSplineFreeFormTransformationTD(const irtkCUBSplineFreeFormTransformationTD &);

  /// Destruct structure and free device memory
  IRTKCU_HOST_API ~irtkCUBSplineFreeFormTransformationTD();

  /// Initialize device FFD given host FFD
  IRTKCU_HOST_API void Initialize(const irtkFreeFormTransformation4D &);

  // ---------------------------------------------------------------------------
  // Coordinate transformation

  /// Convert spatial lattice to world coordinates
  IRTKCU_API void LatticeToWorld(realt &, realt &, realt &) const;

  /// Convert spatial lattice to world coordinates
  IRTKCU_API real3 LatticeToWorld(const real3 &) const;

  /// Convert spatial lattice to world coordinates
  IRTKCU_API real3 LatticeToWorld(const int3 &) const;

  /// Convert spatial world to lattice coordinates
  IRTKCU_API void WorldToLattice(realt &, realt &, realt &) const;

  /// Convert spatial world to lattice coordinates
  IRTKCU_API real3 WorldToLattice(const real3 &) const;


  /// Convert temporal lattice to time coordinate
  IRTKCU_API realt LatticeToTime(realt) const;

  /// Convert temporal lattice to time coordinate
  IRTKCU_API realt LatticeToTime(int) const;

  /// Convert temporal time to lattice coordinate
  IRTKCU_API realt TimeToLattice(realt) const;


  /// Convert lattice to world coordinates
  IRTKCU_API void LatticeToWorld(realt &, realt &, realt &, realt &) const;

  /// Convert lattice to world coordinates
  IRTKCU_API real4 LatticeToWorld(const real4 &) const;

  /// Convert lattice to world coordinates
  IRTKCU_API real4 LatticeToWorld(const int4 &) const;

  /// Convert world to lattice coordinates
  IRTKCU_API void WorldToLattice(realt &, realt &, realt &, realt &) const;

  /// Convert world to lattice coordinates
  IRTKCU_API real4 WorldToLattice(const real4 &) const;


  /// Convert lattice coordinates to index
  IRTKCU_API int LatticeToIndex(int, int, int, int = 0) const;

  /// Convert lattice coordinates to index
  IRTKCU_API int LatticeToIndex(const int3 &, int = 0) const;

  /// Convert lattice coordinates to index
  IRTKCU_API int LatticeToIndex(const int4 &) const;

  /// Convert index to lattice coordinates
  IRTKCU_API int4 IndexToLattice(int) const;

  /// Convert index to lattice coordinates
  IRTKCU_API void IndexToLattice(int, int &, int &, int &, int &) const;

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Evaluate FFD at point in lattice coordinates
  IRTKCU_API real3 EvaluateFFD(realt, realt, realt, realt) const;

  /// Evaluate FFD at point in lattice coordinates
  IRTKCU_API real3 EvaluateFFD(real3 p, realt) const;

  /// Evaluate FFD at point in world coordinates
  IRTKCU_API real3 Evaluate(real3 p, realt) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkCUBSplineFreeFormTransformationTD::irtkCUBSplineFreeFormTransformationTD(const irtkFreeFormTransformation4D &ffd)
:
  _cp(NULL)
{
  Initialize(ffd);
}

// -----------------------------------------------------------------------------
irtkCUBSplineFreeFormTransformationTD::irtkCUBSplineFreeFormTransformationTD(const irtkCUBSplineFreeFormTransformationTD &other)
:
  _origin (other._origin),
  _spacing(other._spacing),
  _size   (other._size),
  _matL2W (other._matL2W),
  _matW2L (other._matW2L),
  _cp     (NULL)
{
  const size_t nbytes = (_size.x+8) * (_size.y+8) * (_size.z+8) * (_size.w+8) * sizeof(real3);
  CudaSafeCall( cudaMalloc(&_cp, nbytes) );
  CudaSafeCall( cudaMemcpy(_cp, other._cp, nbytes, cudaMemcpyDeviceToDevice) );
}

// -----------------------------------------------------------------------------
irtkCUBSplineFreeFormTransformationTD::~irtkCUBSplineFreeFormTransformationTD()
{
  cudaFree(_cp);
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
void irtkCUBSplineFreeFormTransformationTD::Initialize(const irtkFreeFormTransformation4D &ffd)
{
  // Free previous data
  cudaFree(_cp);

  // Initialize FFD attributes
  _origin  = make_real4  (ffd._origin._x, ffd._origin._y, ffd._origin._z, ffd._torigin);
  _spacing = make_real4  (ffd._dx, ffd._dy, ffd._dz, ffd._dt);
  _size    = make_uint4  (ffd._x,  ffd._y,  ffd._z,  ffd._t);
  _matL2W  = make_real3x4(ffd._matL2W);
  _matW2L  = make_real3x4(ffd._matW2L);

  // Allocate device memory for control point data
  const size_t npoints = (_size.x+8) * (_size.y+8) * (_size.z+8) * (_size.w+8);
  const size_t nbytes  = npoints * sizeof(real3);
  CudaSafeCall( cudaMalloc(&_cp, nbytes) );

  // Copy control point data from host to device
  real3 *cp = new real3[npoints];
  if (cp == NULL) {
    cerr << "irtkCUBSplineFreeFormTransformationTD::Initialize: Failed to allocate temporary host memory for data conversion!" << endl;
    exit(1);
  }
  real3 *ptr2cp = cp;
  for (int l = -4; l < _size.w+4; l++) {
    for (int k = -4; k < _size.z+4; k++) {
      for (int j = -4; j < _size.y+4; j++) {
        for (int i = -4; i < _size.x+4; i++) {
          ptr2cp->x = static_cast<realt>(ffd._xdata[l][k][j][i]);
          ptr2cp->y = static_cast<realt>(ffd._ydata[l][k][j][i]);
          ptr2cp->z = static_cast<realt>(ffd._zdata[l][k][j][i]);
          ptr2cp++;
        }
      }
    }
  }
  CudaSafeCall( cudaMemcpy(_cp, cp, nbytes, cudaMemcpyHostToDevice) );
  delete[] cp;
}

// =============================================================================
// Coordinate transformation
// =============================================================================

// -----------------------------------------------------------------------------
real3 irtkCUBSplineFreeFormTransformationTD::LatticeToWorld(const real3 &p) const
{
  return _matL2W * p;
}

// -----------------------------------------------------------------------------
void irtkCUBSplineFreeFormTransformationTD::LatticeToWorld(realt &x, realt &y, realt &z) const
{
  real3 p = make_real3(x, y, z);
  p = LatticeToWorld(p);
  x = p.x, y = p.y, z = p.z;
}

// -----------------------------------------------------------------------------
real3 irtkCUBSplineFreeFormTransformationTD::LatticeToWorld(const int3 &p) const
{
  return LatticeToWorld(make_real3(p));
}

// -----------------------------------------------------------------------------
real3 irtkCUBSplineFreeFormTransformationTD::WorldToLattice(const real3 &p) const
{
  return _matW2L * p;
}

// -----------------------------------------------------------------------------
void irtkCUBSplineFreeFormTransformationTD::WorldToLattice(realt &x, realt &y, realt &z) const
{
  real3 p = make_real3(x, y, z);
  p = WorldToLattice(p);
  x = p.x, y = p.y, z = p.z;
}


// -----------------------------------------------------------------------------
realt irtkCUBSplineFreeFormTransformationTD::LatticeToTime(realt l) const
{
  return _origin.w + l * _spacing.w;
}

// -----------------------------------------------------------------------------
realt irtkCUBSplineFreeFormTransformationTD::LatticeToTime(int l) const
{
  return LatticeToTime(static_cast<realt>(l));
}

// -----------------------------------------------------------------------------
realt irtkCUBSplineFreeFormTransformationTD::TimeToLattice(realt t) const
{
  return (t - _origin.w) / _spacing.w;
}


// -----------------------------------------------------------------------------
void irtkCUBSplineFreeFormTransformationTD::LatticeToWorld(realt &x, realt &y, realt &z, realt &t) const
{
  LatticeToWorld(x, y, z);
  t = LatticeToTime(t);
}

// -----------------------------------------------------------------------------
real4 irtkCUBSplineFreeFormTransformationTD::LatticeToWorld(const real4 &p) const
{
  real4 w = p;
  LatticeToWorld(w.x, w.y, w.z, w.w);
  return w;
}

// -----------------------------------------------------------------------------
real4 irtkCUBSplineFreeFormTransformationTD::LatticeToWorld(const int4 &p) const
{
  return LatticeToWorld(make_real4(p));
}

// -----------------------------------------------------------------------------
void irtkCUBSplineFreeFormTransformationTD::WorldToLattice(realt &x, realt &y, realt &z, realt &t) const
{
  WorldToLattice(x, y, z);
  t = TimeToLattice(t);
}

// -----------------------------------------------------------------------------
real4 irtkCUBSplineFreeFormTransformationTD::WorldToLattice(const real4 &w) const
{
  real4 p = w;
  WorldToLattice(p.x, p.y, p.z, p.w);
  return p;
}


// -----------------------------------------------------------------------------
int irtkCUBSplineFreeFormTransformationTD::LatticeToIndex(int i, int j, int k, int l) const
{
  return ((l * _size.z + k) * _size.y + j) * _size.x + i;
}

// -----------------------------------------------------------------------------
int irtkCUBSplineFreeFormTransformationTD::LatticeToIndex(const int4 &p) const
{
  return LatticeToIndex(p.x, p.y, p.z, p.w);
}

// -----------------------------------------------------------------------------
void irtkCUBSplineFreeFormTransformationTD::IndexToLattice(int idx, int &i, int &j, int &k, int &l) const
{
  int n, m;

  n = _size.x * _size.y * _size.z * _size.w;
  if (idx > n) {
    idx -= n;
    if (idx > n) idx -= n;
  }

  l = idx / (_size.x * _size.y * _size.z);
  n = idx % (_size.x * _size.y * _size.z);
  m = _size.x * _size.y;
  k = n / m;
  j = n % m / _size.x;
  i = n % m % _size.x;
}

// -----------------------------------------------------------------------------
int4 irtkCUBSplineFreeFormTransformationTD::IndexToLattice(int idx) const
{
  int4 p;
  IndexToLattice(idx, p.x, p.y, p.z, p.w);
  return p;
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
real3 irtkCUBSplineFreeFormTransformationTD::EvaluateFFD(realt x, realt y, realt z, realt t) const
{
  realt B_L, B_K, B_J, B_I;

  int4  intp = make_int4 (floorf(x),  floorf(y),  floorf(z),  floorf(t));
  real4 frac = make_real4(x - intp.x, y - intp.y, z - intp.z, t - intp.w);

  const int offset1 =  _size.x + 4;
  const int offset2 = (_size.x + 8) * (_size.y + 4);
  const int offset3 = (_size.x + 8) * (_size.y + 8) * (_size.z + 4);

  const real3 *cp = &_cp[LatticeToIndex(intp - 1)];
  real3        v  = make_real3(.0, .0, .0);

  #pragma unroll
  for (int l = 0; l < 4; l++) {
    B_L = irtkBSpline<realt>::B(l, frac.w);
    #pragma unroll
    for (int k = 0; k < 4; k++) {
      B_K = irtkBSpline<realt>::B(k, frac.z) * B_L;
      #pragma unroll
      for (int j = 0; j < 4; j++) {
        B_J = irtkBSpline<realt>::B(j, frac.y) * B_K;
        #pragma unroll
        for (int i = 0; i < 4; i++) {
          B_I = irtkBSpline<realt>::B(i, frac.x) * B_J;
          v += *(cp++) * B_I;
        }
        cp += offset1;
      }
      cp += offset2;
    }
    cp += offset3;
  }

  return v;
}

// -----------------------------------------------------------------------------
real3 irtkCUBSplineFreeFormTransformationTD::EvaluateFFD(real3 p, realt t) const
{
  return EvaluateFFD(p.x, p.y, p.z, t);
}

// -----------------------------------------------------------------------------
real3 irtkCUBSplineFreeFormTransformationTD::Evaluate(real3 p, realt t) const
{
  return EvaluateFFD(WorldToLattice(p), TimeToLattice(t));
}


#endif
