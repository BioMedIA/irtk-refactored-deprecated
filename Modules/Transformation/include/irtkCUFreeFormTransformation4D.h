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

#ifndef IRTKCUFREEFORMTRANSFORMATION4D_H_
#define IRTKCUFREEFORMTRANSFORMATION4D_H_


/**
 * Base class of spatio-temporal FFD representation for use by CUDA kernels.
 *
 * This structure is intended for use by CUDA device code which is executed
 * within a kernel on the GPGPU.
 *
 * \note No polymorphism, i.e., virtual methods can be used by device code!
 *
 * \note Keep number of data members at a minimum to reduce the memory
 *       footprint of the structure so it can be passed by value to a
 *       CUDA kernel. Current size limit for all kernel arguments is 256 bytes.
 */
class irtkCUFreeFormTransformation4D
{
public:

  // ---------------------------------------------------------------------------
  // Data members

  uint4   _dim;     ///< Dimensions of FFD lattice.
  real3x4 _matW2L;  ///< World to lattice transformation matrix.
  real3x4 _matL2W;  ///< Lattice to world transformation matrix.
  realt   _torigin; ///< Temporal origin of FFD lattice.
  realt   _dt;      ///< Spacing of temporal control points.
  real3  *_data;    ///< Control point data in global device memory.
  bool    _owner;   ///< Whether this copy is responsible for the memory.

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  IRTKCU_HOST_API irtkCUFreeFormTransformation4D();

  /// Constructor
  IRTKCU_HOST_API irtkCUFreeFormTransformation4D(const irtkFreeFormTransformation4D &);

  /// Copy constructor, makes a shallow copy only
  IRTKCU_HOST_API irtkCUFreeFormTransformation4D(const irtkCUFreeFormTransformation4D &);

  /// Destructor
  virtual IRTKCU_HOST_API ~irtkCUFreeFormTransformation4D();

  // ---------------------------------------------------------------------------
  // Initialization

  /// Initialize and copy the control point data to the device memory
  virtual IRTKCU_HOST_API void Initialize(const irtkFreeFormTransformation4D &);

  /// Assignment operator, makes a shallow copy only
  virtual IRTKCU_HOST_API irtkCUFreeFormTransformation4D &operator =(const irtkCUFreeFormTransformation4D &);

  /// Assignment operator, copies control point data to the device memory
  virtual IRTKCU_HOST_API irtkCUFreeFormTransformation4D &operator =(const irtkFreeFormTransformation4D &);

  /// Assignment operator, copies control point data to the device memory
  virtual IRTKCU_HOST_API irtkCUFreeFormTransformation4D &operator =(const irtkFreeFormTransformation4D *);

  /// Reset instance and free (device) memory
  virtual IRTKCU_HOST_API void Reset();

  // ---------------------------------------------------------------------------
  // Attributes

  /// Return number of control points
  IRTKCU_API uint NumberOfControlPoints() const;

  /// Return number of parameters
  IRTKCU_API uint NumberOfDOFs() const;

  // ---------------------------------------------------------------------------
  // Coordinate transformation

  /// Transform world to lattice coordinates
  IRTKCU_API void WorldToLattice(real3 &) const;

  /// Transform time to lattice coordinate
  IRTKCU_API realt TimeToLattice(realt) const;

  /// Transform lattice to world coordinates
  IRTKCU_API void LatticeToWorld(real3 &) const;

  /// Transform lattice to time coordinate
  IRTKCU_API realt LatticeToTime(realt) const;

  /// Convert lattice coordinates to index
  IRTKCU_API int LatticeToIndex(int, int, int = 0, int = 0) const;

  /// Convert lattice coordinates to index
  IRTKCU_API int LatticeToIndex(int3, int = 0) const;

  /// Convert lattice coordinates to index
  IRTKCU_API int LatticeToIndex(int4) const;

  /// Convert lattice coordinates to index
  IRTKCU_API uint LatticeToIndex(uint, uint, uint = 0, uint = 0) const;

  /// Convert control point coordinates to index
  IRTKCU_API uint LatticeToIndex(uint3, uint = 0) const;

  /// Convert control point coordinates to index
  IRTKCU_API uint LatticeToIndex(uint4) const;
  
  /// Convert index to lattice coordinates
  IRTKCU_API uint4 IndexToLattice(uint) const;

  // ---------------------------------------------------------------------------
  // Data access

  /// Get stride of control points data required to advance pointer to data
  /// from one column, one row, one slice, or one frame to another, respectively.
  IRTKCU_API uint4 GetControlPointsStride() const;

  /// Get pointer to control points data
  IRTKCU_API const real3 *GetPointerToControlPoints(int, int, int = 0, int = 0) const;

  /// Get pointer to control points data
  IRTKCU_API const real3 *GetPointerToControlPoints(int3, int = 0) const;

  /// Get pointer to control points data
  IRTKCU_API const real3 *GetPointerToControlPoints(int4) const;

  /// Get pointer to control points data
  IRTKCU_API const real3 *GetPointerToControlPoints(uint, uint, uint = 0, uint = 0) const;

  /// Get pointer to control points data
  IRTKCU_API const real3 *GetPointerToControlPoints(uint3, uint = 0) const;

  /// Get pointer to control points data
  IRTKCU_API const real3 *GetPointerToControlPoints(uint4) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
inline irtkCUFreeFormTransformation4D::irtkCUFreeFormTransformation4D()
:
  _data(NULL),
  _owner(true)
{
}

// -----------------------------------------------------------------------------
inline irtkCUFreeFormTransformation4D
::irtkCUFreeFormTransformation4D(const irtkFreeFormTransformation4D &ffd)
:
  _data (NULL),
  _owner(true)
{
  Initialize(ffd);
}

// -----------------------------------------------------------------------------
inline irtkCUFreeFormTransformation4D
::irtkCUFreeFormTransformation4D(const irtkCUFreeFormTransformation4D &other)
:
  _dim    (other._dim),
  _matW2L (other._matW2L),
  _matL2W (other._matL2W),
  _torigin(other._torigin),
  _dt     (other._dt),
  _data   (other._data),
  _owner  (false)
{
}

// -----------------------------------------------------------------------------
inline irtkCUFreeFormTransformation4D::~irtkCUFreeFormTransformation4D()
{
  Reset();
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
inline void irtkCUFreeFormTransformation4D::Initialize(const irtkFreeFormTransformation4D &ffd)
{
  // Free previously allocated device memory
  Reset();

  // Initialize FFD attributes
  irtkFFDAttributes attr = ffd.GetFFDAttributes();

  _dim     = make_uint4  (attr._x, attr._y, attr._z, attr._t);
  _matW2L  = make_real3x4(attr.GetWorldToLatticeMatrix());
  _matL2W  = make_real3x4(attr.GetLatticeToWorldMatrix());
  _torigin = attr._torigin;
  _dt      = attr._dt;

  // Allocate device memory for control point data
  const size_t numcps = (attr._x+8) * (attr._y+8) * (attr._z+8) * (attr._t+8);
  const size_t nbytes = numcps * sizeof(real3);
  CudaSafeCall( cudaMalloc(&_data, nbytes) );
  _owner = true;

  // Convert control point data in host memory
  real3 *data = new real3[numcps];
  real3 *it   = data;
  if (it == NULL) {
    cerr << "irtkCUFreeFormTransformation4D: Failed to allocate temporary host memory for data conversion!" << endl;
    exit(1);
  }
  double x, y, z;
  for (int l = -4; l < attr._t+4; l++) {
    for (int k = -4; k < attr._z+4; k++) {
      for (int j = -4; j < attr._y+4; j++) {
        for (int i = -4; i < attr._x+4; i++) {
          ffd.Get(i, j, k, l, x, y, z);
          (*it++) = make_real3(x, y, z);
        }
      }
    }
  }

  // Copy control point data from host to device
  CudaSafeCall( cudaMemcpy(_data, data, nbytes, cudaMemcpyHostToDevice) );
  delete[] data;
}

// -----------------------------------------------------------------------------
inline irtkCUFreeFormTransformation4D &
irtkCUFreeFormTransformation4D::operator =(const irtkCUFreeFormTransformation4D &rhs)
{
  this->_dim     = rhs._dim;
  this->_matW2L  = rhs._matW2L;
  this->_matL2W  = rhs._matL2W;
  this->_torigin = rhs._torigin;
  this->_dt      = rhs._dt;
  this->_data    = rhs._data;
  this->_owner   = false;
  return *this;
}

// -----------------------------------------------------------------------------
inline irtkCUFreeFormTransformation4D &
irtkCUFreeFormTransformation4D::operator =(const irtkFreeFormTransformation4D &rhs)
{
  Initialize(rhs);
  return *this;
}

// -----------------------------------------------------------------------------
inline irtkCUFreeFormTransformation4D &
irtkCUFreeFormTransformation4D::operator =(const irtkFreeFormTransformation4D *rhs)
{
  Initialize(*rhs);
  return *this;
}

// -----------------------------------------------------------------------------
inline void irtkCUFreeFormTransformation4D::Reset()
{
  if (_owner) {
    cudaFree(_data);
    _data = NULL;
  }
}

// =============================================================================
// Attributes
// =============================================================================

// -----------------------------------------------------------------------------
inline uint irtkCUFreeFormTransformation4D::NumberOfControlPoints() const
{
  return _dim.x * _dim.y * _dim.z * _dim.w;
}

// -----------------------------------------------------------------------------
inline uint irtkCUFreeFormTransformation4D::NumberOfDOFs() const
{
  return 3 * NumberOfControlPoints();
}

// =============================================================================
// Coordinate transformation
// =============================================================================

// -----------------------------------------------------------------------------
inline void irtkCUFreeFormTransformation4D::WorldToLattice(real3 &p) const
{
  p = _matW2L * p;
}

// -----------------------------------------------------------------------------
inline realt irtkCUFreeFormTransformation4D::TimeToLattice(realt t) const
{
  return (t - _torigin) / _dt;
}

// -----------------------------------------------------------------------------
inline void irtkCUFreeFormTransformation4D::LatticeToWorld(real3 &p) const
{
  p = _matL2W * p;
}

// -----------------------------------------------------------------------------
inline realt irtkCUFreeFormTransformation4D::LatticeToTime(realt t) const
{
  return _torigin + t * _dt;
}

// -----------------------------------------------------------------------------
inline int irtkCUFreeFormTransformation4D::LatticeToIndex(int i, int j, int k, int l) const
{
  return i + _dim.x * (j + _dim.y * (k + _dim.z * l));
}

// -----------------------------------------------------------------------------
inline int irtkCUFreeFormTransformation4D::LatticeToIndex(int3 p, int l) const
{
  return LatticeToIndex(p.x, p.y, p.z, l);
}

// -----------------------------------------------------------------------------
inline int irtkCUFreeFormTransformation4D::LatticeToIndex(int4 p) const
{
  return LatticeToIndex(p.x, p.y, p.z, p.w);
}

// -----------------------------------------------------------------------------
inline uint irtkCUFreeFormTransformation4D::LatticeToIndex(uint i, uint j, uint k, uint l) const
{
  return i + _dim.x * (j + _dim.y * (k + _dim.z * l));
}

// -----------------------------------------------------------------------------
inline uint irtkCUFreeFormTransformation4D::LatticeToIndex(uint3 p, uint l) const
{
  return LatticeToIndex(p.x, p.y, p.z, l);
}

// -----------------------------------------------------------------------------
inline uint irtkCUFreeFormTransformation4D::LatticeToIndex(uint4 p) const
{
  return LatticeToIndex(p.x, p.y, p.z, p.w);
}

// -----------------------------------------------------------------------------
inline uint4 irtkCUFreeFormTransformation4D::IndexToLattice(uint index) const
{
  uint n = NumberOfControlPoints();
  uint m = n / _dim.w;

  if (index > n) {
    index -= n;
    if (index > n) index -= n;
  }

  uint4 cp;
  cp.w = index / m;
  n = index % m;
  m = _dim.x * _dim.y;
  cp.z = n / m;
  cp.y = n % m / _dim.x;
  cp.x = n % m % _dim.x;
  return cp;
}

// =============================================================================
// Data access
// =============================================================================

// -----------------------------------------------------------------------------
inline uint4 irtkCUFreeFormTransformation4D::GetControlPointsStride() const
{
  return make_uint4(1,
                    (_dim.x+4),
                    (_dim.x+8) * (_dim.y+4),
                    (_dim.x+8) * (_dim.y+8) * (_dim.z+4));
}

// -----------------------------------------------------------------------------
inline const real3 *irtkCUFreeFormTransformation4D::GetPointerToControlPoints(int i, int j, int k, int l) const
{
  return &_data[(i+4) + (_dim.x+8) * ((j+4) + (_dim.y+8) * ((k+4) + (_dim.z+8) * (l+4)))];
}

// -----------------------------------------------------------------------------
inline const real3 *irtkCUFreeFormTransformation4D::GetPointerToControlPoints(int3 p, int l) const
{
  return GetPointerToControlPoints(p.x, p.y, p.z, l);
}

// -----------------------------------------------------------------------------
inline const real3 *irtkCUFreeFormTransformation4D::GetPointerToControlPoints(int4 p) const
{
  return GetPointerToControlPoints(p.x, p.y, p.z, p.w);
}

// -----------------------------------------------------------------------------
inline const real3 *irtkCUFreeFormTransformation4D::GetPointerToControlPoints(uint i, uint j, uint k, uint l) const
{
  return &_data[(i+4) + (_dim.x+8) * ((j+4) + (_dim.y+8) * ((k+4) + (_dim.z+8) * (l+4)))];
}

// -----------------------------------------------------------------------------
inline const real3 *irtkCUFreeFormTransformation4D::GetPointerToControlPoints(uint3 p, uint l) const
{
  return GetPointerToControlPoints(p.x, p.y, p.z, l);
}

// -----------------------------------------------------------------------------
inline const real3 *irtkCUFreeFormTransformation4D::GetPointerToControlPoints(uint4 p) const
{
  return GetPointerToControlPoints(p.x, p.y, p.z, p.w);
}




#endif
