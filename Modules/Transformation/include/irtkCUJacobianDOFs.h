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

#ifndef _IRTKCUJACOBIANDOFS_H
#define _IRTKCUJACOBIANDOFS_H

#include <irtkCUMath.h>


/**
 * Sparse matrix for the transformation Jacobian of derivatives w.r.t the parameters.
 *
 * This matrix type only stores the non-zero columns of the Jacobian matrix of
 * the transformation, which contains the derivatives of the transformation
 * w.r.t the transformation parameters. The full Jacobian matrix has dimension
 * 3xN, where N is the number of transformation parameters and the number of
 * rows correspond to the deformation in each spatial dimension (T_x, T_y, T_z).
 */
class irtkCUJacobianDOFs
{
protected:

  uint   _NumberOfColumns; ///< Number of non-zero columns
  uint   _Capacity;        ///< Number of allocated columns
  uint  *_ColumnIndices;   ///< Column indices of non-zero columns
  real3 *_ColumnVectors;   ///< Non-zero column vectors of Jacobian
  bool   _ExternalMemory;  ///< Whether memory is managed externally

public:

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Constructor
  IRTKCU_API irtkCUJacobianDOFs();

  /// Copy constructor
  IRTKCU_API irtkCUJacobianDOFs(const irtkCUJacobianDOFs &);

  /// Destructor  
  IRTKCU_API ~irtkCUJacobianDOFs();

  /// Set external memory to use instead of dynamic allocation of heap memory
  IRTKCU_API void SetExternalMemory(void *, uint);

  /// Set external memory to use instead of dynamic allocation of heap memory
  /// and update input arguments to reflect remaining unused external memory
  IRTKCU_API void UseExternalMemory(void *&, uint &, uint);

  /// Get number of non-zero columns for which storage has been reserved
  IRTKCU_API uint Capacity() const;

  /// Reserve storage for n non-zero columns
  IRTKCU_API void Reserve(uint);

  /// Assignment operator
  IRTKCU_API irtkCUJacobianDOFs &operator =(const irtkCUJacobianDOFs &);

  // ---------------------------------------------------------------------------
  // Element access

  /// Get number of non-zero columns
  IRTKCU_API uint NumberOfNonZeroColumns() const;

  /// Get column index of the i-th non-zero column
  IRTKCU_API uint ColumnIndex(uint) const;

  /// Get i-th non-zero column vector
  IRTKCU_API real3 &ColumnVector(uint);

  /// Get i-th non-zero column vector
  IRTKCU_API const real3 &ColumnVector(uint) const;

  /// Get i-th column vector, inserts new zero vector if necessary
  IRTKCU_API real3 &operator()(uint);

  // ---------------------------------------------------------------------------
  // Operators
  
  /// Add transformation Jacobian to this Jacobian matrix
  IRTKCU_API irtkCUJacobianDOFs &operator +=(const irtkCUJacobianDOFs &);

  /// Multiply this transformation Jacobian by a scalar
  IRTKCU_API irtkCUJacobianDOFs &operator *=(const realt);

  /// Pre-multiply (!) this transformation Jacobian with the given 3x3 matrix
  IRTKCU_API irtkCUJacobianDOFs &operator *=(const real3x3 &);

  // ---------------------------------------------------------------------------
  // Specialized operations

  /// Add transformation Jacobian to this Jacobian matrix
  IRTKCU_API irtkCUJacobianDOFs &add(const irtkCUJacobianDOFs &);

  /// Add fraction of transformation Jacobian to this Jacobian matrix
  IRTKCU_API irtkCUJacobianDOFs &add(const irtkCUJacobianDOFs &, realt);

};


////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
inline irtkCUJacobianDOFs::irtkCUJacobianDOFs()
:
  _NumberOfColumns(0),
  _Capacity(0),
  _ColumnIndices(NULL),
  _ColumnVectors(NULL),
  _ExternalMemory(false)
{
}

// -----------------------------------------------------------------------------
inline irtkCUJacobianDOFs::irtkCUJacobianDOFs(const irtkCUJacobianDOFs &other)
:
  _NumberOfColumns(0),
  _Capacity(0),
  _ColumnIndices(NULL),
  _ColumnVectors(NULL),
  _ExternalMemory(false)
{
  (*this) = other;
}

// -----------------------------------------------------------------------------
inline irtkCUJacobianDOFs::~irtkCUJacobianDOFs()
{
  if (!_ExternalMemory) {
    delete[] _ColumnIndices;
    delete[] _ColumnVectors;
  }
}

// -----------------------------------------------------------------------------
inline void irtkCUJacobianDOFs::SetExternalMemory(void *mem, uint memsz)
{
  _Capacity = memsz / (sizeof(uint) + sizeof(real3));
  if (_Capacity > 0) {
    _ColumnIndices  = reinterpret_cast<uint  *>(mem);
    _ColumnVectors  = reinterpret_cast<real3 *>(_ColumnIndices + _Capacity);
  } else {
    _ColumnIndices = NULL;
    _ColumnVectors = NULL;
  }
  _ExternalMemory = true;
}

// -----------------------------------------------------------------------------
inline void irtkCUJacobianDOFs::UseExternalMemory(void *&mem, uint &memsz, uint c)
{
  const uint chunksz = c * (sizeof(uint) + sizeof(real3));
  if (chunksz > memsz) {
    printf("irtkCUJacobianDOFs::UseExternalMemory: External memory not enough for storing %d column vectors!\n", c);
    #ifndef __CUDACC__
      exit(1);
    #endif
  }
  SetExternalMemory(mem, chunksz);
  mem    = reinterpret_cast<void *>(reinterpret_cast<char *>(mem) + chunksz);
  memsz -= chunksz;
}

// -----------------------------------------------------------------------------
inline uint irtkCUJacobianDOFs::Capacity() const
{
  return _Capacity;
}

// -----------------------------------------------------------------------------
inline void irtkCUJacobianDOFs::Reserve(uint n)
{
  if (n != 0 && n > _Capacity) {
    if (_ExternalMemory) {
      printf("irtkCUJacobianDOFs::Reserve: Run out of external memory (capacity=%d, required=%d)!\n", _Capacity, n);
      #ifndef __CUDACC__
        exit(1);
      #endif
    } else {
      uint  *indices = new uint [n];
      real3 *vectors = new real3[n];
      if (indices == NULL || vectors == NULL) {
        printf("irtkCUJacobianDOFs::Reserve: Failed to allocate memory!\n");
        #ifndef __CUDACC__
          exit(1);
        #endif
      }
      if (_NumberOfColumns > 0) {
        memcpy(indices, _ColumnIndices, _NumberOfColumns * sizeof(uint));
        memcpy(vectors, _ColumnVectors, _NumberOfColumns * sizeof(real3));
      }
      delete[] _ColumnIndices;
      delete[] _ColumnVectors;
      _ColumnIndices = indices;
      _ColumnVectors = vectors;
      _Capacity      = n;      
    }
  }
}

// -----------------------------------------------------------------------------
inline irtkCUJacobianDOFs &irtkCUJacobianDOFs::operator =(const irtkCUJacobianDOFs &rhs)
{
  if (rhs._NumberOfColumns > 0) {
    if (_Capacity < rhs._NumberOfColumns) {
      if (_ExternalMemory) {
        printf("irtkCUJacobianDOFs::operator=: Run out of external memory (capacity=%d, required=%d)!\n", _Capacity, rhs._NumberOfColumns);
        #ifndef __CUDACC__
          exit(1);
        #endif
      } else {
        delete[] _ColumnIndices;
        delete[] _ColumnVectors;
        _ColumnIndices = new uint [rhs._NumberOfColumns];
        _ColumnVectors = new real3[rhs._NumberOfColumns];
        if (_ColumnIndices == NULL || _ColumnVectors == NULL) {
          printf("irtkCUJacobianDOFs::Reserve: Failed to allocate memory!\n");
          #ifndef __CUDACC__
            exit(1);
          #endif
        }
        _Capacity = rhs._NumberOfColumns;
      }
    }
    if (_Capacity >= rhs._NumberOfColumns) {
      memcpy(_ColumnIndices, rhs._ColumnIndices, rhs._NumberOfColumns * sizeof(uint));
      memcpy(_ColumnVectors, rhs._ColumnVectors, rhs._NumberOfColumns * sizeof(real3));
    }
  }
  _NumberOfColumns = rhs._NumberOfColumns;
  return *this;
}

// =============================================================================
// Element access
// =============================================================================

// -----------------------------------------------------------------------------
inline uint irtkCUJacobianDOFs::NumberOfNonZeroColumns() const
{
  return _NumberOfColumns;
}

// -----------------------------------------------------------------------------
inline uint irtkCUJacobianDOFs::ColumnIndex(uint i) const
{
  return _ColumnIndices[i];
}

// -----------------------------------------------------------------------------
inline real3 &irtkCUJacobianDOFs::ColumnVector(uint i)
{
  return _ColumnVectors[i];
}

// -----------------------------------------------------------------------------
inline const real3 &irtkCUJacobianDOFs::ColumnVector(uint i) const
{
  return _ColumnVectors[i];
}

// -----------------------------------------------------------------------------
inline real3 &irtkCUJacobianDOFs::operator()(uint c)
{
  uint i = 0;
  while (i < _NumberOfColumns && _ColumnIndices[i] < c) i++;
  if (i == _NumberOfColumns || _ColumnIndices[i] > c) {
    if (_Capacity == _NumberOfColumns) Reserve(_Capacity + 16);
    for (uint j = _NumberOfColumns; j > i; j--) {
      _ColumnIndices[j] = _ColumnIndices[j-1];
      _ColumnVectors[j] = _ColumnVectors[j-1];
    }
    _ColumnIndices[i] = c;
    _ColumnVectors[i] = make_real3(.0f);
    _NumberOfColumns++;
  }
  return _ColumnVectors[i];
}

// =============================================================================
// Unary operators
// =============================================================================

// -----------------------------------------------------------------------------
inline irtkCUJacobianDOFs &irtkCUJacobianDOFs::operator +=(const irtkCUJacobianDOFs &b)
{
  for (uint i = 0; i < b.NumberOfNonZeroColumns(); i++) {
    this->operator()(b.ColumnIndex(i)) += b.ColumnVector(i);
  }
  return *this;
}

// -----------------------------------------------------------------------------
inline irtkCUJacobianDOFs &irtkCUJacobianDOFs::operator *=(const realt s)
{
  for (uint i = 0; i < _NumberOfColumns; i++) _ColumnVectors[i] *= s;
  return *this;
}

// -----------------------------------------------------------------------------
inline irtkCUJacobianDOFs &irtkCUJacobianDOFs::operator *=(const real3x3 &m)
{
  // Note: Read this operator as: m * (*this)!
  for (uint i = 0; i < _NumberOfColumns; i++) {
    const real3 v = _ColumnVectors[i];
    _ColumnVectors[i].x = m.a.x * v.x + m.a.y * v.y + m.a.z * v.z;
    _ColumnVectors[i].y = m.b.x * v.x + m.b.y * v.y + m.b.z * v.z;
    _ColumnVectors[i].z = m.c.x * v.x + m.c.y * v.y + m.c.z * v.z;
  }
  return *this;
}

// =============================================================================
// Mathematical operations
// =============================================================================

// -----------------------------------------------------------------------------
inline irtkCUJacobianDOFs &irtkCUJacobianDOFs::add(const irtkCUJacobianDOFs &b)
{
  return (*this) += b;
}

// -----------------------------------------------------------------------------
inline irtkCUJacobianDOFs &irtkCUJacobianDOFs::add(const irtkCUJacobianDOFs &b, realt s)
{
  for (uint i = 0; i < b.NumberOfNonZeroColumns(); i++) {
    this->operator()(b.ColumnIndex(i)) += b.ColumnVector(i) * s;
  }
  return *this;
}

// =============================================================================
// Binary operators
// =============================================================================

// -----------------------------------------------------------------------------
/// Calculate column-by-column sum of transformation Jacobian
IRTKCU_API inline irtkCUJacobianDOFs operator +(const irtkCUJacobianDOFs &a, const irtkCUJacobianDOFs &b)
{
  irtkCUJacobianDOFs c = a;
  c += b;
  return c;
}

// -----------------------------------------------------------------------------
/// Multiply transformation Jacobian and scalar
IRTKCU_API inline irtkCUJacobianDOFs operator *(const irtkCUJacobianDOFs &a, double s)
{
  irtkCUJacobianDOFs b = a;
  b *= s;
  return b;
}

// -----------------------------------------------------------------------------
/// Multiply transformation Jacobian and scalar
IRTKCU_API inline irtkCUJacobianDOFs operator *(double s, const irtkCUJacobianDOFs &a)
{
  return a * s;
}

// -----------------------------------------------------------------------------
/// Calculate product of 3x3 matrix and transformation Jacobian
IRTKCU_API inline irtkCUJacobianDOFs operator *(const real3x3 &a, const irtkCUJacobianDOFs &b)
{
  irtkCUJacobianDOFs c = b;
  c *= a;
  return c;
}


#endif
