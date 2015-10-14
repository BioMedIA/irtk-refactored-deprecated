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

#ifndef IRTKIMAGETOIMAGE_CUH_
#define IRTKIMAGETOIMAGE_CUH_


#include <irtkCUCommon.h>
#include <irtkCUImage.h>


namespace irtkCUImageToImage {


// ============================================================================
// Global shared device memory symbols
// ============================================================================

// ----------------------------------------------------------------------------
// Attributes of vector fields
__constant__ unsigned int _X;      ///< Number of voxels in x dimension
__constant__ unsigned int _Y;      ///< Number of voxels in y dimension
__constant__ unsigned int _Z;      ///< Number of voxels in z dimension
__constant__ unsigned int _T;      ///< Number of voxels in t dimension
__constant__ unsigned int _XYZ;    ///< Number of spatial voxels (_X*_Y*_Z)
__constant__ unsigned int _XYZT;   ///< Number of voxels (_X*_Y*_Z*_T)
__constant__ float3x4     _matI2W; ///< Image to world matrix
__constant__ float3x4     _matW2I; ///< World to image matrix

// ============================================================================
// Initialization
// ============================================================================

// ----------------------------------------------------------------------------
/// Copy attributes of image to global shared device memory symbols
///
/// This code is provided as macro because for some reason it has to be in the
/// same .cu file which contains the kernel and the kernel launch host function.
/// At least it only worked for one .cu file, but in another the constants
/// were just zero! See for example the irtkDisplacementToVelocityFieldBCH filter.
/// Here, the irtkVelocityToDisplacementFieldEuler and irtkLieBracketImageFilter
/// are both used. The first had invalid values for _X, _Y, _Z,... but those
/// used by the kernel of the second filter were just fine. As I could not
/// figure a better solution besides copy&paste both __constant__ declarations
/// and the code of this macro into each .cu file, I chose this variant, i.e.,
/// replacing a template function called Initialize which contained the code
/// of this macro by this macro definition. -as12312
#define MemcpyImageAttributesToSymbols(img)                                    \
    do {                                                                       \
      if ((img)->GetX() < 1 || (img)->GetY() < 1 || (img)->GetZ() < 1 || (img)->GetT() < 1) { \
        cerr << "MemcpyImageAttributesToSymbol: Image has invalid size attributes!" << endl; \
        exit(1);                                                               \
      }                                                                        \
                                                                               \
      const unsigned int X    = (img)->GetX();                                 \
      const unsigned int Y    = (img)->GetY();                                 \
      const unsigned int Z    = (img)->GetZ();                                 \
      const unsigned int T    = (img)->GetT();                                 \
      const unsigned int XYZ  = X * Y * Z;                                     \
      const unsigned int XYZT = X * Y * Z * T;                                 \
                                                                               \
      const float3x4 matW2I = make_float3x4((img)->GetWorldToImageMatrix());   \
      const float3x4 matI2W = make_float3x4((img)->GetImageToWorldMatrix());   \
                                                                               \
      CudaSafeCall( cudaMemcpyToSymbol(_X,      &X,      sizeof(X)) );         \
      CudaSafeCall( cudaMemcpyToSymbol(_Y,      &Y,      sizeof(Y)) );         \
      CudaSafeCall( cudaMemcpyToSymbol(_Z,      &Z,      sizeof(Z)) );         \
      CudaSafeCall( cudaMemcpyToSymbol(_T,      &T,      sizeof(T)) );         \
      CudaSafeCall( cudaMemcpyToSymbol(_XYZ,    &XYZ,    sizeof(XYZ)) );       \
      CudaSafeCall( cudaMemcpyToSymbol(_XYZT,   &XYZT,   sizeof(XYZT)) );      \
      CudaSafeCall( cudaMemcpyToSymbol(_matW2I, &matW2I, sizeof(matW2I)) );    \
      CudaSafeCall( cudaMemcpyToSymbol(_matI2W, &matI2W, sizeof(matI2W)) );    \
    } while (false)

// ============================================================================
// Coordinate transformations (device)
// ============================================================================

// ----------------------------------------------------------------------------
/// Map 2D voxel to world coordinates
IRTKCU_DEVICE_API inline float2 image2world(const float2 &p)
{
  return _matI2W * p;
}

// ----------------------------------------------------------------------------
/// Map 3D voxel to world coordinates
IRTKCU_DEVICE_API inline float3 image2world(const float3 &p)
{
  return _matI2W * p;
}

// ----------------------------------------------------------------------------
/// Map 2D world coordinates to voxel
IRTKCU_DEVICE_API inline float2 world2image(const float2 &p)
{
  return _matW2I * p;
}

// ----------------------------------------------------------------------------
/// Map 3D world coordinates to voxel
IRTKCU_DEVICE_API inline float3 world2image(const float3 &p)
{
  return _matW2I * p;
}

// ----------------------------------------------------------------------------
/// Map voxel coordinates to index
IRTKCU_DEVICE_API inline int voxel2index(int i, int j)
{
  return j * _X + i;
}

// ----------------------------------------------------------------------------
/// Map voxel coordinates to index
IRTKCU_DEVICE_API inline int voxel2index(const int2 &v)
{
  return voxel2index(v.x, v.y);
}

// ----------------------------------------------------------------------------
/// Map voxel coordinates to index
IRTKCU_DEVICE_API inline unsigned int voxel2index(unsigned int i, unsigned int j)
{
  return j * _X + i;
}

// ----------------------------------------------------------------------------
/// Map voxel coordinates to index
IRTKCU_DEVICE_API inline unsigned int voxel2index(const uint2 &v)
{
  return voxel2index(v.x, v.y);
}

// ----------------------------------------------------------------------------
/// Map voxel coordinates to index
IRTKCU_DEVICE_API inline int voxel2index(int i, int j, int k)
{
  return (k * _Y + j) * _X + i;
}

// ----------------------------------------------------------------------------
/// Map voxel coordinates to index
IRTKCU_DEVICE_API inline int voxel2index(const int3 &v)
{
  return voxel2index(v.x, v.y, v.z);
}

// ----------------------------------------------------------------------------
/// Map voxel coordinates to index
IRTKCU_DEVICE_API inline unsigned int voxel2index(unsigned int i, unsigned int j, unsigned int k)
{
  return (k * _Y + j) * _X + i;
}

// ----------------------------------------------------------------------------
/// Map voxel coordinates to index
IRTKCU_DEVICE_API inline unsigned int voxel2index(const uint3 &v)
{
  return voxel2index(v.x, v.y, v.z);
}
 
// ----------------------------------------------------------------------------
/// Map index to voxel coordinates
IRTKCU_DEVICE_API inline int3 index2voxel(int idx)
{
  const int n = static_cast<int>(_X * _Y);
  return make_int3((idx % n) % _X, (idx % n) / _X, idx / n);
}

// ----------------------------------------------------------------------------
/// Map index to voxel coordinates
IRTKCU_DEVICE_API inline uint3 index2voxel(unsigned int idx)
{
  const unsigned int n = _X * _Y;
  return make_uint3((idx % n) % _X, (idx % n) / _X, idx / n);
}


} // namespace irtkCUImageToImage


#endif
