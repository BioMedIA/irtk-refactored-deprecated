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

#include <irtkCUCommon.h>
#include <irtkCUImage.h>

#include "irtkImageToImage.cuh"
using namespace irtkCUImageToImage;


namespace irtkCULieBracketImageFilter {


// ---------------------------------------------------------------------------
/// Multiply l with upper 2x2 matrix of r
IRTKCU_API inline float2x2 mult2x2(const float2x2 &l, const float3x4 &r)
{
  float2x2 o;
  o.a.x = l.a.x * r.a.x + l.a.y * r.b.x;
  o.a.y = l.a.x * r.a.y + l.a.y * r.b.y;
  o.b.x = l.b.x * r.a.x + l.b.y * r.b.x;
  o.b.y = l.b.x * r.a.y + l.b.y * r.b.y;
  return o;
}

// ---------------------------------------------------------------------------
/// Multiply l with upper 3x3 matrix of r
IRTKCU_API inline float3x3 mult3x3(const float3x3 &l, const float3x4 &r)
{
  float3x3 o;
  o.a.x = l.a.x * r.a.x + l.a.y * r.b.x + l.a.z * r.c.x;
  o.a.y = l.a.x * r.a.y + l.a.y * r.b.y + l.a.z * r.c.y;
  o.a.z = l.a.x * r.a.z + l.a.y * r.b.z + l.a.z * r.c.z;
  o.b.x = l.b.x * r.a.x + l.b.y * r.b.x + l.b.z * r.c.x;
  o.b.y = l.b.x * r.a.y + l.b.y * r.b.y + l.b.z * r.c.y;
  o.b.z = l.b.x * r.a.z + l.b.y * r.b.z + l.b.z * r.c.z;
  return o;
}

// ---------------------------------------------------------------------------
/// Evaluate Jacobian of 2D vector field at given voxel coordinate
template <class VoxelType>
IRTKCU_DEVICE_API float2x2 Jacobian(const VoxelType *img, int2 p)
{
  float2x2 jac;
  int i, j;

  // Derivatives along x axis
  if (p.x <= 0 || p.x >= _X-1) {
    jac.a.x = 0;
    jac.b.x = 0;
  } else {
    // 1st vector component
    i = voxel2index(make_int2(p.x-1, p.y));
    j = voxel2index(make_int2(p.x+1, p.y));
    jac.a.x = 0.5 * (img[j] - img[i]);
    // 2nd vector component
    i += _XYZ, j += _XYZ;
    jac.b.x = 0.5 * (img[j] - img[i]);
  }

  // Derivatives along y axis
  if (p.y <= 0 || p.y >= _Y-1) {
    jac.a.y = 0;
    jac.b.y = 0;
  } else {
    // 1st vector component
    i = voxel2index(make_int2(p.x, p.y-1));
    j = voxel2index(make_int2(p.x, p.y+1));
    jac.a.y = 0.5 * (img[j] - img[i]);
    // 2nd vector component
    i += _XYZ, j += _XYZ;
    jac.b.y = 0.5 * (img[j] - img[i]);
  }

  return mult2x2(jac, _matW2I);
}

// ---------------------------------------------------------------------------
/// Evaluate Jacobian of 3D vector field at given voxel coordinate
template <class VoxelType>
IRTKCU_DEVICE_API float3x3 Jacobian(const VoxelType *img, int3 p)
{
  float3x3 jac;
  int i, j;

  // Derivatives along x axis
  if (p.x <= 0 || p.x >= _X-1) {
    jac.a.x = 0;
    jac.b.x = 0;
    jac.c.x = 0;
  } else {
    // 1st vector component
    i = voxel2index(make_int3(p.x-1, p.y, p.z));
    j = voxel2index(make_int3(p.x+1, p.y, p.z));
    jac.a.x = 0.5 * (img[j] - img[i]);
    // 2nd vector component
    i += _XYZ, j += _XYZ;
    jac.b.x = 0.5 * (img[j] - img[i]);
    // 3rd vector component
    i += _XYZ, j += _XYZ;
    jac.c.x = 0.5 * (img[j] - img[i]);
  }

  // Derivatives along y axis
  if (p.y <= 0 || p.y >= _Y-1) {
    jac.a.y = 0;
    jac.b.y = 0;
    jac.c.y = 0;
  } else {
    // 1st vector component
    i = voxel2index(make_int3(p.x, p.y-1, p.z));
    j = voxel2index(make_int3(p.x, p.y+1, p.z));
    jac.a.y = 0.5 * (img[j] - img[i]);
    // 2nd vector component
    i += _XYZ, j += _XYZ;
    jac.b.y = 0.5 * (img[j] - img[i]);
    // 3rd vector component
    i += _XYZ, j += _XYZ;
    jac.c.y = 0.5 * (img[j] - img[i]);
  }

  // Derivatives along z axis
  if (p.z <= 0 || p.z >= _Z-1) {
    jac.a.z = 0;
    jac.b.z = 0;
    jac.c.z = 0;
  } else {
    // 1st vector component
    i = voxel2index(make_int3(p.x, p.y, p.z-1));
    j = voxel2index(make_int3(p.x, p.y, p.z+1));
    jac.a.z = 0.5 * (img[j] - img[i]);
    // 2nd vector component
    i += _XYZ, j += _XYZ;
    jac.b.z = 0.5 * (img[j] - img[i]);
    // 3rd vector component
    i += _XYZ, j += _XYZ;
    jac.c.z = 0.5 * (img[j] - img[i]);
  }

  return mult3x3(jac, _matW2I);
}

// ============================================================================
// Kernels
// ============================================================================

// ----------------------------------------------------------------------------
/// Run Lie bracket filter for each 2D voxel
template <class VoxelType>
__global__ void cuLieBracket2D(VoxelType *output, const VoxelType *lv, const VoxelType *rv)
{
  int2 p = make_int2(blockIdx.x * blockDim.x + threadIdx.x,
                     blockIdx.y * blockDim.y + threadIdx.y);
  if (p.x >= _X || p.y >= _Y) return;

  const int x = voxel2index(p); // Index of 1st vector component
  const int y = x + _XYZ;       // Index of 2nd vector component

  const float2x2 lJ = Jacobian(lv, p);
  const float2x2 rJ = Jacobian(rv, p);

  output[x] = (lJ.a.x * rv[x] - lv[x] * rJ.a.x) + (lJ.a.y * rv[y] - lv[y] * rJ.a.y);
  output[y] = (lJ.b.x * rv[x] - lv[x] * rJ.b.x) + (lJ.b.y * rv[y] - lv[y] * rJ.b.y);
}

// ----------------------------------------------------------------------------
/// Run Lie bracket filter for each 3D voxel
template <class VoxelType>
__global__ void cuLieBracket3D(VoxelType *output, const VoxelType *lv, const VoxelType *rv)
{
  int3 p = make_int3(blockIdx.x * blockDim.x + threadIdx.x,
                     blockIdx.y * blockDim.y + threadIdx.y,
                     blockIdx.z * blockDim.z + threadIdx.z);
  if (p.x >= _X || p.y >= _Y || p.z >= _Z) return;

  const int x = voxel2index(p); // Index of 1st vector component
  const int y = x + _XYZ;       // Index of 2nd vector component
  const int z = y + _XYZ;       // Index of 3rd vector component

  const float3x3 lJ = Jacobian(lv, p);
  const float3x3 rJ = Jacobian(rv, p);

  output[x] = (lJ.a.x * rv[x] - lv[x] * rJ.a.x) + (lJ.a.y * rv[y] - lv[y] * rJ.a.y) + (lJ.a.z * rv[z] - lv[z] * rJ.a.z);
  output[y] = (lJ.b.x * rv[x] - lv[x] * rJ.b.x) + (lJ.b.y * rv[y] - lv[y] * rJ.b.y) + (lJ.b.z * rv[z] - lv[z] * rJ.b.z);
  output[z] = (lJ.c.x * rv[x] - lv[x] * rJ.c.x) + (lJ.c.y * rv[y] - lv[y] * rJ.c.y) + (lJ.c.z * rv[z] - lv[z] * rJ.c.z);
}

// ============================================================================
// Host function
// ============================================================================

// ---------------------------------------------------------------------------
/// Compute Lie bracket of two vector fields
template <class VoxelType>
void Run(irtkCUGenericImage<VoxelType> *output, const irtkCUGenericImage<VoxelType> *lv,
                                                const irtkCUGenericImage<VoxelType> *rv)
{
  const int X = output->GetX();
  const int Y = output->GetY();
  const int Z = output->GetZ();
  const int T = output->GetT();

  // Check that input vector fields have same size
  if ((Z == 1 && T != 2) || (Z >  1 && T != 3)) {
    cerr << "irtkCULieBracketImageFilter::Run: Output image is no valid 2D/3D vector field!" << endl;
    exit(1);
  }
  if (lv->GetX() != X || lv->GetY() != Y || lv->GetZ() != Z || lv->GetT() != T) {
    cerr << "irtkCULieBracketImageFilter::Run: Dimensions of left input and output vector fields do not match!" << endl;
    exit(1);
  }
  if (rv->GetX() != X || rv->GetY() != Y || rv->GetZ() != Z || rv->GetT() != T) {
    cerr << "irtkCULieBracketImageFilter::Run: Dimensions of right input and output vector fields do not match!" << endl;
    exit(1);
  }
  if (lv->GetPointerToDevice() == NULL) {
    cerr << "irtkCULieBracketImageFilter::Run: Left input vector field has invalid device pointer!" << endl;
    exit(1);
  }
  if (lv->GetPointerToDevice() == NULL) {
    cerr << "irtkCULieBracketImageFilter::Run: Right input vector field has invalid device pointer!" << endl;
    exit(1);
  }
  if (output->GetPointerToDevice() == NULL) {
    cerr << "irtkCULieBracketImageFilter::Run: Output vector field has invalid device pointer!" << endl;
    exit(1);
  }

  // Initialize constant device memory
  MemcpyImageAttributesToSymbols(output);

  // Output buffer (if necessary)
  VoxelType *d_output = NULL;
  if (lv == output || rv == output) {
    CudaSafeCall( cudaMalloc(&d_output, output->GetNumberOfVoxels() * sizeof(VoxelType)) );
  } else {
    d_output = output->GetPointerToDevice();
  }

  // Execute kernel on device
  IRTKCU_START_TIMING();
  if (T == 2) {
    dim3 blocks, threads;
    DefaultKernelConfigurationND(blocks, threads, X, Y);
    cuLieBracket2D<<<blocks, threads>>>(d_output, lv->GetPointerToDevice(),
                                                  rv->GetPointerToDevice());
    IRTKCU_DEBUG_INTERIM_TIMING(3, "irtkLieBracketImageFilter (cuLieBracket2D)");
  } else {
    dim3 blocks, threads;
    DefaultKernelConfigurationND(blocks, threads, X, Y, Z);
    cuLieBracket3D<<<blocks, threads>>>(d_output, lv->GetPointerToDevice(),
                                                  rv->GetPointerToDevice());
    IRTKCU_DEBUG_INTERIM_TIMING(3, "irtkLieBracketImageFilter (cuLieBracket3D)");
  }
  CudaSyncCheckError();
  IRTKCU_DEBUG_TIMING(2, "irtkLieBracketImageFilter");

  // Copy output buffer
  if (d_output != output->GetPointerToDevice()) {
    CudaSafeCall( cudaMemcpy(output->GetPointerToDevice(), d_output, output->GetNumberOfVoxels() * sizeof(VoxelType), cudaMemcpyDeviceToDevice) );
    CudaSafeCall( cudaFree(d_output) );
    CudaSafeCall( cudaDeviceSynchronize() );
  }
}

// ============================================================================
// Explicit template instantiations
// ============================================================================

#define INSTANTIATE_FOR_(T)                                                   \
    template IRTKCU_HOST_API void Run<T>(                                     \
         irtkCUGenericImage<T>*,                                              \
         const irtkCUGenericImage<T> *,                                       \
         const irtkCUGenericImage<T> *)

INSTANTIATE_FOR_(float);
INSTANTIATE_FOR_(double);


} // namespace irtkCULieBracketImageFilter
