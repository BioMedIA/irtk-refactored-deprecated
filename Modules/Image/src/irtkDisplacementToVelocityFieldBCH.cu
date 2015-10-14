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


namespace irtkCUDisplacementToVelocityFieldBCH {


__constant__ unsigned int _N; ///< Number of voxels


// ============================================================================
// Kernels
// ============================================================================

// ----------------------------------------------------------------------------
template <class VoxelType>
__global__ void cuBCH2(VoxelType *v, const VoxelType *dv)
{
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < _N) v[i] += dv[i];
}

// ----------------------------------------------------------------------------
template <class VoxelType>
__global__ void cuBCH3(VoxelType *v, const VoxelType *dv, const VoxelType *l1)
{
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < _N) v[i] += dv[i] + l1[i] / 2.0;
}

// ----------------------------------------------------------------------------
template <class VoxelType>
__global__ void cuBCH4(VoxelType *v, const VoxelType *dv, const VoxelType *l1, const VoxelType *l2)
{
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < _N) v[i] += dv[i] + l1[i] / 2.0 + l2[i] / 12.0;
}

// ----------------------------------------------------------------------------
template <class VoxelType>
__global__ void cuBCH4(VoxelType *v, const VoxelType *dv, const VoxelType *l1, const VoxelType *l2, const VoxelType *l3)
{
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < _N) v[i] += dv[i] + l1[i] / 2.0 + l2[i] / 12.0 - l3[i] / 12.0;
}

// ----------------------------------------------------------------------------
template <class VoxelType>
__global__ void cuBCH4(VoxelType *v, const VoxelType *dv, const VoxelType *l1, const VoxelType *l2, const VoxelType *l3, const VoxelType *l4)
{
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < _N) v[i] += dv[i] + l1[i] / 2.0 + l2[i] / 12.0 - l3[i] / 12.0 - l4[i] / 24.0;
}

// ============================================================================
// Host function
// ============================================================================

// ----------------------------------------------------------------------------
template <class VoxelType>
IRTKCU_HOST_API void Update(irtkCUGenericImage<VoxelType>       *v,
                            const irtkCUGenericImage<VoxelType> *dv,
                            const irtkCUGenericImage<VoxelType> *l1,
                            const irtkCUGenericImage<VoxelType> *l2,
                            const irtkCUGenericImage<VoxelType> *l3,
                            const irtkCUGenericImage<VoxelType> *l4, int n)
{
  dim3 blocks, threads;
  DefaultKernelConfiguration1D(blocks, threads, v->GetX(), v->GetY(), v->GetZ(), v->GetT());

  const unsigned int N = static_cast<unsigned int>(v->GetNumberOfVoxels());
  CudaSafeCall( cudaMemcpyToSymbol(_N, &N, sizeof(N)) );

  VoxelType       *d_v  = v ->GetPointerToDevice();
  const VoxelType *d_dv = dv->GetPointerToDevice();
  const VoxelType *d_l1 = l1->GetPointerToDevice();
  const VoxelType *d_l2 = l2->GetPointerToDevice();
  const VoxelType *d_l3 = l3->GetPointerToDevice();
  const VoxelType *d_l4 = l4->GetPointerToDevice();

  IRTKCU_START_TIMING();
  if        (n == 2) {
    cuBCH2<VoxelType><<<blocks, threads>>>(d_v, d_dv);
    IRTKCU_DEBUG_INTERIM_TIMING(3, "irtkCUDisplacementToVelocityFieldBCH (cuBCH2)");
  } else if (n == 3) {
    cuBCH3<VoxelType><<<blocks, threads>>>(d_v, d_dv, d_l1);
    IRTKCU_DEBUG_INTERIM_TIMING(3, "irtkCUDisplacementToVelocityFieldBCH (cuBCH3)");
  } else if (n == 4) {
    cuBCH4<VoxelType><<<blocks, threads>>>(d_v, d_dv, d_l1, d_l2);
    IRTKCU_DEBUG_INTERIM_TIMING(3, "irtkCUDisplacementToVelocityFieldBCH (cuBCH4)");
  } else if (n == 5) {
    cuBCH4<VoxelType><<<blocks, threads>>>(d_v, d_dv, d_l1, d_l2, d_l3);
    IRTKCU_DEBUG_INTERIM_TIMING(3, "irtkCUDisplacementToVelocityFieldBCH (cuBCH5)");
  } else if (n == 6) {
    cuBCH4<VoxelType><<<blocks, threads>>>(d_v, d_dv, d_l1, d_l2, d_l3, d_l4);
    IRTKCU_DEBUG_INTERIM_TIMING(3, "irtkCUDisplacementToVelocityFieldBCH (cuBCH6)");
  } else {
    cerr << "irtkCUDisplacementToVelocityFieldBCH::Update: Unsupported number of BCH terms: " << n << endl;
    exit(1);
  }
  CudaSyncCheckError();
  IRTKCU_DEBUG_TIMING(2, "irtkCUDisplacementToVelocityFieldBCH::Update");
}

// =============================================================================
// Explicit template instantiations
// =============================================================================

#define INSTANTIATE_FOR_(T)                                                    \
    template IRTKCU_HOST_API void Update<T>(irtkCUGenericImage<T> *,           \
                                            const irtkCUGenericImage<T> *,     \
                                            const irtkCUGenericImage<T> *,     \
                                            const irtkCUGenericImage<T> *, int n)

INSTANTIATE_FOR_(float);
INSTANTIATE_FOR_(double);


} // namespace irtkCUDisplacementToVelocityFieldBCH
