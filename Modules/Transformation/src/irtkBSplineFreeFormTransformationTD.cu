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
#include <irtkCUTransformation.h>

#include <irtkFreeFormTransformationIntegration.h>


namespace irtkCUBSplineFreeFormTransformationTD {


// =============================================================================
// Displacement field (device)
// =============================================================================

// -----------------------------------------------------------------------------
struct irtkCUDisplacementField
{
  // ---------------------------------------------------------------------------
  uint3   dim; ///< Dimensions of displacement field.
  real3x4 i2w; ///< Image to world transformation matrix.
  real3  *vox; ///< Displacements in global device memory.

  // ---------------------------------------------------------------------------
  /// Constructor
  template <class VoxelType>
  IRTKCU_HOST_API irtkCUDisplacementField(const irtkGenericImage<VoxelType> &disp)
  :
    vox(NULL)
  {
    const size_t n      = disp.GetX() * disp.GetY() * disp.GetZ();
    const size_t nbytes = n * sizeof(real3);

    // Initialize image attributes
    irtkImageAttributes attr = disp.GetImageAttributes();

    dim = make_uint3  (attr._x, attr._y, attr._z);
    i2w = make_real3x4(attr.GetImageToWorldMatrix());

    // Allocate device memory for vector field data
    CudaSafeCall( cudaMalloc(&vox, nbytes) );

    // Convert vector field data in host memory
    real3 *tmp = new real3[n];
    if (tmp == NULL) {
      cerr << "irtkCUDisplacementField: Failed to allocate temporary host memory for data conversion!" << endl;
      exit(1);
    }
    real3 *ptr = tmp;
    for (int k = 0; k < disp.GetZ(); k++) {
      for (int j = 0; j < disp.GetY(); j++) {
        for (int i = 0; i < disp.GetX(); i++) {
          (*ptr++) = make_real3(disp.Get(i, j, k, 0),
                                disp.Get(i, j, k, 1),
                                disp.Get(i, j, k, 2));
        }
      }
    }

    // Copy vector field data from host to device
    CudaSafeCall( cudaMemcpy(vox, tmp, nbytes, cudaMemcpyHostToDevice) );
    delete[] tmp;
  }

  // ---------------------------------------------------------------------------
  /// Destructor
  IRTKCU_HOST_API ~irtkCUDisplacementField()
  {
    // DO NOT free memory here as instance is passed by value
    // (i.e., copied when launching the CUDA kernel).
  }

  // ---------------------------------------------------------------------------
  /// Free device memory
  IRTKCU_HOST_API void FreeDeviceMemory()
  {
    cudaFree(vox);
    vox = NULL;
  }

  // ---------------------------------------------------------------------------
  /// Copy displacements from device memory back to host
  template <class VoxelType>
  IRTKCU_HOST_API void MemcpyToHost(irtkGenericImage<VoxelType> &disp) const
  {
    const size_t n      = dim.x * dim.y * dim.z;
    const size_t nbytes = n * sizeof(real3);

    // Allocate host memory    
    real3 *tmp = new real3[n];
    if (tmp == NULL) {
      cerr << "irtkCUDisplacementField::MemcpyToHost: Failed to allocate temporary host memory for data conversion!" << endl;
      exit(1);
    }

    // Copy vector field data from device to host
    CudaSafeCall( cudaMemcpy(tmp, vox, nbytes, cudaMemcpyDeviceToHost) );
    CudaSafeCall( cudaDeviceSynchronize() );

    // Convert vector field data in host memory
    real3 *ptr = tmp;
    for (int k = 0; k < disp.GetZ(); k++) {
      for (int j = 0; j < disp.GetY(); j++) {
        for (int i = 0; i < disp.GetX(); i++) {
          disp.Put(i, j, k, 0, static_cast<VoxelType>(ptr->x));
          disp.Put(i, j, k, 1, static_cast<VoxelType>(ptr->y));
          disp.Put(i, j, k, 2, static_cast<VoxelType>(ptr->z));
          ptr++;
        }
      }
    }
    delete[] tmp;
  }

  // ---------------------------------------------------------------------------
  /// Transform voxel to world coordinates
  IRTKCU_API real3 ImageToWorld(int3 v) const
  {
    return i2w * make_real3(v);
  }

  // ---------------------------------------------------------------------------
  /// Transform voxel to world coordinates
  IRTKCU_API real3 ImageToWorld(uint3 v) const
  {
    return i2w * make_real3(v);
  }

  // ---------------------------------------------------------------------------
  /// Transform voxel to world coordinates
  IRTKCU_API real3 ImageToWorld(real3 v) const
  {
    return i2w * v;
  }

  // ---------------------------------------------------------------------------
  /// Convert voxel coordinates to index
  IRTKCU_API uint3 IndexToVoxel(unsigned int i) const
  {
    uint3 v;
    const unsigned int n = dim.x * dim.y;
    v.z = i / n;
    v.x = i % n; // intermediate storage only
    v.y = v.x / dim.x;
    v.x = v.x % dim.x;
    return v;
  }

  // ---------------------------------------------------------------------------
  /// Convert voxel coordinates to index
  IRTKCU_API unsigned int VoxelToIndex(uint3 v) const
  {
    return (v.z * dim.y + v.y) * dim.x + v.x;
  }
};

// =============================================================================
// CUDA Kernels
// =============================================================================

// -----------------------------------------------------------------------------
template <unsigned int FFDIM>
__global__ void cuDisplacement(irtkCUDisplacementField               disp,
                               irtkCUBSplineFreeFormTransformation4D ffd,
                               realt t1, realt t2, realt mindt, realt maxdt, realt tol)
{
  const uint3 v = blockIdx * blockDim + threadIdx;
  if (v >= disp.dim) return;
  const unsigned int i = disp.VoxelToIndex(v);
  real3 p1 = disp.ImageToWorld(v);
  real3 p2 = p1 + disp.vox[i];
  p2 = RungeKutta<FFDIM>(ffd, p2, t1, t2, mindt, maxdt, tol);
  disp.vox[i] = p2 - p1;
}

// =============================================================================
// Host functions
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
IRTKCU_HOST_API void Displacement(const irtkBSplineFreeFormTransformationTD &h_ffd,
                                  irtkGenericImage<VoxelType>               &h_disp,
                                  FFDIntegrationMethod                       im,
                                  double t1, double t2, double mindt, double maxdt, double tol)
{
  // Copy data to global device memory
  IRTK_START_TIMING();
  irtkCUBSplineFreeFormTransformation4D d_ffd (h_ffd);
  irtkCUDisplacementField               d_disp(h_disp);
  IRTK_DEBUG_TIMING(2, "irtkBSplineFreeFormTransformationTD::Displacement (H2D)");

  // Evaluate displacements for each voxel
  IRTK_RESET_TIMING();
  {
    dim3 blocks, threads;
    DefaultKernelConfigurationND(blocks, threads, h_disp.GetX(), h_disp.GetY(), h_disp.GetZ());
    IRTKCU_START_TIMING();
    if      (im == FFDIM_RKE1)   cuDisplacement<FFDIM_RKE1>  <<<blocks, threads>>>(d_disp, d_ffd, t1, t2, mindt, maxdt, tol);
    else if (im == FFDIM_RKE2)   cuDisplacement<FFDIM_RKE2>  <<<blocks, threads>>>(d_disp, d_ffd, t1, t2, mindt, maxdt, tol);
    else if (im == FFDIM_RKH2)   cuDisplacement<FFDIM_RKH2>  <<<blocks, threads>>>(d_disp, d_ffd, t1, t2, mindt, maxdt, tol);
    else if (im == FFDIM_RK4)    cuDisplacement<FFDIM_RK4>   <<<blocks, threads>>>(d_disp, d_ffd, t1, t2, mindt, maxdt, tol);
    else if (im == FFDIM_RKEH12) cuDisplacement<FFDIM_RKEH12><<<blocks, threads>>>(d_disp, d_ffd, t1, t2, mindt, maxdt, tol);
    else if (im == FFDIM_RKBS23) cuDisplacement<FFDIM_RKBS23><<<blocks, threads>>>(d_disp, d_ffd, t1, t2, mindt, maxdt, tol);
    else if (im == FFDIM_RKF45)  cuDisplacement<FFDIM_RKF45> <<<blocks, threads>>>(d_disp, d_ffd, t1, t2, mindt, maxdt, tol);
    else if (im == FFDIM_RKCK45) cuDisplacement<FFDIM_RKCK45><<<blocks, threads>>>(d_disp, d_ffd, t1, t2, mindt, maxdt, tol);
    else if (im == FFDIM_RKDP45) cuDisplacement<FFDIM_RKDP45><<<blocks, threads>>>(d_disp, d_ffd, t1, t2, mindt, maxdt, tol);
    else {
      cerr << "irtkCUBSplineFreeFormTransformationTD::Displacement: Invalid integration method: " << im << endl;
      exit(1);
    }
    IRTKCU_DEBUG_TIMING(3, "irtkBSplineFreeFormTransformationTD::Displacement");
  }
  CudaSyncCheckError();
  IRTK_DEBUG_TIMING(2, "irtkBSplineFreeFormTransformationTD::Displacement (main)");

  // Copy displacement field from global device memory
  IRTK_RESET_TIMING();
  d_disp.MemcpyToHost(h_disp);
  IRTK_DEBUG_TIMING(2, "irtkBSplineFreeFormTransformationTD::Displacement (D2H)");

  // Free device memory
  d_disp.FreeDeviceMemory();
}

// =============================================================================
// Explicit template instantiations
// =============================================================================

#define INSTANTIATE_Displacement_FOR_(T)                                       \
    template IRTKCU_HOST_API void Displacement<T>(                             \
         const irtkBSplineFreeFormTransformationTD &h_ffd,                     \
         irtkGenericImage<T>                       &h_disp,                    \
         FFDIntegrationMethod                       im,                        \
         double t1, double t2, double mindt, double maxdt, double tol)

INSTANTIATE_Displacement_FOR_(float);
INSTANTIATE_Displacement_FOR_(double);


} // namespace irtkCUBSplineFreeFormTransformationTD
