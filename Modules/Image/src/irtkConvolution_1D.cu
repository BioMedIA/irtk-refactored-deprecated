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


namespace irtkCUConvolution_1D {


// =============================================================================
// Device functions
// =============================================================================

// -----------------------------------------------------------------------------
// Apply mirror boundary condition along x
IRTKCU_DEVICE_API inline int3 mirror(int3 p)
{
  const int X = static_cast<int>(_X - 1);
  int m, n;

  if (p.x < 0) {
    p.x = -p.x;
    n   = p.x / X;
    m   = p.x - n * X;
    if (n % 2) p.x = X - 1 - m;
    else       p.x = m;
  } else if (p.x > X) {
    p.x = p.x - X;
    n   = static_cast<int>(p.x / X);
    m   = p.x - n * X;
    if (n % 2) p.x = m;
    else       p.x = X - m;
  }

  return p;
}

// =============================================================================
// Kernel
// =============================================================================

// -----------------------------------------------------------------------------
// Perform convolution along x with mirror boundary condition
template <class VoxelType, class KernelType>
__global__ void cuConvolution1D(VoxelType        *output,
                                const VoxelType  *input,
                                const KernelType *kernel, int radius)
{
  const uint3 v = make_uint3(blockIdx.x * blockDim.x + threadIdx.x,
                             blockIdx.y * blockDim.y + threadIdx.y,
                             blockIdx.z * blockDim.z + threadIdx.z);
  if (v.x >= _X || v.y >= _Y || v.z >= _Z) return;

  double val = .0;

  // go to start of input and end of kernel
  int3 p = make_int3(static_cast<int>(v.x) - radius, v.y, v.z);
  int  k = radius + radius;

  // outside left boundary
  while (p.x < 0 && k >= 0) {
    val += static_cast<double>(input[voxel2index(mirror(p))]) * static_cast<double>(kernel[k]);
    p.x++, k--;
  }
  // inside image domain
  while (static_cast<unsigned int>(p.x) < _X && k >= 0) {
    val += static_cast<double>(input[voxel2index(p)]) * static_cast<double>(kernel[k]);
    p.x++, k--;
  }
  // outside right boundary
  while (k >= 0) {
    val += static_cast<double>(input[voxel2index(mirror(p))]) * static_cast<double>(kernel[k]);
    p.x++, k--;
  }

  // output result of convolution
  output[voxel2index(v)] = static_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
// Perform convolution along x with truncation of kernel at boundary
template <class VoxelType, class KernelType>
__global__ void cuConvolution_1D(VoxelType        *output,
                                 const VoxelType  *input,
                                 const KernelType *kernel, int radius)
{
  const uint3 v = make_uint3(blockIdx.x * blockDim.x + threadIdx.x,
                             blockIdx.y * blockDim.y + threadIdx.y,
                             blockIdx.z * blockDim.z + threadIdx.z);
  if (v.x >= _X || v.y >= _Y || v.z >= _Z) return;

  double val = .0;

  // go to start of input and end of kernel
  int3 p = make_int3(static_cast<int>(v.x) - radius, v.y, v.z);
  int  k = radius + radius;

  // outside left boundary
  if (p.x < 0) k += p.x, p.x = 0;
  // inside image domain
  while (static_cast<unsigned int>(p.x) < _X && k >= 0) {
    val += static_cast<double>(input[voxel2index(p)]) * static_cast<double>(kernel[k]);
    p.x++, k--;
  }

  // output normalized result of convolution
  output[voxel2index(v)] = static_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
// Perform normalized convolution along x with truncation of kernel at boundary
template <class VoxelType, class KernelType>
__global__ void cuNormalizedConvolution_1D(VoxelType        *output,
                                           const VoxelType  *input,
                                           const KernelType *kernel, int radius)
{
  const uint3 v = make_uint3(blockIdx.x * blockDim.x + threadIdx.x,
                             blockIdx.y * blockDim.y + threadIdx.y,
                             blockIdx.z * blockDim.z + threadIdx.z);
  if (v.x >= _X || v.y >= _Y || v.z >= _Z) return;

  double val = .0;
  double sum = .0;

  // go to start of input and end of kernel
  int3 p = make_int3(static_cast<int>(v.x) - radius, v.y, v.z);
  int  k = radius + radius;

  // outside left boundary
  if (p.x < 0) k += p.x, p.x = 0;
  // inside image domain
  while (static_cast<unsigned int>(p.x) < _X && k >= 0) {
    val += static_cast<double>(input[voxel2index(p)]) * static_cast<double>(kernel[k]);
    sum += static_cast<double>(kernel[k]);
    p.x++, k--;
  }

  // output normalized result of convolution
  output[voxel2index(v)] = ((sum > .0) ? static_cast<VoxelType>(val / sum) : .0);
}

// -----------------------------------------------------------------------------
// Normalize (i.e., divide) values by given sum of convolution kernel weights
template <class VoxelType>
__global__ void cuNormalize(VoxelType *output, double sum)
{
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < _XYZT) output[i] /= sum;
}

// =============================================================================
// Host function
// =============================================================================

template <class VoxelType, class KernelType>
void Run(irtkCUGenericImage<VoxelType>        *output,
         const irtkCUGenericImage<VoxelType>  *input,
         const irtkCUGenericImage<KernelType> *kernel, bool normalize = false)
{
  // Do not use mirror boundary condition for now to get the same results as
  // with the CPU implementation even though convolution is usually better be
  // done with the mirror boundary condition than simple truncation.
  // -as12312
  const bool mirror_boundary = false;

  dim3 blocks, threads;

  // Check arguments
  if (kernel->GetX() % 2 == 0) {
    cerr << "irtkCUConvolution1D::Run: Size of convolution kernel in x must be odd!" << endl;
    exit(1);
  }
  if (kernel->GetY() != 1 || kernel->GetZ() != 1 || kernel->GetT() != 1) {
    cerr << "irtkCUConvolution1D::Run: Kernel dimensions has to be 1 in y, z, and t!" << endl;
    exit(1);
  }
  if (input->GetPointerToDevice() == NULL) {
    cerr << "irtkCUConvolution1D::Run: Input device pointer invalid!" << endl;
    exit(1);
  }
  if (kernel->GetPointerToDevice() == NULL) {
    cerr << "irtkCUConvolution1D::Run: Kernel device pointer invalid!" << endl;
    exit(1);
  }

  // Initialize
  VoxelType *d_output = NULL;
  if (input == output) {
    CudaSafeCall( cudaMalloc(&d_output, input->GetNumberOfVoxels() * sizeof(VoxelType)) );
  } else {
    output->Initialize(input->GetImageAttributes());
    d_output = output->GetPointerToDevice();
  }
  if (d_output == NULL) {
    cerr << "irtkCUConvolution1D::Run: Output device pointer invalid!" << endl;
    exit(1);
  }

  // Radius of convolution kernel
  const int radius = kernel->GetX() / 2;

  // Number of spatial voxels
  const int N = input->GetX() * input->GetY() * input->GetZ();

  // Copy image attributes to constant device memory
  MemcpyImageAttributesToSymbols(input);

  // Perform convolution...
  IRTKCU_START_TIMING();
  DefaultKernelConfigurationND(blocks, threads, input->GetX(), input->GetY(), input->GetZ());
  // ...with mirror boundary condition
  if (mirror_boundary) {
    double sum = .0;
    if (normalize) {
      KernelType *ptr2kernel = kernel->GetPointerToVoxels();
      for (int x = 0; x < kernel->GetX(); x++) sum += static_cast<double>(ptr2kernel[x]);
    }
    if (normalize && sum == .0) {
      output->CopyFrom(input->GetPointerToVoxels());
    } else {
      for (int t = 0; t < input->GetT(); t++) {
        cuConvolution1D<VoxelType, KernelType><<<blocks, threads>>>(d_output + t * N,
                                                                    input ->GetPointerToDevice(0, 0, 0, t),
                                                                    kernel->GetPointerToDevice(), radius);
        IRTKCU_DEBUG_INTERIM_TIMING(3, "irtkConvolution1D (cuConvolution1D, t=" << t << ")");
      }
      CudaSyncCheckError();
      if (normalize) {
        DefaultKernelConfiguration1D(blocks, threads, input->GetX(), input->GetY(), input->GetZ(), input->GetT());
        cuNormalize<VoxelType><<<blocks, threads>>>(output->GetPointerToDevice(), sum);
        IRTKCU_DEBUG_INTERIM_TIMING(3, "irtkConvolution1D (cuNormalize)");
      }
    }
  // ...with truncation at boundary
  } else {
    if (normalize) {
      for (int t = 0; t < input->GetT(); t++) {
        cuNormalizedConvolution_1D<VoxelType, KernelType><<<blocks, threads>>>(d_output + t * N,
                                                                               input ->GetPointerToDevice(0, 0, 0, t),
                                                                               kernel->GetPointerToDevice(), radius);
        IRTKCU_DEBUG_INTERIM_TIMING(3, "irtkConvolution1D (cuNormalizedConvolution_1D, t=" << t << ")");
      }
    } else {
      for (int t = 0; t < input->GetT(); t++) {
        cuConvolution_1D<VoxelType, KernelType><<<blocks, threads>>>(d_output + t * N,
                                                                     input ->GetPointerToDevice(0, 0, 0, t),
                                                                     kernel->GetPointerToDevice(), radius);
        IRTKCU_DEBUG_INTERIM_TIMING(3, "irtkConvolution1D (cuConvolution_1D, t=" << t << ")");
      }
    }
  }
  CudaSyncCheckError();
  IRTKCU_DEBUG_TIMING(2, "irtkConvolution1D");

  // Finalize
  if (d_output != output->GetPointerToDevice()) {
    CudaSafeCall( cudaMemcpy(output->GetPointerToDevice(), d_output, output->GetT() * N * sizeof(VoxelType), cudaMemcpyDeviceToDevice) );
    CudaSafeCall( cudaFree(d_output) );
    CudaSafeCall( cudaDeviceSynchronize() );
  }
}

// =============================================================================
// Explicit template instantiations
// =============================================================================

#define INSTANTIATE_FOR_(T)                                                    \
    template IRTKCU_HOST_API void Run<T, irtkRealPixel>(                       \
         irtkCUGenericImage<T>*,                                               \
         const irtkCUGenericImage<T> *,                                        \
         const irtkCUGenericImage<irtkRealPixel> *, bool)

INSTANTIATE_FOR_(unsigned char);
INSTANTIATE_FOR_(short);
INSTANTIATE_FOR_(unsigned short);
INSTANTIATE_FOR_(float);
INSTANTIATE_FOR_(double);


} // namespace irtkCUConvolution_1D
