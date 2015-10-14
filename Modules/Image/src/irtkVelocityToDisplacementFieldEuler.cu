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


namespace irtkCUVelocityToDisplacementFieldEuler {


// ============================================================================
// Global variables
// ============================================================================

// ----------------------------------------------------------------------------
// Parameters of numerical integration
__constant__ int   _NumberOfSteps; ///< Number of integration steps
__constant__ float _h;             ///< Length of each integration step

// ----------------------------------------------------------------------------
/// Texture reference to each component of the second vector field
texture<float, cudaTextureType2D, cudaReadModeElementType> _Input2_2D_x;
texture<float, cudaTextureType2D, cudaReadModeElementType> _Input2_2D_y;
texture<float, cudaTextureType3D, cudaReadModeElementType> _Input2_3D_x;
texture<float, cudaTextureType3D, cudaReadModeElementType> _Input2_3D_y;
texture<float, cudaTextureType3D, cudaReadModeElementType> _Input2_3D_z;

// ============================================================================
// Kernels
// ============================================================================

// ----------------------------------------------------------------------------
// Kernel for forward Euler integration of 2D vector field
template <class VoxelType>
__global__ void cuRKE1_2D(VoxelType *output, const VoxelType *input1)
{
  uint2 i0, v0 = make_uint2(blockIdx.x * blockDim.x + threadIdx.x,
                            blockIdx.y * blockDim.y + threadIdx.y);
  if (v0.x >= _X || v0.y >= _Y) return;

  i0.x = voxel2index(v0); // Index of 1st vector component
  i0.y = i0.x + _XYZ;     // Index of 2nd vector component

  const float2 p0 = image2world(make_float2(v0));
  float2       p  = p0;
  if (input1) {
    p.x += static_cast<float>(input1[i0.x]);
    p.y += static_cast<float>(input1[i0.y]);
  }
  for (int s = 0; s < _NumberOfSteps; s++) {
    float2 v = world2image(p);
    p.x += tex2D(_Input2_2D_x, v.x, v.y) * _h;
    p.y += tex2D(_Input2_2D_y, v.x, v.y) * _h;
  }

  output[i0.x] = static_cast<VoxelType>(p.x - p0.x);
  output[i0.y] = static_cast<VoxelType>(p.y - p0.y);
}

// ----------------------------------------------------------------------------
// Kernel for forward Euler integration of 3D vector field
template <class VoxelType>
__global__ void cuRKE1_3D(VoxelType *output, const VoxelType *input1)
{
  uint3 i0, v0 = make_uint3(blockIdx.x * blockDim.x + threadIdx.x,
                            blockIdx.y * blockDim.y + threadIdx.y,
                            blockIdx.z * blockDim.z + threadIdx.z);
  if (v0.x >= _X || v0.y >= _Y || v0.z >= _Z) return;

  i0.x = voxel2index(v0); // Index of 1st vector component
  i0.y = i0.x + _XYZ;     // Index of 2nd vector component
  i0.z = i0.y + _XYZ;     // Index of 3rd vector component

  const float3 p0 = image2world(make_float3(v0));
  float3       p  = p0;
  if (input1) {
    p.x += static_cast<float>(input1[i0.x]);
    p.y += static_cast<float>(input1[i0.y]);
    p.z += static_cast<float>(input1[i0.z]);
  }
  for (int s = 0; s < _NumberOfSteps; s++) {
    float3 v = world2image(p);
    p.x += tex3D(_Input2_3D_x, v.x, v.y, v.z) * _h;
    p.y += tex3D(_Input2_3D_y, v.x, v.y, v.z) * _h;
    p.z += tex3D(_Input2_3D_z, v.x, v.y, v.z) * _h;
  }

  output[i0.x] = static_cast<VoxelType>(p.x - p0.x);
  output[i0.y] = static_cast<VoxelType>(p.y - p0.y);
  output[i0.z] = static_cast<VoxelType>(p.z - p0.z);
}

// ============================================================================
// Host function
// ============================================================================

// ----------------------------------------------------------------------------
// Forward Euler integration of vector field
template <class VoxelType>
IRTKCU_HOST_API void
Run(irtkCUGenericImage<VoxelType>       *output,
    const irtkCUGenericImage<VoxelType> *input1,
    const irtkCUGenericImage<VoxelType> *input2, int n, double dt)
{
  const int X = output->GetX();
  const int Y = output->GetY();
  const int Z = output->GetZ();
  const int T = output->GetT();

  const unsigned int N = static_cast<unsigned int>(X * Y * Z);

  // Check attributes of vector fields
  if ((Z == 1 && T != 2) || (Z >  1 && T != 3)) {
    cerr << "irtkCUVelocityToDisplacementFieldEuler::run: Output image is no valid 2D/3D vector field!" << endl;
    exit(1);
  }
  if (input1 && (input1->GetX() != X || input1->GetY() != Y || input1->GetZ() != Z || input1->GetT() != T)) {
    cerr << "irtkCUVelocityToDisplacementFieldEuler::Run: Dimensions of input1 and output vector fields do not match!" << endl;
    exit(1);
  }
  if (input2->GetX() != X || input2->GetY() != Y || input2->GetZ() != Z || input2->GetT() != T) {
    cerr << "irtkCUVelocityToDisplacementFieldEuler::Run: Dimensions of input2 and output vector fields do not match!" << endl;
    exit(1);
  }
  if (input1 && input1->GetPointerToDevice() == NULL) {
    cerr << "irtkCUVelocityToDisplacementFieldEuler::Run: Invalid input1 device pointer!" << endl;
    exit(1);
  }
  if (input2->GetPointerToDevice() == NULL) {
    cerr << "irtkCUVelocityToDisplacementFieldEuler::Run: Invalid input2 device pointer!" << endl;
    exit(1);
  }
  if (output->GetPointerToDevice() == NULL) {
    cerr << "irtkCUVelocityToDisplacementFieldEuler::Run: Invalid output device pointer!" << endl;
    exit(1);
  }

  // Copy attributes of images to constant device memory
  MemcpyImageAttributesToSymbols(output);

  // Copy parameters to constant device memory
  const float h = static_cast<float>(dt);
  CudaSafeCall( cudaMemcpyToSymbol(_NumberOfSteps, &n, sizeof(n)) );
  CudaSafeCall( cudaMemcpyToSymbol(_h,             &h, sizeof(h)) );

  // Copy components of second vector field to texture memory
  cudaArray *cuInput2[3];
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
  cudaMemcpy3DParms     copyParams  = {0};
  copyParams.extent                 = make_cudaExtent(X, Y, Z);
  copyParams.kind                   = (sizeof(VoxelType) == sizeof(float)) ? cudaMemcpyDeviceToDevice
                                                                           : cudaMemcpyHostToDevice;
  for (int c = 0; c < T; c++) {
    const float *fInput2 = ((copyParams.kind == cudaMemcpyDeviceToDevice)
                                 ? reinterpret_cast<const float *>(         input2->GetPointerToDevice(0, 0, 0, c))
                                 : reinterpret_cast<const float *>(to_float(input2->GetPointerToVoxels(0, 0, 0, c), N)));
    CudaSafeCall( cudaMalloc3DArray(&cuInput2[c], &channelDesc, copyParams.extent) );
    copyParams.srcPtr   = make_cudaPitchedPtr((void *)fInput2, X * sizeof(float), X, Y);
    copyParams.dstArray = cuInput2[c];
    CudaSafeCall( cudaMemcpy3D(&copyParams) );
    if (copyParams.kind == cudaMemcpyHostToDevice) delete[] fInput2;
  }

  // Bind texture to device memory
  if (T == 2) {
    _Input2_2D_x.addressMode[0] = _Input2_2D_y.addressMode[0] = cudaAddressModeClamp;
    _Input2_2D_x.addressMode[1] = _Input2_2D_y.addressMode[1] = cudaAddressModeClamp;
    _Input2_2D_x.filterMode     = _Input2_2D_y.filterMode     = cudaFilterModeLinear;
    _Input2_2D_x.normalized     = _Input2_2D_y.normalized     = false;
    CudaSafeCall( cudaBindTextureToArray( _Input2_2D_x, cuInput2[0], channelDesc) );
    CudaSafeCall( cudaBindTextureToArray( _Input2_2D_y, cuInput2[1], channelDesc) );
  } else {
    _Input2_3D_x.addressMode[0] = _Input2_3D_y.addressMode[0] = _Input2_3D_z.addressMode[0] = cudaAddressModeClamp;
    _Input2_3D_x.addressMode[1] = _Input2_3D_y.addressMode[1] = _Input2_3D_z.addressMode[1] = cudaAddressModeClamp;
    _Input2_3D_x.addressMode[2] = _Input2_3D_y.addressMode[2] = _Input2_3D_z.addressMode[2] = cudaAddressModeClamp;
    _Input2_3D_x.filterMode     = _Input2_3D_y.filterMode     = _Input2_3D_z.filterMode     = cudaFilterModeLinear;
    _Input2_3D_x.normalized     = _Input2_3D_y.normalized     = _Input2_3D_z.normalized     = false;
    CudaSafeCall( cudaBindTextureToArray( _Input2_3D_x, cuInput2[0], channelDesc) );
    CudaSafeCall( cudaBindTextureToArray( _Input2_3D_y, cuInput2[1], channelDesc) );
    CudaSafeCall( cudaBindTextureToArray( _Input2_3D_z, cuInput2[2], channelDesc) );
  }

  // Output buffer (if necessary)
  VoxelType *d_output = NULL;
  if (input2 == output || input1 == output) {
    CudaSafeCall( cudaMalloc(&d_output, output->GetNumberOfVoxels() * sizeof(VoxelType)) );
  } else {
    d_output = output->GetPointerToDevice();
  }

  // Execute kernel on device
  IRTKCU_START_TIMING();
  if (T == 2) {
    dim3 blocks, threads;
    DefaultKernelConfigurationND(blocks, threads, X, Y);
    cuRKE1_2D<VoxelType><<<blocks, threads>>>(d_output, input1 ? input1->GetPointerToDevice() : NULL);
    IRTKCU_DEBUG_INTERIM_TIMING(3, "irtkVelocityToDisplacementFieldEuler (cuRKE1_2D)");
  } else {
    dim3 blocks, threads;
    DefaultKernelConfigurationND(blocks, threads, X, Y, Z);
    cuRKE1_3D<VoxelType><<<blocks, threads>>>(d_output, input1 ? input1->GetPointerToDevice() : NULL);
    IRTKCU_DEBUG_INTERIM_TIMING(3, "irtkVelocityToDisplacementFieldEuler (cuRKE1_3D)");
  }
  CudaSyncCheckError();
  IRTKCU_DEBUG_TIMING(2, "irtkVelocityToDisplacementFieldEuler");

  // Copy output buffer
  if (d_output != output->GetPointerToDevice()) {
    CudaSafeCall( cudaMemcpy(output->GetPointerToDevice(), d_output, output->GetNumberOfVoxels() * sizeof(VoxelType), cudaMemcpyDeviceToDevice) );
    CudaSafeCall( cudaFree(d_output) );
    CudaSafeCall( cudaDeviceSynchronize() );
  }

  // Free memory
  if (T == 2) {
    CudaSafeCall( cudaUnbindTexture(_Input2_2D_x) );
    CudaSafeCall( cudaUnbindTexture(_Input2_2D_y) );
  } else {
    CudaSafeCall( cudaUnbindTexture(_Input2_3D_x) );
    CudaSafeCall( cudaUnbindTexture(_Input2_3D_y) );
    CudaSafeCall( cudaUnbindTexture(_Input2_3D_z) );
  }
  for (int c = 0; c < T; c++) CudaSafeCall( cudaFreeArray(cuInput2[c]) );
}

// ============================================================================
// Explicit template instantiations
// ============================================================================

#define INSTANTIATE_FOR_(T)                                                   \
    template IRTKCU_HOST_API void Run<T>(                                     \
         irtkCUGenericImage<T>*,                                              \
         const irtkCUGenericImage<T> *,                                       \
         const irtkCUGenericImage<T> *, int, double)

INSTANTIATE_FOR_(float);
INSTANTIATE_FOR_(double);


} // irtkCUVelocityToDisplacementFieldEuler
