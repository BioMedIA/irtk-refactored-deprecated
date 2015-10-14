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

#include <irtkImage.h>
#include <irtkDownsampling.h>
#include <irtkUnaryVoxelFunction.h>
#include <irtkConvolutionFunction.h>
#include <irtkScalarFunctionToImage.h>

using namespace irtkUnaryVoxelFunction;
using namespace irtkConvolutionFunction;


// ---------------------------------------------------------------------------
template <class VoxelType>
irtkDownsampling<VoxelType>::irtkDownsampling(int m)
:
  _DownsampleFactorX(m),
  _DownsampleFactorY(m),
  _DownsampleFactorZ(m),
  _Kernel           (NULL),
  _KernelSize       (0),
  _NormalizeKernel  (true)
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
irtkDownsampling<VoxelType>::irtkDownsampling(int mx, int my, int mz)
:
  _DownsampleFactorX(mx),
  _DownsampleFactorY(my),
  _DownsampleFactorZ(mz),
  _Kernel           (NULL),
  _KernelSize       (0),
  _NormalizeKernel  (true)
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void irtkDownsampling<VoxelType>::Kernel(irtkScalarFunction *kernel, int r, bool normalize)
{
  _Kernel          = kernel;
  _KernelSize      = 2 * r + 1;
  _NormalizeKernel = normalize;
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void irtkDownsampling<VoxelType>::Initialize()
{
  // Initialize base class, but not output image except of creating a
  // temporary image instance if set output equals input image instance
  irtkImageToImage<VoxelType>::Initialize(false);
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void irtkDownsampling<VoxelType>::DiscreteKernel(irtkGenericImage<irtkRealPixel> &kernel, double dx)
{
  // Initialize discrete kernel
  kernel.Initialize  (_KernelSize, 1,   1  );
  kernel.PutPixelSize(         dx, 1.0, 1.0);
  // Sample continuous kernel function
  irtkScalarFunctionToImage<irtkRealPixel> generator;
  generator.SetInput (_Kernel);
  generator.SetOutput(&kernel);
  generator.Run();
  // Normalize discrete kernel weights
  if (_NormalizeKernel) {
    irtkRealPixel *v, sum = .0;
    v = kernel.GetPointerToVoxels();
    for (int i = 0; i < _KernelSize; ++i) sum += (*v++);
    v = kernel.GetPointerToVoxels();
    for (int i = 0; i < _KernelSize; ++i) (*v++) /= sum;
  }
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void irtkDownsampling<VoxelType>::Run()
{
  IRTK_START_TIMING();

  // Do the initial set up
  this->Initialize();

  irtkGenericImage<irtkRealPixel> kernel;
  irtkGenericImage<VoxelType>    &input  = *this->GetInput();
  irtkGenericImage<VoxelType>    &output = *this->GetOutput();
  irtkGenericImage<VoxelType>     image;

  if (input.HasBackgroundValue()) {
    image .PutBackgroundValueAsDouble(input.GetBackgroundValueAsDouble());
    output.PutBackgroundValueAsDouble(input.GetBackgroundValueAsDouble());
  }

  // Initial output image attributes
  irtkImageAttributes attr = input.GetImageAttributes();

  // Downsample x dimension
  attr._x  /= _DownsampleFactorX;
  attr._dx *= _DownsampleFactorX;
  output.Initialize(attr);

  if (_Kernel) {
    DiscreteKernel(kernel, input.GetXSize());
    DownsampleConvolvedMirroredForegroundInX<VoxelType, irtkRealPixel> downsample(&input, _DownsampleFactorX, &kernel);
    ParallelForEachVoxel(attr, output, downsample);
  } else {
    DownsampleX<VoxelType> downsample(&input, _DownsampleFactorX);
    ParallelForEachVoxel(attr, output, downsample);
  }

  // Downsample y dimension
  attr._y  /= _DownsampleFactorY;
  attr._dy *= _DownsampleFactorY;
  image.Initialize(attr);

  if (_Kernel) {
    DiscreteKernel(kernel, input.GetYSize());
    DownsampleConvolvedMirroredForegroundInY<VoxelType, irtkRealPixel> downsample(&output, _DownsampleFactorY, &kernel);
    ParallelForEachVoxel(attr, image, downsample);
  } else {
    DownsampleY<VoxelType> downsample(&output, _DownsampleFactorY);
    ParallelForEachVoxel(attr, image, downsample);
  }

  // Downsample z dimension
  attr._z  /= _DownsampleFactorZ;
  attr._dz *= _DownsampleFactorZ;
  output.Initialize(attr);

  if (_Kernel) {
    DiscreteKernel(kernel, input.GetZSize());
    DownsampleConvolvedMirroredForegroundInZ<VoxelType, irtkRealPixel> downsample(&image, _DownsampleFactorZ, &kernel);
    ParallelForEachVoxel(attr, output, downsample);
  } else {
    DownsampleZ<VoxelType> downsample(&image, _DownsampleFactorZ);
    ParallelForEachVoxel(attr, output, downsample);
  }

  // Adjust origin if necessary
  //
  // Note: The first voxel of the output is offset such that the downsampled image
  //       region is centered around the origin of the input image as best as
  //       possible already. Only if the left and right voxel margins are not
  //       identical, the origin of the output must be adjusted by half an input voxel.
  double x = (input.GetX() - 1 - (input.GetX() % _DownsampleFactorX) % 2) / 2.0;
  double y = (input.GetY() - 1 - (input.GetY() % _DownsampleFactorY) % 2) / 2.0;
  double z = (input.GetZ() - 1 - (input.GetZ() % _DownsampleFactorZ) % 2) / 2.0;
  input .ImageToWorld(x, y, z);
  output.PutOrigin   (x, y, z);

  // Do the final cleaning up
  this->Finalize();

  IRTK_DEBUG_TIMING(3, "irtkDownsampling::Run");
}

// ---------------------------------------------------------------------------
// Explicit template instantiations
template class irtkDownsampling<char>;
template class irtkDownsampling<unsigned char>;
template class irtkDownsampling<short>;
template class irtkDownsampling<unsigned short>;
template class irtkDownsampling<int>;
template class irtkDownsampling<float>;
template class irtkDownsampling<double>;
