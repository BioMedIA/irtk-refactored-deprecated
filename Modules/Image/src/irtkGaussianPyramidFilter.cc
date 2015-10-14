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
#include <irtkGaussianPyramidFilter.h>
#include <irtkUnaryVoxelFunction.h>
#include <irtkConvolutionFunction.h>

using namespace irtkUnaryVoxelFunction;
using namespace irtkConvolutionFunction;


// ---------------------------------------------------------------------------
template <class VoxelType>
irtkGaussianPyramidFilter<VoxelType>::irtkGaussianPyramidFilter(int j)
:
  irtkImageToImage<VoxelType>(),
  _InputLevel (0),
  _OutputLevel(j)
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
irtkGaussianPyramidFilter<VoxelType>::irtkGaussianPyramidFilter(int i, int j)
:
  irtkImageToImage<VoxelType>(),
  _InputLevel (i),
  _OutputLevel(j)
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void irtkGaussianPyramidFilter<VoxelType>::Initialize()
{
  // Initialize base class, but not output image except of creating a
  // temporary image instance if set output equals input image instance
  irtkImageToImage<VoxelType>::Initialize(false);
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void irtkGaussianPyramidFilter<VoxelType>::Downsample()
{
  irtkGenericImage<VoxelType> *input  = this->GetInput();
  irtkGenericImage<VoxelType> *temp   = this->GetOutput();
  irtkGenericImage<VoxelType> *output = new irtkGenericImage<VoxelType>();

  if (input->HasBackgroundValue()) {
    temp  ->PutBackgroundValueAsDouble(input->GetBackgroundValueAsDouble());
    output->PutBackgroundValueAsDouble(input->GetBackgroundValueAsDouble());
  }

  // Type of downsampling voxel function
  typedef DownsampleConvolvedExtendedForegroundInX<VoxelType, irtkRealPixel> DownsampleInX;
  typedef DownsampleConvolvedExtendedForegroundInY<VoxelType, irtkRealPixel> DownsampleInY;
  typedef DownsampleConvolvedExtendedForegroundInZ<VoxelType, irtkRealPixel> DownsampleInZ;

  // FIXME: The irtkGaussianPyramidFilter shifts the image origin incorrectly!
  cerr << "WARNING: The irtkGaussianPyramidFilter shifts the image origin" << endl;

  // Smoothing kernel used for downsampling
  const int     size         = 5;
  irtkRealPixel kernel[size] = {1, 4, 6, 4, 1};
  const double  norm         = 1.0 / 16.0;

  for (int i = _InputLevel; i < _OutputLevel; ++i) {
    // Set output of previous iteration as input
    if (i != _InputLevel) {
      input  = output;
      output = new irtkGenericImage<VoxelType>();
    }

    // Input image attributes
    irtkImageAttributes attr = input->GetImageAttributes();
    const int n = attr._dt ? 1 : attr._t;
    attr._dt = 1.0; // s.t. ParallelForEachVoxel downsamples each component

    // Downsample x dimension
    if (attr._x > 1) {
      DownsampleInX downsampleX(input, 2, kernel, size, norm);
      attr._x  /= 2;
      attr._dx *= 2;
      attr._xorigin += attr._xaxis[0] * downsampleX._Offset * attr._dx;
      attr._yorigin += attr._yaxis[0] * downsampleX._Offset * attr._dx;
      attr._zorigin += attr._zaxis[0] * downsampleX._Offset * attr._dx;
      output->Initialize(attr, n);
      ParallelForEachVoxel(attr, output, downsampleX);
    }

    // Downsample y dimension
    if (attr._y > 1) {
      DownsampleInY downsampleY(output, 2, kernel, size, norm);
      attr._y  /= 2;
      attr._dy *= 2;
      attr._xorigin += attr._xaxis[1] * downsampleY._Offset * attr._dy;
      attr._yorigin += attr._yaxis[1] * downsampleY._Offset * attr._dy;
      attr._zorigin += attr._zaxis[1] * downsampleY._Offset * attr._dy;
      temp->Initialize(attr, n);
      ParallelForEachVoxel(attr, temp, downsampleY);
    }

    // Downsample z dimension
    if (attr._z > 1) {
      DownsampleInZ downsampleZ(temp, 2, kernel, size, norm);
      attr._z  /= 2;
      attr._dz *= 2;
      attr._xorigin += attr._xaxis[2] * downsampleZ._Offset * attr._dz;
      attr._yorigin += attr._yaxis[2] * downsampleZ._Offset * attr._dz;
      attr._zorigin += attr._zaxis[2] * downsampleZ._Offset * attr._dz;
      output->Initialize(attr, n);
      ParallelForEachVoxel(attr, output, downsampleZ);
    }

    // Delete intermediate input
    if (input != this->GetInput()) {
      delete input;
      input = NULL;
    }
  }

  // Copy final output to actual output
  *temp = *output;
  delete output;
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void irtkGaussianPyramidFilter<VoxelType>::Upsample()
{
  irtkGenericImage<VoxelType> &input  = *this->GetInput();
  irtkGenericImage<VoxelType> &output = *this->GetOutput();
  irtkGenericImage<VoxelType>  temp;

  // Input image attributes
  irtkImageAttributes attr = input.GetImageAttributes();

  attr._x *= 2; attr._dx /= 2.0;
  attr._y *= 2; attr._dy /= 2.0;
  attr._z *= 2; attr._dz /= 2.0;
  temp  .Initialize(attr);
  output.Initialize(attr);

  cerr << "irtkGaussianPyramidFilter::Upsample: Not implemented" << endl;
  exit(1);
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void irtkGaussianPyramidFilter<VoxelType>::Run()
{
  IRTK_START_TIMING();

  // Do the initial set up
  this->Initialize();

  // Down-/upsample image
  if (_InputLevel == _OutputLevel) {
    this->GetOutput()->Initialize(this->GetInput()->GetImageAttributes());
    this->GetOutput()->CopyFrom  (this->GetInput()->GetPointerToVoxels());
  } else if (_InputLevel < _OutputLevel) {
    this->Downsample();
  } else {
    this->Upsample();
  }

  // Do the final cleaning up
  this->Finalize();

  IRTK_DEBUG_TIMING(5, "irtkGaussianPyramidFilter::Run");
}

// ---------------------------------------------------------------------------
// Explicit template instantiations
template class irtkGaussianPyramidFilter<char>;
template class irtkGaussianPyramidFilter<unsigned char>;
template class irtkGaussianPyramidFilter<short>;
template class irtkGaussianPyramidFilter<unsigned short>;
template class irtkGaussianPyramidFilter<int>;
template class irtkGaussianPyramidFilter<unsigned int>;
template class irtkGaussianPyramidFilter<float>;
template class irtkGaussianPyramidFilter<double>;
