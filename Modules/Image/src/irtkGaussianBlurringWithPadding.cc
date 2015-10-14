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

#include <irtkGaussianBlurringWithPadding.h>

#include <irtkConvolutionFunction.h>
using namespace irtkConvolutionFunction;


template <class VoxelType> irtkGaussianBlurringWithPadding<VoxelType>
::irtkGaussianBlurringWithPadding(double sigma, VoxelType padding_value)
:
  irtkGaussianBlurring<VoxelType>(sigma),
  _PaddingValue(padding_value)
{
}

template <class VoxelType> irtkGaussianBlurringWithPadding<VoxelType>
::irtkGaussianBlurringWithPadding(double xsigma, double ysigma, VoxelType padding_value)
:
  irtkGaussianBlurring<VoxelType>(xsigma, ysigma, .0, .0),
  _PaddingValue(padding_value)
{
}

template <class VoxelType> irtkGaussianBlurringWithPadding<VoxelType>
::irtkGaussianBlurringWithPadding(double xsigma, double ysigma, double zsigma, VoxelType padding_value)
:
  irtkGaussianBlurring<VoxelType>(xsigma, ysigma, zsigma, .0),
  _PaddingValue(padding_value)
{
}

template <class VoxelType> irtkGaussianBlurringWithPadding<VoxelType>
::irtkGaussianBlurringWithPadding(double xsigma, double ysigma, double zsigma, double tsigma, VoxelType padding_value)
:
  irtkGaussianBlurring<VoxelType>(xsigma, ysigma, zsigma, tsigma),
  _PaddingValue(padding_value)
{
}

template <class VoxelType> void irtkGaussianBlurringWithPadding<VoxelType>::Run()
{
  IRTK_START_TIMING();

  // Do the initial set up
  this->Initialize();

  irtkGenericImage<VoxelType> *input  = this->GetInput();
  irtkGenericImage<VoxelType> *output = this->GetOutput();

  const irtkImageAttributes &attr = input->Attributes();
  const int N = ((attr._dt == .0) ? attr._t : 1);

  // Blur along x axis
  if (this->_SigmaX != .0 && input->X() > 1) {
    if (output == this->GetInput()) output = new irtkGenericImage<VoxelType>(attr);
    this->InitializeKernel(this->_SigmaX / input->GetXSize());
    typedef ConvolveTruncatedForegroundInX<irtkRealPixel> Conv;
    for (int n = 0; n < N; ++n) {
      Conv conv(input, _PaddingValue, this->_Kernel->Data(), this->_Kernel->X(), true, n);
      ParallelForEachVoxel(attr, input, output, conv);
    }
    swap(input, output);
  }

  // Blur along y axis
  if (this->_SigmaY != .0 && input->Y() > 1) {
    if (output == this->GetInput()) output = new irtkGenericImage<VoxelType>(attr);
    this->InitializeKernel(this->_SigmaY / input->GetYSize());
    typedef ConvolveTruncatedForegroundInY<irtkRealPixel> Conv;
    for (int n = 0; n < N; ++n) {
      Conv conv(input, _PaddingValue, this->_Kernel->Data(), this->_Kernel->X(), true, n);
      ParallelForEachVoxel(attr, input, output, conv);
    }
    swap(input, output);
  }

  // Blur along z axis
  if (this->_SigmaZ != .0 && input->Z() > 1) {
    if (output == this->GetInput()) output = new irtkGenericImage<VoxelType>(attr);
    this->InitializeKernel(this->_SigmaZ / input->GetZSize());
    typedef ConvolveTruncatedForegroundInZ<irtkRealPixel> Conv;
    for (int n = 0; n < N; ++n) {
      Conv conv(input, _PaddingValue, this->_Kernel->Data(), this->_Kernel->X(), true, n);
      ParallelForEachVoxel(attr, input, output, conv);
    }
    swap(input, output);
  }

  // Blur along t axis
  if (this->_SigmaT != .0 && input->T() > 1) {
    if (output == this->GetInput()) output = new irtkGenericImage<VoxelType>(attr);
    this->InitializeKernel(this->_SigmaT / input->GetTSize());
    typedef ConvolveTruncatedForegroundInT<irtkRealPixel> Conv;
    Conv conv(input, _PaddingValue, this->_Kernel->Data(), this->_Kernel->X(), true);
    ParallelForEachVoxel(attr, input, output, conv);
    swap(input, output);
  }

  // Copy result if last output (i.e., after swap input pointer) was not filter
  // output image and delete possibly additionally allocated image
  if (input != output) {
    if (input == this->GetOutput()) {
      if (output != this->GetInput()) Delete(output);
    } else {
      this->GetOutput()->CopyFrom(input->Data());
      if (input != this->GetInput ()) Delete(input);
    }
  }

  // Do the final cleaning up
  this->Finalize();

  IRTK_DEBUG_TIMING(5, this->NameOfClass());
}

template class irtkGaussianBlurringWithPadding<unsigned char>;
template class irtkGaussianBlurringWithPadding<short>;
template class irtkGaussianBlurringWithPadding<unsigned short>;
template class irtkGaussianBlurringWithPadding<float>;
template class irtkGaussianBlurringWithPadding<double>;
