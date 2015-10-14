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

#include <irtkGaussianBlurring.h>

#include <irtkImage.h>
#include <irtkConvolutionFunction.h>
#include <irtkScalarFunctionToImage.h>

#ifdef USE_CUDA
# include <irtkCUImage.h>
#endif

using namespace irtkConvolutionFunction;


template <class VoxelType> irtkGaussianBlurring<VoxelType>::irtkGaussianBlurring(double sigma)
:
  _SigmaX(sigma),
  _SigmaY(sigma),
  _SigmaZ(sigma),
  _SigmaT(.0),
  _Kernel(NULL)
{
}

template <class VoxelType> irtkGaussianBlurring<VoxelType>::irtkGaussianBlurring(double xsigma, double ysigma, double zsigma, double tsigma)
:
  _SigmaX(xsigma),
  _SigmaY(ysigma),
  _SigmaZ(zsigma),
  _SigmaT(tsigma),
  _Kernel(NULL)
{
}

template <class VoxelType> irtkGaussianBlurring<VoxelType>::~irtkGaussianBlurring()
{
  delete _Kernel;
}

template <class VoxelType> void irtkGaussianBlurring<VoxelType>::SetSigma(double sigma)
{
  _SigmaX = _SigmaY = _SigmaZ = sigma;
  _SigmaT = .0;
}

template <class VoxelType> void irtkGaussianBlurring<VoxelType>::SetSigma(double xsigma, double ysigma, double zsigma, double tsigma)
{
  _SigmaX = xsigma;
  _SigmaY = ysigma;
  _SigmaZ = zsigma;
  _SigmaT = tsigma;
}

template <class VoxelType> void irtkGaussianBlurring<VoxelType>::Initialize()
{
  // Initialize base class
  irtkImageToImage<VoxelType>::Initialize();

  // Instantiate convolution kernel of proper image type
#ifdef USE_CUDA
  if (use_gpu && dynamic_cast<irtkCUGenericImage<VoxelType> *>(this->_input ) != NULL &&
                 dynamic_cast<irtkCUGenericImage<VoxelType> *>(this->_output) != NULL) {
    _Kernel = new irtkCUGenericImage<irtkRealPixel>;
  } else
#endif
  {
    _Kernel = new irtkGenericImage<irtkRealPixel>;
  }
}

template <class VoxelType>
int irtkGaussianBlurring<VoxelType>::KernelSize(double sigma)
{
  return 2 * static_cast<int>(floor(3.0 * sigma)) + 1;
}

template <class VoxelType>
void irtkGaussianBlurring<VoxelType>::InitializeKernel(double sigma)
{
  // Create filter kernel for 1D Gaussian function
  _Kernel->Initialize(KernelSize(sigma), 1, 1);

  // Create scalar function which corresponds to a 1D Gaussian function
  irtkScalarGaussian gaussian(sigma,            .0, .0,  // stddev in x, y, z
                              _Kernel->X() / 2, .0, .0); // center in x, y, z

  // Convert from scalar function to filter kernel
  irtkRealPixel *p = _Kernel->Data();
  for (int i = 0; i < _Kernel->X(); ++i, ++p) {
    (*p) = gaussian.Evaluate(i);
  }
}

template <class VoxelType> void irtkGaussianBlurring<VoxelType>::Run()
{
  IRTK_START_TIMING();

  // Do the initial set up
  this->Initialize();

  irtkGenericImage<VoxelType> *input  = this->GetInput();
  irtkGenericImage<VoxelType> *output = this->GetOutput();

  const irtkImageAttributes &attr = input->Attributes();
  const int N = ((attr._dt == .0) ? attr._t : 1);

  // Blur along x axis
  if (_SigmaX != .0 && input->X() > 1) {
    if (output == this->GetInput()) output = new irtkGenericImage<VoxelType>(attr);
    this->InitializeKernel(_SigmaX / input->GetXSize());
    for (int n = 0; n < N; ++n) {
      ConvolveInX<irtkRealPixel> conv(input, _Kernel->Data(), _Kernel->X(), true, n);
      ParallelForEachVoxel(attr, input, output, conv);
    }
    swap(input, output);
  }

  // Blur along y axis
  if (_SigmaY != .0 && input->Y() > 1) {
    if (output == this->GetInput()) output = new irtkGenericImage<VoxelType>(attr);
    this->InitializeKernel(_SigmaY / input->GetYSize());
    for (int n = 0; n < N; ++n) {
      ConvolveInY<irtkRealPixel> conv(input, _Kernel->Data(), _Kernel->X(), true, n);
      ParallelForEachVoxel(attr, input, output, conv);
    }
    swap(input, output);
  }

  // Blur along z axis
  if (_SigmaZ != .0 && input->Z() > 1) {
    if (output == this->GetInput()) output = new irtkGenericImage<VoxelType>(attr);
    this->InitializeKernel(_SigmaZ / input->GetZSize());
    for (int n = 0; n < N; ++n) {
      ConvolveInZ<irtkRealPixel> conv(input, _Kernel->Data(), _Kernel->X(), true, n);
      ParallelForEachVoxel(attr, input, output, conv);
    }
    swap(input, output);
  }

  // Blur along t axis
  if (_SigmaT != .0 && input->T() > 1) {
    if (output == this->GetInput()) output = new irtkGenericImage<VoxelType>(attr);
    this->InitializeKernel(_SigmaT / input->GetTSize());
    ConvolveInT<irtkRealPixel> conv(input, _Kernel->Data(), _Kernel->X(), true);
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

template <class VoxelType> void irtkGaussianBlurring<VoxelType>::Finalize()
{
  // Finalize base class
  irtkImageToImage<VoxelType>::Finalize();

  // Free convolution kernel
  Delete(_Kernel);
}

template <class VoxelType> void irtkGaussianBlurring<VoxelType>::RunX()
{
  const double ysigma = _SigmaY;
  const double zsigma = _SigmaZ;
  const double tsigma = _SigmaT;

  _SigmaY = _SigmaZ = _SigmaT = .0;
  this->Run();

  _SigmaY = ysigma;
  _SigmaZ = zsigma;
  _SigmaT = tsigma;
}

template <class VoxelType> void irtkGaussianBlurring<VoxelType>::RunY()
{
  const double xsigma = _SigmaX;
  const double zsigma = _SigmaZ;
  const double tsigma = _SigmaT;

  _SigmaX = _SigmaZ = _SigmaT = .0;
  this->Run();

  _SigmaX = xsigma;
  _SigmaZ = zsigma;
  _SigmaT = tsigma;
}

template <class VoxelType> void irtkGaussianBlurring<VoxelType>::RunZ()
{
  const double xsigma = _SigmaX;
  const double ysigma = _SigmaY;
  const double tsigma = _SigmaT;

  _SigmaX = _SigmaY = _SigmaT = .0;
  this->Run();

  _SigmaX = xsigma;
  _SigmaY = ysigma;
  _SigmaT = tsigma;
}

template <class VoxelType> void irtkGaussianBlurring<VoxelType>::RunT()
{
  const double xsigma = _SigmaX;
  const double ysigma = _SigmaY;
  const double zsigma = _SigmaZ;
  const double tsigma = _SigmaT;

  _SigmaX = _SigmaY = _SigmaZ = .0;
  if (_SigmaT == .0) _SigmaT = _SigmaX;
  this->Run();

  _SigmaX = xsigma;
  _SigmaY = ysigma;
  _SigmaZ = zsigma;
  _SigmaT = tsigma;
}

template class irtkGaussianBlurring<unsigned char>;
template class irtkGaussianBlurring<short>;
template class irtkGaussianBlurring<unsigned short>;
template class irtkGaussianBlurring<float>;
template class irtkGaussianBlurring<double>;

