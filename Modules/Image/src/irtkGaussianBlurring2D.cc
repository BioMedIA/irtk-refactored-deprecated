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
#include <irtkConvolution.h>
#include <irtkScalarFunctionToImage.h>

#include <irtkGaussianBlurring2D.h>


#ifdef USE_CUDA
# include <irtkCUImage.h>
#endif


template <class VoxelType> irtkGaussianBlurring2D<VoxelType>::irtkGaussianBlurring2D(double sigma)
:
  irtkGaussianBlurring<VoxelType>(sigma, sigma, .0, .0)
{
}

template <class VoxelType> irtkGaussianBlurring2D<VoxelType>::irtkGaussianBlurring2D(double xsigma, double ysigma)
:
  irtkGaussianBlurring<VoxelType>(xsigma, ysigma, .0, .0)
{
}

template <class VoxelType> irtkGaussianBlurring2D<VoxelType>::~irtkGaussianBlurring2D(void)
{
}


template class irtkGaussianBlurring2D<unsigned char>;
template class irtkGaussianBlurring2D<short>;
template class irtkGaussianBlurring2D<unsigned short>;
template class irtkGaussianBlurring2D<float>;
template class irtkGaussianBlurring2D<double>;

