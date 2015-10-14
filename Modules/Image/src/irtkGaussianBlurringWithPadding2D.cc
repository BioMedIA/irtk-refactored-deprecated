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

#include <irtkGaussianBlurringWithPadding2D.h>


template <class VoxelType>
irtkGaussianBlurringWithPadding2D<VoxelType>
::irtkGaussianBlurringWithPadding2D(double sigma, VoxelType padding_value)
:
  irtkGaussianBlurringWithPadding<VoxelType>(sigma, sigma, .0, .0, padding_value)
{
}

template <class VoxelType>
irtkGaussianBlurringWithPadding2D<VoxelType>
::irtkGaussianBlurringWithPadding2D(double xsigma, double ysigma, VoxelType padding_value)
:
  irtkGaussianBlurringWithPadding<VoxelType>(xsigma, ysigma, .0, .0, padding_value)
{
}

template class irtkGaussianBlurringWithPadding2D<unsigned char>;
template class irtkGaussianBlurringWithPadding2D<short>;
template class irtkGaussianBlurringWithPadding2D<unsigned short>;
template class irtkGaussianBlurringWithPadding2D<float>;
template class irtkGaussianBlurringWithPadding2D<double>;
