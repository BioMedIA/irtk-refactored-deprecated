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

#include <irtkNoise.h>

template <class VoxelType> irtkGaussianNoiseWithPadding<VoxelType>::irtkGaussianNoiseWithPadding() : irtkGaussianNoise<VoxelType>()
{
  _PaddingValue = std::numeric_limits<VoxelType>::min();
}

template <class VoxelType> irtkGaussianNoiseWithPadding<VoxelType>::irtkGaussianNoiseWithPadding(double Mean, double Sigma, VoxelType PaddingValue, VoxelType MinVal, VoxelType MaxVal) : irtkGaussianNoise<VoxelType>(Mean, Sigma, MinVal, MaxVal)
{
  _PaddingValue = PaddingValue;
}

template <class VoxelType> double irtkGaussianNoiseWithPadding<VoxelType>::Run(int x, int y, int z, int t)
{
  if (this->_input->Get(x, y, z, t) > this->_PaddingValue) {
    return this->irtkGaussianNoise<VoxelType>::Run(x, y, z, t);
  } else {
    return this->_PaddingValue;
  }
}

template class irtkGaussianNoiseWithPadding<irtkBytePixel>;
template class irtkGaussianNoiseWithPadding<irtkGreyPixel>;
template class irtkGaussianNoiseWithPadding<irtkRealPixel>;
