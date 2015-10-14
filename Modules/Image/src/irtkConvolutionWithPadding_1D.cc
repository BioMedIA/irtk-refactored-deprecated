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

template <class VoxelType> irtkConvolutionWithPadding_1D<VoxelType>::irtkConvolutionWithPadding_1D(VoxelType padding,
    bool Normalization) : irtkConvolution_1D<VoxelType>(Normalization)
{
  _padding = padding;
}

template <class VoxelType> void irtkConvolutionWithPadding_1D<VoxelType>::PutPaddingValue(VoxelType padding)
{
  _padding = padding;
}

template <class VoxelType> VoxelType irtkConvolutionWithPadding_1D<VoxelType>::GetPaddingValue(void)
{
  return _padding;
}

template <class VoxelType> double irtkConvolutionWithPadding_1D<VoxelType>::Run(int x, int y, int z, int t)
{
  int x1, x2;
  irtkRealPixel *ptr2;
  VoxelType *ptr;
  double val, sum;

  if (this->_input->Get(x, y, z, t) <= this->_padding) return this->_padding;

  // Initialize
  val = 0;
  sum = 0;
  x1 = x - this->_input2->GetX()/2;
  x2 = x + this->_input2->GetX()/2;

  // Check if we use normalization
  if (this->_Normalization == true) {
    // Check whether boundary checking is necessary
    if ((x1 > 0) && (x2 < this->_input->GetX())) {

      // If no, do fast Convolution
      ptr  = this->_input->GetPointerToVoxels(x1, y, z, t);
      ptr2 = this->_input2->GetPointerToVoxels();
      for (x = x1; x <= x2; x++) {
        if (*ptr > this->_padding) {
          val += *ptr2 * *ptr;
          sum += *ptr2;
        }
        ptr++;
        ptr2++;
      }
    } else {

      // If yes, do slow convolution which handles boundaries
      ptr2 = this->_input2->GetPointerToVoxels();
      for (x = x1; x <= x2; x++) {
        if ((x >= 0) && (x < this->_input->GetX())) {
          if (this->_input->Get(x, y, z, t) > this->_padding) {
            val += *ptr2 * this->_input->Get(x, y, z, t);
            sum += *ptr2;
          }
        }
        ptr2++;
      }
    }

    if (sum != 0) {
      return val / sum;
    } else {
      return this->_padding;
    }
  } else {
    // Check whether boundary checking is necessary
    if ((x1 > 0) && (x2 < this->_input->GetX())) {

      // If no, do fast Convolution
      ptr  = this->_input->GetPointerToVoxels(x1, y, z, t);
      ptr2 = this->_input2->GetPointerToVoxels();
      for (x = x1; x <= x2; x++) {
        if (*ptr > this->_padding) {
          val += *ptr2 * *ptr;
        }
        ptr++;
        ptr2++;
      }
    } else {

      // If yes, do slow convolution which handles boundaries
      ptr2 = this->_input2->GetPointerToVoxels();
      for (x = x1; x <= x2; x++) {
        if ((x >= 0) && (x < this->_input->GetX())) {
          if (this->_input->Get(x, y, z, t) > this->_padding) {
            val += *ptr2 * this->_input->Get(x, y, z, t);
          }
        }
        ptr2++;
      }
    }

    return val;
  }
}

template class irtkConvolutionWithPadding_1D<unsigned char>;
template class irtkConvolutionWithPadding_1D<short>;
template class irtkConvolutionWithPadding_1D<unsigned short>;
template class irtkConvolutionWithPadding_1D<float>;
template class irtkConvolutionWithPadding_1D<double>;
