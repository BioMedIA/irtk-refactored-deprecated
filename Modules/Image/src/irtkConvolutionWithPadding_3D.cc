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

template <class VoxelType> irtkConvolutionWithPadding_3D<VoxelType>::irtkConvolutionWithPadding_3D(VoxelType padding,
    bool Normalization) : irtkConvolution_3D<VoxelType>(Normalization)
{
  _padding = padding;
}

template <class VoxelType> void irtkConvolutionWithPadding_3D<VoxelType>::PutPaddingValue(VoxelType padding)
{
  _padding = padding;
}

template <class VoxelType> VoxelType irtkConvolutionWithPadding_3D<VoxelType>::GetPaddingValue(void)
{
  return _padding;
}

template <class VoxelType> double irtkConvolutionWithPadding_3D<VoxelType>::Run(int x, int y, int z, int t)
{
  irtkRealPixel *ptr2;
  irtkRealPixel val, sum;
  int x1, x2, y1, y2, z1, z2;

  if (this->_input->Get(x, y, z, t) <= this->_padding) return this->_padding;

  // Initialize
  val = 0;
  sum = 0;
  x1 = x - this->_input2->GetX()/2;
  x2 = x + this->_input2->GetX()/2;
  y1 = y - this->_input2->GetY()/2;
  y2 = y + this->_input2->GetY()/2;
  z1 = z - this->_input2->GetZ()/2;
  z2 = z + this->_input2->GetZ()/2;

  // Check if we use normalization
  if (this->_Normalization == true) {
    // Check whether boundary checking is necessary
    if ((x1 > 0) && (x2 < this->_input->GetX()) &&
        (y1 > 0) && (y2 < this->_input->GetY()) &&
        (z1 > 0) && (z2 < this->_input->GetZ())) {

      // If no, do fast convolution
      ptr2 = this->_input2->GetPointerToVoxels();
      for (z = z1; z <= z2; z++) {
        for (y = y1; y <= y2; y++) {
          for (x = x1; x <= x2; x++) {
            if (this->_input->Get(x, y, z, t) > this->_padding) {
              val += *ptr2 * this->_input->Get(x, y, z, t);
              sum += *ptr2;
            }
            ptr2++;
          }
        }
      }
    } else {

      // If yes, do slow convolution which handles boundaries
      ptr2 = this->_input2->GetPointerToVoxels();
      for (z = z1; z <= z2; z++) {
        for (y = y1; y <= y2; y++) {
          for (x = x1; x <= x2; x++) {
            if ((x >= 0) && (x < this->_input->GetX()) &&
                (y >= 0) && (y < this->_input->GetY()) &&
                (z >= 0) && (z < this->_input->GetZ())) {
              if (this->_input->Get(x, y, z, t) > this->_padding) {
                val += *ptr2 * this->_input->Get(x, y, z, t);
                sum += *ptr2;
              }
            }
            ptr2++;
          }
        }
      }
    }

    if (sum > 0) {
      return val / sum;
    } else {
      return 0;
    }
  } else {
    // Check whether boundary checking is necessary
    if ((x1 > 0) && (x2 < this->_input->GetX()) &&
        (y1 > 0) && (y2 < this->_input->GetY()) &&
        (z1 > 0) && (z2 < this->_input->GetZ())) {

      // If no, do fast convolution
      ptr2 = this->_input2->GetPointerToVoxels();
      for (z = z1; z <= z2; z++) {
        for (y = y1; y <= y2; y++) {
          for (x = x1; x <= x2; x++) {
            if (this->_input->Get(x, y, z, t) > this->_padding) {
              val += *ptr2 * this->_input->Get(x, y, z, t);
              sum += *ptr2;
            }
            ptr2++;
          }
        }
      }
    } else {

      // If yes, do slow convolution which handles boundaries
      ptr2 = this->_input2->GetPointerToVoxels();
      for (z = z1; z <= z2; z++) {
        for (y = y1; y <= y2; y++) {
          for (x = x1; x <= x2; x++) {
            if ((x >= 0) && (x < this->_input->GetX()) &&
                (y >= 0) && (y < this->_input->GetY()) &&
                (z >= 0) && (z < this->_input->GetZ())) {
              if (this->_input->Get(x, y, z, t) > this->_padding) {
                val += *ptr2 * this->_input->Get(x, y, z, t);
                sum += *ptr2;
              }
            }
            ptr2++;
          }
        }
      }
    }

    return val;
  }
}

template class irtkConvolutionWithPadding_3D<unsigned char>;
template class irtkConvolutionWithPadding_3D<short>;
template class irtkConvolutionWithPadding_3D<unsigned short>;
template class irtkConvolutionWithPadding_3D<float>;
template class irtkConvolutionWithPadding_3D<double>;
