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


#ifdef USE_CUDA
#include <irtkCUImage.h>
namespace irtkCUConvolution_1D {
  template <class VoxelType, class KernelType>
  void Run(irtkCUGenericImage<VoxelType>        *,
           const irtkCUGenericImage<VoxelType>  *,
           const irtkCUGenericImage<KernelType> *, bool);
}
#endif


template <class VoxelType> irtkConvolution_1D<VoxelType>::irtkConvolution_1D(bool Normalization) :
    irtkConvolution<VoxelType>(Normalization)
{}

template <class VoxelType> void irtkConvolution_1D<VoxelType>::SetInput2(irtkGenericImage<irtkRealPixel> *image)
{
  if (image != NULL) {
    _input2 = image;
  } else {
    cerr << this->NameOfClass() << "::SetInput2: Input is not an image\n";
    exit(1);
  }
}


template <class VoxelType> double irtkConvolution_1D<VoxelType>::Run(int x, int y, int z, int t)
{
  int x1, x2;
  irtkRealPixel *ptr2;
  VoxelType *ptr;
  irtkRealPixel val, sum;

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
        val += *ptr2 * *ptr;
        sum += *ptr2;
        ptr++;
        ptr2++;
      }

    } else {

      // If yes, do slow convolution which handles boundaries
      ptr2 = this->_input2->GetPointerToVoxels();
      for (x = x1; x <= x2; x++) {
        if ((x >= 0) && (x < this->_input->GetX())) {
          val += *ptr2 * this->_input->Get(x, y, z, t);
          sum += *ptr2;
        }
        ptr2++;
      }
    }

    //  Normalize filter value by sum of filter elements
    if (sum > 0) {
      return val / sum;
    } else {
      return 0;
    }

  } else {

    // Check whether boundary checking is necessary
    if ((x1 > 0) && (x2 < this->_input->GetX())) {

      // If no, do fast Convolution
      ptr  = this->_input->GetPointerToVoxels(x1, y, z, t);
      ptr2 = this->_input2->GetPointerToVoxels();
      for (x = x1; x <= x2; x++) {
        val += *ptr2 * *ptr;
        ptr++;
        ptr2++;
      }

    } else {

      // If yes, do slow convolution which handles boundaries
      ptr2 = this->_input2->GetPointerToVoxels();
      for (x = x1; x <= x2; x++) {
        if ((x >= 0) && (x < this->_input->GetX())) {
          val += *ptr2 * this->_input->Get(x, y, z, t);
        }
        ptr2++;
      }
    }
    return val;
  }
}

template <class VoxelType> void irtkConvolution_1D<VoxelType>::Initialize()
{
  // Check kernel
  if (this->_input2 == NULL) {
    cerr << this->NameOfClass() << "::Run: Filter has no second input!" << endl;
    exit(1);
  }

  // Check kernel size
  if (this->_input2->GetY() != 1 || this->_input2->GetZ() != 1 || this->_input2->GetT() != 1) {
    cerr << this->NameOfClass() << "::Run: Filter dimensions should be 1 in y, z, and t!" << endl;
    exit(1);
  }

  // Do the initial set up
  this->irtkImageToImage<VoxelType>::Initialize();
}

template <class VoxelType> void irtkConvolution_1D<VoxelType>::Run()
{
#ifdef USE_CUDA
  irtkCUGenericImage<VoxelType>     *output = dynamic_cast<irtkCUGenericImage<VoxelType>     *>(this->_output);
  irtkCUGenericImage<VoxelType>     *input  = dynamic_cast<irtkCUGenericImage<VoxelType>     *>(this->_input);
  irtkCUGenericImage<irtkRealPixel> *kernel = dynamic_cast<irtkCUGenericImage<irtkRealPixel> *>(this->_input2);

  if (use_gpu && input && kernel && output) {
    irtkCUConvolution_1D::Run<VoxelType, irtkRealPixel>(output, input, kernel, this->_Normalization);
  } else
#endif
  {
    IRTK_START_TIMING();
    irtkConvolution<VoxelType>::Run();
    IRTK_DEBUG_TIMING(5, "irtkConvolution1D");
  }
}

template class irtkConvolution_1D<unsigned char>;
template class irtkConvolution_1D<short>;
template class irtkConvolution_1D<unsigned short>;
template class irtkConvolution_1D<float>;
template class irtkConvolution_1D<double>;
