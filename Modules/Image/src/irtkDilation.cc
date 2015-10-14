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

#include <irtkDilation.h>

template <class VoxelType> irtkDilation<VoxelType>::irtkDilation()
{
	// Default connectivity.
	this->_Connectivity = CONNECTIVITY_26;
}

template <class VoxelType> irtkDilation<VoxelType>::~irtkDilation(void)
{
}

template <class VoxelType> void irtkDilation<VoxelType>::Initialize()
{
  // Do the initial set up
  this->irtkImageToImage<VoxelType>::Initialize();

  this->_offsets.Initialize(this->_input, this->_Connectivity);
}

template <class VoxelType> void irtkDilation<VoxelType>::Run()
{
  int i, x, y, z, t, maskSize;
  VoxelType value;
  VoxelType *ptr2current, *ptr2offset;

  // Do the initial set up
  this->Initialize();

  maskSize = this->_offsets.GetSize();

  for (t = 0; t < this->_input->GetT(); t++) {
    for (z = 0; z < this->_input->GetZ(); z++) {
      for (y = 0; y < this->_input->GetY(); y++) {
        for (x = 0; x < this->_input->GetX(); x++) {
          if ((x == 0) || (x == this->_input->GetX()-1) ||
              (y == 0) || (y == this->_input->GetY()-1) ||
              (z == 0) || (z == this->_input->GetZ()-1)) {
            this->_output->Put(x, y, z, t, this->_input->Get(x, y, z, t));
          } else {
            value = this->_input->Get(x, y, z, t);
          	ptr2current = this->_input->GetPointerToVoxels(x, y, z, t);
            for (i = 0; i < maskSize; ++i) {
          		ptr2offset = ptr2current + this->_offsets(i);
            	if (*ptr2offset > value)
            		value = *ptr2offset;
            }
            this->_output->Put(x, y, z, t, value);
          }
        }
      }
    }
  }

  // Do the final cleaning up
  this->Finalize();
}

template class irtkDilation<irtkBytePixel>;
template class irtkDilation<irtkGreyPixel>;
template class irtkDilation<irtkRealPixel>;
