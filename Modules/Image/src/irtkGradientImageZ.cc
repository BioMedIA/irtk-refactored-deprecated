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

#include <irtkGradientImage.h>

template <class VoxelType> double irtkGradientImageZ<VoxelType>::Run(int x, int y, int z, int t)
{

  double previous = this->_input->Get(x, y, z-1, t);
  double next     = this->_input->Get(x, y, z+1, t);

  double gradient = previous - next;

  return gradient;
}


template <class VoxelType> void irtkGradientImageZ<VoxelType>::Run()
{
  int x, y, z, t;

  // Do the initial set up
  this->Initialize();

  // Check image dimensions....
  if (this->_input->GetZ() < 2) {
    cerr<<" irtkGradientImageZ: Dimensions of input image are wrong"<<endl;
    exit(1);
  }

  for ( t = 0; t < this->_input->GetT(); ++t) {
    for ( x = 0; x < this->_input->GetX(); ++x) {
      for ( y = 0; y < this->_input->GetY(); ++y) {
        this->_output->Put(x, y, 0, t, 0);
        for ( z = 1; z < this->_input->GetZ()-1; ++z) {
          this->_output->PutAsDouble(x, y, z, t, this->Run(x, y, z, t));
        }
        this->_output->Put(x, y, z, t, 0);
      }
    }
  }

  // Do the final cleaning up
  this->Finalize();
}


template class  irtkGradientImageZ<irtkBytePixel>;
template class  irtkGradientImageZ<irtkGreyPixel>;
template class  irtkGradientImageZ<irtkRealPixel>;
