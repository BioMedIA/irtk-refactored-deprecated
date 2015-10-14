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

template <class VoxelType> double irtkGradientImageY<VoxelType>::Run(int x, int y, int z, int t)
{

  double previous = this->_input->Get(x, y-1, z, t);
  double next     = this->_input->Get(x, y+1, z, t);

  double gradient = previous - next;

  return gradient;
}


template <class VoxelType> void irtkGradientImageY<VoxelType>::Run()
{
  int x, y, z, t;

  // Do the initial set up
  this->Initialize();

  // Check image dimensions....
  if (this->_input->GetY() < 2) {
    cerr<<" irtkGradientImageY: Dimensions of input image are wrong"<<endl;
    exit(1);
  }

  for ( t = 0; t < this->_input->GetT(); ++t) {
    for ( z = 0; z < this->_input->GetZ(); ++z) {
      for ( x = 0; x < this->_input->GetX(); ++x) {
        this->_output->Put(x, 0, z, t, 0);
        for ( y = 1; y < this->_input->GetY()-1; ++y) {
          this->_output->PutAsDouble(x, y, z, t, this->Run(x, y, z, t));
        }
        this->_output->Put(x, y, z, t, 0);
      }
    }
  }

  // Do the final cleaning up
  this->Finalize();

}


template class  irtkGradientImageY<irtkBytePixel>;
template class  irtkGradientImageY<irtkGreyPixel>;
template class  irtkGradientImageY<irtkRealPixel>;
