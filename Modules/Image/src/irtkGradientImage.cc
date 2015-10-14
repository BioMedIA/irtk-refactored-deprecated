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

template <class VoxelType> irtkGradientImage<VoxelType>::irtkGradientImage()
{
  _Padding = MIN_GREY;
}

template <class VoxelType> double irtkGradientImage<VoxelType>::Run(int x, int y, int z, int t)
{
  double dx, dy, dz;

  if ((x > 0) && (x < this->_input->GetX()-1)
      && this->_input->Get(x-1, y, z, t) > _Padding
      && this->_input->Get(x+1, y, z, t) > _Padding
     ) {
    dx = this->_input->Get(x-1, y, z, t) - this->_input->Get(x+1, y, z, t);
  } else {
    dx = 0;
  }

  if ((y > 0) && (y < this->_input->GetY()-1)
      && this->_input->Get(x, y-1, z, t) > _Padding
      && this->_input->Get(x, y+1, z, t) > _Padding
     ) {
    dy = this->_input->Get(x, y-1, z, t) - this->_input->Get(x, y+1, z, t);
  } else {
    dy = 0;
  }

  if ((z > 0) && (z < this->_input->GetZ()-1)
      && this->_input->Get(x, y, z-1, t) > _Padding
      && this->_input->Get(x, y, z+1, t) > _Padding
     ) {
    dz = this->_input->Get(x, y, z-1, t) - this->_input->Get(x, y, z+1, t);
  } else {
    dz = 0;
  }

  return sqrt(dx*dx + dy*dy + dz*dz);
}

template <class VoxelType> void irtkGradientImage<VoxelType>::Run()
{
  int x, y, z, t;

  // Do the initial set up
  this->Initialize();

  for (t = 0; t < this->_input->GetT(); ++t) {
    for (z = 0; z < this->_input->GetZ(); ++z) {
      for (y = 0; y < this->_input->GetY(); ++y) {
        for (x = 0; x < this->_input->GetX(); ++x) {
          this->_output->PutAsDouble(x, y, z, t, this->Run(x, y, z, t));
        }
      }
    }
  }
  // Do the final cleaning up
  this->Finalize();
}

template class irtkGradientImage<unsigned char>;
template class irtkGradientImage<short>;
template class irtkGradientImage<unsigned short>;
template class irtkGradientImage<float>;
template class irtkGradientImage<double>;
