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

#include <irtkLargestConnectedComponent.h>

template <class VoxelType> irtkLargestConnectedComponent<VoxelType>::irtkLargestConnectedComponent(VoxelType ClusterLabel)
{
  _currentClusterSize = 0;
  _largestClusterSize = 0;
  _Mode2D = false;
  _ClusterLabel = ClusterLabel;
}

template <class VoxelType> irtkLargestConnectedComponent<VoxelType>::~irtkLargestConnectedComponent(void)
{}

template <class VoxelType> void irtkLargestConnectedComponent<VoxelType>::Grow2D(int x, int y, int z, int t)
{
  if ((x >= 0) && (x < this->_input->GetX()) && (y >= 0) && (y < this->_input->GetY()) &&
      (z >= 0) && (z < this->_input->GetZ()) && (t >= 0) && (t < this->_input->GetT())) {
    if ((this->_input->Get(x, y, z, t) == this->_ClusterLabel) && (this->_output->Get(x, y, z, t) == 0)) {
      this->_output->Put(x, y, z, t, 1);
      this->_currentClusterSize++;
      this->Grow2D(x-1, y  , z, t);
      this->Grow2D(x+1, y  , z, t);
      this->Grow2D(x  , y-1, z, t);
      this->Grow2D(x  , y+1, z, t);
    }
  }
}

template <class VoxelType> void irtkLargestConnectedComponent<VoxelType>::Grow3D(int x, int y, int z, int t)
{
  if ((x >= 0) && (x < this->_input->GetX()) && (y >= 0) && (y < this->_input->GetY()) &&
      (z >= 0) && (z < this->_input->GetZ()) && (t >= 0) && (t < this->_input->GetT())) {
    if ((this->_input->Get(x, y, z, t) == this->_ClusterLabel) && (this->_output->Get(x, y, z, t) == 0)) {
      this->_output->Put(x, y, z, t, 1);
      this->_currentClusterSize++;
      this->Grow3D(x-1, y  , z  , t);
      this->Grow3D(x+1, y  , z  , t);
      this->Grow3D(x  , y-1, z  , t);
      this->Grow3D(x  , y+1, z  , t);
      this->Grow3D(x  , y  , z-1, t);
      this->Grow3D(x  , y  , z+1, t);
    }
  }
}

template <class VoxelType> void irtkLargestConnectedComponent<VoxelType>::Run()
{
  int x, y, z, t = 0, largestClusterSeed_x, largestClusterSeed_y, largestClusterSeed_z;

  // Do the initial set up
  this->Initialize();

  // Check if the input is 4D
  if (this->_input->GetT() > 1) {
    cerr << "irtkLargestConnectedComponent<VoxelType>::Run(): 4D images not yet supported\n" << endl;
    exit(1);
  }

  // Do conneted component analysis
  for (z = 0; z < this->_input->GetZ(); z++) {
    for (y = 0; y < this->_input->GetY(); y++) {
      for (x = 0; x < this->_input->GetX(); x++) {
        this->_output->Put(x, y, z, t, 0);
      }
    }
  }

  if (this->_Mode2D == true) {

    for (z = 0; z < this->_input->GetZ(); z++) {
      largestClusterSeed_x = 0;
      largestClusterSeed_y = 0;
      largestClusterSeed_z = z;

      for (y = 0; y < this->_input->GetY(); y++) {
        for (x = 0; x < this->_input->GetX(); x++) {
          this->_currentClusterSize = 0;
          this->Grow2D(x, y, z, t);
          if (this->_currentClusterSize > this->_largestClusterSize) {
            this->_largestClusterSize   = this->_currentClusterSize;
            largestClusterSeed_x = x;
            largestClusterSeed_y = y;
            largestClusterSeed_z = z;
          }
        }
      }

      for (y = 0; y < this->_input->GetY(); y++) {
        for (x = 0; x < this->_input->GetX(); x++) {
          this->_output->Put(x, y, z, t, 0);
        }
      }

      this->Grow2D(largestClusterSeed_x, largestClusterSeed_y, largestClusterSeed_z, t);
    }
  } else {

    largestClusterSeed_x = 0;
    largestClusterSeed_y = 0;
    largestClusterSeed_z = 0;

    for (z = 0; z < this->_input->GetZ(); z++) {
      for (y = 0; y < this->_input->GetY(); y++) {
        for (x = 0; x < this->_input->GetX(); x++) {
          this->_currentClusterSize = 0;
          this->Grow3D(x, y, z, t);
          if (this->_currentClusterSize > this->_largestClusterSize) {
            this->_largestClusterSize   = this->_currentClusterSize;
            largestClusterSeed_x = x;
            largestClusterSeed_y = y;
            largestClusterSeed_z = z;
          }
        }
      }
    }

    for (z = 0; z < this->_input->GetZ(); z++) {
      for (y = 0; y < this->_input->GetY(); y++) {
        for (x = 0; x < this->_input->GetX(); x++) {
          this->_output->Put(x, y, z, t, 0);
        }
      }
    }

    this->Grow3D(largestClusterSeed_x, largestClusterSeed_y, largestClusterSeed_z, t);
  }

  // Do the final cleaning up
  this->Finalize();
}

template class irtkLargestConnectedComponent<irtkBytePixel>;
template class irtkLargestConnectedComponent<irtkGreyPixel>;
template class irtkLargestConnectedComponent<irtkRealPixel>;
