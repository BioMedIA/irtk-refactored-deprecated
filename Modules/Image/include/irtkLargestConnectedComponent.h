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

#ifndef _IRTKLARGESTCONNECTEDCOMPONENT_H

#define _IRTKLARGESTCONNECTEDCOMPONENT_H

#include <irtkImageToImage.h>

/**
 * Class for extracting the largest connected component from a labelled image
 *
 * This class defines and implements the extraction of the largest connected component
 * from a labelled image.
 *
 */

template <class VoxelType>
class irtkLargestConnectedComponent : public irtkImageToImage<VoxelType>
{
  irtkImageFilterMacro(irtkLargestConnectedComponent);

private:

  /// Size of current cluster
  int _currentClusterSize;

  /// Size of largest cluster
  int _largestClusterSize;

  /// Label used to identify labels of interest
  VoxelType _ClusterLabel;

  /// Mode
  bool _Mode2D;

protected:

  /// Recursive region growing
  void Grow2D(int, int, int, int);

  /// Recursive region growing
  void Grow3D(int, int, int, int);

public:

  /// Constructor
  irtkLargestConnectedComponent(VoxelType = 0);

  /// Destructor
  ~irtkLargestConnectedComponent();

  /// Set sigma
  SetMacro(ClusterLabel, VoxelType);

  /// Get sigma
  GetMacro(ClusterLabel, VoxelType);

  /// Set mode
  SetMacro(Mode2D, bool);

  /// Get mode
  GetMacro(Mode2D, bool);

  /// Run filter
  virtual void Run();

};

#endif
