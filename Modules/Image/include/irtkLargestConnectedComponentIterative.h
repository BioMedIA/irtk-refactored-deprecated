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

#ifndef _IRTKLARGESTCONNECTEDCOMPONENTITERATIVE_H

#define _IRTKLARGESTCONNECTEDCOMPONENTITERATIVE_H

#include <irtkImageToImage.h>

/**
 * Class for extracting the largest connected component from a labelled image
 *
 * This class defines and implements the extraction of the largest
 * connected component from a labelled image.  Uses an iterative method so
 * takes longer than the class irtkLargestConnectedComponent.  The stack
 * size needed, however, is lower.
 *
 */

template <class VoxelType>
class irtkLargestConnectedComponentIterative : public irtkImageToImage<VoxelType>
{
  irtkImageFilterMacro(irtkLargestConnectedComponentIterative);

private:

  /// Label used to identify labels of interest
  VoxelType _TargetLabel;

  /// Size of largest cluster with the target label.
  int _largestClusterSize;

  // What label is used for the largest cluster during computation.
  VoxelType _largestClusterLabel;

  // Record of all cluster sizes.
  int *_ClusterSizes;

  // Number of clusters that match the target label.
  int _NumberOfClusters;

  /// Modes
  bool _Mode2D;

  bool _AllClustersMode;

protected:

  /// Run filter
  virtual void Run3D();
  virtual void Run2D();

  // Helper functions
  virtual void ResetMarks();

  virtual void SelectLargestCluster();

  int CheckAdjacency2D(VoxelType& markA, VoxelType& markB);
  int CheckAdjacency3D(VoxelType& markA, VoxelType& markB);

public:

  /// Constructor
  irtkLargestConnectedComponentIterative(VoxelType = 0);

  /// Destructor
  ~irtkLargestConnectedComponentIterative();

  /// Set label sought
  SetMacro(TargetLabel, VoxelType);

  /// Get label sought
  GetMacro(TargetLabel, VoxelType);

  /// Set mode
  SetMacro(Mode2D, bool);

  /// Get mode
  GetMacro(Mode2D, bool);

  /// Set mode
  SetMacro(AllClustersMode, bool);

  /// Get mode
  GetMacro(AllClustersMode, bool);

  /// Run filter
  virtual void Run();
};

#endif
