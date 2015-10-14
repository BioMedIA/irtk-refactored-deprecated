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

#ifndef _IRTKDILATION_H

#define _IRTKDILATION_H

#include <irtkImageToImage.h>

/**
 * Class for dilation of images
 *
 * This class defines and implements the morphological dilation of images.
 *
 */

template <class VoxelType>
class irtkDilation : public irtkImageToImage<VoxelType>
{
  irtkImageFilterMacro(irtkDilation);

protected:

  /// Initialize the filter
  virtual void Initialize();

  /// What connectivity to assume when running the filter.
  irtkConnectivityType _Connectivity;

  // List of voxel offsets of the neighbourhood.
  irtkNeighbourhoodOffsets _offsets;

public:

  /// Constructor
  irtkDilation();

  /// Destructor
  ~irtkDilation();

  /// Run dilation
  virtual void Run();

  SetMacro(Connectivity, irtkConnectivityType);

  GetMacro(Connectivity, irtkConnectivityType);
};

#endif
