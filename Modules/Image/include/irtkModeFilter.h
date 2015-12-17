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

#ifndef _IRTKMODEFILTER_H

#define _IRTKMODEFILTER_H

#include <irtkImageToImage.h>

/**
 * Class for applying mode filter to what should be label images.
 *
 * Assign to each voxel, the modal label of those within a neighbourhood.
 */

template <class VoxelType>
class irtkModeFilter : public irtkImageToImage<VoxelType>
{
  irtkImageFilterMacro(irtkModeFilter);

protected:

  /// Initialize the filter
  virtual void Initialize();

  /// What connectivity to assume when running the filter.
  irtkConnectivityType _Connectivity;

  // List of voxel offsets of the neighbourhood.
  irtkNeighbourhoodOffsets _offsets;

public:

  /// Constructor
  irtkModeFilter();

  /// Destructor
  ~irtkModeFilter();

  /// Run mode filter
  virtual void Run();

  SetMacro(Connectivity, irtkConnectivityType);

  GetMacro(Connectivity, irtkConnectivityType);
};

#endif
