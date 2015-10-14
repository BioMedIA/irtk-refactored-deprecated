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

#ifndef _IRTKRESAMPLINGWITHPADDING_H

#define _IRTKRESAMPLINGWITHPADDING_H

#ifdef HAS_TBB

template <class VoxelType> class irtkMultiThreadedResamplingWithPadding;

#endif

/**
 * Class for resampling of padded images
 *
 * This class defines and implements the resampling of images with arbitrary
 * voxel dimensions.  The new image intensity of the voxels is calculated by
 * interpolation of the old image intensities.  Only linear interpolation is
 * currently supported. If more than 50% of the voxels used for interpolation
 * have intensities smaller or equal to the padding value, the resampled
 * voxel will be filled with the padding value.
 */


template <class VoxelType>
class irtkResamplingWithPadding : public irtkResampling<VoxelType>
{
  irtkImageFilterMacro(irtkResamplingWithPadding);

#ifdef HAS_TBB

  friend class irtkMultiThreadedResamplingWithPadding<VoxelType>;

#endif

protected:

  /// Padding value
  VoxelType _PaddingValue;

  /// Initialize the filter
  virtual void Initialize();

public:

  /// Constructor
  irtkResamplingWithPadding(double, double, double, VoxelType);

  /// Constructor
  irtkResamplingWithPadding(int, int, int, VoxelType);

  /// Constructor
  irtkResamplingWithPadding(int, int, int, double, double, double, VoxelType);

  /// Run the resampling filter
  virtual void Run();

  /// Set padding value
  SetMacro(PaddingValue, VoxelType);

  /// Get padding value
  GetMacro(PaddingValue, VoxelType);

};

#endif
