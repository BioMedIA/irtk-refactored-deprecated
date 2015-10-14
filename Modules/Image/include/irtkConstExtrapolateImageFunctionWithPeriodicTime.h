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

#ifndef _IRTKCONSTEXTRAPOLATEIMAGEFUNCTIONWITHPERIODICTIME_H
#define _IRTKCONSTEXTRAPOLATEIMAGEFUNCTIONWITHPERIODICTIME_H

#include <irtkConstExtrapolateImageFunction.h>
#include <irtkRepeatExtrapolateImageFunction.h>


/**
 * Extrapolation of generic image by padding it with a constant value
 * in the spatial domain, but a periodic extension in the temporal dimension
 */

template <class TImage>
class irtkGenericConstExtrapolateImageFunctionWithPeriodicTime
: public irtkGenericConstExtrapolateImageFunction<TImage>
{
  irtkExtrapolatorMacro(
    irtkGenericConstExtrapolateImageFunctionWithPeriodicTime,
    Extrapolation_ConstWithPeriodicTime
  );

public:

  typedef TImage                        ImageType; ///< Input image type
  typedef typename ImageType::VoxelType VoxelType; ///< Input voxel type
  typedef typename ImageType::RealType  RealType;  ///< Compatible floating-point type

  /// Constructor
  irtkGenericConstExtrapolateImageFunctionWithPeriodicTime(double padding_value = .0)
  :
    irtkGenericConstExtrapolateImageFunction<TImage>(padding_value)
  {}

  /// Destructor
  virtual ~irtkGenericConstExtrapolateImageFunctionWithPeriodicTime() {}

  // Import overloaded non-virtual member functions from base class
  using irtkGenericConstExtrapolateImageFunction<TImage>::Get;

  /// Get image value at an arbitrary discrete image location
  virtual VoxelType Get(int i, int j, int k = 0, int l = 0) const
  {
    irtkRepeatExtrapolateImageFunction::Apply(l, this->T() - 1);
    return irtkGenericConstExtrapolateImageFunction<TImage>::Get(i, j, k, l);
  }

};


/**
 * Extrapolation of any image by padding it with a constant value
 * in the spatial domain, but a periodic extension in the temporal dimension
 */

class irtkConstExtrapolateImageFunctionWithPeriodicTime
: public irtkGenericConstExtrapolateImageFunctionWithPeriodicTime<irtkBaseImage>
{
  irtkObjectMacro(irtkConstExtrapolateImageFunctionWithPeriodicTime);

public:

  /// Constructor
  irtkConstExtrapolateImageFunctionWithPeriodicTime() {}

  /// Destructor
  virtual ~irtkConstExtrapolateImageFunctionWithPeriodicTime() {}

};


#endif
