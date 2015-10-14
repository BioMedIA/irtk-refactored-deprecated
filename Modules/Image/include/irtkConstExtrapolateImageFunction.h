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

#ifndef _IRTKCONSTEXTRAPOLATEIMAGEFUNCTION_H
#define _IRTKCONSTEXTRAPOLATEIMAGEFUNCTION_H


/**
 * Constant extrapolation (i.e., padding) of generic image
 *
 * This discrete image extrapolation function assumes a constant image value
 * for image grid points outside the discrete domain on which the discrete image
 * function is defined.
 */

template <class TImage>
class irtkGenericConstExtrapolateImageFunction
: public irtkGenericExtrapolateImageFunction<TImage>
{
  irtkGenericExtrapolatorMacro(
    irtkGenericConstExtrapolateImageFunction,
    Extrapolation_Const
  );

public:

  /// Default constructor
  irtkGenericConstExtrapolateImageFunction(double padding_value = .0)
  {
    this->_DefaultValue = padding_value;
  }

  /// Destructor
  virtual ~irtkGenericConstExtrapolateImageFunction() {}

  /// Get image value at an arbitrary discrete image location
  virtual VoxelType Get(int i, int j, int k = 0, int l = 0) const
  {
    if (0 <= i && i < this->X() &&
        0 <= j && j < this->Y() &&
        0 <= k && k < this->Z() &&
        0 <= l && l < this->T()) {
      return this->Input()->Get(i, j, k, l);
    } else {
      return voxel_cast<VoxelType>(this->_DefaultValue);
    }
  }

};


/**
 * Constant extrapolation (i.e., padding) of any image
 */

class irtkConstExtrapolateImageFunction
: public irtkGenericConstExtrapolateImageFunction<irtkBaseImage>
{
  irtkObjectMacro(irtkConstExtrapolateImageFunction);

public:

  /// Constructor
  irtkConstExtrapolateImageFunction(double padding_value = .0)
  :
    irtkGenericConstExtrapolateImageFunction<irtkBaseImage>(padding_value)
  {}

  /// Destructor
  virtual ~irtkConstExtrapolateImageFunction() {}

};


#endif
