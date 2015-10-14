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

#ifndef _IRTKREPEATEXTRAPOLATEIMAGEFUNCTION_H
#define _IRTKREPEATEXTRAPOLATEIMAGEFUNCTION_H


/**
 * Extrapolation of generic image by tiling/repeating it
 *
 * This extrapolation implements a periodic extension of the finite image.
 */

template <class TImage>
class irtkGenericRepeatExtrapolateImageFunction
: public irtkIndexExtrapolateImageFunction<TImage>
{
  irtkExtrapolatorMacro(
    irtkGenericRepeatExtrapolateImageFunction,
    Extrapolation_Repeat
  );

public:

  /// Constructor
  irtkGenericRepeatExtrapolateImageFunction() {}

  /// Destructor
  virtual ~irtkGenericRepeatExtrapolateImageFunction() {}

  /// Transform index such that it is inside the range [0, max]
  /// \note Use static function as irtkRepeatExtrapolateImageFunction::Apply.
  static void Apply(int &index, int max)
  {
    if      (index < 0  ) index = max + (index + 1) % (max + 1);
    else if (index > max) index =        index      % (max + 1);
  }

  /// Transform continuous index such that it is inside the range [0, max + 1)
  /// \note Use static function as irtkRepeatExtrapolateImageFunction::Apply.
  static void Apply(double &cindex, int max)
  {
    int    index    = static_cast<int>(floor(cindex));
    double fraction = cindex - index;
    Apply(index, max);
    cindex = index + fraction;
  }

  /// Transform index such that it is inside the range [0, max]
  virtual void TransformIndex(int &index, int max) const
  {
    Apply(index, max);
  }

};


/**
 * Extrapolation of any image by tiling/repeating it
 *
 * This extrapolation implements a periodic extension of the finite image.
 */

class irtkRepeatExtrapolateImageFunction
: public irtkGenericRepeatExtrapolateImageFunction<irtkBaseImage>
{
  irtkObjectMacro(irtkRepeatExtrapolateImageFunction);

public:

  /// Constructor
  irtkRepeatExtrapolateImageFunction() {}

  /// Destructor
  virtual ~irtkRepeatExtrapolateImageFunction() {}

};


#endif
