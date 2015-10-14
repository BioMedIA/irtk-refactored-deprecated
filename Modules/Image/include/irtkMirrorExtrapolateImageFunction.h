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

#ifndef _IRTKMIRROREXTRAPOLATEIMAGEFUNCTION_H
#define _IRTKMIRROREXTRAPOLATEIMAGEFUNCTION_H


/**
 * Extrapolation of generic image by mirroring along the boundaries
 */

template <class TImage>
class irtkGenericMirrorExtrapolateImageFunction
: public irtkIndexExtrapolateImageFunction<TImage>
{
  irtkExtrapolatorMacro(
    irtkGenericMirrorExtrapolateImageFunction,
    Extrapolation_Mirror
  );

public:

  /// Constructor
  irtkGenericMirrorExtrapolateImageFunction() {}

  /// Destructor
  virtual ~irtkGenericMirrorExtrapolateImageFunction() {}

  /// Mirror index at boundary such that it is inside the range [0, max]
  /// \note Use static function as irtkMirrorExtrapolateImageFunction::Apply.
  static void Apply(int &index, int max)
  {
    if (max == 0) {
      index = 0;
    } else if (index < 0) {
      index = -index;
      int n = index / max;
      int m = index - n * max;
      if (n & 1) index = max - m;
      else       index = m;
    } else if (index > max) {
      index -= max;
      int n = index / max;
      int m = index - n * max;
      if (n & 1) index = m;
      else       index = max - m;
    }
  }

  /// Mirror index at boundary such that it is inside the range [0, max]
  virtual void TransformIndex(int &index, int max) const
  {
    Apply(index, max);
  }

};


/**
 * Extrapolation of any image by mirroring along the boundaries
 */

class irtkMirrorExtrapolateImageFunction
: public irtkGenericMirrorExtrapolateImageFunction<irtkBaseImage>
{
  irtkObjectMacro(irtkMirrorExtrapolateImageFunction);

public:

  /// Constructor
  irtkMirrorExtrapolateImageFunction() {}

  /// Destructor
  virtual ~irtkMirrorExtrapolateImageFunction() {}

};


#endif
