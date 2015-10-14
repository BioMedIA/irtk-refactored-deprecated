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

#ifndef _IRTKNEARESTNEIGHBOREXTRAPOLATEIMAGEFUNCTION_H
#define _IRTKNEARESTNEIGHBOREXTRAPOLATEIMAGEFUNCTION_H


/**
 * Nearest neighbor extrapolation of generic image
 *
 * The nearest neighor extrapolation of a discrete image corresponds to a
 * Neumann boundary condition with the derivative value outside the image domain
 * set to zero.
 */

template <class TImage>
class irtkGenericNearestNeighborExtrapolateImageFunction
: public irtkIndexExtrapolateImageFunction<TImage>
{
  irtkExtrapolatorMacro(
    irtkGenericNearestNeighborExtrapolateImageFunction,
    Extrapolation_NN
  );

public:

  /// Constructor
  irtkGenericNearestNeighborExtrapolateImageFunction() {}

  /// Destructor
  virtual ~irtkGenericNearestNeighborExtrapolateImageFunction() {}

  /// Transform index such that it is inside the range [0, max]
  virtual void TransformIndex(int &index, int max) const
  {
    if      (index < 0  ) index = 0;
    else if (index > max) index = max;
  }

};


/**
 * Nearest neighbor extrapolation of any image
 *
 * The nearest neighor extrapolation of a discrete image corresponds to a
 * Neumann boundary condition with the derivative value outside the image domain
 * set to zero.
 */

class irtkNearestNeighborExtrapolateImageFunction
: public irtkGenericNearestNeighborExtrapolateImageFunction<irtkBaseImage>
{
  irtkObjectMacro(irtkNearestNeighborExtrapolateImageFunction);

public:

  /// Constructor
  irtkNearestNeighborExtrapolateImageFunction() {}

  /// Destructor
  virtual ~irtkNearestNeighborExtrapolateImageFunction() {}

};


#endif
