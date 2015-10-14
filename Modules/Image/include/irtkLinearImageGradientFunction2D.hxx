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

#ifndef _IRTKLINEARIMAGEGRADIENTFUNCTION2D_HXX
#define _IRTKLINEARIMAGEGRADIENTFUNCTION2D_HXX

#include <irtkLinearImageGradientFunction2D.h>
#include <irtkLinearImageGradientFunction.hxx>


// -----------------------------------------------------------------------------
template <class TImage>
irtkGenericLinearImageGradientFunction2D<TImage>
::irtkGenericLinearImageGradientFunction2D()
{
  this->NumberOfDimensions(2);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename irtkGenericLinearImageGradientFunction2D<TImage>::GradientType
irtkGenericLinearImageGradientFunction2D<TImage>
::GetInside(double x, double y, double z, double t) const
{
  return this->GetInside2D(x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename irtkGenericLinearImageGradientFunction2D<TImage>::GradientType
irtkGenericLinearImageGradientFunction2D<TImage>
::GetWithPaddingInside(double x, double y, double z, double t) const
{
  return this->GetWithPaddingInside2D(x, y, z, t);
}


#endif
