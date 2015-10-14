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

#ifndef _IRTKFASTCUBICBSPLINEINTERPOLATEIMAGEFUNCTION3D_HXX
#define _IRTKFASTCUBICBSPLINEINTERPOLATEIMAGEFUNCTION3D_HXX

#include <irtkFastCubicBSplineInterpolateImageFunction3D.h>
#include <irtkFastCubicBSplineInterpolateImageFunction.hxx>


// -----------------------------------------------------------------------------
template <class TImage>
irtkGenericFastCubicBSplineInterpolateImageFunction3D<TImage>
::irtkGenericFastCubicBSplineInterpolateImageFunction3D()
{
  this->NumberOfDimensions(3);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename irtkGenericFastCubicBSplineInterpolateImageFunction3D<TImage>::VoxelType
irtkGenericFastCubicBSplineInterpolateImageFunction3D<TImage>
::Get(double x, double y, double z, double t) const
{
  return this->Get3D(x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename irtkGenericFastCubicBSplineInterpolateImageFunction3D<TImage>::VoxelType
irtkGenericFastCubicBSplineInterpolateImageFunction3D<TImage>
::GetWithPadding(double x, double y, double z, double t) const
{
  return this->GetWithPadding3D(x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename TOtherImage::VoxelType
irtkGenericFastCubicBSplineInterpolateImageFunction3D<TImage>
::Get(const TOtherImage *coeff, double x, double y, double z, double t) const
{
  return this->Get3D(coeff, x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage, class TCoefficient>
inline typename TCoefficient::VoxelType
irtkGenericFastCubicBSplineInterpolateImageFunction3D<TImage>
::GetWithPadding(const TOtherImage *input, const TCoefficient *coeff,
                 double x, double y, double z, double t) const
{
  return this->GetWithPadding3D(input, coeff, x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename irtkGenericFastCubicBSplineInterpolateImageFunction3D<TImage>::VoxelType
irtkGenericFastCubicBSplineInterpolateImageFunction3D<TImage>
::GetInside(double x, double y, double z, double t) const
{
  // Use faster coefficient iteration than Get3D(Coefficient(), x, y, z, t)
  return this->GetInside3D(x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename irtkGenericFastCubicBSplineInterpolateImageFunction3D<TImage>::VoxelType
irtkGenericFastCubicBSplineInterpolateImageFunction3D<TImage>
::GetOutside(double x, double y, double z, double t) const
{
  if (this->_InfiniteCoefficient) {
    return Get(this->_InfiniteCoefficient, x, y, z, t);
  } else {
    return Get(x, y, z, t);
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename irtkGenericFastCubicBSplineInterpolateImageFunction3D<TImage>::VoxelType
irtkGenericFastCubicBSplineInterpolateImageFunction3D<TImage>
::GetWithPaddingInside(double x, double y, double z, double t) const
{
  return GetWithPadding(this->Input(), &this->_Coefficient, x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename irtkGenericFastCubicBSplineInterpolateImageFunction3D<TImage>::VoxelType
irtkGenericFastCubicBSplineInterpolateImageFunction3D<TImage>
::GetWithPaddingOutside(double x, double y, double z, double t) const
{
  if (this->Extrapolator() && this->_InfiniteCoefficient) {
    return GetWithPadding(this->Extrapolator(), this->_InfiniteCoefficient, x, y, z, t);
  } else {
    return GetWithPadding(x, y, z, t);
  }
}


#endif
