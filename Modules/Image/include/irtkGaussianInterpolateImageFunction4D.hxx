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

#ifndef _IRTKGAUSSIANINTERPOLATEIMAGEFUNCTION4D_HXX
#define _IRTKGAUSSIANINTERPOLATEIMAGEFUNCTION4D_HXX

#include <irtkGaussianInterpolateImageFunction4D.h>
#include <irtkGaussianInterpolateImageFunction.hxx>


// -----------------------------------------------------------------------------
template <class TImage>
irtkGenericGaussianInterpolateImageFunction4D<TImage>
::irtkGenericGaussianInterpolateImageFunction4D(double sigma)
:
  irtkGenericGaussianInterpolateImageFunction<TImage>(sigma)
{
  this->NumberOfDimensions(4);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename irtkGenericGaussianInterpolateImageFunction4D<TImage>::VoxelType
irtkGenericGaussianInterpolateImageFunction4D<TImage>
::Get(double x, double y, double z, double t) const
{
  return this->Get4D(x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename irtkGenericGaussianInterpolateImageFunction4D<TImage>::VoxelType
irtkGenericGaussianInterpolateImageFunction4D<TImage>
::GetWithPadding(double x, double y, double z, double t) const
{
  return this->GetWithPadding4D(x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename TOtherImage::VoxelType
irtkGenericGaussianInterpolateImageFunction4D<TImage>
::Get(const TOtherImage *input, double x, double y, double z, double t) const
{
  return this->Get4D(input, x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename TOtherImage::VoxelType
irtkGenericGaussianInterpolateImageFunction4D<TImage>
::GetWithPadding(const TOtherImage *input, double x, double y, double z, double t) const
{
  return this->GetWithPadding4D(input, x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename irtkGenericGaussianInterpolateImageFunction4D<TImage>::VoxelType
irtkGenericGaussianInterpolateImageFunction4D<TImage>
::GetInside(double x, double y, double z, double t) const
{
  return Get(this->Input(), x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename irtkGenericGaussianInterpolateImageFunction4D<TImage>::VoxelType
irtkGenericGaussianInterpolateImageFunction4D<TImage>
::GetOutside(double x, double y, double z, double t) const
{
  if (this->Extrapolator()) {
    return Get(this->Extrapolator(), x, y, z, t);
  } else {
    return Get(x, y, z, t);
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename irtkGenericGaussianInterpolateImageFunction4D<TImage>::VoxelType
irtkGenericGaussianInterpolateImageFunction4D<TImage>
::GetWithPaddingInside(double x, double y, double z, double t) const
{
  return GetWithPadding(this->Input(), x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename irtkGenericGaussianInterpolateImageFunction4D<TImage>::VoxelType
irtkGenericGaussianInterpolateImageFunction4D<TImage>
::GetWithPaddingOutside(double x, double y, double z, double t) const
{
  if (this->Extrapolator()) {
    return GetWithPadding(this->Extrapolator(), x, y, z, t);
  } else {
    return GetWithPadding(x, y, z, t);
  }
}


#endif
