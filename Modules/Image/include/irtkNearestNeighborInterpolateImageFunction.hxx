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

#ifndef _IRTKNEARESTNEIGHBORINTERPOLATEIMAGEFUNCTION_HXX
#define _IRTKNEARESTNEIGHBORINTERPOLATEIMAGEFUNCTION_HXX

#include <irtkNearestNeighborInterpolateImageFunction.h>
#include <irtkInterpolateImageFunction.hxx>


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
irtkGenericNearestNeighborInterpolateImageFunction<TImage>
::irtkGenericNearestNeighborInterpolateImageFunction()
{
}

// -----------------------------------------------------------------------------
template <class TImage>
irtkGenericNearestNeighborInterpolateImageFunction<TImage>
::~irtkGenericNearestNeighborInterpolateImageFunction()
{
}

// =============================================================================
// Domain checks
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
void irtkGenericNearestNeighborInterpolateImageFunction<TImage>
::BoundingInterval(double x, int &i, int &I) const
{
  i = I = static_cast<int>(round(x));
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
inline typename irtkGenericNearestNeighborInterpolateImageFunction<TImage>::VoxelType
irtkGenericNearestNeighborInterpolateImageFunction<TImage>
::Get(double x, double y, double z, double t) const
{
  const int i = static_cast<int>(round(x));
  const int j = static_cast<int>(round(y));
  const int k = static_cast<int>(round(z));
  const int l = static_cast<int>(round(t));

  if (this->Input()->IsInside(i, j, k, l)) {
    return this->Input()->Get(i, j, k, l);
  } else {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename irtkGenericNearestNeighborInterpolateImageFunction<TImage>::VoxelType
irtkGenericNearestNeighborInterpolateImageFunction<TImage>
::GetWithPadding(double x, double y, double z, double t) const
{
  const int i = static_cast<int>(round(x));
  const int j = static_cast<int>(round(y));
  const int k = static_cast<int>(round(z));
  const int l = static_cast<int>(round(t));

  if (this->Input()->IsInsideForeground(i, j, k, l)) {
    return this->Input()->Get(i, j, k, l);
  } else {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename TOtherImage::VoxelType
irtkGenericNearestNeighborInterpolateImageFunction<TImage>
::Get(const TOtherImage *input, double x, double y, double z, double t) const
{
  return input->Get(static_cast<int>(round(x)),
                    static_cast<int>(round(y)),
                    static_cast<int>(round(z)),
                    static_cast<int>(round(t)));
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename TOtherImage::VoxelType
irtkGenericNearestNeighborInterpolateImageFunction<TImage>
::GetWithPadding(const TOtherImage *input, double x, double y, double z, double t) const
{
  const int i = static_cast<int>(round(x));
  const int j = static_cast<int>(round(y));
  const int k = static_cast<int>(round(z));
  const int l = static_cast<int>(round(t));

  if (input->IsForeground(i, j, k, l)) {
    return input->Get(i, j, k, l);
  } else {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename irtkGenericNearestNeighborInterpolateImageFunction<TImage>::VoxelType
irtkGenericNearestNeighborInterpolateImageFunction<TImage>
::GetInside(double x, double y, double z, double t) const
{
  return Get(this->Input(), x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename irtkGenericNearestNeighborInterpolateImageFunction<TImage>::VoxelType
irtkGenericNearestNeighborInterpolateImageFunction<TImage>
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
inline typename irtkGenericNearestNeighborInterpolateImageFunction<TImage>::VoxelType
irtkGenericNearestNeighborInterpolateImageFunction<TImage>
::GetWithPaddingInside(double x, double y, double z, double t) const
{
  return GetWithPadding(this->Input(), x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename irtkGenericNearestNeighborInterpolateImageFunction<TImage>::VoxelType
irtkGenericNearestNeighborInterpolateImageFunction<TImage>
::GetWithPaddingOutside(double x, double y, double z, double t) const
{
  if (this->Extrapolator()) {
    return GetWithPadding(this->Extrapolator(), x, y, z, t);
  } else {
    return GetWithPadding(x, y, z, t);
  }
}


#endif
