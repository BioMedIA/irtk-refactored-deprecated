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

#ifndef _IRTKBINARYVOXELFUNCTION_H
#define _IRTKBINARYVOXELFUNCTION_H

#include <irtkVoxelFunction.h>
#include <irtkImageFunction.h>


/**
 * These basic binary voxel functions can be used as VoxelFunc template parameter
 * of the binary ForEachVoxel function templates as follows:
 *
 * \code
 * // Add one image to another in-place
 * irtkGreyImage input (attr);
 * irtkGreyImage output(attr);
 * // Set image in-place to maximum of both images
 * irtkBinaryVoxelFunction::Max max;
 * ForEachVoxel(input, output, max);
 * irtkBinaryVoxelFunction::Add add;
 * ParallelForEachVoxel(input, output, add);
 * // Compute sum-of-squared differences (SSD)
 * irtkGreyImage target(attr);
 * irtkGreyImage source(attr);
 * irtkBinaryVoxelFunction::SSD ssd;
 * ForEachVoxel(target, source, ssd);
 * printf("SSD=%f\n", ssd.value);
 * \endcode
 *
 */

namespace irtkBinaryVoxelFunction {


// =============================================================================
// Copy
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * Copies the voxel value of one image to the corresponding voxel of another
 */
struct Copy : public irtkVoxelFunction
{
  template <class T1, class T2>
  void operator ()(const T1 *in, T2 *out) { *out = *in; }

  template <class TImage, class T1, class T2>
  void operator ()(const TImage&, int, const T1 *in, T2 *out)
  {
    this->operator ()(in, out);
  }

  template <class T1, class T2>
  void operator ()(int, int, int, int, const T1 *in, T2 *out)
  {
    this->operator ()(in, out);
  }
};

// =============================================================================
// Basic mathematical voxel-wise in-place operations
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * Add image to another
 */
struct Add : public irtkVoxelFunction
{
  template <class T1, class T2>
  void operator ()(const T1 *in, T2 *out)
  {
    *out = static_cast<T2>(static_cast<double>(*out) + static_cast<double>(*in));
  }

  template <class TImage, class T1, class T2>
  void operator ()(const TImage&, int, const T1 *in, T2 *out)
  {
    this->operator ()(in, out);
  }

  template <class T1, class T2>
  void operator ()(int, int, int, int, const T1 *in, T2 *out)
  {
    this->operator ()(in, out);
  }
};

// -----------------------------------------------------------------------------
/**
 * Subtract image from another
 */
struct Sub : public irtkVoxelFunction
{
  template <class T1, class T2>
  void operator ()(const T1 *in, T2 *out)
  {
    *out = static_cast<T2>(static_cast<double>(*out) - static_cast<double>(*in));
  }

  template <class TImage, class T1, class T2>
  void operator ()(const TImage&, int, const T1 *in, T2 *out)
  {
    this->operator ()(in, out);
  }

  template <class T1, class T2>
  void operator ()(int, int, int, int, const T1 *in, T2 *out)
  {
    this->operator ()(in, out);
  }
};

// -----------------------------------------------------------------------------
/**
 * Multiplies the voxel value of one image by the value of another
 */
struct Mul : public irtkVoxelFunction
{
  template <class T1, class T2>
  void operator ()(const T1 *in, T2 *out)
  {
    *out = static_cast<T2>(static_cast<double>(*out) * static_cast<double>(*in));
  }

  template <class TImage, class T1, class T2>
  void operator ()(const TImage&, int, const T1 *in, T2 *out)
  {
    this->operator ()(in, out);
  }

  template <class T1, class T2>
  void operator ()(int, int, int, int, const T1 *in, T2 *out)
  {
    this->operator ()(in, out);
  }
};

// -----------------------------------------------------------------------------
/**
 * Divides the voxel value of one image by the value of another
 */
struct Div : public irtkVoxelFunction
{
  template <class T1, class T2>
  void operator ()(const T1 *in, T2 *out)
  {
    double divisor = static_cast<double>(*in);
    if (divisor == .0) {
      *out = static_cast<T2>(0);
    } else {
      *out = static_cast<T2>(static_cast<double>(*out) / divisor);
    }
  }

  template <class TImage, class T1, class T2>
  void operator ()(const TImage&, int, const T1 *in, T2 *out)
  {
    this->operator ()(in, out);
  }

  template <class T1, class T2>
  void operator ()(int, int, int, int, const T1 *in, T2 *out)
  {
    this->operator ()(in, out);
  }
};

// =============================================================================
// Image similarity measures
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * Compute sum-of-squared differences (SSD)
 */
struct SSD : public irtkVoxelReduction
{
  double value;

  SSD()       : value(.0)      {}
  SSD(SSD &o) : value(o.value) {}

  void split(const SSD &)    { value = .0; }
  void join (const SSD &rhs) { value += rhs.value; }

  template <class T1, class T2>
  void operator ()(const T1 *in1, const T2 *in2)
  {
    double diff = static_cast<double>(*in1) - static_cast<double>(*in2);
    value += diff * diff;
  }

  template <class TImage, class T1, class T2>
  void operator ()(const TImage&, int, const T1 *in, T2 *out)
  {
    this->operator ()(in, out);
  }

  template <class T1, class T2>
  void operator ()(int, int, int, int, const T1 *in, T2 *out)
  {
    this->operator ()(in, out);
  }
};

// =============================================================================
// Composition
// =============================================================================
  
// -----------------------------------------------------------------------------
/// Compose two 2D displacement fields: D = D1 ° D2
template <class TReal, class TInterpolator = irtkInterpolateImageFunction>
struct ComposeDisplacementFields2D : irtkVoxelFunction
{
  // ---------------------------------------------------------------------------
  ComposeDisplacementFields2D(TInterpolator *d1, const irtkGenericImage<TReal> *d2)
  :
    _D1(d1->GetInput()),
    _D1Interpolator(d1),
    _D2(d2),
    _y(d2->GetX() * d2->GetY())
  {
  }
  
  // ---------------------------------------------------------------------------
  void operator()(int i, int j, int k, int, const TReal *d2, TReal *dout)
  {
    double d[3]; // 2D displacement field can have constant third component!
    double x1 = i, y1 = j, z1 = k;
    _D2->ImageToWorld(x1, y1, z1);
    double x2 = x1 + d2[_x]; double x = x2;
    double y2 = y1 + d2[_y]; double y = y2;
    _D1            ->WorldToImage(x, y, z1);
    _D1Interpolator->Evaluate (d, x, y, z1);
    x2 = x2 + d[0];
    y2 = y2 + d[1];
    dout[_x] = x2 - x1;
    dout[_y] = y2 - y1;
  }
  
private:
  const irtkBaseImage           *_D1;             ///< First displacement field
  TInterpolator                 *_D1Interpolator; ///< Interpolator of first displacement field
  const irtkGenericImage<TReal> *_D2;             ///< Second displacement field
  
  static const int _x = 0; ///< Offset of x component
  int              _y;     ///< Offset of y component
};

// -----------------------------------------------------------------------------
/// Compose two 3D displacement fields: D = D1 ° D2
template <class TReal, class TInterpolator = irtkInterpolateImageFunction>
struct ComposeDisplacementFields3D : irtkVoxelFunction
{
  // ---------------------------------------------------------------------------
  ComposeDisplacementFields3D(TInterpolator *d1, const irtkGenericImage<TReal> *d2)
  :
    _D1(d1->GetInput()),
    _D1Interpolator(d1),
    _D2(d2),
    _y(d2->GetX() * d2->GetY() * d2->GetZ()),
    _z(2 * _y)
  {
  }
  
  // ---------------------------------------------------------------------------
  void operator()(int i, int j, int k, int, const TReal *d2, TReal *dout)
  {
    double d[3];
    double x1 = i, y1 = j, z1 = k;
    _D2->ImageToWorld(x1, y1, z1);
    double x2 = x1 + d2[_x]; double x = x2;
    double y2 = y1 + d2[_y]; double y = y2;
    double z2 = z1 + d2[_z]; double z = z2;
    _D1            ->WorldToImage(x, y, z);
    _D1Interpolator->Evaluate (d, x, y, z);
    x2 = x2 + d[0];
    y2 = y2 + d[1];
    z2 = z2 + d[2];
    dout[_x] = x2 - x1;
    dout[_y] = y2 - y1;
    dout[_z] = z2 - z1;
  }
  
private:
  const irtkBaseImage           *_D1;             ///< First displacement field
  TInterpolator                 *_D1Interpolator; ///< Interpolator of first displacement field
  const irtkGenericImage<TReal> *_D2;             ///< Second displacement field
  
  static const int _x = 0; ///< Offset of x component
  int              _y;     ///< Offset of y component
  int              _z;     ///< Offset of z component
};


} // namespace irtkBinaryVoxelFunction

#endif
