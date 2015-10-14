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

#ifndef _IRTKINTERPOLATEIMAGEFUNCTION_HXX
#define _IRTKINTERPOLATEIMAGEFUNCTION_HXX

#include <irtkInterpolateImageFunction.h>


////////////////////////////////////////////////////////////////////////////////
// irtkGenericInterpolateImageFunction
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
inline irtkGenericInterpolateImageFunction<TImage>
::irtkGenericInterpolateImageFunction()
{
}

// -----------------------------------------------------------------------------
template <class TImage>
inline irtkGenericInterpolateImageFunction<TImage>
::~irtkGenericInterpolateImageFunction()
{
}

// -----------------------------------------------------------------------------
template <class TImage>
irtkExtrapolateImageFunction *
irtkGenericInterpolateImageFunction<TImage>
::New(irtkExtrapolationMode mode, const irtkBaseImage *image)
{
  const TImage *img = NULL;
  if (image) {
    img = dynamic_cast<const TImage *>(image);
    if (!img) {
      cerr << this->NameOfClass() << "::New(irtkExtrapolationMode): Invalid input image type" << endl;
      exit(1);
    }
  }
  return ExtrapolatorType::New(mode, img);
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
inline void irtkGenericInterpolateImageFunction<TImage>::Input(const irtkBaseImage *input)
{
  irtkInterpolateImageFunction::Input(dynamic_cast<const TImage *>(input));
  if (input && !this->_input) {
    cerr << this->NameOfClass() << "::Input: Invalid input image type" << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
inline void irtkGenericInterpolateImageFunction<TImage>::Initialize(bool coeff)
{
  // Ensure that input has the right type
  if (!this->_input) {
    cerr << this->NameOfClass() << "::Initialize: Missing input image" << endl;
    exit(1);
  } else if (!dynamic_cast<const TImage *>(this->_input)) {
    cerr << this->NameOfClass() << "::Initialize: Invalid input image type" << endl;
    exit(1);
  }
  // Initialize extrapolator, i.e., infinite discrete image
  if (this->_InfiniteInput) {
    if (!dynamic_cast<const ExtrapolatorType *>(this->_InfiniteInput)) {
      cerr << this->NameOfClass() << "::Initialize: Invalid extrapolator type" << endl;
      exit(1);
    }
  }
  // Initialize base class
  irtkInterpolateImageFunction::Initialize(coeff);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline void irtkGenericInterpolateImageFunction<TImage>
::Extrapolator(irtkExtrapolateImageFunction *input, bool owner)
{
  irtkInterpolateImageFunction::Extrapolator(input, owner);
  if (input && !dynamic_cast<const ExtrapolatorType *>(this->_InfiniteInput)) {
    cerr << this->NameOfClass() << "::Extrapolator: Invalid extrapolator type" << endl;
    exit(1);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Instantiation
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
template <class TImage>
irtkGenericInterpolateImageFunction<TImage> *
irtkGenericInterpolateImageFunction<TImage>
::New(irtkInterpolationMode mode, const TImage *image)
{
  mode = InterpolationWithoutPadding(mode);

  irtkGenericInterpolateImageFunction<TImage> *p = NULL;

  int dim = 0;
  if (image) {
    if      (image->Z() == 1)                            dim = 2;
    else if (image->T() == 1 || image->GetTSize() == .0) dim = 3;
    else                                                 dim = 4;
  }

  switch (dim) {
    case 2: {
      switch (mode) {
        case Interpolation_NN:
          p = new irtkGenericNearestNeighborInterpolateImageFunction<TImage>();
          break;
        case Interpolation_Linear:
        case Interpolation_FastLinear:
          p = new irtkGenericLinearInterpolateImageFunction2D<TImage>();
          break;
        case Interpolation_BSpline:
          p = new irtkGenericBSplineInterpolateImageFunction2D<TImage>();
          break;
        case Interpolation_CubicBSpline:
          p = new irtkGenericCubicBSplineInterpolateImageFunction2D<TImage>();
          break;
        case Interpolation_FastCubicBSpline:
          p = new irtkGenericFastCubicBSplineInterpolateImageFunction2D<TImage>();
          break;
        case Interpolation_CSpline:
          p = new irtkGenericCSplineInterpolateImageFunction2D<TImage>();
          break;
        case Interpolation_Gaussian:
          p = new irtkGenericGaussianInterpolateImageFunction2D<TImage>();
          break;
        case Interpolation_Sinc:
          p = new irtkGenericSincInterpolateImageFunction2D<TImage>();
          break;
        default:
          p = NULL;
      }
    }
    case 3: {
      switch (mode) {
        case Interpolation_NN:
          p = new irtkGenericNearestNeighborInterpolateImageFunction<TImage>();
          break;
        case Interpolation_Linear:
        case Interpolation_FastLinear:
          p = new irtkGenericLinearInterpolateImageFunction3D<TImage>();
          break;
        case Interpolation_BSpline:
          p = new irtkGenericBSplineInterpolateImageFunction3D<TImage>();
          break;
        case Interpolation_CubicBSpline:
          p = new irtkGenericCubicBSplineInterpolateImageFunction3D<TImage>();
          break;
        case Interpolation_FastCubicBSpline:
          p = new irtkGenericFastCubicBSplineInterpolateImageFunction3D<TImage>();
          break;
        case Interpolation_CSpline:
          p = new irtkGenericCSplineInterpolateImageFunction3D<TImage>();
          break;
        case Interpolation_Gaussian:
          p = new irtkGenericGaussianInterpolateImageFunction3D<TImage>();
          break;
        case Interpolation_Sinc:
          p = new irtkGenericSincInterpolateImageFunction3D<TImage>();
          break;
        default:
          p = NULL;
      }
    }
    case 4: {
      switch (mode) {
        case Interpolation_NN:
          p = new irtkGenericNearestNeighborInterpolateImageFunction<TImage>();
          break;
        case Interpolation_Linear:
        case Interpolation_FastLinear:
          p = new irtkGenericLinearInterpolateImageFunction4D<TImage>();
          break;
        case Interpolation_BSpline:
          p = new irtkGenericBSplineInterpolateImageFunction4D<TImage>();
          break;
        case Interpolation_CubicBSpline:
          p = new irtkGenericCubicBSplineInterpolateImageFunction4D<TImage>();
          break;
        case Interpolation_FastCubicBSpline:
          p = new irtkGenericFastCubicBSplineInterpolateImageFunction4D<TImage>();
          break;
        case Interpolation_CSpline:
          p = new irtkGenericCSplineInterpolateImageFunction4D<TImage>();
          break;
        case Interpolation_Gaussian:
          p = new irtkGenericGaussianInterpolateImageFunction4D<TImage>();
          break;
        case Interpolation_Sinc:
          p = new irtkGenericSincInterpolateImageFunction4D<TImage>();
          break;
        default:
          p = NULL;
      }
    }
    default: {
      switch (mode) {
        case Interpolation_NN:
          p = new irtkGenericNearestNeighborInterpolateImageFunction<TImage>();
          break;
        case Interpolation_Linear:
        case Interpolation_FastLinear:
          p = new irtkGenericLinearInterpolateImageFunction<TImage>();
          break;
        case Interpolation_BSpline:
          p = new irtkGenericBSplineInterpolateImageFunction<TImage>();
          break;
        case Interpolation_CubicBSpline:
          p = new irtkGenericCubicBSplineInterpolateImageFunction<TImage>();
          break;
        case Interpolation_FastCubicBSpline:
          p = new irtkGenericFastCubicBSplineInterpolateImageFunction<TImage>();
          break;
        case Interpolation_CSpline:
          p = new irtkGenericCSplineInterpolateImageFunction<TImage>();
          break;
        case Interpolation_Gaussian:
          p = new irtkGenericGaussianInterpolateImageFunction<TImage>();
          break;
        case Interpolation_Sinc:
          p = new irtkGenericSincInterpolateImageFunction<TImage>();
          break;
        default:
          p = NULL;
      }
    }
  }

  // Initialize interpolator
  if (p) {
    p->NumberOfDimensions(dim);
    p->Input(image);
  // Throw error if no suitable interpolator available
  } else {
    cerr << "irtkGenericInterpolateImageFunction::New: Interpolation mode (" << mode;
    cerr << ") not supported for " << (dim ? ToString(dim) : "N") << "D images" << endl;
    exit(1);
  }

  return p;
}

// -----------------------------------------------------------------------------
template <class TImage>
irtkGenericInterpolateImageFunction<TImage> *
irtkGenericInterpolateImageFunction<TImage>
::New(irtkInterpolationMode imode, irtkExtrapolationMode emode, const TImage *image)
{
  irtkGenericInterpolateImageFunction<TImage> *p;
  p = irtkGenericInterpolateImageFunction<TImage>::New(imode, image);
  if (emode != Extrapolation_Default) p->Extrapolator(p->New(emode, image), true);
  return p;
}

// -----------------------------------------------------------------------------
// Explicit instantiation macro
#define irtkInterpolatorInstantiations(clsname)                                \
  template class clsname<irtkBaseImage>;                                       \
  template class clsname<irtkGenericImage<irtkGreyPixel> >;                    \
  template class clsname<irtkGenericImage<float> >;                            \
  template class clsname<irtkGenericImage<float2> >;                           \
  template class clsname<irtkGenericImage<float3> >;                           \
  template class clsname<irtkGenericImage<irtkFloat3> >;                       \
  template class clsname<irtkGenericImage<float3x3> >;                         \
  template class clsname<irtkGenericImage<double> >;                           \
  template class clsname<irtkGenericImage<double2> >;                          \
  template class clsname<irtkGenericImage<double3> >;                          \
  template class clsname<irtkGenericImage<irtkDouble3> >;                      \
  template class clsname<irtkGenericImage<double3x3> >


#endif
