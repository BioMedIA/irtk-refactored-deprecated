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

#include <irtkImageFunction.h>
#include <irtkInterpolateImageFunction.hxx>

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkInterpolateImageFunction::irtkInterpolateImageFunction()
:
  _NumberOfDimensions(0),
  _InfiniteInput     (NULL),
  _InfiniteInputOwner(false),
  _x1(.0), _y1(.0), _z1(.0), _t1(.0),
  _x2(.0), _y2(.0), _z2(.0), _t2(.0)
{
}

// -----------------------------------------------------------------------------
irtkInterpolateImageFunction::~irtkInterpolateImageFunction()
{
  if (_InfiniteInputOwner) Delete(_InfiniteInput);
}

// -----------------------------------------------------------------------------
template <class TImage>
irtkInterpolateImageFunction *NewInterpolator(irtkInterpolationMode mode, int dim = 0)
{
  mode = InterpolationWithoutPadding(mode);
  switch (dim) {
    case 2: {
      switch (mode) {
        case Interpolation_NN:               return new irtkGenericNearestNeighborInterpolateImageFunction<TImage>();
        case Interpolation_Linear:           return new irtkGenericLinearInterpolateImageFunction2D<TImage>();
        case Interpolation_FastLinear:       return new irtkGenericLinearInterpolateImageFunction2D<TImage>();
        case Interpolation_BSpline:          return new irtkGenericBSplineInterpolateImageFunction2D<TImage>();
        case Interpolation_CubicBSpline:     return new irtkGenericCubicBSplineInterpolateImageFunction2D<TImage>();
        case Interpolation_FastCubicBSpline: return new irtkGenericFastCubicBSplineInterpolateImageFunction2D<TImage>();
        case Interpolation_CSpline:          return new irtkGenericCSplineInterpolateImageFunction2D<TImage>();
        case Interpolation_Gaussian:         return new irtkGenericGaussianInterpolateImageFunction2D<TImage>();
        case Interpolation_Sinc:             return new irtkGenericSincInterpolateImageFunction2D<TImage>();
        default:                             return NULL;
      }
    }
    case 3: {
      switch (mode) {
        case Interpolation_NN:               return new irtkGenericNearestNeighborInterpolateImageFunction<TImage>();
        case Interpolation_Linear:           return new irtkGenericLinearInterpolateImageFunction3D<TImage>();
        case Interpolation_FastLinear:       return new irtkGenericLinearInterpolateImageFunction3D<TImage>();
        case Interpolation_BSpline:          return new irtkGenericBSplineInterpolateImageFunction3D<TImage>();
        case Interpolation_CubicBSpline:     return new irtkGenericCubicBSplineInterpolateImageFunction3D<TImage>();
        case Interpolation_FastCubicBSpline: return new irtkGenericFastCubicBSplineInterpolateImageFunction3D<TImage>();
        case Interpolation_CSpline:          return new irtkGenericCSplineInterpolateImageFunction3D<TImage>();
        case Interpolation_Gaussian:         return new irtkGenericGaussianInterpolateImageFunction3D<TImage>();
        case Interpolation_Sinc:             return new irtkGenericSincInterpolateImageFunction3D<TImage>();
        default:                             return NULL;
      }
    }
    case 4: {
      switch (mode) {
        case Interpolation_NN:               return new irtkGenericNearestNeighborInterpolateImageFunction<TImage>();
        case Interpolation_Linear:           return new irtkGenericLinearInterpolateImageFunction4D<TImage>();
        case Interpolation_FastLinear:       return new irtkGenericLinearInterpolateImageFunction4D<TImage>();
        case Interpolation_BSpline:          return new irtkGenericBSplineInterpolateImageFunction4D<TImage>();
        case Interpolation_CubicBSpline:     return new irtkGenericCubicBSplineInterpolateImageFunction4D<TImage>();
        case Interpolation_FastCubicBSpline: return new irtkGenericFastCubicBSplineInterpolateImageFunction4D<TImage>();
        case Interpolation_CSpline:          return new irtkGenericCSplineInterpolateImageFunction4D<TImage>();
        case Interpolation_Gaussian:         return new irtkGenericGaussianInterpolateImageFunction4D<TImage>();
        case Interpolation_Sinc:             return new irtkGenericSincInterpolateImageFunction4D<TImage>();
        default:                             return NULL;
      }
    }
    default: {
      switch (mode) {
        case Interpolation_NN:               return new irtkGenericNearestNeighborInterpolateImageFunction<TImage>();
        case Interpolation_Linear:           return new irtkGenericLinearInterpolateImageFunction<TImage>();
        case Interpolation_FastLinear:       return new irtkGenericLinearInterpolateImageFunction<TImage>();
        case Interpolation_BSpline:          return new irtkGenericBSplineInterpolateImageFunction<TImage>();
        case Interpolation_CubicBSpline:     return new irtkGenericCubicBSplineInterpolateImageFunction<TImage>();
        case Interpolation_FastCubicBSpline: return new irtkGenericFastCubicBSplineInterpolateImageFunction<TImage>();
        case Interpolation_CSpline:          return new irtkGenericCSplineInterpolateImageFunction<TImage>();
        case Interpolation_Gaussian:         return new irtkGenericGaussianInterpolateImageFunction<TImage>();
        case Interpolation_Sinc:             return new irtkGenericSincInterpolateImageFunction<TImage>();
        default:                             return NULL;
      }
    }
  }
}

// -----------------------------------------------------------------------------
irtkInterpolateImageFunction *
irtkInterpolateImageFunction::New(irtkInterpolationMode mode, const irtkBaseImage *image)
{
  mode = InterpolationWithoutPadding(mode);
  irtkInterpolateImageFunction *p = NULL;

  // Dimensionality of image to interpolate
  int dim  = 0;
  if (image) {
    if      (image->Z() == 1)                            dim = 2;
    else if (image->T() == 1 || image->GetTSize() == .0) dim = 3;
    else                                                 dim = 4;
  }
  // Instantiate special purpose interpolators
  if (mode == Interpolation_SBased) {
    // Only implemented for 3D scalar images
    if (dim == 3 && (!image || image->N() == 1)) {
      p = new irtkShapeBasedInterpolateImageFunction();
    }
  }
  // Instantiate interpolator for generic image (i.e., instance of irtkGenericImage)
  if (!p && image) {
    typedef irtkGenericImage<char>           CharImage;
    typedef irtkGenericImage<unsigned char>  UCharImage;
    typedef irtkGenericImage<short>          ShortImage;
    typedef irtkGenericImage<unsigned short> UShortImage;
    typedef irtkGenericImage<int>            IntImage;
    typedef irtkGenericImage<unsigned int>   UIntImage;
    typedef irtkGenericImage<float>          FloatImage;
    typedef irtkGenericImage<double>         DoubleImage;

    if      (dynamic_cast<const CharImage   *>(image)) p = NewInterpolator<CharImage>  (mode, dim);
    else if (dynamic_cast<const UCharImage  *>(image)) p = NewInterpolator<UCharImage> (mode, dim);
    else if (dynamic_cast<const ShortImage  *>(image)) p = NewInterpolator<ShortImage> (mode, dim);
    else if (dynamic_cast<const UShortImage *>(image)) p = NewInterpolator<UShortImage>(mode, dim);
    else if (dynamic_cast<const IntImage    *>(image)) p = NewInterpolator<IntImage>   (mode, dim);
    else if (dynamic_cast<const UIntImage   *>(image)) p = NewInterpolator<UIntImage>  (mode, dim);
    else if (dynamic_cast<const FloatImage  *>(image)) p = NewInterpolator<FloatImage> (mode, dim);
    else if (dynamic_cast<const DoubleImage *>(image)) p = NewInterpolator<DoubleImage>(mode, dim);
  }
  // Instantiate interpolator for general image (i.e., subclass of irtkBaseImage)
  if (!p) p = NewInterpolator<irtkBaseImage>(mode, dim);
  // Initialize interpolator
  if (p) {
    p->NumberOfDimensions(dim);
    p->Input(image);
  // Throw error if no suitable interpolator available
  } else {
    cerr << "irtkInterpolateImageFunction::New: Interpolation mode (" << mode;
    cerr << ") not supported for " << (dim ? ToString(dim) : "N") << "D images" << endl;
    exit(1);
  }

  return p;
}

// -----------------------------------------------------------------------------
irtkExtrapolateImageFunction *
irtkInterpolateImageFunction::New(irtkExtrapolationMode mode, const irtkBaseImage *image)
{
  return irtkExtrapolateImageFunction::New(mode, image);
}

// -----------------------------------------------------------------------------
irtkInterpolateImageFunction *
irtkInterpolateImageFunction::New(irtkInterpolationMode imode,
                                  irtkExtrapolationMode emode, const irtkBaseImage *image)
{
  irtkInterpolateImageFunction *p = irtkInterpolateImageFunction::New(imode, image);
  if (emode != Extrapolation_Default) p->Extrapolator(p->New(emode, image), true);
  return p;
}

// -----------------------------------------------------------------------------
void irtkInterpolateImageFunction::Initialize(bool)
{
  // Initialize image function
  irtkImageFunction::Initialize();

  // Check if input is a valid image
  if (Input()->IsEmpty()) {
    cerr << this->NameOfClass() << "::Initialize: Input image has zero extent" << endl;
    exit(1);
  }

  // Determine dimensionality of input (if not specified by subclass/New)
  if (_NumberOfDimensions == 0) {
    if (Input()->Z() > 1) {
      if (Input()->T() > 1 && Input()->GetTSize() != .0) _NumberOfDimensions = 4;
      else                                               _NumberOfDimensions = 3;
    } else                                               _NumberOfDimensions = 2;
  }

  // Default domain within which interpolation can be performed is assumed
  // to be identical to the entire finite image domain
  _x1 = .0;
  _y1 = .0;
  _z1 = .0;
  _t1 = .0;
  _x2 = Input()->X() - 1;
  _y2 = Input()->Y() - 1;
  _z2 = Input()->Z() - 1;
  _t2 = Input()->T() - 1;

  // Initialize extrapolator, i.e., infinite discrete image
  if (_InfiniteInput) {
    _InfiniteInput->Input(this->Input());
    _InfiniteInput->Initialize();
  }
}

////////////////////////////////////////////////////////////////////////////////
// Definitions of interpolate image functions
////////////////////////////////////////////////////////////////////////////////

// ND
#include <irtkNearestNeighborInterpolateImageFunction.hxx>
#include <irtkLinearInterpolateImageFunction.hxx>
#include <irtkBSplineInterpolateImageFunction.hxx>
#include <irtkCubicBSplineInterpolateImageFunction.hxx>
#include <irtkFastCubicBSplineInterpolateImageFunction.hxx>
#include <irtkCSplineInterpolateImageFunction.hxx>
#include <irtkGaussianInterpolateImageFunction.hxx>
#include <irtkSincInterpolateImageFunction.hxx>

// 2D
#include <irtkLinearInterpolateImageFunction2D.hxx>
#include <irtkBSplineInterpolateImageFunction2D.hxx>
#include <irtkCubicBSplineInterpolateImageFunction2D.hxx>
#include <irtkFastCubicBSplineInterpolateImageFunction2D.hxx>
#include <irtkCSplineInterpolateImageFunction2D.hxx>
#include <irtkGaussianInterpolateImageFunction2D.hxx>
#include <irtkSincInterpolateImageFunction2D.hxx>

// 3D
#include <irtkLinearInterpolateImageFunction3D.hxx>
#include <irtkBSplineInterpolateImageFunction3D.hxx>
#include <irtkCubicBSplineInterpolateImageFunction3D.hxx>
#include <irtkFastCubicBSplineInterpolateImageFunction3D.hxx>
#include <irtkCSplineInterpolateImageFunction3D.hxx>
#include <irtkGaussianInterpolateImageFunction3D.hxx>
#include <irtkSincInterpolateImageFunction3D.hxx>

// 4D
#include <irtkLinearInterpolateImageFunction4D.hxx>
#include <irtkBSplineInterpolateImageFunction4D.hxx>
#include <irtkCubicBSplineInterpolateImageFunction4D.hxx>
#include <irtkFastCubicBSplineInterpolateImageFunction4D.hxx>
#include <irtkCSplineInterpolateImageFunction4D.hxx>
#include <irtkGaussianInterpolateImageFunction4D.hxx>
#include <irtkSincInterpolateImageFunction4D.hxx>

////////////////////////////////////////////////////////////////////////////////
// Explicit instantiations
////////////////////////////////////////////////////////////////////////////////

// Base class
irtkInterpolatorInstantiations(irtkGenericInterpolateImageFunction);

// ND
template class irtkGenericNearestNeighborInterpolateImageFunction<irtkByteImage>;
irtkInterpolatorInstantiations(irtkGenericNearestNeighborInterpolateImageFunction);
irtkInterpolatorInstantiations(irtkGenericLinearInterpolateImageFunction);
irtkInterpolatorInstantiations(irtkGenericBSplineInterpolateImageFunction);
irtkInterpolatorInstantiations(irtkGenericCubicBSplineInterpolateImageFunction);
irtkInterpolatorInstantiations(irtkGenericFastCubicBSplineInterpolateImageFunction);
irtkInterpolatorInstantiations(irtkGenericCSplineInterpolateImageFunction);
irtkInterpolatorInstantiations(irtkGenericGaussianInterpolateImageFunction);
irtkInterpolatorInstantiations(irtkGenericSincInterpolateImageFunction);

// 2D
irtkInterpolatorInstantiations(irtkGenericLinearInterpolateImageFunction2D);
irtkInterpolatorInstantiations(irtkGenericBSplineInterpolateImageFunction2D);
irtkInterpolatorInstantiations(irtkGenericCubicBSplineInterpolateImageFunction2D);
irtkInterpolatorInstantiations(irtkGenericFastCubicBSplineInterpolateImageFunction2D);
irtkInterpolatorInstantiations(irtkGenericCSplineInterpolateImageFunction2D);
irtkInterpolatorInstantiations(irtkGenericGaussianInterpolateImageFunction2D);
irtkInterpolatorInstantiations(irtkGenericSincInterpolateImageFunction2D);

// 3D
irtkInterpolatorInstantiations(irtkGenericLinearInterpolateImageFunction3D);
irtkInterpolatorInstantiations(irtkGenericBSplineInterpolateImageFunction3D);
irtkInterpolatorInstantiations(irtkGenericCubicBSplineInterpolateImageFunction3D);
irtkInterpolatorInstantiations(irtkGenericFastCubicBSplineInterpolateImageFunction3D);
irtkInterpolatorInstantiations(irtkGenericCSplineInterpolateImageFunction3D);
irtkInterpolatorInstantiations(irtkGenericGaussianInterpolateImageFunction3D);
irtkInterpolatorInstantiations(irtkGenericSincInterpolateImageFunction3D);

// 4D
irtkInterpolatorInstantiations(irtkGenericLinearInterpolateImageFunction4D);
irtkInterpolatorInstantiations(irtkGenericBSplineInterpolateImageFunction4D);
irtkInterpolatorInstantiations(irtkGenericCubicBSplineInterpolateImageFunction4D);
irtkInterpolatorInstantiations(irtkGenericFastCubicBSplineInterpolateImageFunction4D);
irtkInterpolatorInstantiations(irtkGenericCSplineInterpolateImageFunction4D);
irtkInterpolatorInstantiations(irtkGenericGaussianInterpolateImageFunction4D);
irtkInterpolatorInstantiations(irtkGenericSincInterpolateImageFunction4D);

