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

#include <irtkImageGradientFunction.hxx>


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkImageGradientFunction::irtkImageGradientFunction()
:
  _WrtWorld          (false),
  _DefaultValue      (.0),
  _NumberOfDimensions(0),
  _Input             (NULL),
  _Orientation       (3, 3),
  _InfiniteInput     (NULL),
  _InfiniteInputOwner(false),
  _x1(.0), _y1(.0), _z1(.0), _t1(.0),
  _x2(.0), _y2(.0), _z2(.0), _t2(.0)
{
}

// -----------------------------------------------------------------------------
irtkImageGradientFunction::~irtkImageGradientFunction()
{
  if (_InfiniteInputOwner) Delete(_InfiniteInput);
}

// -----------------------------------------------------------------------------
template <class TImage>
irtkImageGradientFunction *NewInterpolator(irtkInterpolationMode mode, int dim = 0)
{
  mode = InterpolationWithoutPadding(mode);
  switch (dim) {
    case 2: {
      switch (mode) {
        case Interpolation_Linear:     return new irtkGenericLinearImageGradientFunction2D<TImage>();
        case Interpolation_FastLinear: return new irtkGenericFastLinearImageGradientFunction2D<TImage>();
        default:                       return NULL;
      }
    }
    case 3: {
      switch (mode) {
        case Interpolation_Linear:     return new irtkGenericLinearImageGradientFunction3D<TImage>();
        case Interpolation_FastLinear: return new irtkGenericFastLinearImageGradientFunction3D<TImage>();
        default:                       return NULL;
      }
    }
    default: {
      switch (mode) {
        case Interpolation_Linear:     return new irtkGenericLinearImageGradientFunction<TImage>();
        case Interpolation_FastLinear: return new irtkGenericFastLinearImageGradientFunction<TImage>();
        default:                       return NULL;
      }
    }
  }
}

// -----------------------------------------------------------------------------
irtkImageGradientFunction *
irtkImageGradientFunction::New(irtkInterpolationMode mode, const irtkBaseImage *image)
{
  irtkImageGradientFunction *p = NULL;

  // Dimensionality of image to interpolate
  int dim  = 0;
  if (image) {
    if      (image->Z() == 1)                            dim = 2;
    else if (image->T() == 1 || image->GetTSize() == .0) dim = 3;
    else                                                 dim = 4;
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
    cerr << "irtkImageGradientFunction::New: Interpolation mode (" << mode;
    cerr << ") not supported for " << (dim ? ToString(dim) : "N") << "D images" << endl;
    exit(1);
  }

  return p;
}

// -----------------------------------------------------------------------------
irtkExtrapolateImageFunction *
irtkImageGradientFunction::New(irtkExtrapolationMode mode, const irtkBaseImage *image)
{
  return irtkExtrapolateImageFunction::New(mode, image);
}

// -----------------------------------------------------------------------------
irtkImageGradientFunction *
irtkImageGradientFunction::New(irtkInterpolationMode imode,
                               irtkExtrapolationMode emode, const irtkBaseImage *image)
{
  irtkImageGradientFunction *p = irtkImageGradientFunction::New(imode, image);
  if (emode != Extrapolation_Default) p->Extrapolator(p->New(emode, image), true);
  return p;
}

// -----------------------------------------------------------------------------
void irtkImageGradientFunction::Initialize(bool)
{
  // Check if input is a valid image
  if (!_Input) {
    cerr << this->NameOfClass() << "::Initialize: Input image not set" << endl;
    exit(1);
  }
  if (_Input->IsEmpty()) {
    cerr << this->NameOfClass() << "::Initialize: Input image has zero extent" << endl;
    exit(1);
  }

  // Image resolution and orientation matrix
  _VoxelSize._x = (_Input->GetXSize() == .0 ? 1.0 : _Input->GetXSize());
  _VoxelSize._y = (_Input->GetYSize() == .0 ? 1.0 : _Input->GetYSize());
  _VoxelSize._z = (_Input->GetZSize() == .0 ? 1.0 : _Input->GetZSize());
  _Orientation  = _Input->Attributes().GetWorldToImageOrientation();

  // Determine dimensionality of input (if not specified by subclass/New)
  if (_NumberOfDimensions == 0) {
    if (_Input->Z() > 1) {
      if (_Input->T() > 1 && _Input->GetTSize() != .0) _NumberOfDimensions = 4;
      else                                             _NumberOfDimensions = 3;
    } else                                             _NumberOfDimensions = 2;
  }

  // Default domain within which interpolation can be performed is assumed
  // to be identical to the entire finite image domain
  _x1 = .0;
  _y1 = .0;
  _z1 = .0;
  _t1 = .0;
  _x2 = _Input->X() - 1;
  _y2 = _Input->Y() - 1;
  _z2 = _Input->Z() - 1;
  _t2 = _Input->T() - 1;

  // Initialize extrapolator, i.e., infinite discrete image
  if (_InfiniteInput) {
    _InfiniteInput->Input(_Input);
    _InfiniteInput->Initialize();
  }
}

////////////////////////////////////////////////////////////////////////////////
// Definitions of image gradient functions
////////////////////////////////////////////////////////////////////////////////

// ND
#include <irtkLinearImageGradientFunction.hxx>
#include <irtkFastLinearImageGradientFunction.hxx>

// 2D
#include <irtkLinearImageGradientFunction2D.hxx>
#include <irtkFastLinearImageGradientFunction2D.hxx>

// 3D
#include <irtkLinearImageGradientFunction3D.hxx>
#include <irtkFastLinearImageGradientFunction3D.hxx>

////////////////////////////////////////////////////////////////////////////////
// Explicit instantiations
////////////////////////////////////////////////////////////////////////////////

// Base class
irtkGradientInterpolatorInstantiations(irtkGenericImageGradientFunction);

// ND
irtkGradientInterpolatorInstantiations(irtkGenericLinearImageGradientFunction);
irtkGradientInterpolatorInstantiations(irtkGenericFastLinearImageGradientFunction);

// 2D
irtkGradientInterpolatorInstantiations(irtkGenericLinearImageGradientFunction2D);
irtkGradientInterpolatorInstantiations(irtkGenericFastLinearImageGradientFunction2D);

// 3D
irtkGradientInterpolatorInstantiations(irtkGenericLinearImageGradientFunction3D);
irtkGradientInterpolatorInstantiations(irtkGenericFastLinearImageGradientFunction3D);
