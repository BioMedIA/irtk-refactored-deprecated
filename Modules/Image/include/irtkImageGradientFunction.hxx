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

#ifndef _IRTKIMAGEGRADIENTFUNCTION_HXX
#define _IRTKIMAGEGRADIENTFUNCTION_HXX

#include <irtkImageGradientFunction.h>


////////////////////////////////////////////////////////////////////////////////
// irtkGenericInterpolateImageFunction
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
inline irtkGenericImageGradientFunction<TImage>
::irtkGenericImageGradientFunction()
{
}

// -----------------------------------------------------------------------------
template <class TImage>
inline irtkGenericImageGradientFunction<TImage>
::~irtkGenericImageGradientFunction()
{
}

// -----------------------------------------------------------------------------
template <class TImage>
irtkExtrapolateImageFunction *
irtkGenericImageGradientFunction<TImage>
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
inline void irtkGenericImageGradientFunction<TImage>::Input(const irtkBaseImage *input)
{
  irtkImageGradientFunction::Input(dynamic_cast<const TImage *>(input));
  if (input && !this->_Input) {
    cerr << this->NameOfClass() << "::Input: Invalid input image type" << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
inline void irtkGenericImageGradientFunction<TImage>::Initialize(bool coeff)
{
  // Ensure that input has the right type
  if (!this->_Input) {
    cerr << this->NameOfClass() << "::Initialize: Missing input image" << endl;
    exit(1);
  } else if (!dynamic_cast<const TImage *>(this->_Input)) {
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
  irtkImageGradientFunction::Initialize(coeff);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline void irtkGenericImageGradientFunction<TImage>
::Extrapolator(irtkExtrapolateImageFunction *input, bool owner)
{
  irtkImageGradientFunction::Extrapolator(input, owner);
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
irtkGenericImageGradientFunction<TImage> *
irtkGenericImageGradientFunction<TImage>
::New(irtkInterpolationMode mode, const TImage *image)
{
  mode = InterpolationWithoutPadding(mode);

  irtkGenericImageGradientFunction<TImage> *p = NULL;

  int dim = 0;
  if (image) {
    if      (image->Z() == 1)                            dim = 2;
    else if (image->T() == 1 || image->GetTSize() == .0) dim = 3;
    else                                                 dim = 4;
  }

  switch (dim) {
    case 2: {
      switch (mode) {
        case Interpolation_Linear:
          p = new irtkGenericLinearImageGradientFunction2D<TImage>();
          break;
        case Interpolation_FastLinear:
          p = new irtkGenericFastLinearImageGradientFunction2D<TImage>();
          break;
        default:
          p = NULL;
      }
    }
    case 3: {
      switch (mode) {
        case Interpolation_Linear:
          p = new irtkGenericLinearImageGradientFunction3D<TImage>();
          break;
        case Interpolation_FastLinear:
          p = new irtkGenericFastLinearImageGradientFunction3D<TImage>();
          break;
        default:
          p = NULL;
      }
    }
    default: {
      switch (mode) {
        case Interpolation_Linear:
          p = new irtkGenericLinearImageGradientFunction<TImage>();
          break;
        case Interpolation_FastLinear:
          p = new irtkGenericFastLinearImageGradientFunction<TImage>();
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
    cerr << "irtkGenericImageGradientFunction::New: Interpolation mode (" << mode;
    cerr << ") not supported for " << (dim ? ToString(dim) : "N") << "D images" << endl;
    exit(1);
  }

  return p;
}

// -----------------------------------------------------------------------------
template <class TImage>
irtkGenericImageGradientFunction<TImage> *
irtkGenericImageGradientFunction<TImage>
::New(irtkInterpolationMode imode, irtkExtrapolationMode emode, const TImage *image)
{
  irtkGenericImageGradientFunction<TImage> *p;
  p = irtkGenericImageGradientFunction<TImage>::New(imode, image);
  if (emode != Extrapolation_Default) p->Extrapolator(p->New(emode, image), true);
  return p;
}

// -----------------------------------------------------------------------------
// Explicit instantiation macro
#define irtkGradientInterpolatorInstantiations(clsname)                        \
  template class clsname<irtkBaseImage>;                                       \
  template class clsname<irtkGenericImage<irtkGreyPixel> >;                    \
  template class clsname<irtkGenericImage<float> >;                            \
  template class clsname<irtkGenericImage<double> >


#endif
