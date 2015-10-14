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

#ifndef _IRTKLIEBRACKETIMAGEFILTER_H

#define _IRTKLIEBRACKETIMAGEFILTER_H


#include <irtkImage.h>
#include <irtkImageToImage.h>


/**
 * Base class for image filters which compute the Lie bracket of two vector fields.
 */

template <class VoxelType>
class irtkLieBracketImageFilter : public irtkImageToImage<VoxelType>
{
  irtkImageFilterMacro(irtkLieBracketImageFilter);

protected:
  using irtkImageToImage<VoxelType>::SetInput;

  /// Second input vector field
  irtkGenericImage<VoxelType> *_input2;

  /// Constructor
  irtkLieBracketImageFilter();

  /// Initialize filter
  virtual void Initialize();

public:

  /// Construct Lie bracket filter for given image domain
  static irtkLieBracketImageFilter *New(const irtkImageAttributes &, bool = true);

  /// Construct Lie bracket filter for given input vector field
  static irtkLieBracketImageFilter *New(const irtkImage *, bool = true);

  /// Destructor
  virtual ~irtkLieBracketImageFilter();

  /// Set n-th input
  virtual void SetInput(int, irtkGenericImage<VoxelType> *);

  /// Get n-th input
  virtual irtkGenericImage<VoxelType> *GetInput(int);

};

///////////////////////////////////////////////////////////////////////////////
// Inline definitions
///////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
template <class VoxelType>
irtkLieBracketImageFilter<VoxelType>::irtkLieBracketImageFilter()
:
  _input2(NULL)
{
}

// ----------------------------------------------------------------------------
template <class VoxelType>
irtkLieBracketImageFilter<VoxelType>::~irtkLieBracketImageFilter()
{
}

// --------------------------------------------------------------------------
template <class VoxelType>
void irtkLieBracketImageFilter<VoxelType>::SetInput(int i, irtkGenericImage<VoxelType> *image)
{
  if      (i == 0) irtkImageToImage<VoxelType>::SetInput(image);
  else if (i == 1) _input2 = image;
  else {
    cerr << this->NameOfClass() << "::SetInput: Input index out of range: " << i << endl;
    exit(1);
  }
}

// --------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType> *irtkLieBracketImageFilter<VoxelType>::GetInput(int i)
{
  if      (i == 0) return this->_input;
  else if (i == 1) return _input2;
  else {
    cerr << this->NameOfClass() << "::GetInput: Input index out of range: " << i << endl;
    exit(1);
  }
}

// --------------------------------------------------------------------------
template <class VoxelType>
void irtkLieBracketImageFilter<VoxelType>::Initialize()
{
  // Initialize base class
  irtkImageToImage<VoxelType>::Initialize();

  // Check second input
  if (_input2 == NULL) {
    cerr << this->NameOfClass() << "::Initialize: Filter has no second input" << endl;
    exit(1);
  }
  if (!_input2->HasSpatialAttributesOf(this->_input) || _input2->GetT() != this->_input->GetT()) {
    cerr << this->NameOfClass() << "::Initialize: Attributes of input images do not match" << endl;
    exit(1);
  }
}


#endif // _IRTKLIEBRACKETIMAGEFILTER_H

///////////////////////////////////////////////////////////////////////////////
// Filter implementations
///////////////////////////////////////////////////////////////////////////////

#if USE_CUDA
#  include <irtkCUImage.h>
namespace irtkCULieBracketImageFilter {
  template <class VoxelType>
  void Run(irtkCUGenericImage<VoxelType> *, const irtkCUGenericImage<VoxelType> *,
                                            const irtkCUGenericImage<VoxelType> *);
}
#endif

#include <irtkLieBracketImageFilter2D.h>
#include <irtkLieBracketImageFilter3D.h>
#include <irtkDifferenceOfCompositionLieBracketImageFilter3D.h>

///////////////////////////////////////////////////////////////////////////////
// Instantiation of filter implementation
///////////////////////////////////////////////////////////////////////////////

#ifndef _IRTKLIEBRACKETIMAGEFILTERNEW_H
#define _IRTKLIEBRACKETIMAGEFILTERNEW_H


// ----------------------------------------------------------------------------
template <class VoxelType>
irtkLieBracketImageFilter<VoxelType> *
irtkLieBracketImageFilter<VoxelType>::New(const irtkImageAttributes &attr, bool usejac)
{
  if (attr._z > 1) {
    if (usejac) return new irtkLieBracketImageFilter3D<VoxelType>;
    else        return new irtkDifferenceOfCompositionLieBracketImageFilter3D<VoxelType>;
  } else {
    if (usejac) return new irtkLieBracketImageFilter2D<VoxelType>;
    else {
      cerr << NameOfType() << "::New: irtkDifferenceOfCompositionLieBracketImageFilter2D not implemented" << endl;
      exit(1);
    }
  }
}

// ----------------------------------------------------------------------------
template <class VoxelType>
irtkLieBracketImageFilter<VoxelType> *
irtkLieBracketImageFilter<VoxelType>::New(const irtkImage *image, bool usejac)
{
  return irtkLieBracketImageFilter::New(image->GetImageAttributes(), usejac);
}

// ----------------------------------------------------------------------------
template <class VoxelType>
void liebracket(irtkGenericImage<VoxelType> *ov,
                irtkGenericImage<VoxelType> *lv,
                irtkGenericImage<VoxelType> *rv, bool usejac = true)
{
  typedef irtkLieBracketImageFilter<VoxelType> LieBracketFilter;
  auto_ptr<LieBracketFilter> filter(LieBracketFilter::New(ov, usejac));
  filter->SetInput(0, lv);
  filter->SetInput(1, rv);
  filter->SetOutput(ov);
  filter->Run();
}


#endif // _IRTKLIEBRACKETIMAGEFILTERNEW_H
