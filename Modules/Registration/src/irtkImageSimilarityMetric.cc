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

#include <irtkImage.h>
#include <irtkImageFunction.h>
#include <irtkTransformation.h>
#include <irtkGradientImageFilter.h>
#include <irtkHistogram.h>

#include <irtkRegistration.h>           // irtkSimilarityMeasure
#include <irtkImageSimilarityMetric.h>

////////////////////////////////////////////////////////////////////////////////
// New
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
irtkImageSimilarityMetric *irtkImageSimilarityMetric::New(irtkSimilarityMeasure m)
{
  switch (m) {
    //case SSD: return new irtkSumOfSquaredDifferences();
    default:
      cerr << "irtkImageSimilarityMetric::New: Unknown/Missing similarity measure: " << m << endl;
      exit(1);
  }
}

////////////////////////////////////////////////////////////////////////////////
// BaseChannel structure
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
template <class ImageType>
irtkImageSimilarityMetric::BaseChannel<ImageType>::BaseChannel(ImageType *image, bool owner)
:
  _Image     (image),
  _Gradient  (NULL),
  _ImageOwner(owner)
{
}

// -----------------------------------------------------------------------------
template <class ImageType>
irtkImageSimilarityMetric::BaseChannel<ImageType>::BaseChannel(const BaseChannel &other)
:
  _Image     (NULL),
  _Gradient  (NULL),
  _ImageOwner(false)
{
  *this = other;
}

// -----------------------------------------------------------------------------
template <>
irtkImageSimilarityMetric::BaseChannel<irtkImage> &irtkImageSimilarityMetric::BaseChannel<irtkImage>::operator =(const BaseChannel &rhs)
{
  if (_ImageOwner) delete _Image;
  _Image    = (rhs._ImageOwner && rhs._Image) ? irtkImage::New(rhs._Image) : rhs._Image;
  if (_Gradient && rhs._Gradient) {
    *_Gradient = *rhs._Gradient;
  } else {
    delete _Gradient;
    _Gradient = (rhs._Gradient ? new Gradient(*rhs._Gradient) : NULL);
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <class ImageType>
irtkImageSimilarityMetric::BaseChannel<ImageType> &irtkImageSimilarityMetric::BaseChannel<ImageType>::operator =(const BaseChannel &rhs)
{
  if (_ImageOwner) delete _Image;
  _Image    = (rhs._Image && rhs._ImageOwner) ? new ImageType(*rhs._Image) : rhs._Image;
  if (_Gradient && rhs._Gradient) {
    *_Gradient = *rhs._Gradient;
  } else {
    delete _Gradient;
    _Gradient = (rhs._Gradient ? new Gradient(*rhs._Gradient) : NULL);
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <class ImageType>
irtkImageSimilarityMetric::BaseChannel<ImageType>::~BaseChannel()
{
  if (_ImageOwner) delete _Image;    _Image    = NULL;
                   delete _Gradient; _Gradient = NULL;
}

// -----------------------------------------------------------------------------
template <class ImageType>
void irtkImageSimilarityMetric::BaseChannel<ImageType>::ComputeGradient()
{
  typedef irtkGradientImageFilter<Real> GradientImageFilter;
  // Handle case that no image is set gracefully
  if (_Image == NULL) {
    delete _Gradient; _Gradient = NULL;
    return;
  }
  // Initialize gradient image
  if (!_Gradient) _Gradient = new Gradient(_Image->GetImageAttributes(), 3);
  // Compute gradient
  //
  // TODO: Implement gradient filter which does not require input to be
  //       of same type as output image. Moreover, consider convolution
  //       with derivative of Gaussian kernel instead. -as12312
  GradientImageFilter filter(GradientImageFilter::GRADIENT_VECTOR);
  Gradient *image = dynamic_cast<Gradient *>( _Image);
  if (!image)    image = new     Gradient   (*_Image);
  filter.SetInput(image);
  filter.SetOutput(_Gradient);
  if (image->HasBackgroundValue()) {
    filter.SetPadding(image->GetBackgroundValueAsDouble());
  }
  filter.Run();
  if (image != _Image) delete image;
}

////////////////////////////////////////////////////////////////////////////////
// Base class for intensity-based image similarity measures
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkImageSimilarityMetric::irtkImageSimilarityMetric()
:
  _InputImage2World            (NULL),
  _Image2World                 (NULL),
  _Mask                        (NULL),
  _NumberOfFrames              (0),
  _NumberOfChannels            (0),
  _NumberOfVoxels              (0),
  _CurrentFrame                (0),
  _CurrentChannel              (0),
  _CacheWorldCoordinates       (true),
  _CacheGradient               (true),
  _SourceTransformation        (1),
  _NumberOfTargetImagesToUpdate(0),
  _NumberOfSourceImagesToUpdate(0)
{
  _SourceTransformation[0].CacheWorldCoordinatesOn();
  _SourceTransformation[0].CacheDisplacementsOn();
  _TargetTransformation   .CacheWorldCoordinatesOn();
  _TargetTransformation   .CacheDisplacementsOn();
}

// -----------------------------------------------------------------------------
irtkImageSimilarityMetric::irtkImageSimilarityMetric(const irtkImageSimilarityMetric &other)
:
  _Domain                      (other._Domain),
  _InputImage2World            (other._InputImage2World),
  _Image2World                 ((other._Image2World && !other._InputImage2World)
                                  ? new irtkWorldCoordsImage(*other._Image2World)
                                  : other._Image2World),
  _Mask                        (other._Mask),
  _NumberOfFrames              (other._NumberOfFrames),
  _NumberOfChannels            (other._NumberOfChannels),
  _NumberOfVoxels              (other._NumberOfVoxels),
  _CurrentFrame                (other._CurrentFrame),
  _CurrentChannel              (other._CurrentChannel),
  _CacheWorldCoordinates       (other._CacheWorldCoordinates),
  _CacheGradient               (other._CacheGradient),
  _InputTarget                 (other._InputTarget),
  _InputSource                 (other._InputSource),
  _Target                      (other._Target),
  _Source                      (other._Source),
  _TargetTransformation        (other._TargetTransformation),
  _SourceTransformation        (other._SourceTransformation),
  _NumberOfTargetImagesToUpdate(other._NumberOfTargetImagesToUpdate),
  _NumberOfSourceImagesToUpdate(other._NumberOfSourceImagesToUpdate),
  _Gradient                    (other._Gradient)
{
}

// -----------------------------------------------------------------------------
irtkImageSimilarityMetric::~irtkImageSimilarityMetric()
{
  if (_Image2World != _InputImage2World) delete _Image2World;
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
void irtkImageSimilarityMetric::AddTarget(irtkImage *target, bool copy)
{
  if (target->GetT() > 1 && target->GetTSize() > .0) {
    cerr << "irtkImageSimilarityMetric::AddTarget: Target image cannot be image sequence" << endl;
    exit(1);
  }
  if (_InputTarget.size() > 0) {
    if (!target->HasSpatialAttributesOf(_InputTarget[0]._Image)) {
      if (target->GetT() > 0) {
        cerr << "irtkImageSimilarityMetric::AddTarget: Target channels have differing spatial domain" << endl;
      } else {
        cerr << "irtkImageSimilarityMetric::AddTarget: Target channel has differing spatial domain" << endl;
      }
      exit(1);
    }
  }
  if (copy || target->GetT() > 1) {
    for (int l = 0; l < target->GetT(); ++l) {
      irtkImage *image = NULL;
      target->GetFrame(l, image);
      _InputTarget.push_back(InputChannel(image, false));
      // set afterwards to avoid unnecessary copying in InputChannel copy constructor
      _InputTarget.back()._ImageOwner = true;
    }
  } else {
    _InputTarget.push_back(InputChannel(target, false));
  }
}

// -----------------------------------------------------------------------------
void irtkImageSimilarityMetric::SetTarget(irtkImage *target, bool copy)
{
  _InputTarget.clear();
  if (target) AddTarget(target, copy);
}

// -----------------------------------------------------------------------------
void irtkImageSimilarityMetric::SetTarget(irtkImage **target, int num, bool copy)
{
  _InputTarget.clear();
  _InputTarget.reserve(num);
  for (int n = 0; n < num; ++n) AddTarget(target[n], copy);
}

// -----------------------------------------------------------------------------
void irtkImageSimilarityMetric::SetTargetTransformation(const irtkTransformation *current,
                                                        const irtkTransformation *fixed,
                                                        CompositionMode           mode)
{
  _TargetTransformation.SetFixedTransformation(fixed);
  _TargetTransformation.SetTransformation(current, mode);
}

// -----------------------------------------------------------------------------
void irtkImageSimilarityMetric::AddSource(irtkImage *source, bool copy)
{
  if (_InputSource.size() > 0) {
    if (!source->HasSpatialAttributesOf(_InputSource[0][0]._Image)) {
      if (source->GetT() > 0 && source->GetTSize() == .0) {
        cerr << "irtkImageSimilarityMetric::AddSource: Source channels have differing spatial domain" << endl;
      } else {
        cerr << "irtkImageSimilarityMetric::AddSource: Source channel has differing spatial domain" << endl;
      }
      exit(1);
    }
  }
  // Reserve enough entries in vector to ensure that no reallocation
  // takes place during the insert to keep tmp iterator below valid
  _InputSource.reserve(_InputSource.size() + source->GetT());
  // Insert each input frame (dt > 0) or channel (dt == 0)
  // Frames are sorted by increasing time and channels appended
  for (int l = 0; l < source->GetT(); ++l) {
    // Add empty InputChannel to 2D (vector) array at proper position
    const double            time = source->ImageToTime(l);
    InputChannel           *slot = NULL;
    InputSequence::iterator pos  = _InputSource.begin();
    while (pos != _InputSource.end()) {
      const double t = pos->front()._Image->ImageToTime(.0);
      if (t == time) {                // Another channel of existing frame
        pos->push_back(InputChannel());
        slot = &pos->back();
        break;
      } else if (t > time) {          // First channel of another temporal frame
        InputSequence::iterator tmp = pos - 1;
        _InputSource.insert(pos, InputFrame(1));
        slot = &(++tmp)->back();
        break;
      }
      ++pos;
    }
    if (pos == _InputSource.end()) { // First channel of another temporal frame
      _InputSource.push_back(InputFrame(1));
      slot = &_InputSource.back().back();
    }
    // Set image of added InputChannel
    if (copy || source->GetT() > 1) {
      source->GetFrame(l, slot->_Image);
      slot->_ImageOwner = true;
    } else {
      slot->_Image      = source;
      slot->_ImageOwner = false;
      break;
    }
  }
}

// -----------------------------------------------------------------------------
void irtkImageSimilarityMetric::SetSource(irtkImage *source, bool copy)
{
  _InputSource.clear();
  if (source) AddSource(source, copy);
}

// -----------------------------------------------------------------------------
void irtkImageSimilarityMetric::SetSource(irtkImage **source, int num, bool copy)
{
  _InputSource.clear();
  for (int n = 0; n < num; ++n) AddSource(source[n], copy);
}

// -----------------------------------------------------------------------------
void irtkImageSimilarityMetric::SetSourceTransformation(const irtkTransformation *current,
                                                        const irtkTransformation *fixed,
                                                        CompositionMode           mode)
{
  _SourceTransformation[0].SetFixedTransformation(fixed);
  _SourceTransformation[0].SetTransformation(current, mode);
}

// -----------------------------------------------------------------------------
void irtkImageSimilarityMetric::SetInterpolation(irtkInterpolationMode mode)
{
  _SourceTransformation[0].SetInterpolation(mode);
  _TargetTransformation   .SetInterpolation(mode);
}

// -----------------------------------------------------------------------------
void irtkImageSimilarityMetric::Initialize(const irtkImageAttributes &domain)
{
  // Set attributes of domain within which to evaluate similarity
  _Domain           = domain;
  _Mask             = NULL;
  _NumberOfFrames   = _InputSource.size();
  _NumberOfChannels = _InputTarget.size();
  _NumberOfVoxels   = domain._x * domain._y * domain._z;

  // Check inputs
  if (_Domain._x < 1 || _Domain._y < 1 || _Domain._z < 1) {
    cerr << "irtkImageSimilarityMetric::Initialize: Metric evaluation domain is empty" << endl;
    exit(1);
  }
  if (_Domain._t > 1) {
    cerr << "irtkImageSimilarityMetric::Initialize: Metric evaluation domain cannot have four dimensions" << endl;
    exit(1);
  }
  if (_NumberOfChannels == 0) {
    cerr << "irtkImageSimilarityMetric::Initialize: Missing target image" << endl;
    exit(1);
  }
  if (_NumberOfFrames == 0) {
    cerr << "irtkImageSimilarityMetric::Initialize: Missing source image" << endl;
    exit(1);
  }
  for (_CurrentFrame = 0; _CurrentFrame < _NumberOfFrames; ++_CurrentFrame) {
    if (_InputSource[_CurrentFrame].size() != static_cast<size_t>(_NumberOfChannels)) {
      if (_CurrentFrame == 0) {
        cerr << "irtkImageSimilarityMetric::Initialize: Source and target images must have the same number of channels" << endl;
      } else {
        cerr << "irtkImageSimilarityMetric::Initialize: All source frames must have the same number of channels" << endl;
      }
      exit(1);
    }
  }

  // Note: Always convert input images to internally used Image type such that
  //       subclass implementations can rely on a specific scalar type.
  //       This also applies to untransformed target images. Only when the input
  //       image is not transformed, defined on the metric domain grid, and has
  //       the same scalar type as used internally, we can just use the input directly.
  const bool transform_target = (_TargetTransformation.HasFixedTransformation() ||
                                 _TargetTransformation.HasTransformation()      ||
                                 _InputTarget   [0]._Image->GetImageAttributes() != _Domain);
  const bool transform_source = (_SourceTransformation[0].HasFixedTransformation() ||
                                 _SourceTransformation[0].HasTransformation()      ||
                                 _InputSource[0][0]._Image->GetImageAttributes() != _Domain);

  // Allocate memory for transformed target images
  _Target.resize(_NumberOfChannels);
  for (_CurrentChannel = 0; _CurrentChannel < _NumberOfChannels; ++_CurrentChannel) {
    if (!transform_target && InputTarget().GetScalarType() == voxel_info<Scalar>::type()) {
      TargetChannel()._Image      = dynamic_cast<Image *>(&InputTarget());
      TargetChannel()._ImageOwner = false;
    } else {
      TargetChannel()._Image      = new Image(_Domain);
      TargetChannel()._ImageOwner = true;
    }
  }

  // Allocate memory for transformed source images
  _Source.resize(_NumberOfFrames);
  for (_CurrentFrame = 0; _CurrentFrame < _NumberOfFrames; ++_CurrentFrame) {
    SourceFrame().resize(_NumberOfChannels);
    for (_CurrentChannel = 0; _CurrentChannel < _NumberOfChannels; ++_CurrentChannel) {
      if (!transform_source && InputSource().GetScalarType() == voxel_info<Scalar>::type()) {
        SourceChannel()._Image      = dynamic_cast<Image *>(&InputSource());
        SourceChannel()._ImageOwner = false;
      } else {
        SourceChannel()._Image      = new Image(_Domain);
        SourceChannel()._ImageOwner = true;
      }
    }
  }

  // Set background padding value
  _SourceTransformation[0].SetPaddingValue(numeric_limits<double>::quiet_NaN());
  _TargetTransformation   .SetPaddingValue(numeric_limits<double>::quiet_NaN());

  // Pre-compute world coordinates
  if (_InputImage2World) {
    if (_Domain.EqualInSpace(_InputImage2World->GetImageAttributes())) {
      cerr << "irtkImageSimilarityMetric::Initialize: World coordinates map has mismatching spatial attributes" << endl;
      exit(1);
    }
    _Image2World = const_cast<irtkWorldCoordsImage *>(_InputImage2World);
  } else {
    if (_CacheWorldCoordinates) {
      if (!_Image2World) _Image2World = new irtkWorldCoordsImage();
      _Target[0]._Image->ImageToWorld(*_Image2World);
    } else {
      delete _Image2World;
      _Image2World = NULL;
    }
  }
  _SourceTransformation[0].UseWorldCoordinates(_Image2World);
  _TargetTransformation   .UseWorldCoordinates(_Image2World);

  // Set input/output of target transformation helper
  _TargetTransformation.ClearInput();
  _TargetTransformation.ClearOutput();
  for (_CurrentChannel = 0; _CurrentChannel < _NumberOfChannels; ++_CurrentChannel) {
    if (&TargetImage() != &InputTarget()) {
      _TargetTransformation.AddInput (&InputTarget());
      _TargetTransformation.AddOutput(&TargetImage());
    }
  }
  _NumberOfTargetImagesToUpdate = _TargetTransformation.NumberOfImages();

  // Duplicate first source transformation helper with previously set common settings,
  // one copy for each additional input source frame
  _SourceTransformation.resize(1, _SourceTransformation[0]); // Discard previous copies
  _SourceTransformation.resize(_NumberOfFrames, _SourceTransformation[0]);

  // Set input/output of source transformation helpers
  _NumberOfSourceImagesToUpdate.resize(_NumberOfFrames, 0);
  for (_CurrentFrame = 0; _CurrentFrame < _NumberOfFrames; ++_CurrentFrame) {
    _SourceTransformation[_CurrentFrame].ClearInput();
    _SourceTransformation[_CurrentFrame].ClearOutput();
    for (_CurrentChannel = 0; _CurrentChannel < _NumberOfChannels; ++_CurrentChannel) {
      if (&SourceImage() != &InputSource()) {
        _SourceTransformation[_CurrentFrame].AddInput (&InputSource());
        _SourceTransformation[_CurrentFrame].AddOutput(&SourceImage());
      }
    }
    _NumberOfSourceImagesToUpdate[_CurrentFrame] = _SourceTransformation[_CurrentFrame].NumberOfImages();
  }

  // Append target gradient images to inputs/outputs of target transformation helper
  if (_TargetTransformation.HasTransformation()) {
     for (_CurrentChannel = 0; _CurrentChannel < _NumberOfChannels; ++_CurrentChannel) {
       InputTargetChannel().ComputeGradient();
            TargetChannel()._Gradient = new Gradient(_Domain, 3);
       _TargetTransformation.AddInput (InputTargetChannel()._Gradient);
       _TargetTransformation.AddOutput(     TargetChannel()._Gradient);
     }
  }

  // Append source gradient images to inputs/outputs of source transformation helper
  if (_SourceTransformation[0].HasTransformation()) {
    for (_CurrentFrame = 0; _CurrentFrame < _NumberOfFrames; ++_CurrentFrame) {
      for (_CurrentChannel = 0; _CurrentChannel < _NumberOfChannels; ++_CurrentChannel) {
        InputSourceChannel().ComputeGradient();
             SourceChannel()._Gradient = new Gradient(_Domain, 3);
        _SourceTransformation[_CurrentFrame].AddInput (InputSourceChannel()._Gradient);
        _SourceTransformation[_CurrentFrame].AddOutput(     SourceChannel()._Gradient);
      }
    }
  }

  // Initialize transformation helpers
  _TargetTransformation.Initialize();
  for (_CurrentFrame = 0; _CurrentFrame < _NumberOfFrames; ++_CurrentFrame) {
    _SourceTransformation[_CurrentFrame].Initialize();
  }
}

// -----------------------------------------------------------------------------
void irtkImageSimilarityMetric::Initialize(const irtkBinaryImage *mask, const irtkWorldCoordsImage *coords)
{
  if (_Image2World != _InputImage2World) delete _Image2World;
  _InputImage2World = coords;
  if (mask) {
    Initialize(mask->GetImageAttributes());
    _Mask = mask;
  } else if (coords) {
    Initialize(coords->GetImageAttributes());
  } else {
    if (_InputTarget.empty()) {
      cerr << "irtkImageSimilarityMetric::Initialize: Missing target image" << endl;
      exit(1);
    }
    Initialize(_InputTarget[0]._Image->GetImageAttributes());
  }
}

// -----------------------------------------------------------------------------
void irtkImageSimilarityMetric::Initialize(const irtkWorldCoordsImage *coords)
{
  this->Initialize(NULL, coords);
}

// -----------------------------------------------------------------------------
void irtkImageSimilarityMetric::Update(bool image, bool sgradient, bool tgradient)
{
  // Either update only the first n input images excluding the computed gradient
  // images or all images in particular including the transformed image gradients
  sgradient = sgradient && _SourceTransformation[0].HasTransformation();
  tgradient = tgradient && _TargetTransformation   .HasTransformation();

  if (image || sgradient) {
    // Attention: Do not use _CurrentFrame here as loop variable because it may
    //            be used by Update implementation of subclass!
    for (int i = 0; i < _NumberOfFrames; ++i) {
      _SourceTransformation[i].Update(0, (sgradient ? -1 : _NumberOfSourceImagesToUpdate[i]), _Mask);
    }
  }
  if (image || tgradient) {
    _TargetTransformation.Update(0, (tgradient ? -1 : _NumberOfTargetImagesToUpdate), _Mask);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Base class for probabilistc image similarity measures
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Construction/Destructions
// =============================================================================

// -----------------------------------------------------------------------------
irtkProbabilisticSimilarityMetric::irtkProbabilisticSimilarityMetric()
:
  irtkImageSimilarityMetric(),
  _Histogram(NULL)
{
}

// -----------------------------------------------------------------------------
irtkProbabilisticSimilarityMetric::irtkProbabilisticSimilarityMetric(const irtkProbabilisticSimilarityMetric &other)
:
  irtkImageSimilarityMetric(other),
  _Histogram(NULL)
{
}

// -----------------------------------------------------------------------------
irtkProbabilisticSimilarityMetric::~irtkProbabilisticSimilarityMetric()
{
  delete _Histogram;
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
void irtkProbabilisticSimilarityMetric::InitializeBins(int tbins, int sbins)
{
  // Use same number of bins for both images if only one argument provided
  if (sbins < 1) sbins = tbins;

  // Initialize joint histogram
  if (_Histogram) _Histogram->PutNumberOfBins    (tbins, sbins);
  else            _Histogram = new JointHistogram(tbins, sbins);

  // Rescale images for faster filling of joint histogram
  //
  // On either side of the range, 2 extra bins are added to account for
  // the width of the Parzen Window kernel used for histogram smoothing.
  for (_CurrentFrame = 0; _CurrentFrame < _NumberOfFrames; ++_CurrentFrame) {
    for (_CurrentChannel = 0; _CurrentChannel < _NumberOfChannels; ++_CurrentChannel) {
      InputTarget().PutMinMaxAsDouble(2, tbins - 3);
      InputSource().PutMinMaxAsDouble(2, sbins - 3);
    }
  }
}

// -----------------------------------------------------------------------------
void irtkProbabilisticSimilarityMetric::Initialize()
{
  // Initialize base class
  irtkImageSimilarityMetric::Initialize();
  // Initialize this class
  this->InitializeBins(64, 64);
}

// -----------------------------------------------------------------------------
void irtkProbabilisticSimilarityMetric::Initialize(int tbins, int sbins)
{
  // Initialize base class
  irtkImageSimilarityMetric::Initialize();
  // Initialize this class
  this->InitializeBins(tbins, sbins);
}

// -----------------------------------------------------------------------------
void irtkProbabilisticSimilarityMetric::Initialize(const irtkImageAttributes &attr, int tbins, int sbins)
{
  // Initialize base class
  irtkImageSimilarityMetric::Initialize(attr);
  // Initialize this class
  this->InitializeBins(tbins, sbins);
}

// -----------------------------------------------------------------------------
void irtkProbabilisticSimilarityMetric::Initialize(const irtkImageAttributes &attr)
{
  // Initialize base class
  irtkImageSimilarityMetric::Initialize(attr);
  // Initialize this class
  this->InitializeBins(64, 64);
}

// -----------------------------------------------------------------------------
void irtkProbabilisticSimilarityMetric::UpdateHistogram()
{
  // Reset histogram
  _Histogram->Reset();

  // Add histogram samples
  for (_CurrentFrame = 0; _CurrentFrame < _NumberOfFrames; ++_CurrentFrame) {
    for (_CurrentChannel = 0; _CurrentChannel < _NumberOfChannels; ++_CurrentChannel) {
      Scalar *tgt = TargetImage().GetPointerToVoxels();
      Scalar *src = SourceImage().GetPointerToVoxels();
      for (int idx = 0; idx < _NumberOfVoxels; ++idx, ++tgt, ++src) {
        if (EvaluateAt(idx)) _Histogram->Add(round(*tgt), round(*src));
      }
    }
  }

  // Smooth histogram
  _Histogram->Smooth();
}

// -----------------------------------------------------------------------------
void irtkProbabilisticSimilarityMetric::Update(bool image, bool sgradient, bool tgradient)
{
  // Update transformed images
  irtkImageSimilarityMetric::Update(image, sgradient, tgradient);

  // Update joint histogram
  UpdateHistogram();
}
