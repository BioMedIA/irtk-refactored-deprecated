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

#ifndef _IRTKIMAGESEQUENCE_HH
#define _IRTKIMAGESEQUENCE_HH

#include <irtkImageSequence.h>


////////////////////////////////////////////////////////////////////////////////
// Inline definitions of irtkImageChannel
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
template <>
irtkImageChannel<irtkBaseImage>::irtkImageChannel(ImageType *image, bool manage, bool copy)
:
  irtkObject(),
  _Image (image),
  _Manage(manage)
{
  if (image) {
    if (image->GetT() > 0) {
      cerr << "irtkImageChannel::irtkImageChannel: Channel image cannot have fourth dimension" << endl;
      exit(1);
    }
    if (copy && !manage) {
      _Image  = irtkBaseImage::New(image);
      _Manage = true;
    }
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
irtkImageChannel<TImage>::irtkImageChannel(ImageType *image, bool manage, bool copy)
:
  irtkObject(),
  _Image (image),
  _Manage(manage)
{
  if (image) {
    if (image->GetT() > 0) {
      cerr << "irtkImageChannel::irtkImageChannel: Channel image cannot have fourth dimension" << endl;
      exit(1);
    }
    if (copy && !manage) {
      _Image  = new ImageType(image);
      _Manage = true;
    }
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
irtkImageChannel<TImage> &irtkImageChannel<TImage>::operator =(const irtkImageChannel<TImage> &other)
{
  if (_Manage) delete _Image;
  _Image  = other._Image;
  _Manage = other._Manage;
  other._Manage = false; // take over ownership
  return *this;
}

// -----------------------------------------------------------------------------
template <class TImage>
irtkImageChannel<TImage>::irtkImageChannel(const irtkImageChannel<TImage> &other)
:
  irtkObject(other),
  _Image (NULL),
  _Manage(false)
{
  *this = other;
}

// -----------------------------------------------------------------------------
template <> template <class TOtherImage>
irtkImageChannel<irtkBaseImage> &irtkImageChannel<irtkBaseImage>::operator =(const irtkImageChannel<TOtherImage> &other)
{
  if (_Manage) delete _Image;
  _Image  = other._Image;
  _Manage = other._Manage;
  other._Manage = false; // take over ownership
  return *this;
}

// -----------------------------------------------------------------------------
template <> template <class TOtherImage>
irtkImageChannel<irtkBaseImage>::irtkImageChannel(const irtkImageChannel<TOtherImage> &other)
:
  irtkObject(other),
  _Image (NULL),
  _Manage(false)
{
  *this = other;
}

// -----------------------------------------------------------------------------
template <class TImage>
irtkImageChannel<TImage>::~irtkImageChannel()
{
  if (_Manage) delete _Image;
}

// -----------------------------------------------------------------------------
template <>
void irtkImageChannel<irtkBaseImage>::Image(ImageType *image, bool manage, bool copy)
{
  if (_Manage) delete _Image;
  if (image && copy && !manage) {
    _Image  = irtkBaseImage::New(image);
    _Manage = true;
  } else {
    _Image  = image;
    _Manage = manage;
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
void irtkImageChannel<TImage>::Image(ImageType *image, bool manage, bool copy)
{
  if (_Manage) delete _Image;
  if (image && copy && !manage) {
    _Image  = new ImageType(*image);
    _Manage = true;
  } else {
    _Image  = image;
    _Manage = manage;
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
TImage *irtkImageChannel<TImage>::Image() const
{
  return _Image;
}

////////////////////////////////////////////////////////////////////////////////
// Inline definitions of irtkImageFrame
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class TChannel>
inline irtkImageFrame<TChannel>::irtkImageFrame()
:
  irtkObject()
{
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline irtkImageFrame<TChannel>::~irtkImageFrame()
{
}

// =============================================================================
// Frame
// =============================================================================

// -----------------------------------------------------------------------------
template <class TChannel>
inline void irtkImageFrame<TChannel>::Add(ImageType *image, bool manage, bool copy)
{
  if (_Channel.size() > 0) {
    irtkImageAttributes attr = image->GetImageAttributes();
    attr._t  = 1; // may differ
    attr._dt =  _Channel[0].Image().GetTSize();
    if (attr != _Channel[0].Image().GetImageAttributes()) {
      cerr << "irtkImageFrame::Add: Attributes of image do not match those of first channel" << endl;
      exit(1);
    }
  }
  const int T = image->GetT();
  if (T > 0) {
    for (int t = 0; t < T; ++t) {
      irtkBaseImage *channel;
      image.GetFrame(t, channel);
      _Channel.push_back(ChannelType(channel, true, false));
    }
  } else {
    _Channel.push_back(ChannelType(image, manage, copy));
  }
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline void irtkImageFrame<TChannel>::Clear()
{
  _Channel.clear();
}

// =============================================================================
// Channels
// =============================================================================

// -----------------------------------------------------------------------------
template <class TChannel>
inline void irtkImageFrame<TChannel>::NumberOfChannels(int n)
{
  _Channel.resize(n);
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline int irtkImageFrame<TChannel>::NumberOfChannels() const
{
  return static_cast<int>(_Channel.size());
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline typename irtkImageFrame<TChannel>::ChannelType &irtkImageFrame<TChannel>::Channel(int idx)
{
  return _Channel[idx];
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline const typename irtkImageFrame<TChannel>::ChannelType &irtkImageFrame<TChannel>::Channel(int idx) const
{
  return _Channel[idx];
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline void irtkImageFrame<TChannel>::Image(int idx, ImageType *image, bool manage, bool copy)
{
  _Channel[idx].Image(image, manage, copy);
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline typename irtkImageFrame<TChannel>::ImageType *irtkImageFrame<TChannel>::Image(int idx) const
{
  return _Channel[idx].Image();
}

// =============================================================================
// Attributes
// =============================================================================

// -----------------------------------------------------------------------------
template <class TChannel>
inline irtkImageAttributes irtkImageFrame<TChannel>::Attributes() const
{
  const ImageType *image = Image(0);
  return image ? image->GetImageAttributes() : irtkImageAttributes();
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline int irtkImageFrame<TChannel>::NumberOfVoxels() const
{
  const ImageType *image = Image(0);
  return image ? image->GetNumberOfVoxels() : 0;
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline int irtkImageFrame<TChannel>::X() const
{
  const ImageType *image = Image(0);
  return image ? image->GetX() : 0;
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline int irtkImageFrame<TChannel>::Y() const
{
  const ImageType *image = Image(0);
  return image ? image->GetY() : 0;
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline int irtkImageFrame<TChannel>::Z() const
{
  const ImageType *image = Image(0);
  return image ? image->GetZ() : 0;
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline double irtkImageFrame<TChannel>::XSize() const
{
  const ImageType *image = Image(0);
  return image ? image->GetXSize() : .0;
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline double irtkImageFrame<TChannel>::YSize() const
{
  const ImageType *image = Image(0);
  return image ? image->GetYSize() : .0;
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline double irtkImageFrame<TChannel>::ZSize() const
{
  const ImageType *image = Image(0);
  return image ? image->GetZSize() : .0;
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline void irtkImageFrame<TChannel>::GetPixelSize(double *dx, double *dy, double *dz) const
{
  const ImageType *image = Image(0);
  if (image) image->GetPixelSize(dx, dy, dz);
  else {
    if (dx) *dx = .0;
    if (dy) *dy = .0;
    if (dz) *dz = .0;
  }
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline void irtkImageFrame<TChannel>::ImageToWorld(double &x, double &y, double &z) const
{
  const ImageType *image = Image(0);
  if (image) image->ImageToWorld(x, y, z);
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline void irtkImageFrame<TChannel>::WorldToImage(double &x, double &y, double &z) const
{
  const ImageType *image = Image(0);
  if (image) image->WorldToImage(x, y, z);
}

// -----------------------------------------------------------------------------
template <class TChannel>
inline double irtkImageFrame<TChannel>::Time() const
{
  const ImageType *image = Image(0);
  if (image) image->ImageToTime(.0);
}

////////////////////////////////////////////////////////////////////////////////
// Inline definitions of irtkImageSequence
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class TFrame>
inline irtkImageSequence<TFrame>::irtkImageSequence()
:
  irtkObject()
{
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline irtkImageSequence<TFrame>::irtkImageSequence(const irtkImageSequence &other)
:
  irtkObject(other),
  _Frame(other._Frame)
{
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline irtkImageSequence<TFrame> &irtkImageSequence<TFrame>::operator =(const irtkImageSequence &other)
{
  _Frame = other._Frame;
  return *this;
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline irtkImageSequence<TFrame>::~irtkImageSequence()
{
}

// =============================================================================
// Sequence
// =============================================================================

// -----------------------------------------------------------------------------
template <class TFrame>
inline void irtkImageSequence<TFrame>::Add(ImageType *image, bool manage, bool copy)
{
  if (NumberOfFrames() > 0 && NumberOfChannels() > 0) {
    if (!image->HasSpatialAttributesOf(Image(0, 0))) {
      cerr << "irtkImageSequence::Add: Spatial attributes of image differ from those of first frame" << endl;
      exit(1);
    }
  }
  // Reserve enough entries in vector to ensure that no reallocation
  // takes place during the insert to keep tmp iterator below valid
  _Frame.reserve(_Frame.size() + image->GetT());
  // Insert each input frame (dt != 0) or channel (dt == 0)
  // Frames are sorted by increasing time and channels appended
  for (int l = 0; l < image->GetT(); ++l) {
    // Find corresponding frame or add new one if necessary
    const double                         time  = image->ImageToTime(l);
    typename vector<FrameType>::iterator frame = _Frame.begin();
    while (frame != _Frame.end()) {
      const double t = frame->Time();
      if      (t == time) break;
      else if (t  > time) {
        typename vector<FrameType>::iterator tmp = frame - 1;
        _Frame.insert(frame, FrameType());
        frame = ++tmp;
        break;
      }
      ++frame;
    }
    if (frame == _Frame.end()) {
      _Frame.push_back(FrameType());
      frame = _Frame.end() - 1;
    }
    // Add channel to frame
    if (image->GetT() > 0) {
      irtkBaseImage *channel;
      image->GetFrame(l, channel);
      frame->Add(channel, true, false);
    } else {
      frame->Add(image, manage, copy);
    }
  }
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline void irtkImageSequence<TFrame>::Clear()
{
  _Frame.clear();
}

// =============================================================================
// Frames
// =============================================================================

// -----------------------------------------------------------------------------
template <class TFrame>
inline void irtkImageSequence<TFrame>::NumberOfFrames(int n) const
{
  _Frame.resize(n);
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline int irtkImageSequence<TFrame>::NumberOfFrames() const
{
  return static_cast<int>(_Frame.size());
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline typename irtkImageSequence<TFrame>::FrameType &irtkImageSequence<TFrame>::Frame(int f)
{
  return _Frame[f];
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline const typename irtkImageSequence<TFrame>::FrameType &irtkImageSequence<TFrame>::Frame(int f) const
{
  return _Frame[f];
}

// =============================================================================
// Channels
// =============================================================================

// -----------------------------------------------------------------------------
template <class TFrame>
inline void irtkImageSequence<TFrame>::NumberOfChannels(int n) const
{
  for (int f = 0; f < NumberOfChannels(); ++f) Frame(f).NumberOfChannels(n);
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline int irtkImageSequence<TFrame>::NumberOfChannels() const
{
  return NumberOfFrames() > 0 ? Frame(0).NumberOfChannels() : 0;
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline typename irtkImageSequence<TFrame>::ChannelType &irtkImageSequence<TFrame>::Channel(int idx)
{
  const int num = NumberOfChannels();
  return _Frame[idx / num][idx % num];
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline const typename irtkImageSequence<TFrame>::ChannelType &irtkImageSequence<TFrame>::Channel(int idx) const
{
  const int num = NumberOfChannels();
  return _Frame[idx / num][idx % num];
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline typename irtkImageSequence<TFrame>::ChannelType &irtkImageSequence<TFrame>::Channel(int f, int c)
{
  return _Frame[f][c];
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline const typename irtkImageSequence<TFrame>::ChannelType &irtkImageSequence<TFrame>::Channel(int f, int c) const
{
  return _Frame[f][c];
}

// =============================================================================
// Images
// =============================================================================

// -----------------------------------------------------------------------------
template <class TFrame>
inline int irtkImageSequence<TFrame>::NumberOfImages() const
{
  return NumberOfFrames() * NumberOfChannels();
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline typename irtkImageSequence<TFrame>::ImageType *irtkImageSequence<TFrame>::Image(int idx) const
{
  return Channel(idx).Image();
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline typename irtkImageSequence<TFrame>::ImageType *irtkImageSequence<TFrame>::Image(int f, int c) const
{
  return Channel(f, c).Image();
}

// =============================================================================
// Attributes
// =============================================================================

// -----------------------------------------------------------------------------
template <class TFrame>
inline irtkImageAttributes irtkImageSequence<TFrame>::Attributes() const
{
  const ImageType *image = Image(0, 0);
  irtkImageAttributes attr;
  if (image) attr = image->GetImageAttributes();
  attr._t = NumberOfFrames();
  return attr;
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline int irtkImageSequence<TFrame>::NumberOfVoxels() const
{
  const ImageType *image = Image(0, 0);
  return image ? image->NumberOfVoxels() : 0;
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline int irtkImageSequence<TFrame>::X() const
{
  const ImageType *image = Image(0, 0);
  return image ? image->GetX() : 0;
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline int irtkImageSequence<TFrame>::Y() const
{
  const ImageType *image = Image(0, 0);
  return image ? image->GetY() : 0;
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline int irtkImageSequence<TFrame>::Z() const
{
  const ImageType *image = Image(0, 0);
  return image ? image->GetZ() : 0;
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline int irtkImageSequence<TFrame>::T() const
{
  return NumberOfFrames();
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline double irtkImageSequence<TFrame>::XSize() const
{
  const ImageType *image = Image(0, 0);
  return image ? image->GetXSize() : .0;
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline double irtkImageSequence<TFrame>::YSize() const
{
  const ImageType *image = Image(0, 0);
  return image ? image->GetYSize() : .0;
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline double irtkImageSequence<TFrame>::ZSize() const
{
  const ImageType *image = Image(0, 0);
  return image ? image->GetZSize() : .0;
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline void irtkImageSequence<TFrame>::GetPixelSize(double *dx, double *dy, double *dz) const
{
  const ImageType *image = Image(0, 0);
  if (image) image->GetPixelSize(dx, dy, dz);
  else {
    if (dx) *dx = .0;
    if (dy) *dy = .0;
    if (dz) *dz = .0;
  }
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline void irtkImageSequence<TFrame>::ImageToWorld(double &x, double &y, double &z) const
{
  const ImageType *image = Image(0, 0);
  if (image) image->ImageToWorld(x, y, z);
}

// -----------------------------------------------------------------------------
template <class TFrame>
inline void irtkImageSequence<TFrame>::WorldToImage(double &x, double &y, double &z) const
{
  const ImageType *image = Image(0, 0);
  if (image) image->WorldToImage(x, y, z);
}


#endif
