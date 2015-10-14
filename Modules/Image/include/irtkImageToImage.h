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

#ifndef _IRTKIMAGETOIMAGE_H

#define _IRTKIMAGETOIMAGE_H

#include <irtkImage.h>

#ifdef HAS_TBB

template <class VoxelType> class irtkMultiThreadedImageToImage;

#endif

/**
 * Abstract base class for any general image to image filter.
 *
 * This is the abstract base class which defines a common interface for all
 * filters which take an image as input and produce an image as output. Each
 * derived class has to implement all abstract member functions.
 */

template <class VoxelType>
class irtkImageToImage : public irtkObject
{
  irtkObjectMacro(irtkImageToImage);

#ifdef HAS_TBB

  friend class irtkMultiThreadedImageToImage<VoxelType>;

#endif

private:

  /// Debugging flag
  bool _DebugFlag;

protected:

  /// Buffer
  irtkGenericImage<VoxelType> *_tmp;

  /// Input image for filter
  irtkGenericImage<VoxelType> *_input;

  /// Output image for filter
  irtkGenericImage<VoxelType> *_output;

  /** Initialize the filter. This function must be called by any derived
   *  filter class to perform some initialize tasks. */
  void Initialize(bool);

  /** Initialize the filter. This function must be called by any derived
   *  filter class to perform some initialize tasks. */
  virtual void Initialize();

  /** Finalize the filter. This function must be called by any derived
   *  filter class to perform some initialize tasks. */
  virtual void Finalize();

public:

  /// Constructor
  irtkImageToImage();

  /// Destructor
  virtual ~irtkImageToImage();

  /// Set input image for filter
  virtual void SetInput (irtkGenericImage<VoxelType> *);

  /// Get input image of filter
  irtkGenericImage<VoxelType> *GetInput();

  /// Set output image for filter
  virtual void SetOutput(irtkGenericImage<VoxelType> *);

  /// Get output image of filter
  irtkGenericImage<VoxelType> *GetOutput();

  /// Run filter on entire image
  virtual void   Run();

  /// Run filter on single voxel
  virtual double Run(int, int, int, int = 0);

  /** Returns whether the filter requires buffering. Any derived class must
   *  implement this member function to indicate whether the filter should
   *  buffer the input in case that input and output are equal. For example,
   *  filters which only require the voxel value to calculate their output
   *  should return false, otherwise true.
   */
  virtual bool RequiresBuffering() const;

  /// Set debugging flag
  SetMacro(DebugFlag, bool);

  /// Get debugging flag
  GetMacro(DebugFlag, bool);

  /// Print debugging messages if debugging is enabled
  virtual void Debug(const char *);
};

///////////////////////////////////////////////////////////////////////////////
// Auxiliary macro for subclass implementation
///////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
#define irtkImageFilterMacro(clsname)                                          \
  irtkObjectMacro(clsname);                                                    \
  protected:                                                                   \
    /** Indicates that this filter cannot run in-place */                      \
    virtual bool RequiresBuffering() const { return true; }                    \
  private:                                                                     \
    static void _irtkImageFilterMacro_needs_trailing_semicolon()

// ----------------------------------------------------------------------------
#define irtkInPlaceImageFilterMacro(clsname)                                   \
  irtkObjectMacro(clsname);                                                    \
  protected:                                                                   \
    /** Indicates that this filter can run in-place */                         \
    virtual bool RequiresBuffering() const { return false; }                   \
  private:                                                                     \
    static void _irtkInPlaceImageFilterMacro_needs_trailing_semicolon()

///////////////////////////////////////////////////////////////////////////////
// Inline definitions
///////////////////////////////////////////////////////////////////////////////

// ---------------------------------------------------------------------------
template <class VoxelType>
bool irtkImageToImage<VoxelType>::RequiresBuffering() const
{
  return true;
}

// ---------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType> *irtkImageToImage<VoxelType>::GetInput()
{
  return _input;
}

// ---------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType> *irtkImageToImage<VoxelType>::GetOutput()
{
  return _output;
}

#endif
