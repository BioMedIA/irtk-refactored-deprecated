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

#ifndef _IRTKIMAGETOIMAGE2_H
#define _IRTKIMAGETOIMAGE2_H

#include <irtkImageToImage.h> // irtkImageFilterMacro, irtkInPlaceImageFilterMacro


/**
 * Abstract base class for any general image to image filter.
 *
 * This is the abstract base class which defines a common interface for all
 * filters which take an image as input and produce an image as output. Each
 * derived class has to implement all abstract member functions.
 */

class irtkImageToImage2 : public irtkObject
{
  irtkObjectMacro(irtkImageToImage2);

private:

  /// Buffer
  irtkBaseImage *_tmp;

  /// Debugging flag
  bool _DebugFlag;

protected:

  /// Input image for filter
  irtkBaseImage *_input;

  /// Output image for filter
  irtkBaseImage *_output;

  /** Initialize the filter. This function must be called by any derived
   *  filter class to perform some initialize tasks. */
  virtual void Initialize();

  /** Finalize the filter. This function must be called by any derived
   *  filter class to perform some initialize tasks. */
  virtual void Finalize();

public:

  /// Constructor
  irtkImageToImage2();

  /// Deconstuctor
  virtual ~irtkImageToImage2();

  /** Returns whether the filter requires buffering. Any derived class must
   *  implement this member function to indicate whether the filter should
   *  buffer the input in case that input and output are equal. For example,
   *  filters which only require the voxel value to calculate their output
   *  should return false, otherwise true.
   */
  virtual bool RequiresBuffering() const { return true; }

  /// Set input image for filter
  virtual void SetInput (irtkBaseImage *);

  /// Set output image for filter
  virtual void SetOutput(irtkBaseImage *);

  /// Run filter on entire image
  virtual void   Run();

  /// Run filter on single voxel
  virtual double Run(int, int, int, int = 0);

  /// Set debugging flag
  SetMacro(DebugFlag, bool);

  /// Get debugging flag
  GetMacro(DebugFlag, bool);

  /// Print debugging messages if debugging is enabled
  virtual void Debug(const char *);
};


#endif
