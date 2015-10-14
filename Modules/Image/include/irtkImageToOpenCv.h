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

#ifndef _IRTKIMAGETOOPENCV_H

#define _IRTKIMAGETOOPENCV_H

#ifdef HAS_OPENCV
#include <cv.h>

/**
 * Abstract base class for any general image to image filter.
 *
 * This is the abstract base class which defines a common interface for all
 * filters which take an image as input and produce an image as output. Each
 * derived class has to implement all abstract member functions.
 */

template <class VoxelType> class irtkImageToOpenCv : public irtkObject
{
  irtkObjectMacro(irtkImageToOpenCv);

private:

  /// Buffer
  irtkGenericImage<VoxelType> *_tmp;

protected:


  /// Input image for filter
  irtkGenericImage<VoxelType> *_input;

  /// Output image for filter
  IplImage *_output;

  /// min and max value
  VoxelType min,max;

  /** Initialize the filter. This function must be called by any derived
   *  filter class to perform some initialize tasks. */
  virtual void Initialize();

  /** Finalize the filter. This function must be called by any derived
   *  filter class to perform some initialize tasks. */
  virtual void Finalize();

public:

  /// Constructor
  irtkImageToOpenCv();

  /// Deconstuctor
  virtual ~irtkImageToOpenCv();

  /// Set input image for filter
  virtual void SetInput (irtkGenericImage<VoxelType> *);

  /// Set output image for filter
  virtual void SetOutput(IplImage *);

   /// Set input image for filter
  virtual irtkGenericImage<VoxelType> * GetInput ();

  /// Set output image for filter
  virtual IplImage* GetOutput();

  /// Generate output based on input
  virtual void GenOutput ();

  /// Generate input based on output
  virtual void GenInput ();

  /// Run filter on entire image from irtk image to OpenCV image
  virtual void   Run(int = 0);

  /// Run filter on entire image from OpenCV image to irtk image
  virtual void   Invert(int = 0);
};


#endif

#endif
