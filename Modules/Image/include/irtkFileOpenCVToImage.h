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

#ifndef _IRTKFILECVTOIMAGE

#define _IRTKFILECVTOIMAGE

#ifdef HAS_OPENCV
#include <cv.h>
#include <highgui.h>

/**
 * Class for reading images in PGM file format.
 *
 * This is a class which reads images in PGM file format and converts them
 * into images. The PGM (portable graymap) file format is a file format for
 * 2D images and is defined in pgm(1). At the moment only images in PGM raw
 * file format are supported.
 */

class irtkFileOpenCVToImage : public irtkFileToImage
{
  irtkObjectMacro(irtkFileOpenCVToImage);

protected:

  IplImage* _pimage;

  virtual void ReadHeader();

public:

  /// Return whether file has correct header
  static int CheckHeader(const char *);

  /// Get output
  virtual irtkImage *GetOutput();

  /// Set input
  virtual void SetInput (const char *);
};

#endif

#endif
