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

#ifndef _IRTKIMAGETOFILENIFTI_H

#define _IRTKIMAGETOFILENIFTI_H

#ifdef HAS_NIFTI

#include <irtkNIFTI.h>

/**
 * Class for image to NIFTI file writer.
 *
 * This is a class which takes an image as input and produces an image file
 * in NIFTI file format.
 */

class irtkImageToFileNIFTI : public irtkImageToFile
{
  irtkObjectMacro(irtkImageToFileNIFTI);

  /// Filename of header
  char *_headername;

  /// The nifti header
  irtkNIFTIHeader _hdr;

protected:

  /// Initialize filter
  virtual void Initialize();

  /// Finalize filter (empty, overrides parent method).
  virtual void Finalize();

public:

  /// Constructor
  irtkImageToFileNIFTI();

  /// Destructor
  virtual ~irtkImageToFileNIFTI();

  /// Set input
  virtual void SetOutput(const char *);

  /// Write entire image
  virtual void Run();
};

#endif

#endif
