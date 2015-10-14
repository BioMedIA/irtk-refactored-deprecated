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

#ifndef _IRTKIMAGETOFILE_H

#define _IRTKIMAGETOFILE_H

#include <irtkFileToImage.h>

/**
 * Abstract base class for any general image to file filter.
 *
 * This is the abstract base class which defines a common interface for all
 * filters which take an image as input and produce an image file as output.
 * Each derived class has to implement all of the abstract member functions.
 */

class irtkImageToFile : protected irtkCofstream
{
  irtkObjectMacro(irtkImageToFile);

protected:

  /// Pointer to image for input
  irtkImage *_input;

  // Pointer to filename for output
  const char *_output;

  /** Start address of the data in the image file. Should be
   *  initialized by calling Initialize().
   */
  int _start;

  /// Flag whether to reflect X axis
  int _reflectX;

  /// Flag whether to reflect Y axis
  int _reflectY;

  /// Flag whether to reflect Z axis
  int _reflectZ;

  /// Initialize filter
  virtual void Initialize();

  /// Finalize filter
  virtual void Finalize();

public:

  /// Constructor
  irtkImageToFile();

  /// Destructor
  virtual ~irtkImageToFile();

  /** Static constructor. This constructor allocates a derived class which
   *  can be used to read the image file. This involves checking the file
   *  format of the file and dynamically allocating the corresonding derived
   *  class.
   */
  static irtkImageToFile *New(const char *);

  /// Set input
  virtual void SetInput (irtkImage *);

  /// Set output
  virtual void SetOutput(const char *);

  /// Write entire image
  virtual void Run();

};

#include <irtkImageToFilePGM.h>
#ifdef HAS_PNG
#include <irtkImageToFilePNG.h>
#endif
#ifdef HAS_VTK
#include <irtkImageToFileVTK.h>
#endif
#include <irtkImageToFileGIPL.h>
#include <irtkImageToFileANALYZE.h>
#ifdef HAS_NIFTI
#include <irtkImageToFileNIFTI.h>
#endif

#endif
