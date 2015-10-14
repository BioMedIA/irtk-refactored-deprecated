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

#ifndef _IRTKFILENIFTITOIMAGE_H

#define _IRTKFILENIFTITOIMAGE_H

#ifdef HAS_NIFTI

/**
 * Class for reading images in NIFTI file format.
 *
 * This is a class which reads images in NIFTI file format and converts them
 * into images. The NIFTI file format is a file format for 3D and 4D images.
 * More information about the format can be found at http://nifti.nimh.nih.gov/nifti-1/
 */

class irtkFileNIFTIToImage : public irtkFileToImage
{
  irtkObjectMacro(irtkFileNIFTIToImage);

  /// Filename of header
  char *_headername;

protected:

  /// Read header of NIFTI file
  virtual void ReadHeader();

public:

  /// Contructor
  irtkFileNIFTIToImage();

  /// Destructor
  virtual ~irtkFileNIFTIToImage();

  /// Set input
  virtual void SetInput (const char *);

  /// Returns whether file has correct header
  static int CheckHeader(const char *);

  /// Print image file information
  virtual void Print();

};

#endif

#endif
