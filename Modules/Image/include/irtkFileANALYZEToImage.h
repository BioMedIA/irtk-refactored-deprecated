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

#ifndef _IRTKFILEANALYZETOIMAGE_H

#define _IRTKFILEANALYZETOIMAGE_H

/**
 * Class for reading images in ANALYZE file format.
 *
 * This is a class which reads images in ANALYZE file format and converts them
 * into images. The ANALYZE file format is a file format for 3D images. It has
 * been developed by the Mayo Clinic, Rochester, MN. The header information
 * is stored in a file with extension ".hdr" while the image data is stored in
 * a file with with extension ".img".
 */

class irtkFileANALYZEToImage : public irtkFileToImage
{
  irtkObjectMacro(irtkFileANALYZEToImage);

  /// Filename of header
  char *_headername;

protected:

  /// Read header of ANALYZE file
  virtual void ReadHeader();

public:

  /// Contructor
  irtkFileANALYZEToImage();

  /// Destructor
  virtual ~irtkFileANALYZEToImage();

  /// Set input
  virtual void SetInput (const char *);

  /// Returns whether file has correct header
  static int CheckHeader(const char *);
};

#endif
