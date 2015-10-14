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

#ifndef _IRTKFILEGIPLTOIMAGE_H

#define _IRTKFILEGIPLTOIMAGE_H

/**
 * Class for reading images in GIPL file format.
 *
 * This is a class which reads images in GIPL file format and converts them
 * into images. The GIPL file format is a file format for 3D images. It has
 * been used in Guy's Hospital as universal file format for medical images.
 * Supported voxel types are char, unsigned char, short and unsigned short,
 * int, unsigned int and float.
 */

class irtkFileGIPLToImage : public irtkFileToImage
{
  irtkObjectMacro(irtkFileGIPLToImage);

protected:

  /// Read header of GIPL file
  virtual void ReadHeader();

public:

  /// Returns whether file has correct header
  static int CheckHeader(const char *);
};

#endif
