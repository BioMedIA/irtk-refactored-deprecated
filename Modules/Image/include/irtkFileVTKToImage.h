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

#ifndef _IRTKFILEVTKTOIMAGE_H

#define _IRTKFILEVTKTOIMAGE_H

/**
 * Class for reading images in VTK file format.
 *
 * This is a class which reads images in VTK file format and converts them
 * into images. The VTK file format is a file format for 3D images. It is
 * used by the Visualization Toolkit to store structered points. Supported
 * voxel types are char, unsigned char, short and unsigned short.
 */

class irtkFileVTKToImage : public irtkFileToImage
{
  irtkObjectMacro(irtkFileVTKToImage);

protected:

  /// Read header of VTK file
  virtual void ReadHeader();

public:

  /// Returns whether file has correct header
  static int CheckHeader(const char *);
};

#endif
