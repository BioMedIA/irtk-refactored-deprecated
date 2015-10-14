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

#ifndef _IRTKFILETOIMAGE_H

#define _IRTKFILETOIMAGE_H

/**
 * Abstract base class for any general file to image filter.
 *
 * This is the abstract base class which defines a common interface for all
 * filters which take a filename (referrencing an image file) as input and
 * produce an image as output. For each file format a derived class should
 * be created. Each derived class has to implement abstract member functions
 * ReadHeader() and CheckHeader().
 */

class irtkFileToImage : protected irtkCifstream
{
  irtkAbstractMacro(irtkFileToImage);

  /// File name of image file
  const char *_imagename;

protected:

  /// Image attributes
  irtkImageAttributes _attr;

  /// No. of bytes per voxel. Should be initialized by calling ReadHeader().
  int _bytes;

  /// Type of voxels. Should be initialized by calling ReadHeader().
  int _type;

  /// Intensity scaling parameter - slope (default: 1)
  float _slope;

  /// Intensity scaling parameter -  intercept (default: 0)
  float _intercept;

  /// Start of image data
  int _start;

  /// Flag whether to reflect X axis
  int _reflectX;

  /// Flag whether to reflect Y axis
  int _reflectY;

  /// Flag whether to reflect Z axis
  int _reflectZ;

  /// Debug flag
  int _debug;

  /** Read header. This is an abstract function. Each derived class has to
   *  implement this function in order to initialize image dimensions, voxel
   *  dimensions, voxel type and a lookup table which the address for each
   *  slice.
   */
  virtual void ReadHeader() = 0;

public:

  /// Contructor
  irtkFileToImage();

  /// Destructor
  virtual ~irtkFileToImage();

  /** Static constructor. This constructor allocates a derived class which
   *  is then  used to read the image file. This involves checking the file
   *  format of the file and dynamically allocating the corresonding derived
   *  class.
   */
  static irtkFileToImage *New(const char *);

  /// Set input
  virtual void SetInput (const char *);

  /// Get output
  virtual irtkImage *GetOutput();

  /// Get debug flag
  virtual int  GetDebugFlag();

  /// Get slope
  virtual double  GetSlope();

  /// Get intercept
  virtual double  GetIntercept();

  /// Put debug flag
  virtual void PutDebugFlag(int);

  /// Print image file information
  virtual void Print();

  /// Print  debugging info
  virtual void Debug(char *);

  // Return the type of the image's data
  virtual int GetDataType();
};

#include <irtkFilePGMToImage.h>
#ifdef HAS_VTK
#include <irtkFileVTKToImage.h>
#endif
#include <irtkFileGIPLToImage.h>
#ifdef HAS_NIFTI
#include <irtkFileNIFTIToImage.h>
#endif
#ifdef HAS_OPENCV
#include <irtkFileOpenCVToImage.h>
#endif
#include <irtkFileANALYZEToImage.h>

#endif
