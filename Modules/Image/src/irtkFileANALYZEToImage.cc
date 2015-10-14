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

#include <irtkImage.h>

#include <irtkFileToImage.h>

#include <irtkANALYZE.h>

#include <sys/types.h>
#include <sys/stat.h>

irtkFileANALYZEToImage::irtkFileANALYZEToImage()
{
  _headername = NULL;

  // Analyze-specific
  this->_reflectY = true;
}

irtkFileANALYZEToImage::~irtkFileANALYZEToImage()
{
  if (this->_headername != NULL) free(this->_headername);
}

int irtkFileANALYZEToImage::CheckHeader(const char *filename)
{
  return Extension(filename, EXT_LastWithoutGz) == ".hdr";
}

void irtkFileANALYZEToImage::SetInput(const char *headername)
{
  int length;
  char imagename[255];
  struct stat buf;

  // Delete old file name
  if (this->_headername != NULL) free(this->_headername);

  // Copy new file name
  this->_headername = strdup(headername);

  // Parse img file name
  string ext = Extension(headername, EXT_LastWithGz);
  if (ext == ".hdr.gz") {
    length = strlen(headername);
    sprintf(imagename, "%s", headername);
    imagename[length-1] = 'z';
    imagename[length-2] = 'g';
    imagename[length-3] = '.';
    imagename[length-4] = 'g';
    imagename[length-5] = 'm';
    imagename[length-6] = 'i';

    // Check if compressed file exists
    if (stat(imagename, &buf) != 0) {
      cerr << this->NameOfClass() << ": Can't open file " << imagename << endl;
      exit(1);
    }
  } else {
    length = strlen(headername);
    sprintf(imagename, "%s", headername);
    imagename[length-1] = 'g';
    imagename[length-2] = 'm';
    imagename[length-3] = 'i';

    // Check if uncompressed file exists
    if (stat(imagename, &buf) != 0) {
      sprintf(imagename, "%s.gz", headername);
      imagename[length-1] = 'g';
      imagename[length-2] = 'm';
      imagename[length-3] = 'i';

      // Check if gzip compressed file exists
      if (stat(imagename, &buf) != 0) {
        sprintf(imagename, "%s.Z", headername);
        imagename[length-1] = 'g';
        imagename[length-2] = 'm';
        imagename[length-3] = 'i';
        if (stat(imagename, &buf) != 0) {
          cerr << this->NameOfClass() << ": Can't open file " << imagename << endl;
          exit(1);
        }
      }
    }
  }

  this->irtkFileToImage::SetInput(imagename);
}

void irtkFileANALYZEToImage::ReadHeader()
{
  irtkANALYZEHeader hdr;

  // Read header
  hdr.Read(this->_headername);

  // Copy header information
  if (hdr.dims[0] == 3) {
    this->_attr._x  = hdr.dims[1];
    this->_attr._y  = hdr.dims[2];
    this->_attr._z  = hdr.dims[3];
    this->_attr._t  = 1;
    this->_attr._dx = hdr.pixdims[1];
    this->_attr._dy = hdr.pixdims[2];
    this->_attr._dy = hdr.pixdims[3];
    this->_attr._dt = 1;
  } else {
    if (hdr.dims[0] == 4) {
      this->_attr._x  = hdr.dims[1];
      this->_attr._y  = hdr.dims[2];
      this->_attr._z  = hdr.dims[3];
      this->_attr._dx = hdr.pixdims[1];
      this->_attr._dy = hdr.pixdims[2];
      this->_attr._dz = hdr.pixdims[3];
      if (hdr.dims[4] <= 0) {
        this->_attr._t  = 1;
        this->_attr._dt = 1;
      } else {
        this->_attr._t  = hdr.dims[4];
        if (hdr.pixdims[4] < 0.0001) hdr.pixdims[4] = 1;
        this->_attr._dt = hdr.pixdims[4];
      }
    } else {
      cerr << "irtkFileANALYZEToImage::ReadHeader: Invalid no. of dimensions in image" << endl;
      exit(1);
    }
  }

  // If the size of header is not 348 we need to invert swapping
  if (hdr.sizeof_hdr != 348) {
    this->_Swapped = !this->_Swapped;
  }

  switch (hdr.data_type) {
  case ANALYZE_UNSIGNED_CHAR:
    this->_type  = IRTK_VOXEL_UNSIGNED_CHAR;
    this->_bytes = 1;
    break;
  case ANALYZE_SIGNED_SHORT:
    this->_type  = IRTK_VOXEL_SHORT;
    this->_bytes = 2;
    break;
  case ANALYZE_FLOAT:
    this->_type  = IRTK_VOXEL_FLOAT;
    this->_bytes = 4;
    break;
  case ANALYZE_DOUBLE:
    this->_type  = IRTK_VOXEL_DOUBLE;
    this->_bytes = 8;
    break;
  default:
    cerr << this->NameOfClass() << ": Data type " << hdr.data_type << " not supported, trying signed short data type" << endl;
    this->_type  = IRTK_VOXEL_SHORT;
    this->_bytes = 2;
  }

  // Data starts at 0
  this->_start = 0;
}
