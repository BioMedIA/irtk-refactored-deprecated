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

int irtkFilePGMToImage::CheckHeader(const char *filename)
{
  char buffer[255];

  // Create file stream
  ifstream from;

  // Open new file for reading
  from.open(filename);

  // Read header, skip comments
  do {
    from.get(buffer, 255);
    from.seekg(1, ios::cur);
  } while (buffer[0] == '#');

  // Close file
  from.close();

  // Check header
  if (strcmp(buffer, PGM_MAGIC) != 0) {
    return false;
  } else {
    return true;
  }
}

void irtkFilePGMToImage::ReadHeader()
{
  char buffer[255];

  // Read header, skip comments
  do {
    this->ReadAsString(buffer, 255);
  } while (buffer[0] == '#');

  // Check header
  if (strcmp(buffer, PGM_MAGIC) != 0) {
    cerr << this->NameOfClass() << "::Read_Header: Can't read magic number: "
         << buffer << endl;
    exit(1);
  }

  // Read voxel dimensions, skip comments
  do {
    this->ReadAsString(buffer, 255);
  } while (buffer[0] == '#');

  // Parse voxel dimensions
  sscanf(buffer, "%d %d", &this->_attr._x, &this->_attr._y);

  // Ignore maximum greyvalue, skip comments
  do {
    this->ReadAsString(buffer, 255);
  } while (buffer[0] == '#');

  // PGM files support only 2D images, so set z and t to 1
  this->_attr._z = 1;
  this->_attr._t = 1;

  // PGM files do not have voxel dimensions, so set them to default values
  this->_attr._dx = 1;
  this->_attr._dy = 1;
  this->_attr._dz = 1;
  this->_attr._dt = 1;

  // PGM files have voxels which are unsigned char
  this->_type  = IRTK_VOXEL_UNSIGNED_CHAR;
  this->_bytes = 1;

  // Data starts here
  this->_start = this->Tell();
}

