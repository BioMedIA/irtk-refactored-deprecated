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

#include <irtkImageToFile.h>

#ifdef HAS_PNG

#include <png.h>


void irtkImageToFilePNG::Initialize()
{
  // Initialize base class
  this->irtkImageToFile::Initialize();

  if (this->_input->GetZ() != 3) {
    cerr << this->NameOfClass() << " supports only images with (z = 3, e.g. R-G-B)" << endl;
    exit(1);
  }
  
  if (dynamic_cast<irtkGenericImage<unsigned char> *>(_input) == NULL){
    cerr << this->NameOfClass() << " supports only images of voxel type unsigned char" << endl;
    exit(1);
  }
}

void irtkImageToFilePNG::Run()
{
  int x, y;
  irtkGenericImage<unsigned char> *image;
  
  // Initialize filter
  this->Initialize();

  png_structp png_ptr = png_create_write_struct
                        (PNG_LIBPNG_VER_STRING, (png_voidp)NULL, NULL, NULL);

  if (!png_ptr) {
    cerr << "irtkImageToFilePNG::Run: Unable to write PNG file!" << endl;
    exit(1);
  }

  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) {
    png_destroy_write_struct(&png_ptr,
                             (png_infopp)NULL);
    cerr << "irtkImageToFilePNG::Run: Unable to write PNG file!" << endl;;
    exit(1);
  }

  if (setjmp(png_jmpbuf(png_ptr))) {
    png_destroy_write_struct(&png_ptr, &info_ptr);
    exit(1);
  }

  // Open file
  FILE *fp = fopen(_output, "wb");

  // Initialize PNG I/O
  png_init_io(png_ptr, fp);

  // Initialize header
  png_set_IHDR(png_ptr, info_ptr, _input->GetX(), _input->GetY(),
               8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
               PNG_COMPRESSION_TYPE_DEFAULT,
               PNG_FILTER_TYPE_DEFAULT);

  // Write header
  png_write_info(png_ptr, info_ptr);

  // Copy image
  png_byte *data = new png_byte  [3*_input->GetX()*_input->GetY()];
  png_byte **ptr = new png_byte *[  _input->GetY()];
  png_byte *ptr2data = data;

  image = dynamic_cast<irtkGenericImage<unsigned char> *>(_input);
  for (y = 0; y < image->GetY(); y++) {
    for (x = 0; x < image->GetX(); x++) {
      *ptr2data = image->Get(x, y, 0, 0);
      ptr2data++;
      *ptr2data = image->Get(x, y, 1, 0);
      ptr2data++;
      *ptr2data = image->Get(x, y, 2, 0);
      ptr2data++;
    }
    ptr[_input->GetY() - y - 1] = &(data[_input->GetX()*y*3]);
  }
  png_write_image(png_ptr, ptr);
  png_write_end(png_ptr, info_ptr);

  // Delete pointers
  delete [] data;
  delete [] ptr;

  // Destroy PNG data
  png_destroy_write_struct(&png_ptr, &info_ptr);

  // Close file
  fclose(fp);

  // Finalize filter
  this->Finalize();
}

#endif
