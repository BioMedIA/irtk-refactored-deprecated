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

irtkImageToFile::irtkImageToFile()
{
  _input  = NULL;
  _output = NULL;
  _start   = 0;
  _reflectX = false;
  _reflectY = false;
  _reflectZ = false;
}

irtkImageToFile::~irtkImageToFile()
{
  _input  = NULL;
  if (_output != NULL)
    _output = NULL;
  _start = 0;
}

irtkImageToFile *irtkImageToFile::New(const char *imagename)
{
  irtkImageToFile *writer = NULL;

  // Check format for PNG
  if (strstr(imagename, ".png") != NULL) {
#ifdef HAS_PNG
    writer = new irtkImageToFilePNG;
    writer->SetOutput(imagename);
#else
    cerr << "irtkImageToFile::New: Cannot write PNG files. Make sure IRTK is built with PNG support." << endl;
    exit(1);
#endif
  }
  
  // Check format for GIPL
  if (strstr(imagename, ".gipl") != NULL) {
    writer = new irtkImageToFileGIPL;
    writer->SetOutput(imagename);
  }

  // Check format for VTK
  if (strstr(imagename, ".vtk") != NULL) {
#ifdef HAS_VTK
    writer = new irtkImageToFileVTK;
    writer->SetOutput(imagename);
#else
    cerr << "irtkImageToFile::New: Cannot write VTK files. Make sure IRTK is built with VTK support." << endl;
    exit(1);
#endif
  }

  // Check format for PGM
  if (strstr(imagename, ".pgm") != NULL) {
    writer = new irtkImageToFilePGM;
    writer->SetOutput(imagename);
  }

  // Check format for ANALYZE
  if (strstr(imagename, ".hdr") != NULL) {
    writer = new irtkImageToFileANALYZE;
    writer->SetOutput(imagename);
  }

  // Check format for NIFTI
  if (strstr(imagename, ".nii") != NULL) {
#ifdef HAS_NIFTI
    writer = new irtkImageToFileNIFTI;
    writer->SetOutput(imagename);
#else
    cerr << "irtkImageToFile::New: Cannot write NifTI files. Make sure IRTK is built with NifTI support." << endl;
    exit(1);
#endif
  }

  // Check for default file format
  if (writer == NULL) {
    writer = new irtkImageToFileGIPL;
    writer->SetOutput(imagename);
  }

  return writer;
}

void irtkImageToFile::SetInput(irtkImage *image)
{
  if (image != NULL) {
    _input = image;
  } else {
    cerr << this->NameOfClass() << "::SetInput: Output is not an image\n";
    exit(1);
  }
}

void irtkImageToFile::SetOutput(const char *name)
{
  if (name != NULL) {
    _output = name;
  } else {
    cerr << this->NameOfClass() << "::SetInput: Output is not a filename\n";
    exit(1);
  }
}

void irtkImageToFile::Initialize()
{
  // Check inputs and outputs
  if (_input == NULL) {
    cerr << this->NameOfClass() << "::Run: Filter has no input" << endl;
    exit(1);
  }
  if (_output == NULL) {
    cerr << this->NameOfClass() << "::Run: Filter has no output" << endl;
    exit(1);
  }

  // Open file for writing
  this->Open(_output);

  // Reflect if necessary
  if (_reflectX == true) _input->ReflectX();
  if (_reflectY == true) _input->ReflectY();
  if (_reflectZ == true) _input->ReflectZ();
}

void irtkImageToFile::Finalize()
{
  // Close file
  this->Close();

  // Reflect back if necessary
  if (_reflectX == true) _input->ReflectX();
  if (_reflectY == true) _input->ReflectY();
  if (_reflectZ == true) _input->ReflectZ();
}

void irtkImageToFile::Run()
{
  // Initialize filter
  this->Initialize();

  // Write data
  switch (this->_input->GetScalarType()) {
  case IRTK_VOXEL_CHAR: {
      this->WriteAsChar((char *)_input->GetScalarPointer(), _input->GetNumberOfVoxels(), _start);
      break;
    }
  case IRTK_VOXEL_UNSIGNED_CHAR: {
      this->WriteAsUChar((unsigned char *)_input->GetScalarPointer(), _input->GetNumberOfVoxels(), _start);
      break;
    }
  case IRTK_VOXEL_SHORT: {
      this->WriteAsShort((short *)_input->GetScalarPointer(), _input->GetNumberOfVoxels(), _start);
      break;
    }
  case IRTK_VOXEL_UNSIGNED_SHORT: {
      this->WriteAsUShort((unsigned short *)_input->GetScalarPointer(), _input->GetNumberOfVoxels(), _start);
      break;
    }
  case IRTK_VOXEL_FLOAT: {
      this->WriteAsFloat((float *)_input->GetScalarPointer(), _input->GetNumberOfVoxels(), _start);
      break;
    }
  case IRTK_VOXEL_DOUBLE: {
      this->WriteAsDouble((double *)_input->GetScalarPointer(), _input->GetNumberOfVoxels(), _start);
      break;
    }
  default:
  	cerr << "irtkImageToFile::Run(): Unknown voxel type" << endl;
  	exit(1);
  }
  
  // Finalize filter
  this->Finalize();
}
