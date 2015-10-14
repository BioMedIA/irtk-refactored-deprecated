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

#include <irtkScalarFunctionToImage.h>

template <class VoxelType>
irtkScalarFunctionToImage<VoxelType>::irtkScalarFunctionToImage(bool UseWorldCorrdinates)
{
  // Set in- and outputs
  _input  = NULL;
  _output = NULL;

  // Default parameters
  _DebugFlag = false;
  _UseWorldCoordinates = UseWorldCorrdinates;
}

template <class VoxelType>
irtkScalarFunctionToImage<VoxelType>::~irtkScalarFunctionToImage()
{
  // Set in- and outputs
  _input  = NULL;
  _output = NULL;
}

template <class VoxelType>
void irtkScalarFunctionToImage<VoxelType>::SetInput(irtkScalarFunction *irtkScalarfunction)
{
  if (irtkScalarfunction != NULL) {
    _input = irtkScalarfunction;
  } else {
    cerr << "irtkScalarFunctionToImage::SetInput: Input is not a scalar function\n";
    exit(1);
  }
}

template <class VoxelType>
void irtkScalarFunctionToImage<VoxelType>::SetOutput(irtkGenericImage<VoxelType> *image)
{
  if (image != NULL) {
    _output = image;
  } else {
    cerr << "irtkScalarFunctionToImage::SetOutput: Output is not an image\n";
    exit(1);
  }
}

template <class VoxelType>
void irtkScalarFunctionToImage<VoxelType>::Debug(char *message)
{
  if (_DebugFlag == true)
    cout << message << endl;
}

template <class VoxelType>
void irtkScalarFunctionToImage<VoxelType>::Run()
{
  int i, j, k, l;
  double x, y, z;

  if (_UseWorldCoordinates == true) {
    // Calculate scalar function using world coordinates
    for (l = 0; l < _output->GetT(); l++) {
      for (k = 0; k < _output->GetZ(); k++) {
        for (j = 0; j < _output->GetY(); j++) {
          for (i = 0; i < _output->GetX(); i++) {
            x = i;
            y = j;
            z = k;
            _output->ImageToWorld(x, y, z);
            _output->PutAsDouble(i, j, k, l, _input->Evaluate(x, y, z));
            if (fabs(static_cast<double>(_output->Get(i, j, k, l)))
                < FLT_MIN) {
              _output->Put(i, j, k, l, 0);
            }
          }
        }
      }
    }
  } else {
    // Calculate scalar function using image coordinates
    for (l = 0; l < _output->GetT(); l++) {
      for (k = 0; k < _output->GetZ(); k++) {
        for (j = 0; j < _output->GetY(); j++) {
          for (i = 0; i < _output->GetX(); i++) {
            _output->PutAsDouble(i, j, k, l, _input->Evaluate(i, j, k));
            if (fabs(static_cast<double>(_output->Get(i, j, k, l)))
                < FLT_MIN) {
              _output->Put(i, j, k, l, 0);
            }
          }
        }
      }
    }
  }
}

template class irtkScalarFunctionToImage<irtkBytePixel>;
template class irtkScalarFunctionToImage<irtkGreyPixel>;
template class irtkScalarFunctionToImage<irtkRealPixel>;
