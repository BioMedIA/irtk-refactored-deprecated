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

#include <irtkHessianImageFilter.h>

template <class VoxelType> irtkHessianImageFilter<VoxelType>::irtkHessianImageFilter(int type)
{
  _type           = type;
  _UseVoxelSize   = true;
  _UseOrientation = false;
  _Padding        = MIN_GREY;
}

template <class VoxelType> void irtkHessianImageFilter<VoxelType>::Initialize()
{
  // Print debugging information
  this->Debug("irtkImageToImage::Initialize");

  // Check inputs and outputs
  if (this->_input == NULL) {
    cerr << this->NameOfClass() << "::Run: Filter has no input" << endl;
    exit(1);
  }

  if (this->_output == NULL) {
    cerr << this->NameOfClass() << "::Run: Filter has no output" << endl;
    exit(1);
  }

  if (this->_input->IsEmpty() == true) {
    cerr << this->NameOfClass() << "::Run: Input is empty" << endl;
    exit(1);
  }

  if (this->_input->GetT() > 1) {
    cerr << this->NameOfClass() << "::Run: Only implemented for images with t = 1" << endl;
    exit(1);
  }

  // Check whether filter requires buffering
  if (this->RequiresBuffering()) {
    this->Debug("irtkHessianImageFilter::Initialize: Filter requires buffering");

    // Check whether filter has external buffer
    if (this->_input == this->_output) {
      this->Debug("irtkHessianImageFilter::Initialize: Filter has internal buffer");
      this->_tmp    = this->_output;
      this->_output = new irtkGenericImage<VoxelType>;
    } else {
      this->Debug("irtkHessianImageFilter::Initialize: Filter has external buffer");
      this->_tmp    = NULL;
    }
  } else {
    this->Debug("irtkHessianImageFilter::Initialize: Filter requires no buffering");
  }

  // Make sure that output has the correct dimensions
  if (_type == HESSIAN_MATRIX) {
    irtkImageAttributes attr = this->_input->GetImageAttributes();
    attr._t = 9;
    this->_output->Initialize(attr);
  } else if (_type == HESSIAN_VECTOR) {
    irtkImageAttributes attr = this->_input->GetImageAttributes();
    attr._t = 6;
    this->_output->Initialize(attr);
  } else {
    if (this->_input != this->_output) this->_output->Initialize(this->_input->GetImageAttributes());
  }
}

template <class VoxelType> void irtkHessianImageFilter<VoxelType>::Run()
{
  double dxx, dxy, dxz, dyy, dyz, dzz, dii, dij, dik, djj, djk, dkk;
  int x, y, z, x1, y1, z1, x2, y2, z2;

  // Do the initial set up
  this->Initialize();

  const irtkImageAttributes &attr = this->_input->GetImageAttributes();
  irtkMatrix                 R    = attr.GetWorldToImageOrientation();

  for (z = 0; z < this->_input->GetZ(); ++z) {
    z1 = z - 1;
    if (z1 < 0) z1 = 0;
    z2 = z + 1;
    if (z2 > this->_input->GetZ()-1) z2 = this->_input->GetZ()-1;

    for (y = 0; y < this->_input->GetY(); ++y) {
      y1 = y - 1;
      if (y1 < 0) y1 = 0;
      y2 = y + 1;
      if (y2 > this->_input->GetY()-1) y2 = this->_input->GetY()-1;

      for (x = 0; x < this->_input->GetX(); ++x) {
        x1 = x - 1;
        if (x1 < 0) x1 = 0;
        x2 = x + 1;
        if (x2 > this->_input->GetX()-1) x2 = this->_input->GetX()-1;

        // Compute derivatives
        if (x1 != x2 &&
            this->_input->Get(x, y, z)  > _Padding &&
            this->_input->Get(x1, y, z) > _Padding &&
            this->_input->Get(x2, y, z) > _Padding) {
          dxx = (this->_input->Get(x2, y, z) - 2.0 * this->_input->Get(x, y, z) + this->_input->Get(x1, y, z));
        } else {
          dxx = .0;
        }
        if (x1 != x2 &&
            y1 != y2 &&
            this->_input->Get(x1, y1, z) > _Padding &&
            this->_input->Get(x1, y2, z) > _Padding &&
            this->_input->Get(x2, y1, z) > _Padding &&
            this->_input->Get(x2, y2, z) > _Padding) {
          dxy = (this->_input->Get(x2, y2, z) - this->_input->Get(x2, y1, z) - this->_input->Get(x1, y2, z) + this->_input->Get(x1, y1, z)) / ((x2 - x1) * (y2 - y1));
        } else {
          dxy = .0;
        }
        if (x1 != x2 &&
            z1 != z2 &&
            this->_input->Get(x1, y, z1) > _Padding &&
            this->_input->Get(x1, y, z2) > _Padding &&
            this->_input->Get(x2, y, z1) > _Padding &&
            this->_input->Get(x2, y, z2) > _Padding) {
          dxz = (this->_input->Get(x2, y, z2) - this->_input->Get(x2, y, z1) - this->_input->Get(x1, y, z2) + this->_input->Get(x1, y, z1)) / ((x2 - x1) * (z2 - z1));
        } else {
          dxz = .0;
        }

        if (y1 != y2 &&
            this->_input->Get(x, y, z)  > _Padding &&
            this->_input->Get(x, y1, z) > _Padding &&
            this->_input->Get(x, y2, z) > _Padding) {
          dyy = (this->_input->Get(x, y2, z) - 2.0 * this->_input->Get(x, y, z) + this->_input->Get(x, y1, z));
        } else {
          dyy = .0;
        }

        if (y1 != y2 &&
            z1 != z2 &&
            this->_input->Get(x, y1, z1) > _Padding &&
            this->_input->Get(x, y1, z2) > _Padding &&
            this->_input->Get(x, y2, z1) > _Padding &&
            this->_input->Get(x, y2, z2) > _Padding) {
          dyz = (this->_input->Get(x, y2, z2) - this->_input->Get(x, y2, z1) - this->_input->Get(x, y1, z2) + this->_input->Get(x, y1, z1)) / ((y2 - y1) * (z2 - z1));
        } else {
          dyz = .0;
        }

        if (z1 != z2 &&
            this->_input->Get(x, y, z)  > _Padding &&
            this->_input->Get(x, y, z1) > _Padding &&
            this->_input->Get(x, y, z2) > _Padding) {
          dzz = (this->_input->Get(x, y, z2) - 2.0 * this->_input->Get(x, y, z) + this->_input->Get(x, y, z1));
        } else {
          dzz = .0;
        }

        if (_UseVoxelSize) {
          dxx /= (this->_input->GetXSize() * this->_input->GetXSize());
          dxy /= (this->_input->GetXSize() * this->_input->GetYSize());
          dxz /= (this->_input->GetXSize() * this->_input->GetZSize());
          dyy /= (this->_input->GetYSize() * this->_input->GetYSize());
          dyz /= (this->_input->GetYSize() * this->_input->GetZSize());
          dzz /= (this->_input->GetZSize() * this->_input->GetZSize());
        }
        if (_UseOrientation) {
          // Using numerator-layout for matrix calculus.
          // http://en.wikipedia.org/wiki/Matrix_calculus#Numerator-layout_notation
          //
          // Expression computed here is transpose(R) * Hessian * R = transpose(Hessian * R) * R
          dii = dxx, dij = dxy, dik = dxz, djj = dyy, djk = dyz, dkk = dzz;
          dxx = R(0, 0) * (R(0, 0) * dii + R(1, 0) * dij + R(2, 0) * dik) + R(1, 0) * (R(0, 0) * dij + R(1, 0) * djj + R(2, 0) * djk) + R(2, 0) * (R(0, 0) * dik + R(1, 0) * djk + R(2, 0) * dkk);
          dxy = R(0, 1) * (R(0, 0) * dii + R(1, 0) * dij + R(2, 0) * dik) + R(1, 1) * (R(0, 0) * dij + R(1, 0) * djj + R(2, 0) * djk) + R(2, 1) * (R(0, 0) * dik + R(1, 0) * djk + R(2, 0) * dkk);
          dxz = R(0, 2) * (R(0, 0) * dii + R(1, 0) * dij + R(2, 0) * dik) + R(1, 2) * (R(0, 0) * dij + R(1, 0) * djj + R(2, 0) * djk) + R(2, 2) * (R(0, 0) * dik + R(1, 0) * djk + R(2, 0) * dkk);
          dyy = R(0, 1) * (R(0, 1) * dii + R(1, 1) * dij + R(2, 1) * dik) + R(1, 1) * (R(0, 1) * dij + R(1, 1) * djj + R(2, 1) * djk) + R(2, 1) * (R(0, 1) * dik + R(1, 1) * djk + R(2, 1) * dkk);
          dyz = R(0, 2) * (R(0, 1) * dii + R(1, 1) * dij + R(2, 1) * dik) + R(1, 2) * (R(0, 1) * dij + R(1, 1) * djj + R(2, 1) * djk) + R(2, 2) * (R(0, 1) * dik + R(1, 1) * djk + R(2, 1) * dkk);
          dzz = R(0, 2) * (R(0, 2) * dii + R(1, 2) * dij + R(2, 2) * dik) + R(1, 2) * (R(0, 2) * dij + R(1, 2) * djj + R(2, 2) * djk) + R(2, 2) * (R(0, 2) * dik + R(1, 2) * djk + R(2, 2) * dkk);
        }

        switch (_type) {
          case HESSIAN_XX:
            this->_output->PutAsDouble(x, y, z, 0, dxx);
            break;
          case HESSIAN_XY:
            this->_output->PutAsDouble(x, y, z, 0, dxy);
            break;
          case HESSIAN_XZ:
            this->_output->PutAsDouble(x, y, z, 0, dxz);
            break;
          case HESSIAN_YY:
            this->_output->PutAsDouble(x, y, z, 0, dyy);
            break;
          case HESSIAN_YZ:
            this->_output->PutAsDouble(x, y, z, 0, dyz);
            break;
          case HESSIAN_ZZ:
            this->_output->PutAsDouble(x, y, z, 0, dzz);
            break;
          case HESSIAN_VECTOR:
            this->_output->PutAsDouble(x, y, z, 0, dxx);
            this->_output->PutAsDouble(x, y, z, 1, dxy);
            this->_output->PutAsDouble(x, y, z, 2, dxz);
            this->_output->PutAsDouble(x, y, z, 3, dyy);
            this->_output->PutAsDouble(x, y, z, 4, dyz);
            this->_output->PutAsDouble(x, y, z, 5, dzz);
            break;
          case HESSIAN_MATRIX:
            this->_output->PutAsDouble(x, y, z, 0, dxx);
            this->_output->PutAsDouble(x, y, z, 1, dxy);
            this->_output->PutAsDouble(x, y, z, 2, dxz);
            this->_output->PutAsDouble(x, y, z, 3, dxy);
            this->_output->PutAsDouble(x, y, z, 4, dyy);
            this->_output->PutAsDouble(x, y, z, 5, dyz);
            this->_output->PutAsDouble(x, y, z, 6, dxz);
            this->_output->PutAsDouble(x, y, z, 7, dyz);
            this->_output->PutAsDouble(x, y, z, 8, dzz);
            break;
          default:
            cerr << this->NameOfClass() << "::Run: Unknown gradient computation" << endl;
            exit(1);
        }
      }
    }
  }

  // Do the final cleaning up
  this->Finalize();
}

template class irtkHessianImageFilter<unsigned char>;
template class irtkHessianImageFilter<short>;
template class irtkHessianImageFilter<float>;
template class irtkHessianImageFilter<double>;
