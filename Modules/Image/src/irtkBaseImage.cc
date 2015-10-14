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
#include <irtkNIFTI.h>

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkBaseImage::irtkBaseImage()
:
  _mask     (NULL),
  _maskOwner(false),
  _bg       (-numeric_limits<double>::infinity()),
  _bgSet    (false)
{
  UpdateMatrix();
}

// -----------------------------------------------------------------------------
irtkBaseImage::irtkBaseImage(const irtkImageAttributes &attr, int n)
:
  _mask     (NULL),
  _maskOwner(false),
  _bg       (-numeric_limits<double>::infinity()),
  _bgSet    (false)
{
  _attr = attr;
  if (n > 1) _attr._t = n, _attr._dt = .0; // i.e., vector image with n components
  PutAttributes(_attr);
}

// -----------------------------------------------------------------------------
irtkBaseImage::irtkBaseImage(const irtkBaseImage &image)
:
  irtkObject     (image),
  _attr          (image._attr),
  _NumberOfVoxels(image._NumberOfVoxels), 
  _matI2W        (image._matI2W),
  _matW2I        (image._matW2I),
  _mask          (image._maskOwner ? new irtkBinaryImage(*image._mask) : image._mask),
  _maskOwner     (image._maskOwner),
  _bg            (image._bg),
  _bgSet         (image._bgSet)
{
}

// -----------------------------------------------------------------------------
irtkBaseImage::~irtkBaseImage()
{
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
irtkBaseImage *irtkBaseImage::New(const char *fname)
{
  irtkFileToImage *reader = irtkFileToImage::New(fname);
  irtkBaseImage   *image  = reader->GetOutput();
  delete reader;
  return image;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline irtkBaseImage *NewImage(const irtkBaseImage *base)
{
  const irtkGenericImage<VoxelType> *image;
  if ((image = dynamic_cast<const irtkGenericImage<VoxelType> *>(base))) {
    return new irtkGenericImage<VoxelType>(*image);
  } else {
    cerr << "irtkBaseImage::New: Input image is not of expected type" << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
irtkBaseImage *irtkBaseImage::New(const irtkBaseImage *image)
{
  switch (image->GetDataType()) {
    case IRTK_VOXEL_CHAR:           return NewImage<char>          (image);
    case IRTK_VOXEL_UNSIGNED_CHAR:  return NewImage<unsigned char> (image);
    case IRTK_VOXEL_SHORT:          return NewImage<short>         (image);
    case IRTK_VOXEL_UNSIGNED_SHORT: return NewImage<unsigned short>(image);
    case IRTK_VOXEL_INT:            return NewImage<int>           (image);
    case IRTK_VOXEL_UNSIGNED_INT:   return NewImage<unsigned int>  (image);
    case IRTK_VOXEL_FLOAT:          return NewImage<float>         (image);
    case IRTK_VOXEL_DOUBLE:         return NewImage<double>        (image);
    default:
      cerr << "irtkBaseImage::New: Cannot allocate image of unknown type: "
           << image->GetDataType() << endl;
      exit(1);
  }
}

// -----------------------------------------------------------------------------
irtkBaseImage *irtkBaseImage::New(int dtype)
{
  switch (dtype) {
    case IRTK_VOXEL_CHAR:           return new irtkGenericImage<char>;
    case IRTK_VOXEL_UNSIGNED_CHAR:  return new irtkGenericImage<unsigned char>;
    case IRTK_VOXEL_SHORT:          return new irtkGenericImage<short>;
    case IRTK_VOXEL_UNSIGNED_SHORT: return new irtkGenericImage<unsigned short>;
    case IRTK_VOXEL_INT:            return new irtkGenericImage<int>;
    case IRTK_VOXEL_UNSIGNED_INT:   return new irtkGenericImage<unsigned int>;
    case IRTK_VOXEL_FLOAT:          return new irtkGenericImage<float>;
    case IRTK_VOXEL_DOUBLE:         return new irtkGenericImage<double>;
    default:
      cerr << "irtkBaseImage::New: Cannot allocate image of unknown type: " << dtype << endl;
      exit(1);
  }
}

// -----------------------------------------------------------------------------
irtkBaseImage *irtkBaseImage::Copy() const
{
  return irtkBaseImage::New(this);
}

// -----------------------------------------------------------------------------
void irtkBaseImage::UpdateMatrix()
{
  _matI2W = _attr.GetImageToWorldMatrix();
  _matW2I = _attr.GetWorldToImageMatrix();
}

// -----------------------------------------------------------------------------
void irtkBaseImage::PutAttributes(const irtkImageAttributes &attr)
{
  _attr           = attr;
  _NumberOfVoxels = attr.NumberOfLatticePoints();
  UpdateMatrix();
}

// -----------------------------------------------------------------------------
irtkBaseImage& irtkBaseImage::operator =(const irtkBaseImage &image)
{
  // Copy image attributes
  this->Initialize(image.Attributes());
  // Copy/cast image data
  for (int idx = 0; idx < _NumberOfVoxels; idx++) {
    this->PutAsDouble(idx, image.GetAsDouble(idx));
  }
  // Copy foreground region
  if (image.OwnsMask()) _mask = new irtkBinaryImage(*image.GetMask());
  else                  _mask = const_cast<irtkBinaryImage *>(image.GetMask());
  _maskOwner = image.OwnsMask();
  if (image.HasBackgroundValue()) {
    this->PutBackgroundValueAsDouble(image.GetBackgroundValueAsDouble());
  }
  return *this;
}

// =============================================================================
// Lattice
// =============================================================================

// -----------------------------------------------------------------------------
void irtkBaseImage::Orientation(int &i, int &j, int &k) const
{
#ifdef HAS_NIFTI
  mat44 m;
  for (int r = 0; r < 4; ++r) {
    for (int c = 0; c < 4; ++c) {
      m.m[r][c] = _matI2W(r, c); // FIXME: With or without affine transformation ?
    }
  }
  nifti_mat44_to_orientation(m, &i, &j, &k);
#else
  cerr << "irtkBaseImage::Orientation: Requires NIFTI support!"
          " Please recompile with NIFTI enabled" << endl;
  exit(1);
#endif
}

// -----------------------------------------------------------------------------
void irtkBaseImage::ImageToWorld(irtkWorldCoordsImage &i2w, bool _3D) const
{
  if (_attr._z == 1 && !_3D) {
    const int nvoxs = _attr._x * _attr._y;
    i2w.Initialize(_attr, 2);
    double *wx = i2w.GetPointerToVoxels();
    double *wy = wx + nvoxs;
    for (int j = 0; j < _attr._y; j++) {
      for (int i = 0; i < _attr._x; i++) {
        (*wx++) = _matI2W(0, 0) * i + _matI2W(0, 1) * j + _matI2W(0, 3);
        (*wy++) = _matI2W(1, 0) * i + _matI2W(1, 1) * j + _matI2W(1, 3);
      }
    }
  } else {
    const int nvoxs = _attr._x * _attr._y * _attr._z;
    i2w.Initialize(_attr, 3);
    double *wx = i2w.GetPointerToVoxels();
    double *wy = wx + nvoxs;
    double *wz = wy + nvoxs;
    for (int k = 0; k < _attr._z; k++) {
      for (int j = 0; j < _attr._y; j++) {
        for (int i = 0; i < _attr._x; i++) {
          (*wx++) = _matI2W(0, 0) * i + _matI2W(0, 1) * j + _matI2W(0, 2) * k + _matI2W(0, 3);
          (*wy++) = _matI2W(1, 0) * i + _matI2W(1, 1) * j + _matI2W(1, 2) * k + _matI2W(1, 3);
          (*wz++) = _matI2W(2, 0) * i + _matI2W(2, 1) * j + _matI2W(2, 2) * k + _matI2W(2, 3);
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
void irtkBaseImage::ImageToWorld(double *i2w, bool _3D) const
{
  if (_attr._z == 1 && !_3D) {
    for (int j = 0; j < _attr._y; j++) {
      for (int i = 0; i < _attr._x; i++) {
        (*i2w++) = _matI2W(0, 0) * i + _matI2W(0, 1) * j + _matI2W(0, 3);
        (*i2w++) = _matI2W(1, 0) * i + _matI2W(1, 1) * j + _matI2W(1, 3);
      }
    }
  } else {
    for (int k = 0; k < _attr._z; k++) {
      for (int j = 0; j < _attr._y; j++) {
        for (int i = 0; i < _attr._x; i++) {
          (*i2w++) = _matI2W(0, 0) * i + _matI2W(0, 1) * j + _matI2W(0, 2) * k + _matI2W(0, 3);
          (*i2w++) = _matI2W(1, 0) * i + _matI2W(1, 1) * j + _matI2W(1, 2) * k + _matI2W(1, 3);
          (*i2w++) = _matI2W(2, 0) * i + _matI2W(2, 1) * j + _matI2W(2, 2) * k + _matI2W(2, 3);
        }
      }
    }
  }
}

// =============================================================================
// Common image statistics
// =============================================================================

// -----------------------------------------------------------------------------
void irtkBaseImage::GetMinMaxAsDouble(double &min, double &max) const
{
  min = max = .0;

  bool   first = true;
  double value;

  for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
    if (IsForeground(idx)) {
      value = this->GetAsDouble(idx);
      if (first) {
        min   = max = value;
        first = false;
      } else {
        if (value < min) min = value;
        if (value > max) max = value;
      }
    }
  }
}

// -----------------------------------------------------------------------------
void irtkBaseImage::PutMinMaxAsDouble(double min, double max)
{
  // Check if full range can be represented by data type
  if (min < this->GetDataTypeMin()) {
    cerr << this->NameOfClass() << "::PutMinMaxAsDouble: "
         << "Requested minimum value is out of range for the voxel type: "
         << this->GetDataType() << endl;
    exit(1);
  }
  if (max > this->GetDataTypeMax()) {
    cerr << this->NameOfClass() << "::PutMinMaxAsDouble: "
         << "Requested maximum value is out of range for the voxel type: "
         << this->GetDataType() << endl;
    exit(1);
  }

  // Get current min and max values
  double min_val, max_val;
  this->GetMinMaxAsDouble(min_val, max_val);

  // Rescale foreground to desired [min, max] range
  const double slope = (max - min) / (max_val - min_val);
  const double inter = min - slope * min_val;
  for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
    if (IsForeground(idx)) {
      this->PutAsDouble(idx, inter + slope * this->GetAsDouble(idx));
    }
  }
}

// =============================================================================
// Foreground region
// =============================================================================

// -----------------------------------------------------------------------------
void irtkBaseImage::PutMask(irtkBinaryImage *mask, bool owner)
{
  if (_maskOwner) delete _mask;
  _mask      = mask;
  _maskOwner = owner;
}

// -----------------------------------------------------------------------------
void irtkBaseImage::InitializeMask(int t, bool force)
{
  if (!_mask) {
    _mask      = new irtkBinaryImage(_attr, t);
    _maskOwner = 2;
    force      = true;
  }
  if (force) {
    if (!_bgSet) {
      cerr << "irtkBaseImage::InitializeMask: Background value not set!" << endl;
      exit(1);
    }
    if (_mask->GetT() == _attr._t) {
      irtkBinaryPixel *ptr2msk = _mask->GetPointerToVoxels();
      for (int l = 0; l < _attr._t; l++) {
        for (int k = 0; k < _attr._z; k++) {
          for (int j = 0; j < _attr._y; j++) {
            for (int i = 0; i < _attr._x; i++) {
              (*ptr2msk++) = (this->GetAsDouble(i, j, k, l) == _bg ? false : true);
            }
          }
        }
      }
    } else {
      irtkBinaryPixel *ptr2msk = _mask->GetPointerToVoxels();
      for (int k = 0; k < _attr._z; k++) {
        for (int j = 0; j < _attr._y; j++) {
          for (int i = 0; i < _attr._x; i++) {
            (*ptr2msk++) = (this->GetAsDouble(i, j, k, 0) == _bg ? false : true);
          }
        }
      }
      for (int l = 1; l < _attr._t; l++) {
        irtkBinaryPixel *ptr2msk = _mask->GetPointerToVoxels();
        for (int k = 0; k < _attr._z; k++) {
          for (int j = 0; j < _attr._y; j++) {
            for (int i = 0; i < _attr._x; i++) {
              if (*ptr2msk == true) {
                (*ptr2msk) = (this->GetAsDouble(i, j, k, l) == _bg ? false : true);
              }
              ++ptr2msk;
            }
          }
        }
      }
    }
  }
  if (_maskOwner > 1) {
    if (_maskOwner == numeric_limits<int>::max()) {
      cerr << "irtkBaseImage::InitializeMask: Too many nested calls! Do also not forget to call ClearMask at some point." << endl;
      exit(1);
    }
    _maskOwner++;
  }
}

// -----------------------------------------------------------------------------
void irtkBaseImage::ClearMask(bool force)
{
  if (_maskOwner > 2) {
    _maskOwner = (force ? 2 : _maskOwner-1);
    if (_maskOwner == 2) delete _mask, _mask = NULL;
  }
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void irtkBaseImage::Print(irtkIndent indent) const
{
  cout << indent << "Image lattice:" << endl;
  _attr.Print(indent + 1);
  cout << indent << "Foreground mask:  " << ToString(_mask != NULL) << endl;
  cout << indent << "Background value: ";
  if (_bgSet) cout << _bg;
  else        cout << "n/a";
  cout << endl;
}
