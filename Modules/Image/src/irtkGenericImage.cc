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
#include <irtkImageToFile.h>
#include <irtkMatrix3x3.h>

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
// Note: Base class irtkBaseImage must be initialized before calling this function!
template <class VoxelType>
void irtkGenericImage<VoxelType>::AllocateImage(VoxelType *data)
{
  // Delete existing mask (if any)
  if (_maskOwner) Delete(_mask);
  // Free previously allocated memory
  Deallocate(_matrix, _data);
  if (_dataOwner) Deallocate(_data);
  _dataOwner = false;
  // Initialize memory
  const int nvox = _attr.NumberOfLatticePoints();
  if (nvox > 0) {
    if (data) {
      _data      = data;
      _dataOwner = false;
    } else {
      _data      = CAllocate<VoxelType>(nvox);
      _dataOwner = true;
    }
    Allocate(_matrix, _attr._x, _attr._y, _attr._z, _attr._t, _data);
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType>::irtkGenericImage()
:
  _matrix   (NULL),
  _data     (NULL),
  _dataOwner(false)
{
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType>::irtkGenericImage(const char *fname)
:
  _matrix   (NULL),
  _data     (NULL),
  _dataOwner(false)
{
  Read(fname);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType>::irtkGenericImage(int x, int y, int z, int t, VoxelType *data)
:
  _matrix   (NULL),
  _data     (NULL),
  _dataOwner(false)
{
  irtkImageAttributes attr;
  attr._x = x;
  attr._y = y;
  attr._z = z;
  attr._t = t;
  PutAttributes(attr);
  AllocateImage(data);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType>::irtkGenericImage(int x, int y, int z, int t, int n, VoxelType *data)
:
  _matrix   (NULL),
  _data     (NULL),
  _dataOwner(false)
{
  if (t > 1 && n > 1) {
    cerr << "irtkGenericImage::irtkGenericImage: 5D images not supported! Use 4D image with vector voxel type instead." << endl;
    exit(1);
  }
  irtkImageAttributes attr;
  if (n > 1) t = n, attr._dt = .0; // i.e., vector image with n components
  attr._x = x;
  attr._y = y;
  attr._z = z;
  attr._t = t;
  PutAttributes(attr);
  AllocateImage(data);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType>::irtkGenericImage(const irtkImageAttributes &attr, VoxelType *data)
:
  irtkBaseImage(attr),
  _matrix   (NULL),
  _data     (NULL),
  _dataOwner(false)
{
  AllocateImage(data);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType>::irtkGenericImage(const irtkImageAttributes &attr, int n, VoxelType *data)
:
  irtkBaseImage(attr, n),
  _matrix   (NULL),
  _data     (NULL),
  _dataOwner(false)
{
  AllocateImage(data);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType>::irtkGenericImage(const irtkBaseImage &image)
:
  irtkBaseImage(image),
  _matrix   (NULL),
  _data     (NULL),
  _dataOwner(false)
{
  // Initialize image
  AllocateImage();
  // Copy/cast data
  VoxelType *ptr = _data;
  for (int idx = 0; idx < _NumberOfVoxels; ++idx, ++ptr) {
    (*ptr) = voxel_cast<VoxelType>(image.GetAsVector(idx));
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType>::irtkGenericImage(const irtkGenericImage &image)
:
  irtkBaseImage(image),
  _matrix   (NULL),
  _data     (NULL),
  _dataOwner(false)
{
  if (image._dataOwner) {
    AllocateImage();
    memcpy(_data, image._data, _NumberOfVoxels * sizeof(VoxelType));
  } else {
    AllocateImage(const_cast<VoxelType *>(image.Data()));
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType> template <class VoxelType2>
irtkGenericImage<VoxelType>::irtkGenericImage(const irtkGenericImage<VoxelType2> &image)
:
  irtkBaseImage(image),
  _matrix   (NULL),
  _data     (NULL),
  _dataOwner(false)
{
  AllocateImage();
  VoxelType        *ptr1 = this->Data();
  const VoxelType2 *ptr2 = image.Data();
  for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
    ptr1[idx] = voxel_cast<VoxelType>(ptr2[idx]);
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType>::~irtkGenericImage()
{
  Deallocate(_matrix, _data);
  if (_dataOwner) Deallocate(_data);
  if (_maskOwner) Delete(_mask);
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkBaseImage *irtkGenericImage<VoxelType>::Copy() const
{
  return new irtkGenericImage<VoxelType>(*this);
}

// -----------------------------------------------------------------------------
template <class VoxelType> void irtkGenericImage<VoxelType>::Initialize()
{
  if (_matrix) *this = VoxelType();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>::Initialize(const irtkImageAttributes &a, int n, VoxelType *data)
{
  // Initialize attributes
  irtkImageAttributes attr(a);
  if (n >= 1) attr._t = n, attr._dt = .0; // i.e., vector image with n components
  // Initialize memory
  if (_attr._x != attr._x || _attr._y != attr._y || _attr._z != attr._z || _attr._t != attr._t) {
    PutAttributes(attr);
    AllocateImage(data);
  } else {
    PutAttributes(attr);
    if (_dataOwner) *this = VoxelType();
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>::Initialize(const irtkImageAttributes &attr, int n)
{
  this->Initialize(attr, n, NULL);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>::Initialize(const irtkImageAttributes &attr, VoxelType *data)
{
  this->Initialize(attr, -1, data);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>::Initialize(int x, int y, int z, int t, int n, VoxelType *data)
{
  irtkImageAttributes attr(_attr);
  attr._x = x;
  attr._y = y;
  attr._z = z;
  attr._t = t;
  this->Initialize(attr, n, data);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>::Initialize(int x, int y, int z, int t, VoxelType *data)
{
  this->Initialize(x, y, z, t, 1, data);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>::CopyFrom(const VoxelType *data)
{
  if (_data != data) {
    memcpy(_data, data, _NumberOfVoxels * sizeof(VoxelType));
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>::CopyFrom(const irtkBaseImage &image)
{
  for (int l = 0; l < _attr._t; l++) {
    for (int k = 0; k < _attr._z; k++) {
      for (int j = 0; j < _attr._y; j++) {
        for (int i = 0; i < _attr._x; i++) {
          _matrix[l][k][j][i] = voxel_cast<VoxelType>(image.GetAsVector(i, j, k, l));
        }
      }
    }
  }
  if (_maskOwner) delete _mask;
  if (image.OwnsMask()) {
    _mask      = new irtkBinaryImage(*image.GetMask());
    _maskOwner = true;
  } else {
    _mask      = const_cast<irtkBinaryImage *>(image.GetMask());
    _maskOwner = false;
  }
  if (image.HasBackgroundValue()) {
    this->PutBackgroundValueAsDouble(image.GetBackgroundValueAsDouble());
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>::CopyFrom(const irtkGenericImage &image)
{
  CopyFrom(image.Data());
  if (_maskOwner) delete _mask;
  if (image.OwnsMask()) {
    _mask      = new irtkBinaryImage(*image.GetMask());
    _maskOwner = true;
  } else {
    _mask      = const_cast<irtkBinaryImage *>(image.GetMask());
    _maskOwner = false;
  }
  if (image.HasBackgroundValue()) {
    this->PutBackgroundValueAsDouble(image.GetBackgroundValueAsDouble());
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType>& irtkGenericImage<VoxelType>::operator=(VoxelType scalar)
{
  VoxelType *ptr = this->Data();
  for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
    ptr[idx] = scalar;
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType>& irtkGenericImage<VoxelType>::operator=(const irtkBaseImage &image)
{
  if (this != &image) {
    this->Initialize(image.Attributes());
    this->CopyFrom(image);
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType>& irtkGenericImage<VoxelType>::operator=(const irtkGenericImage &image)
{
  if (this != &image) {
    this->Initialize(image.Attributes());
    this->CopyFrom(image);
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType> void irtkGenericImage<VoxelType>::Clear()
{
  Deallocate(_matrix, _data);
  if (_dataOwner) Deallocate(_data);
  if (_maskOwner) Delete(_mask);
  _attr = irtkImageAttributes();
}

// =============================================================================
// Region-of-interest extraction
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>
::GetRegion(irtkGenericImage<VoxelType> &image, int k, int m) const
{
  int i, j;
  double x1, y1, z1, t1, x2, y2, z2, t2;

  if ((k < 0) || (k >= _attr._z) || (m < 0) || (m >= _attr._t)) {
    cerr << "irtkGenericImage<VoxelType>::GetRegion: Parameter out of range" << endl;
    exit(1);
  }

  // Initialize
  irtkImageAttributes attr = this->Attributes();
  attr._z = 1;
  attr._t = 1;
  attr._xorigin = 0;
  attr._yorigin = 0;
  attr._zorigin = 0;
  image.Initialize(attr);

  // Calculate position of first voxel in roi in original image
  x1 = 0;
  y1 = 0;
  z1 = k;
  this->ImageToWorld(x1, y1, z1);
  t1 = this->ImageToTime(m);

  // Calculate position of first voxel in roi in new image
  x2 = 0;
  y2 = 0;
  z2 = 0;
  t2 = 0;
  image.ImageToWorld(x2, y2, z2);
  t2 = image.ImageToTime(0);

  // Shift origin of new image accordingly
  image.PutOrigin(x1 - x2, y1 - y2, z1 - z2, t1 - t2);

  // Copy region
  for (j = 0; j < _attr._y; j++) {
    for (i = 0; i < _attr._x; i++) {
      image._matrix[0][0][j][i] = _matrix[m][k][j][i];
    }
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>
::GetRegion(int k, int m) const
{
  irtkGenericImage<VoxelType> image;
  this->GetRegion(image, k, m);
  return image;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>
::GetRegion(irtkBaseImage *&base, int k, int m) const
{
  irtkGenericImage<VoxelType> *image = dynamic_cast<irtkGenericImage<VoxelType> *>(base);
  if (image == NULL) {
    delete base;
    image = new irtkGenericImage<VoxelType>();
    base  = image;
  }
  this->GetRegion(*image, k, m);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>
::GetRegion(irtkGenericImage<VoxelType> &image, int i1, int j1, int k1,
                                                int i2, int j2, int k2) const
{
  int i, j, k, l;
  double x1, y1, z1, x2, y2, z2;

  if ((i1 < 0) || (i1 >= i2) ||
      (j1 < 0) || (j1 >= j2) ||
      (k1 < 0) || (k1 >= k2) ||
      (i2 > _attr._x) || (j2 > _attr._y) || (k2 > _attr._z)) {
    cerr << "irtkGenericImage<VoxelType>::GetRegion: Parameter out of range\n";
    exit(1);
  }

  // Initialize
  irtkImageAttributes attr = this->Attributes();
  attr._x = i2 - i1;
  attr._y = j2 - j1;
  attr._z = k2 - k1;
  attr._xorigin = 0;
  attr._yorigin = 0;
  attr._zorigin = 0;
  image.Initialize(attr);

  // Calculate position of first voxel in roi in original image
  x1 = i1;
  y1 = j1;
  z1 = k1;
  this->ImageToWorld(x1, y1, z1);

  // Calculate position of first voxel in roi in new image
  x2 = 0;
  y2 = 0;
  z2 = 0;
  image.ImageToWorld(x2, y2, z2);

  // Shift origin of new image accordingly
  image.PutOrigin(x1 - x2, y1 - y2, z1 - z2);

  // Copy region
  for (l = 0; l < _attr._t; l++) {
    for (k = k1; k < k2; k++) {
      for (j = j1; j < j2; j++) {
        for (i = i1; i < i2; i++) {
          image._matrix[l][k-k1][j-j1][i-i1] = _matrix[l][k][j][i];
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>
::GetRegion(int i1, int j1, int k1, int i2, int j2, int k2) const
{
  irtkGenericImage<VoxelType> image;
  this->GetRegion(image, i1, j1, k1, i2, j2, k2);
  return image;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>
::GetRegion(irtkBaseImage *&base, int i1, int j1, int k1, int i2, int j2, int k2) const
{
  irtkGenericImage<VoxelType> *image = dynamic_cast<irtkGenericImage<VoxelType> *>(base);
  if (image == NULL) {
    delete base;
    image = new irtkGenericImage<VoxelType>();
    base  = image;
  }
  this->GetRegion(*image, i1, j1, k1, i2, j2, k2);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>
::GetRegion(irtkGenericImage<VoxelType> &image, int i1, int j1, int k1, int l1,
                                                int i2, int j2, int k2, int l2) const
{
  int i, j, k, l;
  double x1, y1, z1, x2, y2, z2;

  if ((i1 < 0) || (i1 >= i2) ||
      (j1 < 0) || (j1 >= j2) ||
      (k1 < 0) || (k1 >= k2) ||
      (l1 < 0) || (l1 >= l2) ||
      (i2 > _attr._x) || (j2 > _attr._y) || (k2 > _attr._z) || (l2 > _attr._t)) {
    cerr << "irtkGenericImage<VoxelType>::GetRegion: Parameter out of range\n";
    exit(1);
  }

  // Initialize
  irtkImageAttributes attr = this->Attributes();
  attr._x = i2 - i1;
  attr._y = j2 - j1;
  attr._z = k2 - k1;
  attr._t = l2 - l1;
  attr._xorigin = 0;
  attr._yorigin = 0;
  attr._zorigin = 0;
  image.Initialize(attr);

  // Calculate position of first voxel in roi in original image
  x1 = i1;
  y1 = j1;
  z1 = k1;
  this->ImageToWorld(x1, y1, z1);

  // Calculate position of first voxel in roi in new image
  x2 = 0;
  y2 = 0;
  z2 = 0;
  image.ImageToWorld(x2, y2, z2);

  // Shift origin of new image accordingly
  image.PutOrigin(x1 - x2, y1 - y2, z1 - z2);

  // Copy region
  for (l = l1; l < l2; l++) {
    for (k = k1; k < k2; k++) {
      for (j = j1; j < j2; j++) {
        for (i = i1; i < i2; i++) {
          image._matrix[l-l1][k-k1][j-j1][i-i1] = _matrix[l][k][j][i];
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>
::GetRegion(int i1, int j1, int k1, int l1, int i2, int j2, int k2, int l2) const
{
  irtkGenericImage<VoxelType> image;
  this->GetRegion(image, i1, j1, k1, l1, i2, j2, k2, l2);
  return image;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>
::GetRegion(irtkBaseImage *&base, int i1, int j1, int k1, int l1, int i2, int j2, int k2, int l2) const
{
  irtkGenericImage<VoxelType> *image = dynamic_cast<irtkGenericImage<VoxelType> *>(base);
  if (image == NULL) {
    delete base;
    image = new irtkGenericImage<VoxelType>();
    base  = image;
  }
  this->GetRegion(*image, i1, j1, k1, l1, i2, j2, k2, l2);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>::GetFrame(irtkGenericImage<VoxelType> &image, int l1, int l2) const
{
  if (l2 < 0) l2 = l1;

  if ((l2 < 0) || (l1 >= _attr._t)) {
    cerr << "irtkGenericImage<VoxelType>::GetFrame: Parameter out of range\n";
    exit(1);
  }

  if (l1 < 0) l1 = 0;
  if (l2 >= _attr._t) l2 = _attr._t - 1;

  // Initialize
  irtkImageAttributes attr = this->Attributes();
  attr._t       = l2 - l1 + 1;
  attr._torigin = this->ImageToTime(l1);
  image.Initialize(attr);

  // Copy region
  for (int l = l1; l <= l2; l++) {
    for (int k = 0; k < _attr._z; k++) {
      for (int j = 0; j < _attr._y; j++) {
        for (int i = 0; i < _attr._x; i++) {
          image._matrix[l-l1][k][j][i] = _matrix[l][k][j][i];
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>::GetFrame(int l1, int l2) const
{
  irtkGenericImage<VoxelType> image;
  this->GetFrame(image, l1, l2);
  return image;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>::GetFrame(irtkBaseImage *&base, int l1, int l2) const
{
  irtkGenericImage<VoxelType> *image = dynamic_cast<irtkGenericImage<VoxelType> *>(base);
  if (image == NULL) {
    delete base;
    image = new irtkGenericImage<VoxelType>();
    base  = image;
  }
  this->GetFrame(*image, l1, l2);
}

// =============================================================================
// Image arithmetic
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType>& irtkGenericImage<VoxelType>::operator+=(const irtkGenericImage &image)
{
  if (image.Attributes() != this->Attributes()) {
    cerr << "irtkGenericImage<VoxelType>::operator+=: Size mismatch in images" << endl;
    this->Attributes().Print();
    image.Attributes().Print();
    exit(1);
  }
  VoxelType       *ptr1 = this->Data();
  const VoxelType *ptr2 = image.Data();
  for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
    if (IsForeground(idx) && image.IsForeground(idx)) ptr1[idx] += ptr2[idx];
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType>& irtkGenericImage<VoxelType>::operator-=(const irtkGenericImage &image)
{
  if (image.Attributes() != this->Attributes()) {
    cerr << "irtkGenericImage<VoxelType>::operator-=: Size mismatch in images" << endl;
    this->Attributes().Print();
    image.Attributes().Print();
    exit(1);
  }
  VoxelType       *ptr1 = this->Data();
  const VoxelType *ptr2 = image.Data();
  for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
    if (IsForeground(idx) && image.IsForeground(idx)) ptr1[idx] -= ptr2[idx];
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType>& irtkGenericImage<VoxelType>::operator*=(const irtkGenericImage &image)
{
  if (image.Attributes() != this->Attributes()) {
    cerr << "irtkGenericImage<VoxelType>::operator*=: Size mismatch in images" << endl;
    this->Attributes().Print();
    image.Attributes().Print();
    exit(1);
  }
  VoxelType       *ptr1 = this->Data();
  const VoxelType *ptr2 = image.Data();
  for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
    if (IsForeground(idx) && image.IsForeground(idx)) ptr1[idx] *= ptr2[idx];
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType>& irtkGenericImage<VoxelType>::operator/=(const irtkGenericImage &image)
{
  if (image.Attributes() != this->Attributes()) {
    cerr << "irtkGenericImage<VoxelType>::operator/=: Size mismatch in images" << endl;
    this->Attributes().Print();
    image.Attributes().Print();
    exit(1);
  }
  VoxelType       *ptr1 = this->Data();
  const VoxelType *ptr2 = image.Data();
  for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
    if (IsForeground(idx) && image.IsForeground(idx)) {
      if (ptr2[idx] == VoxelType()) {
        if (HasBackgroundValue()) {
          ptr1[idx] = voxel_cast<VoxelType>(GetBackgroundValueAsDouble());
        } else {
          ptr1[idx] = VoxelType();
        }
      } else {
        ptr1[idx] /= ptr2[idx];
      }
    }
  }
  return *this;
}

template <> irtkGenericImage<float3x3 > &irtkGenericImage<float3x3 >::operator/=(const irtkGenericImage &)
{
  cerr << "irtkGenericImage<float3x3>::operator /=: Not implemented" << endl;
  exit(1);
}

template <> irtkGenericImage<double3x3> &irtkGenericImage<double3x3>::operator/=(const irtkGenericImage &)
{
  cerr << "irtkGenericImage<double3x3>::operator /=: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType>& irtkGenericImage<VoxelType>::operator+=(ScalarType scalar)
{
  VoxelType *ptr = this->Data();
  for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
    if (IsForeground(idx)) ptr[idx] += scalar;
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType>& irtkGenericImage<VoxelType>::operator-=(ScalarType scalar)
{
  VoxelType *ptr = this->Data();
  for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
    if (IsForeground(idx)) ptr[idx] -= scalar;
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType>& irtkGenericImage<VoxelType>::operator*=(ScalarType scalar)
{
  VoxelType *ptr = this->Data();
  for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
    if (IsForeground(idx)) ptr[idx] *= scalar;
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType>& irtkGenericImage<VoxelType>::operator/=(ScalarType scalar)
{
  if (scalar) {
    VoxelType *ptr = this->Data();
    for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
      if (IsForeground(idx)) ptr[idx] /= scalar;
    }
  } else {
    cerr << "irtkGenericImage<VoxelType>::operator/=: Division by zero" << endl;
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>::operator+(const irtkGenericImage &image) const
{
  irtkGenericImage<VoxelType> tmp(*this);
  tmp += image;
  return tmp;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>::operator-(const irtkGenericImage &image) const
{
  irtkGenericImage<VoxelType> tmp(*this);
  tmp -= image;
  return tmp;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>::operator*(const irtkGenericImage &image) const
{
  irtkGenericImage<VoxelType> tmp(*this);
  tmp *= image;
  return tmp;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>::operator/(const irtkGenericImage &image) const
{
  irtkGenericImage<VoxelType> tmp(*this);
  tmp /= image;
  return tmp;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>::operator+(ScalarType scalar) const
{
  irtkGenericImage<VoxelType> tmp(*this);
  tmp += scalar;
  return tmp;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>::operator-(ScalarType scalar) const
{
  irtkGenericImage<VoxelType> tmp(*this);
  tmp -= scalar;
  return tmp;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>::operator*(ScalarType scalar) const
{
  irtkGenericImage<VoxelType> tmp(*this);
  tmp *= scalar;
  return tmp;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>::operator/(ScalarType scalar) const
{
  irtkGenericImage<VoxelType> tmp(*this);
  tmp /= scalar;
  return tmp;
}

// =============================================================================
// Thresholding
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>::PutBackgroundValueAsDouble(double value, bool threshold)
{
  irtkBaseImage::PutBackgroundValueAsDouble(value);
  if (threshold) {
    const VoxelType bg = voxel_cast<VoxelType>(this->_bg);
    VoxelType *ptr = this->GetPointerToVoxels();
    for (int idx = 0; idx < _NumberOfVoxels; ++idx, ++ptr) {
      if (*ptr < bg) *ptr = bg;
    }
  }
}

template <> void irtkGenericImage<float3x3>::PutBackgroundValueAsDouble(double value, bool threshold)
{
  irtkBaseImage::PutBackgroundValueAsDouble(value, threshold);
}

template <> void irtkGenericImage<double3x3>::PutBackgroundValueAsDouble(double value, bool threshold)
{
  irtkBaseImage::PutBackgroundValueAsDouble(value, threshold);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType> &irtkGenericImage<VoxelType>::operator>=(VoxelType pixel)
{
  VoxelType *ptr = this->GetPointerToVoxels();
  for (int idx = 0; idx < _NumberOfVoxels; idx++) {
    if (IsForeground(idx) && ptr[idx] > pixel) ptr[idx] = pixel;
  }
  return *this;
}

template <> irtkGenericImage<float3x3> &irtkGenericImage<float3x3>::operator>=(float3x3)
{
  cerr << "irtkGenericImage<float3x3 >::operator >=: Not implemented" << endl;
  exit(1);
}

template <> irtkGenericImage<double3x3> &irtkGenericImage<double3x3>::operator>=(double3x3)
{
  cerr << "irtkGenericImage<double3x3 >::operator >=: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType>& irtkGenericImage<VoxelType>::operator<=(VoxelType pixel)
{
  VoxelType *ptr = this->GetPointerToVoxels();
  for (int idx = 0; idx < _NumberOfVoxels; idx++) {
    if (IsForeground(idx) && ptr[idx] < pixel) ptr[idx] = pixel;
  }
  return *this;
}

template <> irtkGenericImage<float3x3> &irtkGenericImage<float3x3>::operator<=(float3x3)
{
  cerr << "irtkGenericImage<float3x3 >::operator <=: Not implemented" << endl;
  exit(1);
}

template <> irtkGenericImage<double3x3> &irtkGenericImage<double3x3>::operator<=(double3x3)
{
  cerr << "irtkGenericImage<double3x3 >::operator <=: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>::operator>(VoxelType pixel) const
{
  irtkGenericImage<VoxelType> image(*this);
  image >= pixel;
  return image;
}

template <> irtkGenericImage<float3x3> irtkGenericImage<float3x3>::operator>(float3x3) const
{
  cerr << "irtkGenericImage<float3x3>::operator >: Not implemented" << endl;
  exit(1);
}

template <> irtkGenericImage<double3x3> irtkGenericImage<double3x3>::operator>(double3x3) const
{
  cerr << "irtkGenericImage<double3x3>::operator >: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType> irtkGenericImage<VoxelType>::operator<(VoxelType pixel) const
{
  irtkGenericImage<VoxelType> image(*this);
  image <= pixel;
  return image;
}

template <> irtkGenericImage<float3x3> irtkGenericImage<float3x3>::operator<(float3x3 ) const
{
  cerr << "irtkGenericImage<float3x3>::operator <: Not implemented" << endl;
  exit(1);
}

template <> irtkGenericImage<double3x3> irtkGenericImage<double3x3>::operator<(double3x3) const
{
  cerr << "irtkGenericImage<double3x3>::operator <: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkBinaryImage irtkGenericImage<VoxelType>::operator!=(VoxelType pixel) const
{
  irtkBinaryImage mask(_attr);
  const VoxelType *ptr1 = this->GetPointerToVoxels();
  irtkBinaryPixel *ptr2 = mask .GetPointerToVoxels();
  for (int idx = 0; idx < _NumberOfVoxels; ++idx, ++ptr1, ++ptr2) {
    *ptr2 = (*ptr1 != pixel);
  }
  return mask;
}

// =============================================================================
// Common image statistics
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>::GetMinMax(VoxelType &min, VoxelType &max) const
{
  min = max = VoxelType();
  
  const VoxelType *ptr   = this->Data();
  bool             first = true;
  
  for (int idx = 0; idx < _NumberOfVoxels; ++idx, ++ptr) {
    if (IsForeground(idx)) {
      if (first) {
        min   = max = *ptr;
        first = false;
      } else {
        if (*ptr < min) min = *ptr;
        if (*ptr > max) max = *ptr;
      }
    }
  }
}

template <> void irtkGenericImage<float3x3 >::GetMinMax(VoxelType &, VoxelType &) const
{
  cerr << "irtkGenericImage<float3x3>::GetMinMax: Not implemented" << endl;
  exit(1);
}

template <> void irtkGenericImage<double3x3>::GetMinMax(VoxelType &, VoxelType &) const
{
  cerr << "irtkGenericImage<double3x3>::GetMinMax: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>::GetMinMax(VoxelType &min, VoxelType &max, VoxelType pad) const
{
  min = max = VoxelType();

  const VoxelType *ptr   = this->Data();
  bool             first = true;
  
  for (int idx = 0; idx < _NumberOfVoxels; ++idx, ++ptr) {
    if (*ptr != pad) {
      if (first) {
        min   = max = *ptr;
        first = false;
      } else {
        if (*ptr < min) min = *ptr;
        if (*ptr > max) max = *ptr;
      }
    }
  }
}

template <> void irtkGenericImage<float3x3 >::GetMinMax(VoxelType &, VoxelType &, VoxelType) const
{
  cerr << "irtkGenericImage<float3x3>::GetMinMax: Not implemented" << endl;
  exit(1);
}

template <> void irtkGenericImage<double3x3>::GetMinMax(VoxelType &, VoxelType &, VoxelType) const
{
  cerr << "irtkGenericImage<double3x3>::GetMinMax: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>::PutMinMax(VoxelType min, VoxelType max)
{
  VoxelType min_val, max_val;
  this->GetMinMax(min_val, max_val);
  VoxelType *ptr = this->Data();
  RealType slope = voxel_cast<RealType>(max  - min)  / voxel_cast<RealType>(max_val - min_val);
  RealType inter = voxel_cast<RealType>(min) - slope * voxel_cast<RealType>(min_val);
  for (int idx = 0; idx < _NumberOfVoxels; ++idx, ++ptr) {
    if (IsForeground(idx)) *ptr = static_cast<VoxelType>(inter + slope * static_cast<RealType>(*ptr));
  }
}

template <> void irtkGenericImage<float3x3 >::PutMinMax(VoxelType, VoxelType)
{
  cerr << "irtkGenericImage<float3x3>::PutMinMax: Not implemented" << endl;
  exit(1);
}

template <> void irtkGenericImage<double3x3>::PutMinMax(VoxelType, VoxelType)
{
  cerr << "irtkGenericImage<double3x3>::PutMinMax: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
typename irtkGenericImage<VoxelType>::RealType
irtkGenericImage<VoxelType>::GetAverage(int toggle) const
{
  const VoxelType  zero = voxel_cast<VoxelType>(0);
  RealType         avg  = voxel_cast<RealType >(0);
  const VoxelType *ptr;

  if (toggle) {
    int n = 0;
    ptr = this->Data();
    for (int i = 0; i < _NumberOfVoxels; i++) {
      if (IsForeground(i) && (*ptr) > zero) n++;
      ++ptr;
    }
    ptr = this->Data();
    for (int i = 0; i < _NumberOfVoxels; i++) {
      if (IsForeground(i) && (*ptr) > zero) {
        avg += voxel_cast<RealType>(*ptr) / static_cast<double>(n);
      }
      ++ptr;
    }
  } else {
    ptr = this->Data();
    for (int i = 0; i < _NumberOfVoxels; i++) {
      if (IsForeground(i)) {
        avg += voxel_cast<RealType>(*ptr) / static_cast<double>(_NumberOfVoxels);
      }
      ++ptr;
    }
  }

  return avg;
}

template <> typename irtkGenericImage<float3x3 >::RealType irtkGenericImage<float3x3 >::GetAverage(int) const
{
  cerr << "irtkGenericImage<float3x3>::GetAverage: Not implemented" << endl;
  exit(1);
}

template <> typename irtkGenericImage<double3x3>::RealType irtkGenericImage<double3x3>::GetAverage(int) const
{
  cerr << "irtkGenericImage<double3x3>::GetAverage: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
typename irtkGenericImage<VoxelType>::RealType
irtkGenericImage<VoxelType>::GetSD(int toggle) const
{
  const VoxelType  zero = voxel_cast<VoxelType>(0);
  RealType         std  = voxel_cast<RealType >(0);
  const RealType   avg  = this->GetAverage(toggle);
  const VoxelType *ptr;

  if (toggle) {
    int n = 0;
    ptr = this->Data();
    for (int i = 0; i < _NumberOfVoxels; i++) {
      if (IsForeground(i) && (*ptr) > zero) n++;
      ++ptr;
    }
    ptr = this->Data();
    for (int i = 0; i < _NumberOfVoxels; i++) {
      if (IsForeground(i) && (*ptr) > zero) {
        std += pow(voxel_cast<RealType>(*ptr) - avg, 2) / static_cast<double>(n);
      }
      ++ptr;
    }
  } else {
    ptr = this->Data();
    for (int i = 0; i < _NumberOfVoxels; i++) {
      if (IsForeground(i)) {
        std += pow(voxel_cast<RealType>(*ptr) - avg, 2) / static_cast<double>(_NumberOfVoxels);
      }
      ++ptr;
    }
  }

  return sqrt(std);
}

template <> typename irtkGenericImage<float3x3 >::RealType irtkGenericImage<float3x3 >::GetSD(int) const
{
  cerr << "irtkGenericImage<float3x3>::GetSD: Not implemented" << endl;
  exit(1);
}

template <> typename irtkGenericImage<double3x3>::RealType irtkGenericImage<double3x3>::GetSD(int) const
{
  cerr << "irtkGenericImage<double3x3>::GetSD: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>::GetMaxPosition(irtkPoint& p, int ds, int t) const
{
  if (ds < 0) ds = -ds;
  
  this->WorldToImage(p);
  const int x = static_cast<int>(round(p._x));
  const int y = static_cast<int>(round(p._y));
  const int z = static_cast<int>(round(p._z));
  
  const VoxelType *ptr = this->Data(0, 0, 0, t);
  VoxelType        max = this->Get(x, y, z, t);
  
  p._z = static_cast<double>(z);
  for (int j = y - ds; j <= y + ds; j++) {
    for (int i = x - ds; i <= x + ds; i++) {
      if (IsForeground(i, j, z, t) && *ptr > max) {
        p._x = static_cast<double>(i);
        p._y = static_cast<double>(j);
        max  = *ptr;
      }
      ++ptr;
    }
  }
  
  this->ImageToWorld(p);
}

template <> void irtkGenericImage<float3x3 >::GetMaxPosition(irtkPoint &, int, int) const
{
  cerr << "irtkGenericImage<float3x3>::GetMaxPosition: Not implemented" << endl;
  exit(1);
}

template <> void irtkGenericImage<double3x3>::GetMaxPosition(irtkPoint &, int, int) const
{
  cerr << "irtkGenericImage<double3x3>::GetMaxPosition: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>::GravityCenter(irtkPoint& p, int ds, int t) const
{
  if (ds < 0) ds = -ds;

  this->WorldToImage(p);
  const int x = static_cast<int>(round(p._x));
  const int y = static_cast<int>(round(p._y));
  const int z = static_cast<int>(round(p._z));

  double si = .0, sj = .0, sk = .0, sw = .0, ssw = .0;
  for (int k = z - ds; k <= z + ds; k++) {
    for (int j = y - ds; j <= y + ds; j++) {
      for (int i = x - ds; i <= x + ds; i++) {
        if (IsForeground(i, j, k, t)) {
          sw   = this->GetAsDouble(i, j, k, t);
          si  += i * sw;
          sj  += j * sw;
          sk  += k * sw;
          ssw +=     sw;
        }
      }
    }
  }

  if (ssw) p._x = si/ssw, p._y = sj/ssw, p._z = sk/ssw;
  this->ImageToWorld(p);
}

template <> void irtkGenericImage<float3x3 >::GravityCenter(irtkPoint &, int, int) const
{
  cerr << "irtkGenericImage<float3x3>::GravityCenter: Not implemented" << endl;
  exit(1);
}

template <> void irtkGenericImage<double3x3>::GravityCenter(irtkPoint &, int, int) const
{
  cerr << "irtkGenericImage<double3x3>::GravityCenter: Not implemented" << endl;
  exit(1);
}

// =============================================================================
// Common image manipulations
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>::ReflectX()
{
  for (int t = 0; t < _attr._t; ++t) {
    for (int z = 0; z < _attr._z; ++z) {
      for (int y = 0; y < _attr._y; ++y) {
        for (int x = 0; x < _attr._x / 2; ++x) {
          swap(_matrix[t][z][y][x], _matrix[t][z][y][_attr._x-(x+1)]);
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>::ReflectY()
{
  for (int t = 0; t < _attr._t; ++t) {
    for (int z = 0; z < _attr._z; ++z) {
      for (int y = 0; y < _attr._y / 2; ++y) {
        for (int x = 0; x < _attr._x; ++x) {
          swap(_matrix[t][z][y][x], _matrix[t][z][_attr._y-(y+1)][x]);
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>::ReflectZ()
{
  for (int t = 0; t < _attr._t; ++t) {
    for (int z = 0; z < _attr._z / 2; ++z) {
      for (int y = 0; y < _attr._y; ++y) {
        for (int x = 0; x < _attr._x; ++x) {
          swap(_matrix[t][z][y][x], _matrix[t][_attr._z-(z+1)][y][x]);
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>::FlipXY(int modifyOrigin)
{
  // TODO: Implement irtkBaseImage::FlipXY which flips the foreground mask (if any),
  //       adjusts the attributes, and updates the coordinate transformation matrices.
  //       The subclass then only needs to reshape the image _matrix data itself.

  // Allocate memory
  VoxelType ****matrix = Allocate<VoxelType>(_attr._y, _attr._x, _attr._z, _attr._t);

  // Flip image in memory
  for (int l = 0; l < _attr._t; l++) {
    for (int k = 0; k < _attr._z; k++) {
      for (int j = 0; j < _attr._y; j++) {
        for (int i = 0; i < _attr._x; i++) {
          matrix[l][k][i][j] = _matrix[l][k][j][i];
        }
      }
    }
  }

  // Swap image dimensions
  swap(_attr._x, _attr._y);

  // Swap voxel dimensions
  swap(_attr._dx, _attr._dy);

  // Swap origin coordinates
  if (modifyOrigin > 0) swap(_attr._xorigin, _attr._yorigin);

  // Reshape image matrix
  //
  // Attention: DO NOT just swap the pointers to the data elements as this
  //            changes the memory location of the image data. This is not
  //            predictable by users of the class which may still hold a pointer
  //            to the old memory and in particular complicates the synchronization
  //            of host and device memory in irtkCUGenericImage used by CUDA code.
  _matrix = Reshape(_matrix, _attr._x, _attr._y, _attr._z, _attr._t);

  // Copy flipped image
  CopyFrom(matrix[0][0][0]);

  // Deallocate memory
  Deallocate(matrix);

  // Update coordinate transformation
  UpdateMatrix();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>::FlipXZ(int modifyOrigin)
{
  // TODO: Implement irtkBaseImage::FlipXZ which flips the foreground mask (if any),
  //       adjusts the attributes, and updates the coordinate transformation matrices.
  //       The subclass then only needs to reshape the image _matrix data itself.

  // Allocate memory
  VoxelType ****matrix = Allocate<VoxelType>(_attr._z, _attr._y, _attr._x, _attr._t);

  // Flip image in memory
  for (int l = 0; l < _attr._t; l++) {
    for (int k = 0; k < _attr._z; k++) {
      for (int j = 0; j < _attr._y; j++) {
        for (int i = 0; i < _attr._x; i++) {
          matrix[l][i][j][k] = _matrix[l][k][j][i];
        }
      }
    }
  }
 
  // Swap image dimensions
  swap(_attr._x, _attr._z);

  // Swap voxel dimensions
  swap(_attr._dx, _attr._dz);

  // Swap origin coordinates
  if (modifyOrigin > 0) swap(_attr._xorigin, _attr._zorigin);

  // Reshape image data
  //
  // Attention: DO NOT just swap the pointers to the data elements as this
  //            changes the memory location of the image data. This is not
  //            predictable by users of the class which may still hold a pointer
  //            to the old memory and in particular complicates the synchronization
  //            of host and device memory in irtkCUGenericImage used by CUDA code.
  _matrix = Reshape(_matrix, _attr._x, _attr._y, _attr._z, _attr._t);
  
  // Copy flipped image
  CopyFrom(matrix[0][0][0]);

  // Deallocate memory
  Deallocate(matrix);

  // Update coordinate transformation
  UpdateMatrix();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>::FlipYZ(int modifyOrigin)
{
  // TODO: Implement irtkBaseImage::FlipYZ which flips the foreground mask (if any),
  //       adjusts the attributes, and updates the coordinate transformation matrices.
  //       The subclass then only needs to reshape the image _matrix data itself.

  // Allocate memory for flipped image
  VoxelType ****matrix = Allocate<VoxelType>(_attr._x, _attr._z, _attr._y, _attr._t);

  // Flip image in memory  
  for (int l = 0; l < _attr._t; l++) {
    for (int k = 0; k < _attr._z; k++) {
      for (int j = 0; j < _attr._y; j++) {
        for (int i = 0; i < _attr._x; i++) {
          matrix[l][j][k][i] = _matrix[l][k][j][i];
        }
      }
    }
  }

  // Swap image dimensions
  swap(_attr._y, _attr._z);

  // Swap voxel dimensions
  swap(_attr._dy, _attr._dz);

  // Swap origin coordinates
  if (modifyOrigin > 0) swap(_attr._yorigin, _attr._zorigin);

  // Reshape image matrix
  //
  // Attention: DO NOT just swap the pointers to the data elements as this
  //            changes the memory location of the image data. This is not
  //            predictable by users of the class which may still hold a pointer
  //            to the old memory and in particular complicates the synchronization
  //            of host and device memory in irtkCUGenericImage used by CUDA code.
  _matrix = Reshape(_matrix, _attr._x, _attr._y, _attr._z, _attr._t);
  
  // Copy flipped image
  CopyFrom(matrix[0][0][0]);

  // Deallocate memory
  Deallocate(matrix);

  // Update coordinate transformation
  UpdateMatrix();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>::FlipXT(int modifyOrigin)
{
  // TODO: Implement irtkBaseImage::FlipXT which flips the foreground mask (if any),
  //       adjusts the attributes, and updates the coordinate transformation matrices.
  //       The subclass then only needs to reshape the image _matrix data itself.

  // Allocate memory
  VoxelType ****matrix = Allocate<VoxelType>(_attr._t, _attr._y, _attr._z, _attr._x);

  // Flip image in memory
  for (int l = 0; l < _attr._t; l++) {
    for (int k = 0; k < _attr._z; k++) {
      for (int j = 0; j < _attr._y; j++) {
        for (int i = 0; i < _attr._x; i++) {
          matrix[i][k][j][l] = _matrix[l][k][j][i];
        }
      }
    }
  }

  // Swap image dimensions
  swap(_attr._x, _attr._t);

  // Swap voxel dimensions
  swap(_attr._dx, _attr._dt);

  // Swap origin coordinates
  if (modifyOrigin > 0) swap(_attr._xorigin, _attr._torigin);

  // Reshape image matrix
  //
  // Attention: DO NOT just swap the pointers to the data elements as this
  //            changes the memory location of the image data. This is not
  //            predictable by users of the class which may still hold a pointer
  //            to the old memory and in particular complicates the synchronization
  //            of host and device memory in irtkCUGenericImage used by CUDA code.
  _matrix = Reshape(_matrix, _attr._x, _attr._y, _attr._z, _attr._t);
  
  // Copy flipped image
  CopyFrom(matrix[0][0][0]);

  // Deallocate memory
  Deallocate(matrix);

  // Update coordinate transformation
  UpdateMatrix();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>::FlipYT(int modifyOrigin)
{
  // TODO: Implement irtkBaseImage::FlipYT which flips the foreground mask (if any),
  //       adjusts the attributes, and updates the coordinate transformation matrices.
  //       The subclass then only needs to reshape the image _matrix data itself.

  // Allocate memory
  VoxelType ****matrix = Allocate<VoxelType>(_attr._x, _attr._t, _attr._z, _attr._y);

  for (int l = 0; l < _attr._t; l++) {
    for (int k = 0; k < _attr._z; k++) {
      for (int j = 0; j < _attr._y; j++) {
        for (int i = 0; i < _attr._x; i++) {
          matrix[j][k][l][i] = _matrix[l][k][j][i];
        }
      }
    }
  }

  // Swap image dimensions
  swap(_attr._y, _attr._t);

  // Swap voxel dimensions
  swap(_attr._dy, _attr._dt);

  // Swap origin coordinates
  if (modifyOrigin > 0) swap(_attr._yorigin, _attr._torigin);

  // Reshape image matrix
  //
  // Attention: DO NOT just swap the pointers to the data elements as this
  //            changes the memory location of the image data. This is not
  //            predictable by users of the class which may still hold a pointer
  //            to the old memory and in particular complicates the synchronization
  //            of host and device memory in irtkCUGenericImage used by CUDA code.
  _matrix = Reshape(_matrix, _attr._x, _attr._y, _attr._z, _attr._t);
  
  // Copy flipped image
  CopyFrom(matrix[0][0][0]);

  // Deallocate memory
  Deallocate(matrix);

  // Update coordinate transformation
  UpdateMatrix();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>::FlipZT(int modifyOrigin)
{
  // TODO: Implement irtkBaseImage::FlipZT which flips the foreground mask (if any),
  //       adjusts the attributes, and updates the coordinate transformation matrices.
  //       The subclass then only needs to reshape the image _matrix data itself.

  // Allocate memory
  VoxelType ****matrix = Allocate<VoxelType>(_attr._x, _attr._y, _attr._t, _attr._z);

  for (int l = 0; l < _attr._t; l++) {
    for (int k = 0; k < _attr._z; k++) {
      for (int j = 0; j < _attr._y; j++) {
        for (int i = 0; i < _attr._x; i++) {
          matrix[k][l][j][i] = _matrix[l][k][j][i];
        }
      }
    }
  }

  // Swap image dimensions
  swap(_attr._z, _attr._t);

  // Swap voxel dimensions
  swap(_attr._dz, _attr._dt);

  // Swap origin coordinates
  if (modifyOrigin > 0) swap(_attr._zorigin, _attr._torigin);

  // Reshape image matrix
  //
  // Attention: DO NOT just swap the pointers to the data elements as this
  //            changes the memory location of the image data. This is not
  //            predictable by users of the class which may still hold a pointer
  //            to the old memory and in particular complicates the synchronization
  //            of host and device memory in irtkCUGenericImage used by CUDA code.
  _matrix = Reshape(_matrix, _attr._x, _attr._y, _attr._z, _attr._t);
  
  // Copy flipped image
  CopyFrom(matrix[0][0][0]);

  // Deallocate memory
  Deallocate(matrix);

  // Update coordinate transformation
  UpdateMatrix();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
bool irtkGenericImage<VoxelType>::CropPad(int margin)
{
  irtkImageAttributes attr = _attr;
  // Determine lower bound along x axis: i1
  int i1 = _attr._x;
  for (int l = 0; l < _attr._t; ++l) {
    for (int k = 0; k < _attr._z; ++k) {
      for (int j = 0; j < _attr._y; ++j) {
        for (int i = 0; i < _attr._x; ++i) {
          if (IsForeground(i, j, k, l)) {
            if (i < i1) i1 = i;
            break;
          }
        }
      }
    }
  }
  // Determine upper bound along x axis: i2
  int i2 = -1;
  for (int l = 0; l < _attr._t; ++l) {
    for (int k = 0; k < _attr._z; ++k) {
      for (int j = 0; j < _attr._y; ++j) {
        for (int i = _attr._x - 1; i >= i1; --i) {
          if (IsForeground(i, j, k, l)) {
            if (i > i2) i2 = i;
            break;
          }
        }
      }
    }
  }
  // Determine lower bound along y axis: j1
  int j1 = _attr._y;
  for (int l = 0; l < _attr._t; ++l) {
    for (int k = 0; k < _attr._z; ++k) {
      for (int i = i1; i <= i2; ++i) {
        for (int j = 0; j < _attr._y; ++j) {
          if (IsForeground(i, j, k, l)) {
            if (j < j1) j1 = j;
            break;
          }
        }
      }
    }
  }
  // Determine upper bound along y axis: j2
  int j2 = -1;
  for (int l = 0; l < _attr._t; ++l) {
    for (int k = 0; k < _attr._z; ++k) {
      for (int i = i1; i <= i2; ++i) {
        for (int j = _attr._y - 1; j >= j1; --j) {
          if (IsForeground(i, j, k, l)) {
            if (j > j2) j2 = j;
            break;
          }
        }
      }
    }
  }
  // Determine lower bound along z axis: k1
  int k1 = _attr._z;
  for (int l = 0; l < _attr._t; ++l) {
    for (int j = j1; j <= j2; ++j) {
      for (int i = i1; i <= i2; ++i) {
        for (int k = 0; k < _attr._z; ++k) {
          if (IsForeground(i, j, k, l)) {
            if (k < k1) k1 = k;
            break;
          }
        }
      }
    }
  }
  // Determine upper bound along z axis: k2
  int k2 = -1;
  for (int l = 0; l < _attr._t; ++l) {
    for (int j = j1; j <= j2; ++j) {
      for (int i = i1; i <= i2; ++i) {
        for (int k = _attr._z - 1; k >= k1; --k) {
          if (IsForeground(i, j, k, l)) {
            if (k > k2) k2 = k;
            break;
          }
        }
      }
    }
  }
  // Determine lower bound along t axis: l1
  int l1 = _attr._t;
  for (int k = k1; k <= k2; ++k) {
    for (int j = j1; j <= j2; ++j) {
      for (int i = i1; i <= i2; ++i) {
        for (int l = 0; l < _attr._t; ++l) {
          if (IsForeground(i, j, k, l)) {
            if (l < l1) l1 = l;
            break;
          }
        }
      }
    }
  }
  // Determine upper bound along t axis: l2
  int l2 = -1;
  for (int k = k1; k <= k2; ++k) {
    for (int j = j1; j <= j2; ++j) {
      for (int i = i1; i <= i2; ++i) {
        for (int l = _attr._t - 1; l >= l1; --l) {
          if (IsForeground(i, j, k, l)) {
            if (l > l2) l2 = l;
            break;
          }
        }
      }
    }
  }
  // Do nothing if all voxels are background, but report it
  if (i1 > i2 || j1 > j2 || k1 > k2 || l1 > l2) return false;
  // Convert upper index bounds to margin widths
  i2 = (_attr._x - 1) - i2;
  j2 = (_attr._y - 1) - j2;
  k2 = (_attr._z - 1) - k2;
  l2 = (_attr._t - 1) - l2;
  // Leave a margin of background voxels with specified width
  // Note: Negative value gives the number of voxels to add.
  if (_attr._x > 1) i1 -= margin, i2 -= margin;
  if (_attr._y > 1) j1 -= margin, j2 -= margin;
  if (_attr._z > 1) k1 -= margin, k2 -= margin;
  if (_attr._t > 1) l1 -= margin, l2 -= margin;
  // Do nothing, if nothing to be done
  if (i1 == 0 && i2 == 0 && j1 == 0 && j2 == 0 &&
      k1 == 0 && k2 == 0 && l1 == 0 && l2 == 0) return true;
  // Adjust image lattice
  attr._x -= i1 + i2;
  attr._y -= j1 + j2;
  attr._z -= k1 + k2;
  attr._t -= l1 + l2;
  attr._xorigin = 0.5 * ((_attr._x - 1) + (i1 - i2));
  attr._yorigin = 0.5 * ((_attr._y - 1) + (j1 - j2));
  attr._zorigin = 0.5 * ((_attr._z - 1) + (k1 - k2));
  _attr.LatticeToWorld(attr._xorigin, attr._yorigin, attr._zorigin);
  attr._torigin = _attr.LatticeToTime(l1);
  // Convert upper margin widths to index bounds
  i2 = (_attr._x - 1) - i2;
  j2 = (_attr._y - 1) - j2;
  k2 = (_attr._z - 1) - k2;
  l2 = (_attr._t - 1) - l2;
  // Copy remaining voxels and pad lattice where needed
  const int nvoxels = attr.NumberOfLatticePoints();
  VoxelType *data       = Allocate<VoxelType>(nvoxels);
  VoxelType *data_iter  = data;
  for (int l = l1; l <= l2; ++l) {
    for (int k = k1; k <= k2; ++k) {
      for (int j = j1; j <= j2; ++j) {
        for (int i = i1; i <= i2; ++i, ++data_iter) {
          if (0 <= i && i < _attr._x &&
              0 <= j && j < _attr._y &&
              0 <= k && k < _attr._z &&
              0 <= l && l < _attr._t) {
            (*data_iter) = Get(i, j, k, l);
          } else {
            // Padded voxel to extend margin
            (*data_iter) = voxel_cast<VoxelType>(_bg);
          }
        }
      }
    }
  }
  // Initialize new image lattice
  this->Initialize(attr);
  data_iter = data;
  for (int l = 0; l < _attr._t; ++l) {
    for (int k = 0; k < _attr._z; ++k) {
      for (int j = 0; j < _attr._y; ++j) {
        for (int i = 0; i < _attr._x; ++i, ++data_iter) {
          Put(i, j, k, l, (*data_iter));
        }
      }
    }
  }
  // Free temporary allocated memory
  Deallocate(data);
  return true;
}

// =============================================================================
// VTK interface
// =============================================================================
#ifdef HAS_VTK

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>::ImageToVTK(vtkStructuredPoints *vtk) const
{
  if (this->ImageToVTKScalarType() == VTK_VOID) {
    cerr << "irtkGenericImage::ImageToVTK: Cannot convert image to VTK structured points" << endl;
    exit(1);
  }
  double x = 0, y = 0, z = 0;
  this->ImageToWorld(x, y, z);
  vtk->SetOrigin    (x, y, z);
  vtk->SetDimensions(_attr._x,  _attr._y,  _attr._z);
  vtk->SetSpacing   (_attr._dx, _attr._dy, _attr._dz);
#if VTK_MAJOR_VERSION >= 6
  vtk->AllocateScalars(this->ImageToVTKScalarType(), 1);
#else
  vtk->SetScalarType(this->ImageToVTKScalarType());
  vtk->AllocateScalars();
#endif
  const int        nvox = _attr._x * _attr._y * _attr._z;
  const VoxelType *ptr1 = this->Data();
  VoxelType       *ptr2 = reinterpret_cast<VoxelType *>(vtk->GetScalarPointer());
  for (int i = 0; i < nvox; ++i, ++ptr1) {
    for (int l = 0; l < _attr._t; ++l, ++ptr2) *ptr2 = ptr1[l * nvox];
  }
}
template <>
void irtkGenericImage<irtkMatrix3x3>::ImageToVTK(vtkStructuredPoints *) const
{
  cerr << "irtkGenericImage<irtkMatrix3x3>::VTKToImage: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class Type>
void irtkGenericImage<Type>::VTKToImage(vtkStructuredPoints *)
{
  cerr << this->NameOfClass() << "::VTKToImage: Not implemented" << endl;
  exit(1);
}

#endif
// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>::Read(const char *fname)
{
  // Read image
  irtkFileToImage *reader = irtkFileToImage::New(fname);
  irtkBaseImage   *image  = reader->GetOutput();
  // Convert image
  switch (image->GetDataType()) {
    case IRTK_VOXEL_CHAR:           { *this = *(dynamic_cast<irtkGenericImage<char>           *>(image)); } break;
    case IRTK_VOXEL_UNSIGNED_CHAR:  { *this = *(dynamic_cast<irtkGenericImage<unsigned char>  *>(image)); } break;
    case IRTK_VOXEL_SHORT:          { *this = *(dynamic_cast<irtkGenericImage<short>          *>(image)); } break;
    case IRTK_VOXEL_UNSIGNED_SHORT: { *this = *(dynamic_cast<irtkGenericImage<unsigned short> *>(image)); } break;
    case IRTK_VOXEL_INT:            { *this = *(dynamic_cast<irtkGenericImage<int>            *>(image)); } break;
    case IRTK_VOXEL_FLOAT:          { *this = *(dynamic_cast<irtkGenericImage<float>          *>(image)); } break;
    case IRTK_VOXEL_DOUBLE:         { *this = *(dynamic_cast<irtkGenericImage<double>         *>(image)); } break;
    default:
      cerr << "irtkFileToImage::GetOutput: Unknown data type: " << image->GetDataType() << endl;
      exit(1);
  }
  // Apply rescaling function
  if (reader->GetSlope() != .0) {
    switch (this->GetScalarType()) {
      case IRTK_VOXEL_FLOAT: {
        *this *= static_cast<float>(reader->GetSlope());
        *this += static_cast<float>(reader->GetIntercept());
        break;
      }
      case IRTK_VOXEL_DOUBLE: {
        *this *= static_cast<double>(reader->GetSlope());
        *this += static_cast<double>(reader->GetIntercept());
        break;
      }
      default: {
        if ((reader->GetSlope() != 1.0) || (reader->GetIntercept() != .0)) {
          cerr << this->NameOfClass() << "::Read: WARNING: Ignoring slope and intercept, use real data type instead" << endl;
        }
      }
    }
  }
  // Clean up
  delete reader;
  delete image;
}

template <> void irtkGenericImage<float3x3 >::Read(const char *)
{
  cerr << "irtkGenericImage<float3x3>::Read: Not implemented" << endl;
  exit(1);
}

template <> void irtkGenericImage<double3x3>::Read(const char *)
{
  cerr << "irtkGenericImage<double3x3>::Read: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkGenericImage<VoxelType>::Write(const char *fname) const
{
  irtkImageToFile *writer = irtkImageToFile::New(fname);
  writer->SetInput(const_cast<irtkGenericImage<VoxelType> *>(this));
  writer->Run();
  delete writer;
}

template <> void irtkGenericImage<float3x3 >::Write(const char *) const
{
  cerr << "irtkGenericImage<float3x3>::Write: Not implemented" << endl;
  exit(1);
}

template <> void irtkGenericImage<double3x3>::Write(const char *) const
{
  cerr << "irtkGenericImage<double3x3>::Read: Not implemented" << endl;
  exit(1);
}

// =============================================================================
// Explicit template instantiations
// =============================================================================

template class irtkGenericImage<char>;
template class irtkGenericImage<unsigned char>;
template class irtkGenericImage<short>;
template class irtkGenericImage<unsigned short>;
template class irtkGenericImage<int>;
template class irtkGenericImage<unsigned int>;
template class irtkGenericImage<float>;
template class irtkGenericImage<float2>;
template class irtkGenericImage<float3>;
template class irtkGenericImage<float4>;
template class irtkGenericImage<double>;
template class irtkGenericImage<double2>;
template class irtkGenericImage<double3>;
template class irtkGenericImage<double4>;
template class irtkGenericImage<float3x3>;
template class irtkGenericImage<double3x3>;

template irtkGenericImage<char>::irtkGenericImage(const irtkGenericImage<unsigned char> &);
template irtkGenericImage<char>::irtkGenericImage(const irtkGenericImage<short> &);
template irtkGenericImage<char>::irtkGenericImage(const irtkGenericImage<unsigned short> &);
template irtkGenericImage<char>::irtkGenericImage(const irtkGenericImage<int> &);
template irtkGenericImage<char>::irtkGenericImage(const irtkGenericImage<float> &);
template irtkGenericImage<char>::irtkGenericImage(const irtkGenericImage<double> &);

template irtkGenericImage<unsigned char>::irtkGenericImage(const irtkGenericImage<char> &);
template irtkGenericImage<unsigned char>::irtkGenericImage(const irtkGenericImage<short> &);
template irtkGenericImage<unsigned char>::irtkGenericImage(const irtkGenericImage<unsigned short> &);
template irtkGenericImage<unsigned char>::irtkGenericImage(const irtkGenericImage<int> &);
template irtkGenericImage<unsigned char>::irtkGenericImage(const irtkGenericImage<float> &);
template irtkGenericImage<unsigned char>::irtkGenericImage(const irtkGenericImage<double> &);

template irtkGenericImage<short>::irtkGenericImage(const irtkGenericImage<char> &);
template irtkGenericImage<short>::irtkGenericImage(const irtkGenericImage<unsigned char> &);
template irtkGenericImage<short>::irtkGenericImage(const irtkGenericImage<unsigned short> &);
template irtkGenericImage<short>::irtkGenericImage(const irtkGenericImage<int> &);
template irtkGenericImage<short>::irtkGenericImage(const irtkGenericImage<float> &);
template irtkGenericImage<short>::irtkGenericImage(const irtkGenericImage<double> &);

template irtkGenericImage<unsigned short>::irtkGenericImage(const irtkGenericImage<char> &);
template irtkGenericImage<unsigned short>::irtkGenericImage(const irtkGenericImage<unsigned char> &);
template irtkGenericImage<unsigned short>::irtkGenericImage(const irtkGenericImage<short> &);
template irtkGenericImage<unsigned short>::irtkGenericImage(const irtkGenericImage<int> &);
template irtkGenericImage<unsigned short>::irtkGenericImage(const irtkGenericImage<float> &);
template irtkGenericImage<unsigned short>::irtkGenericImage(const irtkGenericImage<double> &);

template irtkGenericImage<int>::irtkGenericImage(const irtkGenericImage<char> &);
template irtkGenericImage<int>::irtkGenericImage(const irtkGenericImage<unsigned char> &);
template irtkGenericImage<int>::irtkGenericImage(const irtkGenericImage<short> &);
template irtkGenericImage<int>::irtkGenericImage(const irtkGenericImage<unsigned short> &);
template irtkGenericImage<int>::irtkGenericImage(const irtkGenericImage<float> &);
template irtkGenericImage<int>::irtkGenericImage(const irtkGenericImage<double> &);

template irtkGenericImage<float>::irtkGenericImage(const irtkGenericImage<char> &);
template irtkGenericImage<float>::irtkGenericImage(const irtkGenericImage<unsigned char> &);
template irtkGenericImage<float>::irtkGenericImage(const irtkGenericImage<short> &);
template irtkGenericImage<float>::irtkGenericImage(const irtkGenericImage<unsigned short> &);
template irtkGenericImage<float>::irtkGenericImage(const irtkGenericImage<int> &);
template irtkGenericImage<float>::irtkGenericImage(const irtkGenericImage<double> &);

template irtkGenericImage<double>::irtkGenericImage(const irtkGenericImage<char> &);
template irtkGenericImage<double>::irtkGenericImage(const irtkGenericImage<unsigned char> &);
template irtkGenericImage<double>::irtkGenericImage(const irtkGenericImage<short> &);
template irtkGenericImage<double>::irtkGenericImage(const irtkGenericImage<unsigned short> &);
template irtkGenericImage<double>::irtkGenericImage(const irtkGenericImage<int> &);
template irtkGenericImage<double>::irtkGenericImage(const irtkGenericImage<float> &);

template irtkGenericImage<float2  >::irtkGenericImage(const irtkGenericImage<double2  > &);
template irtkGenericImage<float3  >::irtkGenericImage(const irtkGenericImage<double3  > &);
template irtkGenericImage<float4  >::irtkGenericImage(const irtkGenericImage<double4  > &);
template irtkGenericImage<float3x3>::irtkGenericImage(const irtkGenericImage<double3x3> &);

template irtkGenericImage<double2  >::irtkGenericImage(const irtkGenericImage<float2  > &);
template irtkGenericImage<double3  >::irtkGenericImage(const irtkGenericImage<float3  > &);
template irtkGenericImage<double4  >::irtkGenericImage(const irtkGenericImage<float4  > &);
template irtkGenericImage<double3x3>::irtkGenericImage(const irtkGenericImage<float3x3> &);

// TODO: Remove deprecated template instantiations below
template class irtkGenericImage<irtkVector3D<float>  >;
template class irtkGenericImage<irtkVector3D<double> >;
template irtkGenericImage<irtkVector3D<double> >::irtkGenericImage(const irtkGenericImage<irtkVector3D<float> > &);
