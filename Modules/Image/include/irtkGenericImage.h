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

#ifndef _IRTKGENERICIMAGE_H
#define _IRTKGENERICIMAGE_H

#include <irtkBaseImage.h>
#include <irtkVoxelCast.h>


/**
 * Generic class for 2D or 3D images
 *
 * This class implements generic 2D and 3D images. It provides functions
 * for accessing, reading, writing and manipulating images. This class can
 * be used for images with arbitrary voxel types using templates.
 */

template <class TVoxel>
class irtkGenericImage : public irtkBaseImage
{
  irtkObjectMacro(irtkGenericImage);

  // ---------------------------------------------------------------------------
  // Types

public:

  /// Voxel type
  typedef TVoxel VoxelType;

  /// Floating point type corresponding to voxel type
  /// \note The VoxelType as well as the RealType may be a matrix/vector type!
  typedef typename voxel_info<VoxelType>::RealType RealType;

  /// Scalar type corresponding to voxel type
  typedef typename voxel_info<VoxelType>::ScalarType ScalarType;

  // ---------------------------------------------------------------------------
  // Data members

protected:

  /// Pointer array for access to image data
  ///
  /// \note The image data is stored in a contiguous memory block which can
  ///       be alternatively referred to as 1D data array using \c _data.
  VoxelType ****_matrix;

  /// Pointer to image data
  VoxelType *_data;

  /// Whether image data memory itself is owned by this instance
  bool _dataOwner;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Allocate image memory
  void AllocateImage(VoxelType * = NULL);

public:

  /// Default constructor
  irtkGenericImage();

  /// Constructor from image file
  irtkGenericImage(const char *);

  /// Constructor for given image size
  explicit irtkGenericImage(int, int, int = 1, int = 1, VoxelType *data = NULL);

  /// Constructor for given image size
  explicit irtkGenericImage(int, int, int, int, int, VoxelType *data = NULL);

  /// Constructor for given image attributes
  explicit irtkGenericImage(const irtkImageAttributes &, VoxelType *data = NULL);

  /// Constructor for given image attributes
  explicit irtkGenericImage(const irtkImageAttributes &, int, VoxelType *data = NULL);

  /// Copy constructor for image
  explicit irtkGenericImage(const irtkBaseImage &);

  /// Copy constructor for image
  irtkGenericImage(const irtkGenericImage &);

  /// Copy constructor for image of different type
  template <class TVoxel2>
  irtkGenericImage(const irtkGenericImage<TVoxel2> &);

  /// Destructor
  ~irtkGenericImage();

  // ---------------------------------------------------------------------------
  // Initialization

  /// Create copy of this image
  virtual irtkBaseImage *Copy() const;

  /// Initialize a previously allocated image
  virtual void Initialize();

  /// Initialize an image
  virtual void Initialize(const irtkImageAttributes &, int, VoxelType *data);

  /// Initialize an image
  void Initialize(const irtkImageAttributes &, int);

  /// Initialize an image
  void Initialize(const irtkImageAttributes &, VoxelType *data = NULL);

  /// Initialize an image
  void Initialize(int, int, int, int, int, VoxelType *data = NULL);

  /// Initialize an image
  void Initialize(int, int, int = 1, int = 1, VoxelType *data = NULL);

  /// Copy image data from 1D array
  void CopyFrom(const VoxelType *);

  /// Copy image data from other image of same size
  void CopyFrom(const irtkBaseImage &);

  /// Copy image data from other image of same size
  void CopyFrom(const irtkGenericImage &);

  /// Assign constant value to each voxel
  irtkGenericImage& operator= (VoxelType);

  /// Assignment operator with implicit cast to double and then VoxelType
  irtkGenericImage<VoxelType>& operator= (const irtkBaseImage &);

  /// Assignment operator
  irtkGenericImage<VoxelType>& operator= (const irtkGenericImage &);

  /// Assignment operator with implicit cast
  template <class TVoxel2>
  irtkGenericImage<VoxelType>& operator= (const irtkGenericImage<TVoxel2> &);

  /// Clear an image
  void Clear();

  // ---------------------------------------------------------------------------
  // Lattice

  /// Number of vector components per voxel
  int N() const;

  /// Function to convert pixel to index
  /// more efficient than overwritten base class implementation
  int VoxelToIndex(int, int, int = 0, int = 0) const;

  // ---------------------------------------------------------------------------
  // Image data access

  /// Function for pixel get access
  VoxelType Get(int) const;

  /// Function for pixel get access
  VoxelType Get(int, int, int = 0, int = 0) const;

  /// Function for pixel put access
  void Put(int, VoxelType);

  /// Function for pixel put access
  void Put(int, int, VoxelType);

  /// Function for pixel put access
  void Put(int, int, int, VoxelType);

  /// Function for pixel put access
  void Put(int, int, int, int, VoxelType);

  /// Function for pixel access from via operators
  VoxelType& operator()(int);

  /// Function for pixel access from via operators
  const VoxelType& operator()(int) const;

  /// Function for pixel access from via operators
  VoxelType& operator()(int, int, int = 0, int = 0);

  /// Function for pixel access from via operators
  const VoxelType& operator()(int, int, int = 0, int = 0) const;

  // ---------------------------------------------------------------------------
  // Type independent access to scalar image data

  /// Function for pixel get access as double
  virtual double GetAsDouble(int) const;

  /// Function for pixel get access as double
  virtual double GetAsDouble(int, int, int = 0, int = 0) const;

  /// Function for pixel put access
  virtual void PutAsDouble(int, double);
  
  /// Function for pixel put access
  virtual void PutAsDouble(int, int, double);
  
  /// Function for pixel put access
  virtual void PutAsDouble(int, int, int, double);
  
  /// Function for pixel put access
  virtual void PutAsDouble(int, int, int, int, double);

  /// Function for pixel get access as double
  virtual void GetAsVector(irtkVector &, int) const;

  /// Function for pixel get access as double
  virtual void GetAsVector(irtkVector &, int, int, int = 0, int = 0) const;

  /// Function for pixel get access as double
  virtual irtkVector GetAsVector(int) const;

  /// Function for pixel get access as double
  virtual irtkVector GetAsVector(int, int, int = 0, int = 0) const;

  /// Function for pixel put access
  virtual void PutAsVector(int, const irtkVector &);

  /// Function for pixel put access
  virtual void PutAsVector(int, int, const irtkVector &);

  /// Function for pixel put access
  virtual void PutAsVector(int, int, int, const irtkVector &);

  /// Function for pixel put access
  virtual void PutAsVector(int, int, int, int, const irtkVector &);

  // ---------------------------------------------------------------------------
  // Access to raw image data
  
  /// Get raw pointer to contiguous image data
  VoxelType *Data(int = 0);

  /// Get raw pointer to contiguous image data
  VoxelType *Data(int, int, int = 0, int = 0);

  /// Get raw pointer to contiguous image data
  const VoxelType *Data(int = 0) const;

  /// Get raw pointer to contiguous image data
  const VoxelType *Data(int, int, int = 0, int = 0) const;

  /// Get raw pointer to contiguous image data
  virtual void *GetDataPointer(int = 0);

  /// Get raw pointer to contiguous image data
  virtual const void *GetDataPointer(int = 0) const;

  /// Get raw pointer to contiguous image data
  virtual void *GetDataPointer(int, int, int = 0, int = 0);

  /// Get raw pointer to contiguous image data
  virtual const void *GetDataPointer(int, int, int = 0, int = 0) const;

  /// Get enumeration value corresponding to voxel type
  virtual int GetDataType() const;

  /// Get size of each voxel in bytes
  virtual int GetDataTypeSize() const;

  /// Minimum value a pixel can hold without overflowing
  virtual double GetDataTypeMin() const;
  
  /// Maximum value a pixel can hold without overflowing
  virtual double GetDataTypeMax() const;

  // ---------------------------------------------------------------------------
  // Region-of-interest extraction

  /// Get image consisting of specified 2D slice
  irtkGenericImage GetRegion(int, int) const;

  /// Get image consisting of specified 2D slice
  void GetRegion(irtkGenericImage &, int, int) const;

  /// Get image consisting of specified 2D slice
  virtual void GetRegion(irtkBaseImage *&, int, int) const;

  /// Get image consisting of specified 3D subregion
  irtkGenericImage GetRegion(int, int, int,
                             int, int, int) const;

  /// Get image consisting of specified 3D subregion
  void GetRegion(irtkGenericImage &, int, int, int,
                                     int, int, int) const;

  /// Get image consisting of specified 3D subregion
  virtual void GetRegion(irtkBaseImage *&, int, int, int,
                                           int, int, int) const;

  /// Get image consisting of specified 4D subregion
  irtkGenericImage GetRegion(int, int, int, int,
                             int, int, int, int) const;

  /// Get image consisting of specified 4D subregion
  void GetRegion(irtkGenericImage &, int, int, int, int,
                                     int, int, int, int) const;

  /// Get image consisting of specified 4D subregion
  virtual void GetRegion(irtkBaseImage *&, int, int, int, int,
                                           int, int, int, int) const;

  /// Get time instance (i.e., frame) or channel of image
  irtkGenericImage GetFrame(int, int = -1) const;

  /// Get time instance (i.e., frame) or channel of image
  void GetFrame(irtkGenericImage &, int, int = -1) const;

  /// Get time instance (i.e., frame) or channel of image
  virtual void GetFrame(irtkBaseImage *&, int, int = -1) const;

  // ---------------------------------------------------------------------------
  // Image arithmetic

  /// Equality operator
  /// \note Use explicit negation for inequality comparison.
  ///       The overloaded != operator is used for binarization of the image.
  template <class TVoxel2>
  bool operator== (const irtkGenericImage<TVoxel2> &) const;

  irtkGenericImage& operator+=(const irtkGenericImage &); ///< Add image
  irtkGenericImage& operator-=(const irtkGenericImage &); ///< Subtract image
  irtkGenericImage& operator*=(const irtkGenericImage &); ///< Multipy voxels
  irtkGenericImage& operator/=(const irtkGenericImage &); ///< Divide voxels

  irtkGenericImage& operator+=(ScalarType); ///< Add scalar
  irtkGenericImage& operator-=(ScalarType); ///< Subtract scalar
  irtkGenericImage& operator*=(ScalarType); ///< Multiply by scalar
  irtkGenericImage& operator/=(ScalarType); ///< Divide by scalar

  irtkGenericImage  operator+ (const irtkGenericImage &) const; ///< Add images
  irtkGenericImage  operator- (const irtkGenericImage &) const; ///< Subtract images
  irtkGenericImage  operator* (const irtkGenericImage &) const; ///< Multiply images voxel-wise
  irtkGenericImage  operator/ (const irtkGenericImage &) const; ///< Divide images voxel-wise

  irtkGenericImage  operator+ (ScalarType) const; ///< Add scalar to image
  irtkGenericImage  operator- (ScalarType) const; ///< Subtract scalar from image
  irtkGenericImage  operator* (ScalarType) const; ///< Multiply image by scalar
  irtkGenericImage  operator/ (ScalarType) const; ///< Divide image by scalar

  // ---------------------------------------------------------------------------
  // Thresholding

  // Import other overload
  using irtkBaseImage::PutBackgroundValueAsDouble;

  /// Put background value
  virtual void PutBackgroundValueAsDouble(double, bool);

  irtkGenericImage& operator>=(VoxelType);       ///< Clamp image given upper threshold
  irtkGenericImage& operator<=(VoxelType);       ///< Clamp image given lower threshold

  irtkGenericImage  operator> (VoxelType) const; ///< Clamp image given upper threshold
  irtkGenericImage  operator< (VoxelType) const; ///< Clamp image given lower threshold

  /// Get binary mask for voxels which are not equal the scalar
  irtkBinaryImage operator!=(VoxelType) const;

  // ---------------------------------------------------------------------------
  // Common image statistics

  /// Minimum and maximum pixel values get accessor
  void GetMinMax(VoxelType &, VoxelType &) const;

  /// Minimum and maximum pixel values get accessor with padding
  void GetMinMax(VoxelType &, VoxelType &, VoxelType) const;

  /// Linearly rescale intensities
  void PutMinMax(VoxelType, VoxelType);

  /// Average pixel values get accessor
  RealType GetAverage(int = 1) const;

  /// Standard Deviation of the pixels
  RealType GetSD(int = 1) const;

  /// Get Max Intensity position around the point
  void GetMaxPosition(irtkPoint &, int = 1, int = 0) const;
  
  /// Get Gravity center position of a given window
  void GravityCenter(irtkPoint &, int = 1, int = 0) const;

  // ---------------------------------------------------------------------------
  // Common image manipulations

  virtual void ReflectX();  ///< Reflect image along x
  virtual void ReflectY();  ///< Reflect image along y
  virtual void ReflectZ();  ///< Reflect image along z

  virtual void FlipXY(int); ///< Flip x and y axis
  virtual void FlipXZ(int); ///< Flip x and z axis
  virtual void FlipYZ(int); ///< Flip y and z axis
  virtual void FlipXT(int); ///< Flip x and t axis
  virtual void FlipYT(int); ///< Flip y and t axis
  virtual void FlipZT(int); ///< Flip z and t axis

  bool CropPad(int margin = 0); ///< Crop/pad image background

  // ---------------------------------------------------------------------------
  // VTK interface
#ifdef HAS_VTK

  /// Convert image to VTK structured points
  virtual void ImageToVTK(vtkStructuredPoints *) const;

  /// Convert VTK structured points to image
  virtual void VTKToImage(vtkStructuredPoints *);

#endif

  // ---------------------------------------------------------------------------
  // I/O

  /// Read image from file
  virtual void Read(const char *);
  
  /// Write image to file
  virtual void Write(const char *) const;

  // ---------------------------------------------------------------------------
  // Deprecated

  /// Minimum and maximum pixel values get accessor
  /// \deprecated Use respective overloaded method of GetMinMax instead.
  void GetMinMax(VoxelType *, VoxelType *) const;
  
  /// Minimum and maximum pixel values get accessor with padding
  /// \deprecated Use respective overloaded method of GetMinMax instead.
  void GetMinMax(VoxelType *, VoxelType *, VoxelType) const;

  /// Minimum and maximum pixel values get accessor with padding
  /// \deprecated Use respective overloaded method of GetMinMax instead.
  void GetMinMaxPad(VoxelType *, VoxelType *, VoxelType) const;

  /// \returns Raw pointer to contiguous image data.
  /// \deprecated Use Data instead.
  VoxelType *GetPointerToVoxels(int = 0, int = 0, int = 0, int = 0);
  
  /// \returns Raw pointer to contiguous image data.
  /// \deprecated Use Data instead.
  const VoxelType *GetPointerToVoxels(int = 0, int = 0, int = 0, int = 0) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType> template <class VoxelType2>
irtkGenericImage<VoxelType>& irtkGenericImage<VoxelType>::operator=(const irtkGenericImage<VoxelType2> &image)
{
  this->Initialize(image.GetImageAttributes());
  VoxelType        *ptr1 = this->GetPointerToVoxels();
  const VoxelType2 *ptr2 = image.GetPointerToVoxels();
  for (int idx = 0; idx < _NumberOfVoxels; idx++) {
    ptr1[idx] = voxel_cast<VoxelType>(ptr2[idx]);
  }
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
  return *this;
}

// =============================================================================
// Lattice
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
inline int irtkGenericImage<VoxelType>::N() const
{
  return voxel_info<VoxelType>::vector_size();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline int irtkGenericImage<VoxelType>::VoxelToIndex(int x, int y, int z, int t) const
{
  return (&_matrix[t][z][y][x] - &_matrix[0][0][0][0]);
}

// =============================================================================
// Image data access
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void irtkGenericImage<VoxelType>::Put(int index, VoxelType val)
{
  _data[index] = val;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void irtkGenericImage<VoxelType>::Put(int x, int y, VoxelType val)
{
  _matrix[0][0][y][x] = val;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void irtkGenericImage<VoxelType>::Put(int x, int y, int z, VoxelType val)
{
  _matrix[0][z][y][x] = val;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void irtkGenericImage<VoxelType>::Put(int x, int y, int z, int t, VoxelType val)
{
  _matrix[t][z][y][x] = val;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType irtkGenericImage<VoxelType>::Get(int index) const
{
  return _data[index];
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType irtkGenericImage<VoxelType>::Get(int x, int y, int z, int t) const
{
  return _matrix[t][z][y][x];
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType &irtkGenericImage<VoxelType>::operator ()(int index)
{
  return _data[index];
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType &irtkGenericImage<VoxelType>::operator ()(int index) const
{
  return _data[index];
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType& irtkGenericImage<VoxelType>::operator()(int x, int y, int z, int t)
{
  return _matrix[t][z][y][x];
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType& irtkGenericImage<VoxelType>::operator()(int x, int y, int z, int t) const
{
  return _matrix[t][z][y][x];
}

// =============================================================================
// Image arithmetics
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType> template <class VoxelType2>
bool irtkGenericImage<VoxelType>::operator==(const irtkGenericImage<VoxelType2> &image) const
{
  if (this->GetImageAttributes() != image.GetImageAttributes()) return false;
  const VoxelType  *ptr1 = this->GetPointerToVoxels();
  const VoxelType2 *ptr2 = image.GetPointerToVoxels();
  for (int idx = 0; idx < image; ++idx) {
    if (IsForeground(idx) && image.IsForeground(idx) && ptr1[idx] != ptr2[idx]) {
      return false;
    }
  }
  return true;
}

// =============================================================================
// Type independent access to scalar image data
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void irtkGenericImage<VoxelType>::PutAsDouble(int index, double val)
{
  _data[index] = voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void irtkGenericImage<VoxelType>::PutAsDouble(int x, int y, double val)
{
  _matrix[0][0][y][x] = voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void irtkGenericImage<VoxelType>::PutAsDouble(int x, int y, int z, double val)
{
  _matrix[0][z][y][x] = voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void irtkGenericImage<VoxelType>::PutAsDouble(int x, int y, int z, int t, double val)
{
  _matrix[t][z][y][x] = voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline double irtkGenericImage<VoxelType>::GetAsDouble(int index) const
{
  return voxel_cast<double>(_data[index]);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline double irtkGenericImage<VoxelType>::GetAsDouble(int x, int y, int z, int t) const
{
  return voxel_cast<double>(_matrix[t][z][y][x]);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void irtkGenericImage<VoxelType>::PutAsVector(int index, const irtkVector &value)
{
  _data[index] = voxel_cast<VoxelType>(value);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void irtkGenericImage<VoxelType>::PutAsVector(int x, int y, const irtkVector &value)
{
  _matrix[0][0][y][x] = voxel_cast<VoxelType>(value);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void irtkGenericImage<VoxelType>::PutAsVector(int x, int y, int z, const irtkVector &value)
{
  _matrix[0][z][y][x] = voxel_cast<VoxelType>(value);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void irtkGenericImage<VoxelType>::PutAsVector(int x, int y, int z, int t, const irtkVector &value)
{
  _matrix[t][z][y][x] = voxel_cast<VoxelType>(value);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void irtkGenericImage<VoxelType>::GetAsVector(irtkVector &value, int index) const
{
  value = voxel_cast<irtkVector>(_data[index]);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void irtkGenericImage<VoxelType>::GetAsVector(irtkVector &value, int x, int y, int z, int t) const
{
  value = voxel_cast<irtkVector>(_matrix[t][z][y][x]);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline irtkVector irtkGenericImage<VoxelType>::GetAsVector(int index) const
{
  return voxel_cast<irtkVector>(_data[index]);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline irtkVector irtkGenericImage<VoxelType>::GetAsVector(int x, int y, int z, int t) const
{
  return voxel_cast<irtkVector>(_matrix[t][z][y][x]);
}

// =============================================================================
// Access to raw image data
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType *irtkGenericImage<VoxelType>::Data(int i)
{
  return _data + i;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType *irtkGenericImage<VoxelType>::Data(int i) const
{
  return _data + i;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType *irtkGenericImage<VoxelType>::Data(int x, int y, int z, int t)
{
  return &_matrix[t][z][y][x];
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType *irtkGenericImage<VoxelType>::Data(int x, int y, int z, int t) const
{
  return &_matrix[t][z][y][x];
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void *irtkGenericImage<VoxelType>::GetDataPointer(int i)
{
  return _data + i;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline const void *irtkGenericImage<VoxelType>::GetDataPointer(int i) const
{
  return _data + i;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void *irtkGenericImage<VoxelType>::GetDataPointer(int x, int y, int z, int t)
{
  return &_matrix[t][z][y][x];
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline const void *irtkGenericImage<VoxelType>::GetDataPointer(int x, int y, int z, int t) const
{
  return &_matrix[t][z][y][x];
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline int irtkGenericImage<VoxelType>::GetDataType() const
{
  return voxel_info<VoxelType>::type();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline int irtkGenericImage<VoxelType>::GetDataTypeSize() const
{
  return static_cast<int>(sizeof(VoxelType));
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline double irtkGenericImage<VoxelType>::GetDataTypeMin() const
{
  return voxel_limits<VoxelType>::min();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline double irtkGenericImage<VoxelType>::GetDataTypeMax() const
{
  return voxel_limits<VoxelType>::max();
}

// =============================================================================
// Deprecated
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void irtkGenericImage<VoxelType>::GetMinMax(VoxelType *min, VoxelType *max) const
{
  this->GetMinMax(*min, *max);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void irtkGenericImage<VoxelType>::GetMinMax(VoxelType *min, VoxelType *max, VoxelType pad) const
{
  this->GetMinMax(*min, *max, pad);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void irtkGenericImage<VoxelType>::GetMinMaxPad(VoxelType *min, VoxelType *max, VoxelType pad) const
{
  this->GetMinMax(*min, *max, pad);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType *irtkGenericImage<VoxelType>::GetPointerToVoxels(int x, int y, int z, int t)
{
  return &_matrix[t][z][y][x];
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType *irtkGenericImage<VoxelType>::GetPointerToVoxels(int x, int y, int z, int t) const
{
  return &_matrix[t][z][y][x];
}

////////////////////////////////////////////////////////////////////////////////
// Common specializations
////////////////////////////////////////////////////////////////////////////////

typedef irtkGenericImage<irtkBytePixel> irtkByteImage; ///< Unsigned char image
typedef irtkGenericImage<irtkGreyPixel> irtkGreyImage; ///< Short image
typedef irtkGenericImage<irtkRealPixel> irtkRealImage; ///< Float image


#endif

