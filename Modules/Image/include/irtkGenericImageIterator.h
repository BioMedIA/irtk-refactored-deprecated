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

#ifndef _irtkGenericImageIterator_H

#define _irtkGenericImageIterator_H

#include <irtkImageIterator.h>


/**
 * Non-const image iterator
 */
template <class VoxelType>
class irtkGenericImageIterator : public irtkImageIterator
{
public:

  /// Constructor
  irtkGenericImageIterator(const irtkImageAttributes &, VoxelType * = NULL);

  /// Constructor
  irtkGenericImageIterator(irtkGenericImage<VoxelType> &);

  /// Constructor
  irtkGenericImageIterator(irtkGenericImage<VoxelType> *);

  /// Copy constructor
  irtkGenericImageIterator(const irtkConstImageIterator &);

  /// Assignment operator
  irtkGenericImageIterator &operator =(const irtkGenericImageIterator &);

  /// Destructor
  ~irtkGenericImageIterator();

  /// Get pointer to current iterator position
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  VoxelType *Current() const;

  /// Get pointer to current iterator position
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  ///
  /// \param[in] t Channel/Component/Frame offset relative to current iterator position.
  ///              For example, set iterator region to only the first channel/frame and
  ///              then access other channels/vector components/frames using this method.
  VoxelType *Current(int) const;

  /// Get pointer to current iterator position and post-increment iterator
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  VoxelType *Next();

  /// Get pointer to current iterator position and post-increment iterator
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  ///
  /// \param[in] t Channel/Component/Frame offset relative to current iterator position.
  ///              For example, set iterator region to only the first channel/frame and
  ///              then access other channels/vector components/frames using this method.
  VoxelType *Next(int);

  /// Get reference to voxel value at current iterator position
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  VoxelType &Value() const;

  /// Get reference to voxel value at current iterator position
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  ///
  /// \param[in] t Channel/Component/Frame offset relative to current iterator position.
  ///              For example, set iterator region to only the first channel/frame and
  ///              then access other channels/vector components/frames using this method.
  VoxelType &Value(int t) const;

  /// Get current voxel value casted to double
  double ValueAsDouble() const;

  /// Get current voxel value casted to double
  ///
  /// \param[in] t Channel/Component/Frame offset relative to current iterator position.
  ///              For example, set iterator region to only the first channel/frame and
  ///              then access other channels/vector components/frames using this method.
  double ValueAsDouble(int) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// ---------------------------------------------------------------------------
template <class VoxelType>
inline irtkGenericImageIterator<VoxelType>::irtkGenericImageIterator(const irtkImageAttributes &attr, VoxelType *data)
:
  irtkImageIterator(attr, data, voxel_info<VoxelType>::type())
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline irtkGenericImageIterator<VoxelType>::irtkGenericImageIterator(irtkGenericImage<VoxelType> &image)
:
  irtkImageIterator(image)
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline irtkGenericImageIterator<VoxelType>::irtkGenericImageIterator(irtkGenericImage<VoxelType> *image)
:
  irtkImageIterator(image)
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline irtkGenericImageIterator<VoxelType>::irtkGenericImageIterator(const irtkConstImageIterator &other)
:
  irtkImageIterator(other)
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline irtkGenericImageIterator<VoxelType> &irtkGenericImageIterator<VoxelType>::operator =(const irtkGenericImageIterator &rhs)
{
  irtkImageIterator::operator =(rhs);
  return *this;
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline irtkGenericImageIterator<VoxelType>::~irtkGenericImageIterator()
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType *irtkGenericImageIterator<VoxelType>::Current() const
{
  return irtkImageIterator::Current<VoxelType>();
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType *irtkGenericImageIterator<VoxelType>::Current(int t) const
{
  return irtkImageIterator::Current<VoxelType>(t);
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType *irtkGenericImageIterator<VoxelType>::Next()
{
  return irtkImageIterator::Next<VoxelType>();
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType *irtkGenericImageIterator<VoxelType>::Next(int t)
{
  return irtkImageIterator::Next<VoxelType>(t);
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType &irtkGenericImageIterator<VoxelType>::Value() const
{
  return irtkImageIterator::Value<VoxelType>();
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType &irtkGenericImageIterator<VoxelType>::Value(int t) const
{
  return irtkImageIterator::Value<VoxelType>(t);
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline double irtkGenericImageIterator<VoxelType>::ValueAsDouble() const
{
  return static_cast<double>(Value());
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline double irtkGenericImageIterator<VoxelType>::ValueAsDouble(int t) const
{
  return static_cast<double>(Value(t));
}


#endif
