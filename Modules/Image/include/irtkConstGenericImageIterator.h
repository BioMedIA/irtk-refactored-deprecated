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

#ifndef _irtkConstGenericImageIterator_H

#define _irtkConstGenericImageIterator_H


/**
 * Const image iterator
 */
template <class VoxelType>
class irtkConstGenericImageIterator : public irtkConstImageIterator
{
public:

  /// Constructor
  irtkConstGenericImageIterator(const irtkImageAttributes &, const VoxelType * = NULL);

  /// Constructor
  irtkConstGenericImageIterator(irtkGenericImage<VoxelType> &);

  /// Constructor
  irtkConstGenericImageIterator(irtkGenericImage<VoxelType> *);

  /// Copy constructor
  irtkConstGenericImageIterator(const irtkConstImageIterator &);

  /// Assignment operator
  irtkConstGenericImageIterator &operator =(const irtkConstGenericImageIterator &);

  /// Destructor
  ~irtkConstGenericImageIterator();

  /// Get pointer to current iterator position
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  const VoxelType *Current() const;

  /// Get pointer to current iterator position
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  ///
  /// \param[in] t Channel/Component/Frame offset relative to current iterator position.
  ///              For example, set iterator region to only the first channel/frame and
  ///              then access other channels/vector components/frames using this method.
  const VoxelType *Current(int) const;

  /// Get pointer to current iterator position and post-increment iterator
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  const VoxelType *Next();

  /// Get pointer to current iterator position and post-increment iterator
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  ///
  /// \param[in] t Channel/Component/Frame offset relative to current iterator position.
  ///              For example, set iterator region to only the first channel/frame and
  ///              then access other channels/vector components/frames using this method.
  const VoxelType *Next(int);

  /// Get reference to voxel value at current iterator position
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  const VoxelType &Value() const;

  /// Get reference to voxel value at current iterator position
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  ///
  /// \param[in] t Channel/Component/Frame offset relative to current iterator position.
  ///              For example, set iterator region to only the first channel/frame and
  ///              then access other channels/vector components/frames using this method.
  const VoxelType &Value(int t) const;

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
inline irtkConstGenericImageIterator<VoxelType>::irtkConstGenericImageIterator(const irtkImageAttributes &attr, const VoxelType *data)
:
  irtkConstImageIterator(attr, data, voxel_info<VoxelType>::type())
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline irtkConstGenericImageIterator<VoxelType>::irtkConstGenericImageIterator(irtkGenericImage<VoxelType> &image)
:
  irtkConstImageIterator(image)
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline irtkConstGenericImageIterator<VoxelType>::irtkConstGenericImageIterator(irtkGenericImage<VoxelType> *image)
:
  irtkConstImageIterator(image)
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline irtkConstGenericImageIterator<VoxelType>::irtkConstGenericImageIterator(const irtkConstImageIterator &other)
:
  irtkConstImageIterator(other)
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline irtkConstGenericImageIterator<VoxelType> &irtkConstGenericImageIterator<VoxelType>
::operator =(const irtkConstGenericImageIterator &rhs)
{
  irtkConstImageIterator::operator =(rhs);
  return *this;
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline irtkConstGenericImageIterator<VoxelType>::~irtkConstGenericImageIterator()
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType *irtkConstGenericImageIterator<VoxelType>::Current() const
{
  return irtkConstImageIterator::Current<VoxelType>();
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType *irtkConstGenericImageIterator<VoxelType>::Current(int t) const
{
  return irtkConstImageIterator::Current<VoxelType>(t);
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType *irtkConstGenericImageIterator<VoxelType>::Next()
{
  return irtkConstImageIterator::Next<VoxelType>();
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType *irtkConstGenericImageIterator<VoxelType>::Next(int t)
{
  return irtkConstImageIterator::Next<VoxelType>(t);
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType &irtkConstGenericImageIterator<VoxelType>::Value() const
{
  return irtkConstImageIterator::Value<VoxelType>();
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType &irtkConstGenericImageIterator<VoxelType>::Value(int t) const
{
  return irtkConstImageIterator::Value<VoxelType>(t);
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline double irtkConstGenericImageIterator<VoxelType>::ValueAsDouble() const
{
  return static_cast<double>(Value());
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline double irtkConstGenericImageIterator<VoxelType>::ValueAsDouble(int t) const
{
  return static_cast<double>(Value(t));
}


#endif
