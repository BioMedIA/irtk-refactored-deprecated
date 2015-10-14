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

#ifndef _irtkImageIterator_H

#define _irtkImageIterator_H

#include <irtkConstImageIterator.h>


/**
 * Base class of non-const image iterator
 */
class irtkImageIterator : public irtkConstImageIterator
{
public:

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Constructor
  irtkImageIterator(const irtkImageAttributes &, int);

  /// Constructor
  irtkImageIterator(const irtkImageAttributes &, void * = NULL, int = IRTK_VOXEL_UNKNOWN);

  /// Constructor
  irtkImageIterator(irtkBaseImage &);

  /// Constructor
  irtkImageIterator(irtkBaseImage *);

  /// Copy constructor
  irtkImageIterator(const irtkConstImageIterator &);

  /// Assignment operator
  irtkImageIterator &operator =(const irtkImageIterator &);

  /// Destructor
  virtual ~irtkImageIterator();

  // ---------------------------------------------------------------------------
  // Iteration

  /// Get pointer to current iterator position
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  template <class VoxelType>
  VoxelType *Current() const;

  /// Get pointer to current iterator position
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  ///
  /// \param[in] t Channel/Component/Frame offset relative to current iterator position.
  ///              For example, set iterator region to only the first channel/frame and
  ///              then access other channels/vector components/frames using this method.
  template <class VoxelType>
  VoxelType *Current(int) const;

  /// Get pointer to current iterator position and post-increment iterator
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  template <class VoxelType>
  VoxelType *Next();

  /// Get pointer to current iterator position and post-increment iterator
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  ///
  /// \param[in] t Channel/Component/Frame offset relative to current iterator position.
  ///              For example, set iterator region to only the first channel/frame and
  ///              then access other channels/vector components/frames using this method.
  template <class VoxelType>
  VoxelType *Next(int);

  /// Get reference to voxel value at current iterator position
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  template <class VoxelType>
  VoxelType &Value() const;

  /// Get reference to voxel value at current iterator position
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  ///
  /// \param[in] t Channel/Component/Frame offset relative to current iterator position.
  ///              For example, set iterator region to only the first channel/frame and
  ///              then access other channels/vector components/frames using this method.
  template <class VoxelType>
  VoxelType &Value(int t) const;

};

//////////////////////////////////////////////////////////////////////////////
// Inline definitions
//////////////////////////////////////////////////////////////////////////////

// ===========================================================================
// Construction/Destruction
// ===========================================================================

// ---------------------------------------------------------------------------
inline irtkImageIterator::irtkImageIterator(const irtkImageAttributes &attr, int type)
:
  irtkConstImageIterator(attr, type)
{
}

// ---------------------------------------------------------------------------
inline irtkImageIterator::irtkImageIterator(const irtkImageAttributes &attr, void *data, int type)
:
  irtkConstImageIterator(attr, data, type)
{
}

// ---------------------------------------------------------------------------
inline irtkImageIterator::irtkImageIterator(irtkBaseImage &image)
:
  irtkConstImageIterator(image)
{
}

// ---------------------------------------------------------------------------
inline irtkImageIterator::irtkImageIterator(irtkBaseImage *image)
:
  irtkConstImageIterator(image)
{
}

// ---------------------------------------------------------------------------
inline irtkImageIterator::irtkImageIterator(const irtkConstImageIterator &other)
:
  irtkConstImageIterator(other)
{
}

// ---------------------------------------------------------------------------
inline irtkImageIterator &irtkImageIterator::operator =(const irtkImageIterator &rhs)
{
  irtkConstImageIterator::operator =(rhs);
  return *this;
}

// ---------------------------------------------------------------------------
inline irtkImageIterator::~irtkImageIterator()
{
}

// ===========================================================================
// Iteration
// ===========================================================================

// ---------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType *irtkImageIterator::Current() const
{
  return reinterpret_cast<VoxelType *>(const_cast<char *>(_Next));
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType *irtkImageIterator::Current(int t) const
{
  return reinterpret_cast<VoxelType *>(const_cast<char *>(_Next) + t * _XYZ);
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType *irtkImageIterator::Next()
{
  VoxelType *current = reinterpret_cast<VoxelType *>(const_cast<char *>(_Next));
  this->operator ++();
  return current;
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType *irtkImageIterator::Next(int t)
{
  VoxelType *current = reinterpret_cast<VoxelType *>(const_cast<char *>(_Next) + t * _XYZ);
  this->operator ++();
  return current;
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType &irtkImageIterator::Value() const
{
  return *reinterpret_cast<VoxelType *>(const_cast<char *>(_Next));
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType &irtkImageIterator::Value(int t) const
{
  return *reinterpret_cast<VoxelType *>(const_cast<char *>(_Next) + t * _XYZ);
}


#endif
