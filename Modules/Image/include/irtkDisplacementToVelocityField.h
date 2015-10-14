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

#ifndef _IRTKDISPLACEMENTTOVELOCITYFIELD_H

#define _IRTKDISPLACEMENTTOVELOCITYFIELD_H


#include <irtkImageToImage.h>


/**
 * Base class of image filters which compute a stationary velocity field
 * given a displacement field.
 */

template <class VoxelType>
class irtkDisplacementToVelocityField : public irtkImageToImage<VoxelType>
{
  irtkImageFilterMacro(irtkDisplacementToVelocityField);

protected:

  /// Constructor
  irtkDisplacementToVelocityField();

public:  

  /// Destructor
  virtual ~irtkDisplacementToVelocityField();

};

///////////////////////////////////////////////////////////////////////////////
// Inline definitions
///////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
template <class VoxelType>
irtkDisplacementToVelocityField<VoxelType>::irtkDisplacementToVelocityField()
:
  irtkImageToImage<VoxelType>()
{
}

// ----------------------------------------------------------------------------
template <class VoxelType>
irtkDisplacementToVelocityField<VoxelType>::~irtkDisplacementToVelocityField()
{
}

///////////////////////////////////////////////////////////////////////////////
// Actual implementations
///////////////////////////////////////////////////////////////////////////////

#include <irtkDisplacementToVelocityFieldBCH.h>

///////////////////////////////////////////////////////////////////////////////
// log function for vector fields
///////////////////////////////////////////////////////////////////////////////

template <class VoxelType>
void log(irtkGenericImage<VoxelType> *d)
{
  irtkDisplacementToVelocityFieldBCH<VoxelType> dtov;
  dtov.SetInput (d);
  dtov.SetOutput(d);
  dtov.Run();
}


#endif
