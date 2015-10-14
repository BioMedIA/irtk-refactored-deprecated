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

#ifndef _IRTKPOINTTOIMAGE_H

#define _IRTKPOINTTOIMAGE_H

/**
 * Class for point to image filter.
 *
 * This class uses a point set to produce an image as output. The filter
 * loops through point of the input point set and sets the nearest voxel
 * in the output image to 1. All voxels of the output image are initially
 * set to 0.
 */

template <class VoxelType> class irtkPointToImage : public irtkObject
{
  irtkObjectMacro(irtkPointToImage);

protected:

  /// Input for the image to point filter
  irtkPointSet *_input;

  /// Output for the image to point filter
  irtkGenericImage<VoxelType> *_output;

  /// Debugging flag
  bool _DebugFlag;

  /// Flag to use world or image coordinates for point locations
  bool _UseWorldCoordinates;

public:

  /// Constructor (using world coordinates by default)
  irtkPointToImage(bool = true);

  // Deconstuctor
  virtual ~irtkPointToImage();

  /// Set input
  virtual void SetInput (irtkPointSet  *);

  /// Set output
  virtual void SetOutput(irtkGenericImage<VoxelType> *);

  /// Run image to point filter
  virtual void Run();

  /// Set debugging flag
  SetMacro(DebugFlag, bool);

  /// Get debugging flag
  GetMacro(DebugFlag, bool);

  /// Print debugging messages if debugging is enabled
  virtual void Debug(char *);
};

#include <irtkFill.h>

#endif
