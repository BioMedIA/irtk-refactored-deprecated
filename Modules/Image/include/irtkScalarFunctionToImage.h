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

#ifndef _IRTKSCALARFUNCTIONTOIMAGE_H

#define _IRTKSCALARFUNCTIONTOIMAGE_H

/**
 * Class for scalar function to image filter.
 *
 * This class uses a scalar function to produce an image as output. The
 * filter loops through each voxel of the output image and calculates its
 * intensity as the value of the scalar function as a function of spatial
 * location.
 */

template <class VoxelType> class irtkScalarFunctionToImage : public irtkObject
{
  irtkObjectMacro(irtkScalarFunctionToImage);

private:

  /// Debugging flag
  bool _DebugFlag;

  /// Flag to use world or image coordinates for scalar function evaluation
  bool _UseWorldCoordinates;

protected:

  /// Input for the filter
  irtkScalarFunction *_input;

  /// Output for the filter
  irtkGenericImage<VoxelType> *_output;

public:

  /// Constructor (using world coordinates by default)
  irtkScalarFunctionToImage(bool = true);

  /// Deconstuctor
  virtual ~irtkScalarFunctionToImage();

  /// Set input
  virtual void SetInput (irtkScalarFunction *);

  /// Set output
  virtual void SetOutput(irtkGenericImage<VoxelType> *);

  /// Run the filter on entire image
  virtual void   Run();

  /// Set debugging flag
  SetMacro(DebugFlag, bool);

  /// Get debugging flag
  GetMacro(DebugFlag, bool);

  /// Print debugging messages if debugging is enabled
  virtual void Debug(char *);

};

#endif
