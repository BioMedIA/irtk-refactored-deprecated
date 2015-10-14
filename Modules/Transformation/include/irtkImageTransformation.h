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

#ifndef _IRTKIMAGETRANSFORMATION_H

#define _IRTKIMAGETRANSFORMATION_H

#include <irtkImage.h>

#include <irtkImageFunction.h>

#include <irtkTransformation.h>


class irtkMultiThreadedImageHomogeneousTransformation;

/// Image used to cache displacements
///
/// \todo Use irtkVector3D<double> as voxel type to enable the caching of
///       a 4D displacement field which is generated by a 4D transformation
///       or a 3D SV FFD integrated over multiple time intervals.
class irtkImageTransformationCache : public irtkRealImage
{
  irtkObjectMacro(irtkImageTransformationCache);
 
  /// Whether the cache requires an update
  irtkPublicAttributeMacro(bool, Modified);

public:

  /// Constructor
  irtkImageTransformationCache() : _Modified(true) {}

  /// Destructor
  virtual ~irtkImageTransformationCache() {}
};


/**
 * Filter for image transformations.
 *
 * This class implements an image filter which takes an input image and a
 * transformation and computes the corresponding output image. The output
 * image is calculated by looping over the voxel locations and calculating
 * the corresponding voxel locations in the input image. The intensity of
 * the voxels of the output image is by interpolation from the input image.
 * Note, that the transformation is assumed to map the input image into the
 * output image and is therefore inverted during the execution of the filter.
 * All calculations are carried out using world coordinates rather than image
 * coordinates.
 *
 */

class irtkImageTransformation : public irtkObject
{
  irtkObjectMacro(irtkImageTransformation);
 
  friend class irtkMultiThreadedImageHomogeneousTransformation;

  /// Cached displacements
  irtkAggregateMacro(irtkImageTransformationCache, Cache);

  /// Whether the cache was allocated by this instance
  irtkAttributeMacro(bool, CacheOwner);

  /// Interpolation mode to use for interpolating cached displacements
  irtkPublicAttributeMacro(irtkInterpolationMode, CacheInterpolation);

  /// Extrapolation mode to use when interpolating cached displacements
  irtkPublicAttributeMacro(irtkExtrapolationMode, CacheExtrapolation);

  /// Interpolator for cached displacements
  irtkComponentMacro(irtkInterpolateImageFunction, DisplacementField);

  /// Number of points for which the required inverse transformation was invalid
  irtkReadOnlyAttributeMacro(int, NumberOfSingularPoints);

protected:

  /// Input for the image to image filter
  irtkImage *_input;

  /// Output for the image to image filter
  irtkImage *_output;

  /// Transformation
  irtkTransformation *_transformation;

  /// Interpolation
  irtkImageFunction *_interpolator;

  /// Padding value in target (voxels in the target image with this
  /// value will be ignored)
  double _TargetPaddingValue;

  /// Padding value in source (voxels outside the source image will
  /// be set to this value)
  double _SourcePaddingValue;

  /// Scale factor for intensities in transformed image
  double _ScaleFactor;

  /// Offset for intensities in transformed image
  double _Offset;

  /// Temporal offset of input image (in addition to input _torigin)
  double _InputTimeOffset;

  /// Temporal offset of output image (in addition to output _torigin)
  double _OutputTimeOffset;

  /// Flag whether to invert transformation
  int _Invert;

  /// Flag whether input is 2D
  int _2D;

  /// Initialize filter
  virtual void Initialize();

public:

  /** Constructor. This constructs an transformation filter with a given
   *  interpolation mode and padding value. By default the interpolation
   *  mode is set to trilinear and invert is off
   */
  irtkImageTransformation();

  /** Static constructor. This constructs an transformation filter for a
   *  given transformation. If there is a fast transformation filter, it
   *  will use this fast transformation filter, otherwise the standard
   *  transformation filter is used.
   */
  static irtkImageTransformation *New(irtkTransformation *);

  /// Destructor
  virtual ~irtkImageTransformation();

  /// Sets input image
  virtual void SetInput(irtkImage *);

  /// Sets input image and transformation
  virtual void SetInput(irtkImage *, irtkTransformation *);

  /// Sets output image
  virtual void SetOutput(irtkImage *);

  /// Sets transformation
  virtual void SetTransformation(irtkTransformation *);

  /// Set displacement field cache
  virtual void SetCache(irtkImageTransformationCache *);

  /// Puts the target padding value
  virtual void PutTargetPaddingValue(double);

  /// Gets the target padding value
  virtual double GetTargetPaddingValue();

  /// Puts the source padding value
  virtual void PutSourcePaddingValue(double);

  /// Gets the source padding value
  virtual double GetSourcePaddingValue();

  /// Gets the interpolator
  virtual irtkImageFunction *GetInterpolator(void);

  /// Sets the interpolator
  virtual void PutInterpolator(irtkImageFunction *);

  /// Sets scale and offset
  virtual void PutScaleFactorAndOffset(double, double);

  /// Sets temporal offset for input image
  virtual void PutInputTimeOffset(double);

  /// Sets temporal offset for output image
  virtual void PutOutputTimeOffset(double);

  /// Invert on
  virtual void InvertOn(void);

  /// Invert off
  virtual void InvertOff(void);

  /// 2D mode on
  virtual void TwoDOn(void);

  /// 2D mode off
  virtual void TwoDOff(void);

  /// Runs the filter
  virtual void Run();

};


#endif
