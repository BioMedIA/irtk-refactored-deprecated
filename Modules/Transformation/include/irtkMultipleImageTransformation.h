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

#ifndef _IRTKMULTIPLEIMAGETRANSFORMATION_H
#define _IRTKMULTIPLEIMAGETRANSFORMATION_H

#include <irtkObject.h>
#include <irtkEnums.h>

// Forward declarations
class irtkBinaryImage;
class irtkInterpolateImageFunction;
class irtkImage;
class irtkTransformation;

/**
 * Transforms multiple images with persistently changing transformation
 *
 * The different images may represent, for example, different channels or
 * additional images computed from a given input image such as the image gradient.
 * It is more efficient to transform all these images together because the
 * mapping from voxel indices to world coordinates, transformation displacements,
 * and interpolation coefficients need only be computed once.
 *
 * This image transformation helper consideres furthermore two transformations:
 * one fixed transformation and one changing transformation. During deformable
 * registration of a source image to a target image, the fixed transformation
 * usually corresponds either to a global transformation obtained by rigid or affine
 * pre-registration or to a fixed deformable transformation of already registered
 * levels which is currently not being optimized for. This transformation is thus
 * assumed to not change whereas the changing transformation may persistently change.
 * This distiction is important when updating the transformed images, in which case
 * the during initialization pre-computed displacements of the fixed transformation
 * can be re-used and only those of the changing transformation have to be re-evaluated.
 *
 * If caching of displacements is enabled, the transformations are evaluated for all
 * voxels at once and the displacements stored in dense displacement field. This is
 * more efficient in case of multiple images or for some transformation models, i.e.,
 * those parameterized by velocities.
 *
 */

class irtkMultipleImageTransformation : public irtkObject
{
  irtkObjectMacro(irtkMultipleImageTransformation);

  // ---------------------------------------------------------------------------
  // Types
public:

  /// Modes of how given transformations are composed with each other
  enum CompositionMode
  {
    COMPOSITION_DEFAULT,
    COMPOSITION_ADDITIVE, ///< Additive transformation model, i.e., T1 + T2
    COMPOSITION_FLUID     ///< Fluid    transformation model, i.e., T1 ° T2 = T2(T1(x))
  };

  /// Type of world coordinates map
  typedef irtkWorldCoordsImage          Image2WorldMap;

  /// Type of dense displacement field
  typedef irtkGenericImage<double>      DisplacementField;

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  irtkMultipleImageTransformation();

  /// Copy constructor
  irtkMultipleImageTransformation(const irtkMultipleImageTransformation &);

  /// Assignment operator
  irtkMultipleImageTransformation &operator =(const irtkMultipleImageTransformation &);

  /// Destructor
  ~irtkMultipleImageTransformation();

  // ---------------------------------------------------------------------------
  // Initialization
private:

  /// Initialize image interpolators
  void InitializeInterpolator();

  /// Clear interpolator
  void ClearInterpolator();

public:

  /// Number of input images
  int NumberOfInputImages() const;

  /// Number of output images
  int NumberOfOutputImages() const;

  /// Number of input(/output) images
  int NumberOfImages() const;

  /// Reset all input images
  void ClearInput();

  /// Set single/first input image
  void SetInput(const irtkImage *);

  /// Add input image
  void AddInput(const irtkImage *);

  /// Set multiple input images at once
  void SetInput(const irtkImage **, int);

  /// Set single/first input image and fixed transformation
  /// For backwards compatibility with irtkImageTransformation
  void SetInput(const irtkImage *, const irtkTransformation *);

  /// Set fixed transformation
  void SetFixedTransformation(const irtkTransformation *);

  /// Get fixed transformation
  const irtkTransformation *GetFixedTransformation() const;

  /// Whether a fixed transformation was set
  bool HasFixedTransformation() const;

  /// Set changing transformation
  void SetTransformation(const irtkTransformation *, CompositionMode = COMPOSITION_DEFAULT);

  /// Whether to apply the inverse or forward transformation
  void Invert(bool);

  /// Apply inverse transformation
  void InvertOn();

  /// Apply forward transformation
  void InvertOff();

  /// Get changine transformation
  const irtkTransformation *GetTransformation() const;

  /// Whether a changing transformation was set
  bool HasTransformation() const;

  /// Set interpolation mode
  void SetInterpolation(irtkInterpolationMode);

  /// Get interpolation mode
  irtkInterpolationMode GetInterpolation() const;

  /// Set intensity rescaling function
  void SetScaleFactorAndOffset(double, double);

  /// Set padding value for voxels outside the overlap region
  ///
  /// By default, if no padding value was set using this method (or it was set to NaN),
  /// the background value of the input is used as padding value. If the input has
  /// no background value, the background value of the output is used instead. If the
  /// output has also no background value set, the padding value defaults to zero.
  ///
  /// It is recommended to set either the background value of the input or output
  /// images instead of using this method.
  void SetPaddingValue(double);

  /// Get padding value for voxels outside the overlap region (NaN if not set)
  ///
  /// The actual padding value used for each input/output image pair, is set as new
  /// background value of the respective output by Initialize. Thus, use
  /// irtkImageBase::GetBackgroundValueAsDouble to retrieve this value.
  double GetPaddingValue() const;

  /// Enable/disable foreground checks
  ///
  /// If foreground checks are enabled, the irtkInterpolateImageFunction::IsForeground
  /// method is used to determine whether an input position is within the foreground
  /// of the (first) input image. Only if all voxels required for the interpolation
  /// are within the foreground region, the output is set to the respective
  /// interpolated value. Otherwise, it is set to the outside value.
  void ForegroundCheck(bool);

  /// Enable foreground checks \sa ForegroundCheck
  void ForegroundCheckOn();

  /// Enable foreground checks \sa ForegroundCheck
  void ForegroundCheckOff();

  /// Reset all output images
  void ClearOutput();

  /// Set single/first output image
  void SetOutput(irtkImage *);

  /// Add output image
  void AddOutput(irtkImage *);

  /// Set multiple output images at once
  void SetOutput(irtkImage **, int);

  /// Initialize filter
  void Initialize();

  /// Clear filter
  void Clear();

  // ---------------------------------------------------------------------------
  // Cache
public:

  /// Use given pre-computed output image to world coordinates map
  void UseWorldCoordinates(const Image2WorldMap *);

  /// Enable/disable pre-computation of world coordinates map
  void CacheWorldCoordinates(bool);

  /// Whether pre-computation of world coordinates map is enabled
  bool CacheWorldCoordinates() const;

  /// Enable pre-computation of world coordinates map
  void CacheWorldCoordinatesOn();

  /// Disable pre-computation of world coordinates map
  void CacheWorldCoordinatesOff();

  /// Enable/disable pre-computation of displacement maps
  void CacheDisplacements(bool);

  /// Whether pre-computation of displacement maps is enabled
  bool CacheDisplacements() const;

  /// Enable pre-computation of displacement maps
  void CacheDisplacementsOn();

  /// Disable pre-computation of displacement maps
  void CacheDisplacementsOff();

protected:

  /// Pre-compute world coordinates of output voxels
  void ComputeWorldCoordinates(bool = false);

  /// Apply fixed transformation to pre-computed world coordinates of output voxels
  void ApplyFixedTransformation();

  /// Pre-compute fixed displacements of output voxels
  void ComputeFixedDisplacement();

  /// Pre-compute changing displacements of output voxels
  void ComputeDisplacement();

  /// Clear cached data
  void ClearCache();

  // ---------------------------------------------------------------------------
  // Transformation
protected:

  /// Update transformed images with known unique input voxel type
  template <class VoxelType>
  void RunWithScalarType(int = -1, int = 1, const irtkBinaryImage * = NULL);

  /// Update transformed images with known interpolator type
  template <class InterpolateImageFunction, class VoxelType>
  void RunWithInterpolator(int = -1, int = 1, const irtkBinaryImage * = NULL);

public:

  /// Update transformed images the first time or force update
  void Run(int = -1, int = 1, const irtkBinaryImage * = NULL);

  /// Update transformed images if required, i.e., changing transformation set
  void Update(int = -1, int = 1, const irtkBinaryImage * = NULL);

  // ---------------------------------------------------------------------------
  // Members
protected:

  // Source images
  vector<const irtkImage *>              _Input;         ///< Input images
  irtkInterpolationMode                  _Interpolation; ///< How to interpolate images
  vector<irtkInterpolateImageFunction *> _Interpolator;  ///< Image interpolator

  int _UniqueScalarType; ///< Unique scalar type of input images or IRTK_VOXEL_UNKNOWN

  // Transformations
  const irtkTransformation *_FixedTransformation; ///< Fixed transformation
  const irtkTransformation *_Transformation;      ///< Changing transformation
  CompositionMode           _Composition;         ///< How to compose transformations
  bool                      _Invert;              ///< Whether to invert the transformations

  // Transformed images
  vector<irtkImage *>   _Output;           ///< Transformed output images
  irtkBinaryImage      *_Mask;             ///< Foreground mask of transformed images
  const Image2WorldMap *_InputImage2World; ///< User supplied world coordinates map

  // Settings
  double           _ScaleFactor;           ///< Rescale slope
  double           _Offset;                ///< Rescale intercept
  double           _DefaultPaddingValue;   ///< User set padding value
  double          *_InputPaddingValue;     ///< Value used for regions outside of overlap
  double           _OutputPaddingValue;    ///< Existing output padding value
  irtkBinaryImage *_RegionOfInterest;      ///< Default output region-of-interest
  bool             _CheckIfForeground;     ///< Only interpolate foreground values
  bool             _CacheWorldCoordinates; ///< Whether to use lookup table of world coordinates
  bool             _CacheDisplacements;    ///< Whether to use dense displacement fields
  bool             _Update;                ///< Whether output needs update

  // Cache - use protected member methods to modify
private:
  Image2WorldMap    *_Image2World;                ///< Lookup table of world coordinates
  DisplacementField *_FixedDisplacement;          ///< Dense fixed displacement field
  DisplacementField *_Displacement;               ///< Dense changing displacement field
  bool               _FixedTransformationApplied; ///< Whether fixed transformation was applied
                                                  ///< to lookup table of world coordinates

  // ---------------------------------------------------------------------------
  // Deprecated - for compatibility with irtkImageTransformation
protected:
  irtkImageFunction *_UnusedInterpolator; ///< Set via deprecated PutInterpolator

public:
  void PutTargetPaddingValue(double);           // Use mask argument of Run/Update
  double GetTargetPaddingValue() const;         // Not supported
  void PutSourcePaddingValue(double);           // Use SetPaddingValue instead
  double GetSourcePaddingValue() const;         // Use GetPaddingValue instead
  irtkImageFunction *GetInterpolator() const;   // Not supported
  void PutInterpolator(irtkImageFunction *);    // Use SetInterpolation instead
  void PutScaleFactorAndOffset(double, double); // Use SetScaleFactorAndOffset instead
  void PutInputTimeOffset(double);              // Not implemented
  void PutOutputTimeOffset(double);             // Not implemented
  void TwoDOn();                                // Not implemented
  void TwoDOff();                               // Not implemented

  // ---------------------------------------------------------------------------
  // Tests
  FRIEND_TEST(irtkImageSimilarityMetric, Initialize);
  FRIEND_TEST(irtkImageSimilarityMetric, Update);
};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline int irtkMultipleImageTransformation::NumberOfInputImages() const
{
  return _Input.size();
}

// -----------------------------------------------------------------------------
inline int irtkMultipleImageTransformation::NumberOfOutputImages() const
{
  return _Output.size();
}

// -----------------------------------------------------------------------------
inline int irtkMultipleImageTransformation::NumberOfImages() const
{
  return NumberOfInputImages();
}

// -----------------------------------------------------------------------------
inline void irtkMultipleImageTransformation::ClearInput()
{
  _Input.resize(0);
}

// -----------------------------------------------------------------------------
inline void irtkMultipleImageTransformation::SetInput(const irtkImage *input)
{
  if (input) {
    _Input.resize(1);
    _Input[0] = input;
  } else {
    _Input.resize(0);
  }
}

// -----------------------------------------------------------------------------
inline void irtkMultipleImageTransformation::AddInput(const irtkImage *input)
{
  if (input) _Input.push_back(input);
}

// -----------------------------------------------------------------------------
inline void irtkMultipleImageTransformation::SetInput(const irtkImage **input, int n)
{
  _Input.resize(0);
  _Input.reserve(n);
  if (n > 0) {
    SetInput(input[0]);
    for (int i = 1; i < n; i++) AddInput(input[i]);
  }
}

// -----------------------------------------------------------------------------
inline void irtkMultipleImageTransformation::SetFixedTransformation(const irtkTransformation *transformation)
{
  _FixedTransformation = transformation;
  _Update              = true;
}

// -----------------------------------------------------------------------------
inline void irtkMultipleImageTransformation::SetInput(const irtkImage *input, const irtkTransformation *transformation)
{
  SetInput(input);
  SetFixedTransformation(transformation);
}

// -----------------------------------------------------------------------------
inline const irtkTransformation *irtkMultipleImageTransformation::GetFixedTransformation() const
{
  return _FixedTransformation;
}

// -----------------------------------------------------------------------------
inline bool irtkMultipleImageTransformation::HasFixedTransformation() const
{
  return static_cast<bool>(GetFixedTransformation());
}

// -----------------------------------------------------------------------------
inline void irtkMultipleImageTransformation::SetTransformation(const irtkTransformation *transformation, CompositionMode m)
{
  _Transformation = transformation;
  _Composition    = m;
  _Update         = true;
}

// -----------------------------------------------------------------------------
inline const irtkTransformation *irtkMultipleImageTransformation::GetTransformation() const
{
  return _Transformation;
}

// -----------------------------------------------------------------------------
inline bool irtkMultipleImageTransformation::HasTransformation() const
{
  return static_cast<bool>(GetTransformation());
}

// -----------------------------------------------------------------------------
inline void irtkMultipleImageTransformation::SetInterpolation(irtkInterpolationMode m)
{
  _Interpolation = m;
}

// -----------------------------------------------------------------------------
inline irtkInterpolationMode irtkMultipleImageTransformation::GetInterpolation() const
{
  return _Interpolation;
}

// -----------------------------------------------------------------------------
inline void irtkMultipleImageTransformation::SetScaleFactorAndOffset(double scale, double offset)
{
  _ScaleFactor = scale;
  _Offset      = offset;
}

// -----------------------------------------------------------------------------
inline void irtkMultipleImageTransformation::ClearOutput()
{
  _Output.resize(0);
}

// -----------------------------------------------------------------------------
inline void irtkMultipleImageTransformation::SetOutput(irtkImage *output)
{
  if (output) {
    _Output.resize(1);
    _Output[0] = output;
  } else {
    _Output.resize(0);
  }
}

// -----------------------------------------------------------------------------
inline void irtkMultipleImageTransformation::AddOutput(irtkImage *output)
{
  if (output) _Output.push_back(output);
}

// -----------------------------------------------------------------------------
inline void irtkMultipleImageTransformation::SetOutput(irtkImage **output, int n)
{
  _Output.reserve(n);
  _Output.resize(0);
  if (n > 0) {
    SetOutput(output[0]);
    for (int i = 1; i < n; i++) AddOutput(output[i]);
  }
}

// -----------------------------------------------------------------------------
inline void irtkMultipleImageTransformation::SetPaddingValue(double value)
{
  _DefaultPaddingValue = value;
}

// -----------------------------------------------------------------------------
inline double irtkMultipleImageTransformation::GetPaddingValue() const
{
  return _DefaultPaddingValue;
}

// -----------------------------------------------------------------------------
inline void irtkMultipleImageTransformation::ForegroundCheck(bool b)
{
  _CheckIfForeground = b;
}

// -----------------------------------------------------------------------------
inline void irtkMultipleImageTransformation::ForegroundCheckOn()
{
  _CheckIfForeground = true;
}

// -----------------------------------------------------------------------------
inline void irtkMultipleImageTransformation::ForegroundCheckOff()
{
  _CheckIfForeground = false;
}

// -----------------------------------------------------------------------------
inline void irtkMultipleImageTransformation::Invert(bool invert)
{
  _Invert = invert;
}

// -----------------------------------------------------------------------------
inline void irtkMultipleImageTransformation::InvertOn()
{
  _Invert = true;
}

// -----------------------------------------------------------------------------
inline void irtkMultipleImageTransformation::InvertOff()
{
  _Invert = false;
}

// -----------------------------------------------------------------------------
inline void irtkMultipleImageTransformation::UseWorldCoordinates(const Image2WorldMap *i2w)
{
  _InputImage2World = i2w;
}

// -----------------------------------------------------------------------------
inline void irtkMultipleImageTransformation::CacheWorldCoordinates(bool cache)
{
  _CacheWorldCoordinates = cache;
}

// -----------------------------------------------------------------------------
inline bool irtkMultipleImageTransformation::CacheWorldCoordinates() const
{
  return _CacheWorldCoordinates;
}

// -----------------------------------------------------------------------------
inline void irtkMultipleImageTransformation::CacheWorldCoordinatesOn()
{
  _CacheWorldCoordinates = true;
}

// -----------------------------------------------------------------------------
inline void irtkMultipleImageTransformation::CacheWorldCoordinatesOff()
{
  _CacheWorldCoordinates = false;
}

// -----------------------------------------------------------------------------
inline void irtkMultipleImageTransformation::CacheDisplacements(bool cache)
{
  _CacheDisplacements = cache;
}

// -----------------------------------------------------------------------------
inline bool irtkMultipleImageTransformation::CacheDisplacements() const
{
  return _CacheDisplacements;
}

// -----------------------------------------------------------------------------
inline void irtkMultipleImageTransformation::CacheDisplacementsOn()
{
  _CacheDisplacements = true;
}

// -----------------------------------------------------------------------------
inline void irtkMultipleImageTransformation::CacheDisplacementsOff()
{
  _CacheDisplacements = false;
}


#endif
