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

#include <irtkMultipleImageTransformation.h>
#include <irtkImage.h>
#include <irtkImageFunction.h>
#include <irtkTransformation.h>

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkMultipleImageTransformation::irtkMultipleImageTransformation()
:
  _Interpolation             (Interpolation_Linear),
  _FixedTransformation       (NULL),
  _Transformation            (NULL),
  _Composition               (COMPOSITION_DEFAULT),
  _Invert                    (false),
  _Mask                      (NULL),
  _InputImage2World          (NULL),
  _ScaleFactor               (1.0),
  _Offset                    (0.0),
  _DefaultPaddingValue       (numeric_limits<double>::quiet_NaN()),
  _InputPaddingValue         (NULL),
  _OutputPaddingValue        (numeric_limits<double>::quiet_NaN()),
  _RegionOfInterest          (NULL),
  _CheckIfForeground         (false),
  _CacheWorldCoordinates     (false),
  _CacheDisplacements        (false),
  _Update                    (true),
  _Image2World               (NULL),
  _FixedDisplacement         (NULL),
  _Displacement              (NULL),
  _FixedTransformationApplied(false),
  _UnusedInterpolator        (NULL)
{
}

// -----------------------------------------------------------------------------
irtkMultipleImageTransformation::irtkMultipleImageTransformation(const irtkMultipleImageTransformation &other)
:
  _FixedTransformation(NULL),
  _Transformation     (NULL),
  _Mask               (NULL),
  _InputImage2World   (NULL),
  _InputPaddingValue  (NULL),
  _RegionOfInterest   (NULL),
  _Image2World        (NULL),
  _FixedDisplacement  (NULL),
  _Displacement       (NULL),
  _UnusedInterpolator (NULL)
{
  *this = other;
}

// -----------------------------------------------------------------------------
irtkMultipleImageTransformation &irtkMultipleImageTransformation::operator =(const irtkMultipleImageTransformation &rhs)
{
  // Copy input/output
  _Input                 = rhs._Input;
  _Interpolation         = rhs._Interpolation;
  _FixedTransformation   = rhs._FixedTransformation;
  _Transformation        = rhs._Transformation;
  _Composition           = rhs._Composition;
  _Invert                = rhs._Invert;
  _Output                = rhs._Output;
  _Mask                  = rhs._Mask;
  _InputImage2World      = rhs._InputImage2World;

  // Copy settings
  _ScaleFactor           = rhs._ScaleFactor;
  _Offset                = rhs._Offset;
  _DefaultPaddingValue   = rhs._DefaultPaddingValue;
  if (rhs._InputPaddingValue) {
    _InputPaddingValue = new double[_Output.size()];
    memcpy(_InputPaddingValue, rhs._InputPaddingValue, _Output.size() * sizeof(double));
  } else {
    _InputPaddingValue = NULL;
  }
  _OutputPaddingValue    = rhs._OutputPaddingValue;
  _RegionOfInterest      = rhs._RegionOfInterest ? new irtkBinaryImage(*rhs._RegionOfInterest) : NULL;
  _CheckIfForeground     = rhs._CheckIfForeground;
  _CacheWorldCoordinates = rhs._CacheWorldCoordinates;
  _CacheDisplacements    = rhs._CacheDisplacements;
  _Update                = rhs._Update;

  // Copy cache
  ClearCache();
  _Image2World                = rhs._Image2World       ? new Image2WorldMap   (*rhs._Image2World)       : NULL;
  _FixedDisplacement          = rhs._FixedDisplacement ? new DisplacementField(*rhs._FixedDisplacement) : NULL;
  _Displacement               = rhs._Displacement      ? new DisplacementField(*rhs._Displacement)      : NULL;
  _FixedTransformationApplied = rhs._FixedTransformationApplied;

  // Copy deprecated members
  _UnusedInterpolator = rhs._UnusedInterpolator;

  // Initialize interpolators if rhs.Initialize has been called already
  if (rhs._Interpolator.size() > 0) InitializeInterpolator();
  else                              ClearInterpolator();
  return *this;
}

// -----------------------------------------------------------------------------
irtkMultipleImageTransformation::~irtkMultipleImageTransformation()
{
  Clear();
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
void irtkMultipleImageTransformation::InitializeInterpolator()
{
  ClearInterpolator();
  _Interpolator.resize(_Input.size());
  for (size_t n = 0; n < _Input.size(); n++) {
    _Interpolator[n] = irtkInterpolateImageFunction::New(_Interpolation, const_cast<irtkImage *>(_Input[n]));
    _Interpolator[n]->SetInput(const_cast<irtkImage *>(_Input[n]));
    _Interpolator[n]->Initialize();
  }
}

// -----------------------------------------------------------------------------
void irtkMultipleImageTransformation::ClearInterpolator()
{
  for (size_t n = 0; n < _Interpolator.size(); n++) delete _Interpolator[n];
  _Interpolator.clear();
}

// -----------------------------------------------------------------------------
void irtkMultipleImageTransformation
::Initialize()
{
  // Clear previously initialized data
  Clear();

  // Check inputs
  if (_Input.size() == 0) {
    cerr << "irtkMultipleImageTransformation::Initialize: Missing input image" << endl;
    exit(1);
  }
  _UniqueScalarType = _Input[0]->GetScalarType();
  for (size_t n = 1; n < _Input.size(); ++n) {
    // Image attributes must be identical except for the number of components
    //
    // This is required because the Update methods check only once if a given
    // world coordinate is inside the image volume using the attributes of the
    // first input image for faster execution. As transformations usually are
    // applied to images defined in the same coordinate space, this is no real
    // limitation of this filter. Otherwise, the user can always use multiple
    // instances of this filter with different inputs.
    if (!_Input[n]->HasSpatialAttributesOf(_Input[0])) {
      cerr << "irtkMultipleImageTransformation::Initialize: Mismatch of spatial domain attributes of input images" << endl;
      exit(1);
    }
    // Invalidate unique scalar type if input data type differs
    if (_Input[n]->GetScalarType() != _UniqueScalarType) _UniqueScalarType = IRTK_VOXEL_UNKNOWN;
  }
  for (size_t n = 0; n < _Input.size(); ++n) {
    if (_Input[n]->GetT() > 1 && _Input[n]->GetTSize() != .0) {
      cerr << "irtkMultipleImageTransformation::Initialize: Input images must be channels/vector images corresponding to the same temporal frame" << endl;
      exit(1);
    }
  }

  // Check outputs
  if (_Output.size() == 0) {
    cerr << "irtkMultipleImageTransformation::Initialize: Missing output image" << endl;
    exit(1);
  }
  if (_Output.size() != _Input.size()) {
    cerr << "irtkMultipleImageTransformation::Initialize: Number of outputs does not match number of inputs" << endl;
    exit(1);
  }
  for (size_t n = 1; n < _Output.size(); ++n) {
    if (!_Output[n]->HasSpatialAttributesOf(_Output[0])) {
      cerr << "irtkMultipleImageTransformation::Initialize: Mismatch of spatial domain attributes of output images" << endl;
      exit(1);
    }
  }
  for (size_t n = 0; n < _Output.size(); ++n) {
    if (_Output[n]->GetT() > 1 && _Output[n]->GetTSize() != .0) {
      cerr << "irtkMultipleImageTransformation::Initialize: Output images must be channels/vector images corresponding to the same temporal frame" << endl;
      exit(1);
    }
  }

  // Set default composition mode
  if (_Composition == COMPOSITION_DEFAULT) {
    _Composition = (HasFixedTransformation() && strstr(GetFixedTransformation()->NameOfClass(), "Fluid") != NULL)
                   ? COMPOSITION_FLUID : COMPOSITION_ADDITIVE;
  }

  // Initialize interpolators
  InitializeInterpolator();

  // Set background value of outputs
  _InputPaddingValue = new double[_Output.size()];
  _Mask              = NULL;
  for (size_t n = 0; n < _Output.size(); ++n) {
    if (IsNaN(_DefaultPaddingValue)) {
      if      (_Input [n]->HasBackgroundValue()) _InputPaddingValue[n] = _Input [n]->GetBackgroundValueAsDouble();
      else if (_Output[n]->HasBackgroundValue()) _InputPaddingValue[n] = _Output[n]->GetBackgroundValueAsDouble();
      else                                       _InputPaddingValue[n] = .0;
    } else {
      _InputPaddingValue[n] = _DefaultPaddingValue;
    }
    _Output[n]->PutBackgroundValueAsDouble(_InputPaddingValue[n]);
    if (!_Mask && _Output[n]->HasMask()) _Mask = _Output[n]->GetMask();
  }

  // Initialize default region-of-interest
  if (!IsNaN(_OutputPaddingValue)) {
    _RegionOfInterest = new irtkBinaryImage(_Output[0]->GetImageAttributes(), 1);
    for (int idx = 0; idx < _RegionOfInterest->GetNumberOfVoxels(); ++idx) {
      _RegionOfInterest->Put(idx, static_cast<irtkBinaryPixel>(_Output[0]->GetAsDouble(idx) > _OutputPaddingValue));
    }
  }

  _Update = true;
}

// -----------------------------------------------------------------------------
void irtkMultipleImageTransformation::Clear()
{
  delete[] _InputPaddingValue; _InputPaddingValue = NULL;
  delete   _RegionOfInterest;  _RegionOfInterest  = NULL;
  ClearCache();
  ClearInterpolator();
}

// =============================================================================
// Cache
// =============================================================================

// -----------------------------------------------------------------------------
void irtkMultipleImageTransformation::ClearCache()
{
  if (_Image2World != _InputImage2World) delete _Image2World;
  _Image2World = NULL;
  delete _FixedDisplacement; _FixedDisplacement = NULL;
  delete _Displacement;      _Displacement      = NULL;
  _FixedTransformationApplied = false;
}

// -----------------------------------------------------------------------------
void irtkMultipleImageTransformation::ComputeWorldCoordinates(bool overwrite)
{
  if (_InputImage2World && (!_Image2World || overwrite)) {
    _Image2World = const_cast<Image2WorldMap *>(_InputImage2World);
    _FixedTransformationApplied = false;
  } else if (_CacheWorldCoordinates) {
    if (!_Image2World || overwrite) {
      if (!_Image2World) _Image2World = new Image2WorldMap();
      _Output[0]->ImageToWorld(*_Image2World);
      _FixedTransformationApplied = false;
    }
  } else {
    delete _Image2World;
    _Image2World = NULL;
  }
}

// -----------------------------------------------------------------------------
void irtkMultipleImageTransformation::ApplyFixedTransformation()
{
  if (!_FixedTransformationApplied) {
    // If input world coordinates map given, copy it now
    if (_InputImage2World && (!_Image2World || _Image2World == _InputImage2World)) {
      _Image2World = new Image2WorldMap(*_InputImage2World);
    // Force pre-computation of world coordinates if not done yet
    } else {
      bool cache = _CacheWorldCoordinates;
      _CacheWorldCoordinates = true;
      ComputeWorldCoordinates();
      _CacheWorldCoordinates = cache;
    }
    // Apply fixed transformation to world coordinates
    const double t  = _Input [0]->ImageToTime(.0);
    const double t0 = _Output[0]->ImageToTime(.0);
    Image2WorldMap::VoxelType *wx = _Image2World->GetPointerToVoxels(0, 0, 0, 0);
    Image2WorldMap::VoxelType *wy = _Image2World->GetPointerToVoxels(0, 0, 0, 1);
    Image2WorldMap::VoxelType *wz = _Image2World->GetPointerToVoxels(0, 0, 0, 2);
    for (int k = 0; k < _Output[0]->GetZ(); k++) {
      for (int j = 0; j < _Output[0]->GetY(); j++) {
        for (int i = 0; i < _Output[0]->GetX(); i++) {
          if (_Invert) _FixedTransformation->Inverse  (*wx, *wy, *wz, t, t0);
          else         _FixedTransformation->Transform(*wx, *wy, *wz, t, t0);
          ++wx, ++wy, ++wz;
        }
      }
    }
    // Clear fixed displacements
    delete _FixedDisplacement;
    _FixedDisplacement          = NULL;
    _FixedTransformationApplied = true;
  }
}

// -----------------------------------------------------------------------------
void irtkMultipleImageTransformation::ComputeFixedDisplacement()
{
  if (_CacheDisplacements && !_FixedTransformationApplied) {
    const double t  = _Input [0]->ImageToTime(.0);
    const double t0 = _Output[0]->ImageToTime(.0);
    if (!_FixedDisplacement) _FixedDisplacement = new DisplacementField();
    _FixedDisplacement->Initialize(_Output[0]->GetImageAttributes(), 3);
    if (_Invert) _FixedTransformation->InverseDisplacement(*_FixedDisplacement, t, t0, _Image2World);
    else         _FixedTransformation->Displacement       (*_FixedDisplacement, t, t0, _Image2World);
  } else {
    delete _FixedDisplacement;
    _FixedDisplacement = NULL;
  }
}

// -----------------------------------------------------------------------------
void irtkMultipleImageTransformation::ComputeDisplacement()
{
  if (_CacheDisplacements) {
    const double t  = _Input [0]->ImageToTime(.0);
    const double t0 = _Output[0]->ImageToTime(.0);
    if (!_Displacement) _Displacement = new DisplacementField();
    _Displacement->Initialize(_Output[0]->GetImageAttributes(), 3);
    if (_Composition == COMPOSITION_FLUID && _Invert) {
      ComputeWorldCoordinates(true);
      _Transformation     ->InverseDisplacement(*_Displacement, t, t0, _Image2World);
      _FixedTransformation->InverseDisplacement(*_Displacement, t, t0, _Image2World);
    } else {
      irtkWorldCoordsImage *wc = _FixedTransformationApplied ? NULL : _Image2World;
      if (_Invert) _Transformation->InverseDisplacement(*_Displacement, t, t0, wc);
      else         _Transformation->Displacement       (*_Displacement, t, t0, wc);
    }
  } else {
    delete _Displacement;
    _Displacement = NULL;
  }
}

// =============================================================================
// Auxiliary macros for irtkMultipleImageTransformationUpdateBody implementation
// =============================================================================

// -----------------------------------------------------------------------------
#define _VARIABLES_AND_CONSTANTS1                                              \
  const irtkImage          *ri = (*_Input) [0];                                \
  InterpolateImageFunction *rf                                                 \
      = reinterpret_cast<InterpolateImageFunction *>(                          \
          const_cast<irtkInterpolateImageFunction *>((*_Interpolator)[0])      \
        );                                                                     \
                                                                               \
  const int bi = re.cols ().begin();                                           \
  const int bj = re.rows ().begin();                                           \
  const int bk = re.pages().begin();                                           \
  const int ei = re.cols ().end();                                             \
  const int ej = re.rows ().end();                                             \
  const int ek = re.pages().end()
#define _VARIABLES_AND_CONSTANTS2                                              \
  const int s1 =  _X - (ei - bi);                                              \
  const int s2 = (_Y - (ej - bj)) * _X - (ei - bi)
#define _VARIABLES_AND_CONSTANTS                                               \
  _VARIABLES_AND_CONSTANTS1;                                                   \
  _VARIABLES_AND_CONSTANTS2

// -----------------------------------------------------------------------------
#define _INIT_REGION_ITERATORS(type, image)                                    \
    const type::VoxelType *image##_x = image->GetPointerToVoxels(bi, bj, bk);  \
    const type::VoxelType *image##_y = image##_x + _NumberOfVoxels;            \
    const type::VoxelType *image##_z = image##_y + _NumberOfVoxels
// -----------------------------------------------------------------------------
#define _INIT_NON_CONST_REGION_ITERATORS(type, image)                          \
    type::VoxelType *image##_x = image->GetPointerToVoxels(bi, bj, bk);        \
    type::VoxelType *image##_y = image##_x + _NumberOfVoxels;                  \
    type::VoxelType *image##_z = image##_y + _NumberOfVoxels

// -----------------------------------------------------------------------------
#define _MOVE_REGION_ITERATORS(image, inc)                                     \
  image##_x += inc, image##_y += inc, image##_z += inc
#define _MOVE_REGION_ITERATORS1(image)                                         \
  _MOVE_REGION_ITERATORS(image,  1); } /* end of loop over i */                \
  _MOVE_REGION_ITERATORS(image, s1); } /* end of loop over j */                \
  _MOVE_REGION_ITERATORS(image, s2); } /* end of loop over k */
#define _MOVE_REGION_ITERATORS2(image1, image2)                                \
  _MOVE_REGION_ITERATORS(image1,  1);                                          \
  _MOVE_REGION_ITERATORS(image2,  1); } /* end of loop over i */               \
  _MOVE_REGION_ITERATORS(image1, s1);                                          \
  _MOVE_REGION_ITERATORS(image2, s1); } /* end of loop over j */               \
  _MOVE_REGION_ITERATORS(image1, s2);                                          \
  _MOVE_REGION_ITERATORS(image2, s2); } /* end of loop over k */
#define _MOVE_REGION_ITERATORS3(image1, image2, image3)                        \
  _MOVE_REGION_ITERATORS(image1,  1);                                          \
  _MOVE_REGION_ITERATORS(image2,  1);                                          \
  _MOVE_REGION_ITERATORS(image3,  1); } /* end of loop over i */               \
  _MOVE_REGION_ITERATORS(image1, s1);                                          \
  _MOVE_REGION_ITERATORS(image2, s1);                                          \
  _MOVE_REGION_ITERATORS(image3, s1); } /* end of loop over j */               \
  _MOVE_REGION_ITERATORS(image1, s2);                                          \
  _MOVE_REGION_ITERATORS(image2, s2);                                          \
  _MOVE_REGION_ITERATORS(image3, s2); } /* end of loop over k */
#define _MOVE_REGION_ITERATORS4(image1, image2, image3, image4)                \
  _MOVE_REGION_ITERATORS(image1,  1);                                          \
  _MOVE_REGION_ITERATORS(image2,  1);                                          \
  _MOVE_REGION_ITERATORS(image3,  1);                                          \
  _MOVE_REGION_ITERATORS(image4,  1); } /* end of loop over i */               \
  _MOVE_REGION_ITERATORS(image1, s1);                                          \
  _MOVE_REGION_ITERATORS(image2, s1);                                          \
  _MOVE_REGION_ITERATORS(image3, s1);                                          \
  _MOVE_REGION_ITERATORS(image4, s1); } /* end of loop over j */               \
  _MOVE_REGION_ITERATORS(image1, s2);                                          \
  _MOVE_REGION_ITERATORS(image2, s2);                                          \
  _MOVE_REGION_ITERATORS(image3, s2);                                          \
  _MOVE_REGION_ITERATORS(image4, s2); } /* end of loop over k */

// -----------------------------------------------------------------------------
#define _INVALID_COORD numeric_limits<Image2WorldMap::VoxelType>::quiet_NaN()
#define _HEADER()                                                              \
  double x, y, z;                                                              \
  for (int k = bk; k < ek; k++) {                                              \
    for (int j = bj; j < ej;  j++) {                                           \
      for (int i = bi; i < ei; i++) {                                          \
        if (!_RegionOfInterest || _RegionOfInterest->Get(i, j, k)) {           \
          do {} while(false)
//          index to world coordinate transformation code
#define _TRANSFORM_FOOTER(IsInsideImpl)                                        \
          ri->WorldToImage(x, y, z);                                           \
          if (!rf->IsInsideImpl(x, y, z)) x = _INVALID_COORD;                  \
        } else x = _INVALID_COORD;                                             \
        *_Output2Input_x = x, *_Output2Input_y = y, *_Output2Input_z = z
// or
#define _UPDATE_FOOTER(IsInsideImpl, ...)                                      \
          ri->WorldToImage(x, y, z);                                           \
          if (rf->IsInsideImpl(x, y, z)) {                                     \
            __VA_ARGS__                                                        \
          } else PutOutsideValue(i, j, k);                                     \
        } else PutOutsideValue(i, j, k)
//} } } /* end of loop over i, j, k */ - see _MOVE_REGION_ITERATORS

// -----------------------------------------------------------------------------
#define _TRANSFORM_FOOTER0(IsInsideImpl)                                       \
  _TRANSFORM_FOOTER(IsInsideImpl);                                             \
  _MOVE_REGION_ITERATORS1(_Output2Input)
// or
#define _TRANSFORM_FOOTER1(IsInsideImpl, image1)                               \
  _TRANSFORM_FOOTER(IsInsideImpl);                                             \
  _MOVE_REGION_ITERATORS2(_Output2Input, image1)
// or
#define _TRANSFORM_FOOTER2(IsInsideImpl, image1, image2)                       \
  _TRANSFORM_FOOTER(IsInsideImpl);                                             \
  _MOVE_REGION_ITERATORS3(_Output2Input, image1, image2)
// or
#define _TRANSFORM_FOOTER3(IsInsideImpl, image1, image2, image3)               \
  _TRANSFORM_FOOTER(IsInsideImpl);                                             \
  _MOVE_REGION_ITERATORS4(_Output2Input, image1, image2, image3)

// -----------------------------------------------------------------------------
#define _UPDATE_FOOTER0(IsInsideImpl, ...)                                     \
  _UPDATE_FOOTER(IsInsideImpl, __VA_ARGS__);                                   \
  } } } /* end of loop over i, j, k */ do {} while(false)
// or
#define _UPDATE_FOOTER1(IsInsideImpl, image1, ...)                             \
  _UPDATE_FOOTER(IsInsideImpl, __VA_ARGS__);                                   \
  _MOVE_REGION_ITERATORS1(image1)
// or
#define _UPDATE_FOOTER2(IsInsideImpl, image1, image2, ...)                     \
  _UPDATE_FOOTER(IsInsideImpl, __VA_ARGS__);                                   \
  _MOVE_REGION_ITERATORS2(image1, image2)
// or
#define _UPDATE_FOOTER3(IsInsideImpl, image1, image2, image3, ...)             \
  _UPDATE_FOOTER(IsInsideImpl, __VA_ARGS__);                                   \
  _MOVE_REGION_ITERATORS3(image1, image2, image3)

// -----------------------------------------------------------------------------
//void TransformUsingXXXXX(const blocked_range3d<int> &re)
//{
#define _TRANSFORM0(...)                                                       \
  _VARIABLES_AND_CONSTANTS;                                                    \
  _INIT_NON_CONST_REGION_ITERATORS(Image2WorldMap, _Output2Input);             \
  if (_CheckIfForeground) {                                                    \
    _HEADER();                               /* start of loop over i, j, k */  \
      __VA_ARGS__                            /* coordinate tranformation   */  \
    _TRANSFORM_FOOTER0(IsForeground);        /* end of loop over i, j, k   */  \
  } else {                                                                     \
    _HEADER();                               /* start of loop over i, j, k */  \
      __VA_ARGS__                            /* coordinate tranformation   */  \
    _TRANSFORM_FOOTER0(IsInside);            /* end of loop over i, j, k   */  \
  } do {} while(false)
//}

// -----------------------------------------------------------------------------
//void TransformUsingXXXXX(const blocked_range3d<int> &re)
//{
#define _TRANSFORM1(type1, image1, ...)                                        \
  _VARIABLES_AND_CONSTANTS;                                                    \
  _INIT_NON_CONST_REGION_ITERATORS(Image2WorldMap, _Output2Input);             \
  _INIT_REGION_ITERATORS(type1, image1);                                       \
  if (_CheckIfForeground) {                                                    \
    _HEADER();                               /* start of loop over i, j, k */  \
      __VA_ARGS__                            /* coordinate tranformation   */  \
    _TRANSFORM_FOOTER1(IsForeground, image1);/* end of loop over i, j, k   */  \
  } else {                                                                     \
    _HEADER();                               /* start of loop over i, j, k */  \
      __VA_ARGS__                            /* coordinate tranformation   */  \
    _TRANSFORM_FOOTER1(IsInside, image1);    /* end of loop over i, j, k   */  \
  } do {} while(false)
//}

// -----------------------------------------------------------------------------
//void TransformUsingXXXXX(const blocked_range3d<int> &re)
//{
#define _TRANSFORM2(type1, image1, type2, image2, ...)                         \
  _VARIABLES_AND_CONSTANTS;                                                    \
  _INIT_NON_CONST_REGION_ITERATORS(Image2WorldMap, _Output2Input);             \
  _INIT_REGION_ITERATORS(type1, image1);                                       \
  _INIT_REGION_ITERATORS(type2, image2);                                       \
  if (_CheckIfForeground) {                                                    \
    _HEADER();                                /* start of loop over i, j, k */ \
      __VA_ARGS__                             /* coordinate tranformation   */ \
    _TRANSFORM_FOOTER2(IsForeground, image1, image2); /* end of loop        */ \
  } else {                                                                     \
    _HEADER();                                /* start of loop over i, j, k */ \
      __VA_ARGS__                             /* coordinate tranformation   */ \
    _TRANSFORM_FOOTER2(IsInside, image1, image2); /* end of loop over i,j,k */ \
  } do {} while(false)
//}

// -----------------------------------------------------------------------------
//void TransformUsingXXXXX(const blocked_range3d<int> &re)
//{
#define _TRANSFORM3(type1, image1, type2, image2, type3, image3, ...)          \
  _VARIABLES_AND_CONSTANTS;                                                    \
  _INIT_NON_CONST_REGION_ITERATORS(Image2WorldMap, _Output2Input);             \
  _INIT_REGION_ITERATORS(type1, image1);                                       \
  _INIT_REGION_ITERATORS(type2, image2);                                       \
  _INIT_REGION_ITERATORS(type3, image3);                                       \
  if (_CheckIfForeground) {                                                    \
    _HEADER();                                /* start of loop over i, j, k */ \
      __VA_ARGS__                             /* coordinate tranformation   */ \
    _TRANSFORM_FOOTER3(IsForeground, image1, image2, image3); /* end of loop*/ \
  } else {                                                                     \
    _HEADER();                                /* start of loop over i, j, k */ \
      __VA_ARGS__                             /* coordinate tranformation   */ \
    _TRANSFORM_FOOTER3(IsInside, image1, image2, image3); /* end of loop    */ \
  } do {} while(false)
//}

// -----------------------------------------------------------------------------
//void UpdateUsingXXXXX(const blocked_range3d<int> &re)
//{
#define _UPDATE0(...)                                                          \
  _VARIABLES_AND_CONSTANTS1;                                                   \
  if (_CheckIfForeground) {                                                    \
    if (_ScaleFactor == 1.0 && _Offset == 0.0) {                               \
      if (typeid(OutputVoxelType1) == typeid(double)) {                        \
        _HEADER();                            /* start of loop over i, j, k */ \
          __VA_ARGS__                         /* coordinate tranformation   */ \
        _UPDATE_FOOTER0(IsForeground,                                          \
          UpdateDoubleWithoutRescaling(i, j, k, x, y, z);                      \
        );                                    /* end of loop over i, j, k   */ \
      } else {                                                                 \
        double *v = new double[_MaxNumberOfComponents];                        \
        _HEADER();                            /* start of loop over i, j, k */ \
          __VA_ARGS__                         /* coordinate tranformation   */ \
        _UPDATE_FOOTER0(IsForeground,                                          \
          UpdateWithoutRescaling(i, j, k, x, y, z, v);                         \
        );                                    /* end of loop over i, j, k   */ \
        delete[] v;                                                            \
      }                                                                        \
    } else {                                                                   \
      double *v = new double[_MaxNumberOfComponents];                          \
      _HEADER();                              /* start of loop over i, j, k */ \
        __VA_ARGS__                           /* coordinate tranformation   */ \
      _UPDATE_FOOTER0(IsForeground,                                            \
        UpdateWithRescaling(i, j, k, x, y, z, v);                              \
      );                                      /* end of loop over i, j, k   */ \
      delete[] v;                                                              \
    }                                                                          \
  } else {                                                                     \
    if (_ScaleFactor == 1.0 && _Offset == 0.0) {                               \
      if (typeid(OutputVoxelType1) == typeid(double)) {                        \
        _HEADER();                            /* start of loop over i, j, k */ \
          __VA_ARGS__                         /* coordinate tranformation   */ \
        _UPDATE_FOOTER0(IsInside,                                              \
          UpdateDoubleWithoutRescaling(i, j, k, x, y, z);                      \
        );                                    /* end of loop over i, j, k   */ \
      } else {                                                                 \
        double *v = new double[_MaxNumberOfComponents];                        \
        _HEADER();                            /* start of loop over i, j, k */ \
          __VA_ARGS__                         /* coordinate tranformation   */ \
        _UPDATE_FOOTER0(IsInside,                                              \
          UpdateWithoutRescaling(i, j, k, x, y, z, v);                         \
        );                                    /* end of loop over i, j, k   */ \
        delete[] v;                                                            \
      }                                                                        \
    } else {                                                                   \
      double *v = new double[_MaxNumberOfComponents];                          \
      _HEADER();                              /* start of loop over i, j, k */ \
        __VA_ARGS__                           /* coordinate tranformation   */ \
      _UPDATE_FOOTER0(IsInside,                                                \
        UpdateWithRescaling(i, j, k, x, y, z, v);                              \
      );                                      /* end of loop over i, j, k   */ \
      delete[] v;                                                              \
    }                                                                          \
  } do {} while(false)
//}

// -----------------------------------------------------------------------------
//void UpdateUsingXXXXX(const blocked_range3d<int> &re)
//{
#define _UPDATE1(type1, image1, ...)                                           \
  _VARIABLES_AND_CONSTANTS;                                                    \
  _INIT_REGION_ITERATORS(type1, image1);                                       \
  if (_CheckIfForeground) {                                                    \
    if (_ScaleFactor == 1.0 && _Offset == 0.0) {                               \
      if (typeid(OutputVoxelType1) == typeid(double)) {                        \
        _HEADER();                            /* start of loop over i, j, k */ \
          __VA_ARGS__                         /* coordinate tranformation   */ \
        _UPDATE_FOOTER1(IsForeground, image1,                                  \
          UpdateDoubleWithoutRescaling(i, j, k, x, y, z);                      \
        );                                    /* end of loop over i, j, k   */ \
      } else {                                                                 \
        double *v = new double[_MaxNumberOfComponents];                        \
        _HEADER();                            /* start of loop over i, j, k */ \
          __VA_ARGS__                         /* coordinate tranformation   */ \
        _UPDATE_FOOTER1(IsForeground, image1,                                  \
          UpdateWithoutRescaling(i, j, k, x, y, z, v);                         \
        );                                    /* end of loop over i, j, k   */ \
        delete[] v;                                                            \
      }                                                                        \
    } else {                                                                   \
      double *v = new double[_MaxNumberOfComponents];                          \
      _HEADER();                              /* start of loop over i, j, k */ \
        __VA_ARGS__                           /* coordinate tranformation   */ \
      _UPDATE_FOOTER1(IsForeground, image1,                                    \
        UpdateWithRescaling(i, j, k, x, y, z, v);                              \
      );                                      /* end of loop over i, j, k   */ \
      delete[] v;                                                              \
    }                                                                          \
  } else {                                                                     \
    if (_ScaleFactor == 1.0 && _Offset == 0.0) {                               \
      if (typeid(OutputVoxelType1) == typeid(double)) {                        \
        _HEADER();                            /* start of loop over i, j, k */ \
          __VA_ARGS__                         /* coordinate tranformation   */ \
        _UPDATE_FOOTER1(IsInside, image1,                                      \
          UpdateDoubleWithoutRescaling(i, j, k, x, y, z);                      \
        );                                    /* end of loop over i, j, k   */ \
      } else {                                                                 \
        double *v = new double[_MaxNumberOfComponents];                        \
        _HEADER();                            /* start of loop over i, j, k */ \
          __VA_ARGS__                         /* coordinate tranformation   */ \
        _UPDATE_FOOTER1(IsInside, image1,                                      \
          UpdateWithoutRescaling(i, j, k, x, y, z, v);                         \
        );                                    /* end of loop over i, j, k   */ \
        delete[] v;                                                            \
      }                                                                        \
    } else {                                                                   \
      double *v = new double[_MaxNumberOfComponents];                          \
      _HEADER();                              /* start of loop over i, j, k */ \
        __VA_ARGS__                           /* coordinate tranformation   */ \
      _UPDATE_FOOTER1(IsInside, image1,                                        \
        UpdateWithRescaling(i, j, k, x, y, z, v);                              \
      );                                      /* end of loop over i, j, k   */ \
      delete[] v;                                                              \
    }                                                                          \
  } do {} while(false)
//}

// -----------------------------------------------------------------------------
//void UpdateUsingXXXXX(const blocked_range3d<int> &re)
//{
#define _UPDATE2(type1, image1, type2, image2, ...)                            \
  _VARIABLES_AND_CONSTANTS;                                                    \
  _INIT_REGION_ITERATORS(type1, image1);                                       \
  _INIT_REGION_ITERATORS(type2, image2);                                       \
  if (_CheckIfForeground) {                                                    \
    if (_ScaleFactor == 1.0 && _Offset == 0.0) {                               \
      if (typeid(OutputVoxelType1) == typeid(double)) {                        \
        _HEADER();                            /* start of loop over i, j, k */ \
          __VA_ARGS__                         /* coordinate tranformation   */ \
        _UPDATE_FOOTER2(IsForeground, image1, image2,                          \
          UpdateDoubleWithoutRescaling(i, j, k, x, y, z);                      \
        );                                    /* end of loop over i, j, k   */ \
      } else {                                                                 \
        double *v = new double[_MaxNumberOfComponents];                        \
        _HEADER();                            /* start of loop over i, j, k */ \
          __VA_ARGS__                         /* coordinate tranformation   */ \
        _UPDATE_FOOTER2(IsForeground, image1, image2,                          \
          UpdateWithoutRescaling(i, j, k, x, y, z, v);                         \
        );                                    /* end of loop over i, j, k   */ \
        delete[] v;                                                            \
      }                                                                        \
    } else {                                                                   \
      double *v = new double[_MaxNumberOfComponents];                          \
      _HEADER();                              /* start of loop over i, j, k */ \
        __VA_ARGS__                           /* coordinate tranformation   */ \
      _UPDATE_FOOTER2(IsForeground, image1, image2,                            \
        UpdateWithRescaling(i, j, k, x, y, z, v);                              \
      );                                      /* end of loop over i, j, k   */ \
      delete[] v;                                                              \
    }                                                                          \
  } else {                                                                     \
    if (_ScaleFactor == 1.0 && _Offset == 0.0) {                               \
      if (typeid(OutputVoxelType1) == typeid(double)) {                        \
        _HEADER();                            /* start of loop over i, j, k */ \
          __VA_ARGS__                         /* coordinate tranformation   */ \
        _UPDATE_FOOTER2(IsInside, image1, image2,                              \
          UpdateDoubleWithoutRescaling(i, j, k, x, y, z);                      \
        );                                    /* end of loop over i, j, k   */ \
      } else {                                                                 \
        double *v = new double[_MaxNumberOfComponents];                        \
        _HEADER();                            /* start of loop over i, j, k */ \
          __VA_ARGS__                         /* coordinate tranformation   */ \
        _UPDATE_FOOTER2(IsInside, image1, image2,                              \
          UpdateWithoutRescaling(i, j, k, x, y, z, v);                         \
        );                                    /* end of loop over i, j, k   */ \
        delete[] v;                                                            \
      }                                                                        \
    } else {                                                                   \
      double *v = new double[_MaxNumberOfComponents];                          \
      _HEADER();                              /* start of loop over i, j, k */ \
        __VA_ARGS__                           /* coordinate tranformation   */ \
      _UPDATE_FOOTER2(IsInside, image1, image2,                                \
        UpdateWithRescaling(i, j, k, x, y, z, v);                              \
      );                                      /* end of loop over i, j, k   */ \
      delete[] v;                                                              \
    }                                                                          \
  } do {} while(false)
//}

// -----------------------------------------------------------------------------
//void UpdateUsingXXXXX(const blocked_range3d<int> &re)
//{
#define _UPDATE3(type1, image1, type2, image2, type3, image3, ...)             \
  _VARIABLES_AND_CONSTANTS;                                                    \
  _INIT_REGION_ITERATORS(type1, image1);                                       \
  _INIT_REGION_ITERATORS(type2, image2);                                       \
  _INIT_REGION_ITERATORS(type3, image3);                                       \
  if (_CheckIfForeground) {                                                    \
    if (_ScaleFactor == 1.0 && _Offset == 0.0) {                               \
      if (typeid(OutputVoxelType1) == typeid(double)) {                        \
        _HEADER();                            /* start of loop over i, j, k */ \
          __VA_ARGS__                         /* coordinate tranformation   */ \
        _UPDATE_FOOTER3(IsForeground, image1, image2, image3,                  \
          UpdateDoubleWithoutRescaling(i, j, k, x, y, z);                      \
        );                                    /* end of loop over i, j, k   */ \
      } else {                                                                 \
        double *v = new double[_MaxNumberOfComponents];                        \
        _HEADER();                            /* start of loop over i, j, k */ \
          __VA_ARGS__                         /* coordinate tranformation   */ \
        _UPDATE_FOOTER3(IsForeground, image1, image2, image3,                  \
          UpdateWithoutRescaling(i, j, k, x, y, z, v);                         \
        );                                    /* end of loop over i, j, k   */ \
        delete[] v;                                                            \
      }                                                                        \
    } else {                                                                   \
      double *v = new double[_MaxNumberOfComponents];                          \
      _HEADER();                              /* start of loop over i, j, k */ \
        __VA_ARGS__                           /* coordinate tranformation   */ \
      _UPDATE_FOOTER3(IsForeground, image1, image2, image3,                    \
        UpdateWithRescaling(i, j, k, x, y, z, v);                              \
      );                                      /* end of loop over i, j, k   */ \
      delete[] v;                                                              \
    }                                                                          \
  } else {                                                                     \
    if (_ScaleFactor == 1.0 && _Offset == 0.0) {                               \
      if (typeid(OutputVoxelType1) == typeid(double)) {                        \
        _HEADER();                            /* start of loop over i, j, k */ \
          __VA_ARGS__                         /* coordinate tranformation   */ \
        _UPDATE_FOOTER3(IsInside, image1, image2, image3,                      \
          UpdateDoubleWithoutRescaling(i, j, k, x, y, z);                      \
        );                                    /* end of loop over i, j, k   */ \
      } else {                                                                 \
        double *v = new double[_MaxNumberOfComponents];                        \
        _HEADER();                            /* start of loop over i, j, k */ \
          __VA_ARGS__                         /* coordinate tranformation   */ \
        _UPDATE_FOOTER3(IsInside, image1, image2, image3,                      \
          UpdateWithoutRescaling(i, j, k, x, y, z, v);                         \
        );                                    /* end of loop over i, j, k   */ \
        delete[] v;                                                            \
      }                                                                        \
    } else {                                                                   \
      double *v = new double[_MaxNumberOfComponents];                          \
      _HEADER();                              /* start of loop over i, j, k */ \
        __VA_ARGS__                           /* coordinate tranformation   */ \
      _UPDATE_FOOTER3(IsInside, image1, image2, image3,                        \
        UpdateWithRescaling(i, j, k, x, y, z, v);                              \
      );                                      /* end of loop over i, j, k   */ \
      delete[] v;                                                              \
    }                                                                          \
  } do {} while(false)
//}

// =============================================================================
// irtkMultipleImageTransformationUpdateBody
// =============================================================================

/**
 * Single-/Multi-threaded body of irtkMultipleImageTransformation::Update.
 *
 * Explicitly casting the interpolate image function to the particular type
 * of the interpolator instances tells the compiler at compile time which subclass
 * is being used. It therefore can deduce at compile time which virtual methods
 * will be called and thus inline the code for better performance.
 *
 * \note The integer prefix of the Transform_ and Update_ methods corresponds
 *       to the _Implementation which is encoded by a 4-bit code. Each bit
 *       specifies which of the pre-computed data is available. The entry
 *       operator() selects the appropriate optimized implementation according
 *       to this bit code. Meaning of the bits from most to least significant:
 *       "<_Displacement1><_Displacement2><_Transformation2><_Image2World>".
 */
template <class InterpolateImageFunction = irtkInterpolateImageFunction,
          class OutputVoxelType1         = double>
class irtkMultipleImageTransformationUpdateBody
{
public:

  typedef irtkMultipleImageTransformation::Image2WorldMap    Image2WorldMap;
  typedef irtkMultipleImageTransformation::DisplacementField DisplacementField;

  // -------------------------------------------------------------------------
  // Data members
  vector<const irtkImage *>              *_Input;
  Image2WorldMap                         *_Image2World;
  vector<irtkInterpolateImageFunction *> *_Interpolator;
  const irtkTransformation               *_Transformation1;
  DisplacementField                      *_Displacement1;
  const irtkTransformation               *_Transformation2;
  DisplacementField                      *_Displacement2;
  Image2WorldMap                         *_Output2Input;
  bool                                    _Output2InputOwner;
  bool                                    _Invert;
  double                                  _ScaleFactor;
  double                                  _Offset;
  vector<irtkImage *>                    *_Output;
  const irtkBinaryImage                  *_RegionOfInterest;
  irtkBinaryImage                        *_Mask;
  int                                     _Begin;
  int                                     _End;
  int                                     _X, _Y;
  int                                     _NumberOfVoxels;
  int                                     _MaxNumberOfComponents;
  int                                     _UniqueScalarType;
  double                                 *_InputPaddingValue;
  bool                                    _CheckIfForeground;

  double _t, _t0;

  // -------------------------------------------------------------------------
  /// Default constructor
  irtkMultipleImageTransformationUpdateBody()
  :
    _Input                (NULL),
    _Image2World          (NULL),
    _Interpolator         (NULL),
    _Transformation1      (NULL),
    _Displacement1        (NULL),
    _Transformation2      (NULL),
    _Displacement2        (NULL),
    _Output2Input         (NULL),
    _Output2InputOwner    (false),
    _Invert               (false),
    _ScaleFactor          (1.0),
    _Offset               (0.0),
    _Output               (NULL),
    _RegionOfInterest     (NULL),
    _Mask                 (NULL),
    _Begin                (0),
    _End                  (1),
    _X                    (0),
    _Y                    (0),
    _NumberOfVoxels       (0),
    _MaxNumberOfComponents(0),
    _UniqueScalarType     (voxel_info<OutputVoxelType1>::type()),
    _InputPaddingValue    (NULL),
    _CheckIfForeground    (false),
    _t                    ( 0.0),
    _t0                   (-1.0)
  {
  }

  // -------------------------------------------------------------------------
  /// Copy constructor
  irtkMultipleImageTransformationUpdateBody(const irtkMultipleImageTransformationUpdateBody &other)
  :
    _Input                (other._Input),
    _Image2World          (other._Image2World),
    _Interpolator         (other._Interpolator),
    _Transformation1      (other._Transformation1),
    _Displacement1        (other._Displacement1),
    _Transformation2      (other._Transformation2),
    _Displacement2        (other._Displacement2),
    _Output2Input         (other._Output2Input),
    _Output2InputOwner    (other._Output2InputOwner),
    _Invert               (other._Invert),
    _ScaleFactor          (other._ScaleFactor),
    _Offset               (other._Offset),
    _Output               (other._Output),
    _RegionOfInterest     (other._RegionOfInterest),
    _Mask                 (other._Mask),
    _Begin                (other._Begin),
    _End                  (other._End),
    _X                    (other._X),
    _Y                    (other._Y),
    _NumberOfVoxels       (other._NumberOfVoxels),
    _MaxNumberOfComponents(other._MaxNumberOfComponents),
    _UniqueScalarType     (voxel_info<OutputVoxelType1>::type()),
    _InputPaddingValue    (other._InputPaddingValue),
    _CheckIfForeground    (other._CheckIfForeground),
    _t                    (other._t),
    _t0                   (other._t0)
  {
    if (_Output2InputOwner && _Output2Input) (*_Output2Input) =  (*other._Output2Input);
  }

private:

  // ===========================================================================
  // Two step update with pre-computed output to input index map
  // ===========================================================================

  // ---------------------------------------------------------------------------
  /// Transform indices using lookup table of world coordinates
  void TransformUsingWorldCoordinates(const blocked_range3d<int> &re) const
  {
    _TRANSFORM1(Image2WorldMap, _Image2World,
      x = *_Image2World_x;
      y = *_Image2World_y;
      z = *_Image2World_z;
    );
  }

  // ---------------------------------------------------------------------------
  /// Transform indices using single transformation
  void TransformUsingTransformation(const blocked_range3d<int> &re) const
  {
    _TRANSFORM0(
      x = i, y = j, z = k;
      _Output->front()->ImageToWorld(x, y, z);
      if (_Invert) _Transformation1->Inverse  (x, y, z, _t, _t0);
      else         _Transformation1->Transform(x, y, z, _t, _t0);
    );
  }

  // ---------------------------------------------------------------------------
  /// Transform indices using lookup table of world coordinates and single transformation
  void TransformUsingWorldCoordinatesAndTransformation(const blocked_range3d<int> &re) const
  {
    _TRANSFORM1(Image2WorldMap, _Image2World,
      x = *_Image2World_x;
      y = *_Image2World_y;
      z = *_Image2World_z;
      if (_Invert) _Transformation1->Inverse  (x, y, z, _t, _t0);
      else         _Transformation1->Transform(x, y, z, _t, _t0);
    );
  }

  // ---------------------------------------------------------------------------
  /// Transform indices using single dense displacement field
  void TransformUsingDisplacement(const blocked_range3d<int> &re) const
  {
    _TRANSFORM1(DisplacementField, _Displacement1,
      x = i, y = j, z = k;
      _Output->front()->ImageToWorld(x, y, z);
      x += *_Displacement1_x;
      y += *_Displacement1_y;
      z += *_Displacement1_z;
    );
  }

  // ---------------------------------------------------------------------------
  /// Transform indices using lookup table of world coordinates and single dense
  /// displacement field
  void TransformUsingWorldCoordinatesAndDisplacement(const blocked_range3d<int> &re) const
  {
    _TRANSFORM2(Image2WorldMap, _Image2World, DisplacementField, _Displacement1,
      x = *_Image2World_x + *_Displacement1_x;
      y = *_Image2World_y + *_Displacement1_y;
      z = *_Image2World_z + *_Displacement1_z;
    );
  }

  // ---------------------------------------------------------------------------
  /// Transform using fluid composition of transformations
  void TransformUsingCompositionOfTransformations(const blocked_range3d<int> &re) const
  {
    _TRANSFORM0(
      x = i, y = j, z = k;
      _Output->front()->ImageToWorld(x, y, z);
      if (_Invert) {
        _Transformation1->Inverse  (x, y, z, _t, _t0);
        _Transformation2->Inverse  (x, y, z, _t, _t0);
      } else {
        _Transformation1->Transform(x, y, z, _t, _t0);
        _Transformation2->Transform(x, y, z, _t, _t0);
      }
    );
  }

  // ---------------------------------------------------------------------------
  /// Transform using lookup table of world coordinates and fluid composition of transformations
  void TransformUsingWorldCoordinatesAndCompositionOfTransformations(const blocked_range3d<int> &re) const
  {
    _TRANSFORM1(Image2WorldMap, _Image2World,
      x = *_Image2World_x;
      y = *_Image2World_y;
      z = *_Image2World_z;
      if (_Invert) {
        _Transformation1->Inverse  (x, y, z, _t, _t0);
        _Transformation2->Inverse  (x, y, z, _t, _t0);
      } else {
        _Transformation1->Transform(x, y, z, _t, _t0);
        _Transformation2->Transform(x, y, z, _t, _t0);
      }
    );
  }

  // ---------------------------------------------------------------------------
  /// Transform indices using fluid composition of dense displacement field and transformation
  void TransformUsingCompositionOfDisplacementAndTransformation(const blocked_range3d<int> &re) const
  {
    _TRANSFORM1(DisplacementField, _Displacement1,
      x = i, y = j, z = k;
      _Output->front()->ImageToWorld(x, y, z);
      x += *_Displacement1_x;
      y += *_Displacement1_y;
      z += *_Displacement1_z;
      if (_Invert) _Transformation2->Inverse  (x, y, z, _t, _t0);
      else         _Transformation2->Transform(x, y, z, _t, _t0);
    );
  }

  // ---------------------------------------------------------------------------
  /// Transform indices using lookup table of world coordinates and fluid
  /// composition of dense displacement field and transformation
  void TransformUsingWorldCoordinatesAndCompositionOfDisplacementAndTransformation(const blocked_range3d<int> &re) const
  {
    _TRANSFORM2(Image2WorldMap, _Image2World, DisplacementField, _Displacement1,
      x = *_Image2World_x + *_Displacement1_x;
      y = *_Image2World_y + *_Displacement1_y;
      z = *_Image2World_z + *_Displacement1_z;
      if (_Invert) _Transformation2->Inverse  (x, y, z, _t, _t0);
      else         _Transformation2->Transform(x, y, z, _t, _t0);
    );
  }

  // ---------------------------------------------------------------------------
  /// Transform indices using additive composition of dense displacement fields
  void TransformUsingSumOfDisplacements(const blocked_range3d<int> &re) const
  {
    _TRANSFORM2(DisplacementField, _Displacement1, DisplacementField, _Displacement2,
      x = i, y = j, z = k;
      _Output->front()->ImageToWorld(x, y, z);
      x = *_Displacement1_x + *_Displacement2_x;
      y = *_Displacement1_y + *_Displacement2_y;
      z = *_Displacement1_z + *_Displacement2_z;
    );
  }

  // ---------------------------------------------------------------------------
  /// Transform indices using lookup table of world coordinates and additive
  /// composition of dense displacement fields
  void TransformUsingWorldCoordinatesAndSumOfDisplacements(const blocked_range3d<int> &re) const
  {
    _TRANSFORM3(Image2WorldMap, _Image2World, DisplacementField, _Displacement1, DisplacementField, _Displacement2,
      x = *_Image2World_x + *_Displacement1_x + *_Displacement2_x;
      y = *_Image2World_y + *_Displacement1_y + *_Displacement2_y;
      z = *_Image2World_z + *_Displacement1_z + *_Displacement2_z;
    );
  }

  // -------------------------------------------------------------------------
  /// Update output of scalar type possibly differing from OutputVoxelType1 template
  /// argument given a pre-computed map from output voxel indices to input indices
  template <class OutputVoxelType>
  void Update(int n, const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  _X - (ei - bi);
    const int s2 = (_Y - (ej - bj)) * _X - (ei - bi);

    irtkInterpolateImageFunction    *&input  = (*_Interpolator)[n];
    const Image2WorldMap::VoxelType  *x      = _Output2Input->GetPointerToVoxels(bi, bj, bk, 0);
    const Image2WorldMap::VoxelType  *y      = _Output2Input->GetPointerToVoxels(bi, bj, bk, 1);
    const Image2WorldMap::VoxelType  *z      = _Output2Input->GetPointerToVoxels(bi, bj, bk, 2);
    irtkImage                       *&output = (*_Output)[n];
    const int                         T      = (*_Output)[n]->GetT();

    // Update double output without rescaling
    if (typeid(OutputVoxelType) == typeid(double) && _ScaleFactor == 1.0 && _Offset == .0) {
      double *o;
      for (int k = bk; k != ek; k++) {
        for (int j = bj; j != ej;  j++) {
          for (int i = bi; i != ei; i++) {
            if (IsNaN(*x)) {
              o = reinterpret_cast<double *>(output->GetScalarPointer(i, j, k));
              for (int l = 0; l < T; l++, o += _NumberOfVoxels) *o = _InputPaddingValue[n];
            } else {
              reinterpret_cast<InterpolateImageFunction *>(input)->EvaluateInside(
                  reinterpret_cast<double *>(output->GetScalarPointer(i, j, k)), *x, *y, *z, _NumberOfVoxels);
            }
            ++x, ++y, ++z;
          }
          x += s1, y += s1, z += s1;
        }
        x += s2, y += s2, z += s2;
      }
    // Update output with rescaling
    } else {
      OutputVoxelType *o;
      double *v = new double[_MaxNumberOfComponents];
      for (int k = bk; k != ek; k++) {
        for (int j = bj; j != ej;  j++) {
          for (int i = bi; i != ei; i++) {
            if (IsNaN(*x)) {
              o = reinterpret_cast<OutputVoxelType *>(output->GetScalarPointer(i, j, k));
              for (int l = 0; l < T; l++, o += _NumberOfVoxels) {
                *o = static_cast<OutputVoxelType>(_InputPaddingValue[n]);
              }
            } else {
              reinterpret_cast<InterpolateImageFunction *>(input)->EvaluateInside(v, *x, *y, *z);
              o = reinterpret_cast<OutputVoxelType *>(output->GetScalarPointer(i, j, k));
              for (int l = 0; l < T; l++, o += _NumberOfVoxels) {
                *o = static_cast<OutputVoxelType>(_ScaleFactor * v[l] + _Offset);
              }
            }
            ++x, ++y, ++z;
          }
          x += s1, y += s1, z += s1;
        }
        x += s2, y += s2, z += s2;
      }
      delete[] v;
    }
  }

  // -------------------------------------------------------------------------
  /// Update specified output of unknown scalar type given a pre-computed map
  /// from output voxel indices to input indices
  void DefaultUpdate(int n, const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  _X - (ei - bi);
    const int s2 = (_Y - (ej - bj)) * _X - (ei - bi);

    irtkInterpolateImageFunction    *&input  = (*_Interpolator)[n];
    const Image2WorldMap::VoxelType  *x      = _Output2Input->GetPointerToVoxels(bi, bj, bk, 0);
    const Image2WorldMap::VoxelType  *y      = _Output2Input->GetPointerToVoxels(bi, bj, bk, 1);
    const Image2WorldMap::VoxelType  *z      = _Output2Input->GetPointerToVoxels(bi, bj, bk, 2);
    irtkImage                       *&output = (*_Output)[n];
    const int                         T      = (*_Output)[n]->GetT();

    double *v = new double[_MaxNumberOfComponents];
    for (int k = bk; k != ek; k++) {
      for (int j = bj; j != ej;  j++) {
        for (int i = bi; i != ei; i++) {
          if (IsNaN(*x)) {
            for (int l = 0; l < T; l++) output->PutAsDouble(i, j, k, l, _InputPaddingValue[n]);
          } else {
            reinterpret_cast<InterpolateImageFunction *>(input)->EvaluateInside(v, *x, *y, *z);
            for (int l = 0; l < T; l++) {
              output->PutAsDouble(i, j, k, l, _ScaleFactor * v[l] + _Offset);
            }
          }
          ++x, ++y, ++z;
        }
        x += s1, y += s1, z += s1;
      }
      x += s2, y += s2, z += s2;
    }
    delete[] v;
  }

  // ---------------------------------------------------------------------------
  /// Update foreground mask given a pre-computed map from output voxel indices to input indices
  void UpdateMask(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  _X - (ei - bi);
    const int s2 = (_Y - (ej - bj)) * _X - (ei - bi);

    const Image2WorldMap::VoxelType *x = _Output2Input->GetPointerToVoxels(bi, bj, bk);
    for (int k = bk; k != ek; k++) {
      for (int j = bj; j != ej;  j++) {
        for (int i = bi; i != ei; i++) {
          _Mask->Put(i, j, k, static_cast<irtkBinaryPixel>(IsNaN(*x++)));
        }
        x += s1;
      }
      x += s2;
    }
  }

  // ===========================================================================
  // One step update without explicit output to input index map
  // (all output images must have scalar type OutputVoxelType1)
  // ===========================================================================

  // ---------------------------------------------------------------------------
  /// Set components of all outputs to outside value
  void PutOutsideValue(int i, int j, int k) const
  {
    // Set output voxel to outside value
    OutputVoxelType1 *po;
    for (int n = _Begin; n != _End; n++) {
      irtkImage *&output = (*_Output)[n];
      po = reinterpret_cast<OutputVoxelType1 *>(output->GetScalarPointer(i, j, k));
      for (int l = 0; l < output->GetT(); l++, po += _NumberOfVoxels) {
        *po = static_cast<OutputVoxelType1>(_InputPaddingValue[n]);
      }
    }
    // Update foreground mask
    if (_Mask) _Mask->Put(i, j, k, static_cast<irtkBinaryPixel>(false));
  }

  // ---------------------------------------------------------------------------
  /// Update components of all outputs at a given voxel known the index mapping
  /// \attention Use only when OutputVoxelType1 is double!
  void UpdateDoubleWithoutRescaling(int    i, int    j, int    k,       // output voxel indices
                                    double x, double y, double z) const // input  voxel indices
  {
    // Set output voxel to interpolated value
    for (int n = _Begin; n != _End; n++) {
      irtkImage *&output = (*_Output)[n];
      reinterpret_cast<InterpolateImageFunction *>((*_Interpolator)[n])->EvaluateInside(
          reinterpret_cast<double *>(output->GetScalarPointer(i, j, k)), x, y, z, _NumberOfVoxels);
    }
    // Update foreground mask
    if (_Mask) _Mask->Put(i, j, k, static_cast<irtkBinaryPixel>(true));
  }

  // ---------------------------------------------------------------------------
  /// Update components of all outputs at a given voxel known the index mapping
  void UpdateWithoutRescaling(int    i, int    j, int    k, // output voxel indices
                              double x, double y, double z, // input voxel indices
                              double *v) const              // pre-allocated memory for interpolation
  {
    // Set output voxel to interpolated value
    OutputVoxelType1 *po;
    for (int n = _Begin; n != _End; n++) {
      irtkImage *&output = (*_Output)[n];
      reinterpret_cast<InterpolateImageFunction *>((*_Interpolator)[n])->EvaluateInside(v, x, y, z);
      po = reinterpret_cast<OutputVoxelType1 *>(output->GetScalarPointer(i, j, k));
      for (int l = 0; l < output->GetT(); l++, po += _NumberOfVoxels) {
        *po = static_cast<OutputVoxelType1>(v[l]);
      }
    }
    // Update foreground mask
    if (_Mask) _Mask->Put(i, j, k, static_cast<irtkBinaryPixel>(true));
  }

  // ---------------------------------------------------------------------------
  /// Update components of all outputs at a given voxel known the index mapping
  /// with rescaling of the interpolated value using _ScaleFactor and _Offset
  void UpdateWithRescaling(int    i, int    j, int    k, // output voxel indices
                           double x, double y, double z, // input voxel indices
                           double *v) const              // pre-allocated memory for interpolation
  {
    // Set output voxel to rescaled interpolated value
    OutputVoxelType1 *po;
    for (int n = _Begin; n != _End; n++) {
      irtkImage *&output = (*_Output)[n];
      reinterpret_cast<InterpolateImageFunction *>((*_Interpolator)[n])->EvaluateInside(v, x, y, z);
      po = reinterpret_cast<OutputVoxelType1 *>(output->GetScalarPointer(i, j, k));
      for (int l = 0; l < output->GetT(); l++, po += _NumberOfVoxels) {
        *po = static_cast<OutputVoxelType1>(_ScaleFactor * v[l] + _Offset);
      }
    }
    // Update foreground mask
    if (_Mask) _Mask->Put(i, j, k, static_cast<irtkBinaryPixel>(true));
  }

  // ---------------------------------------------------------------------------
  /// Update using lookup table of world coordinates
  void UpdateUsingWorldCoordinates(const blocked_range3d<int> &re) const
  {
    _UPDATE1(Image2WorldMap, _Image2World,
      x = *_Image2World_x;
      y = *_Image2World_y;
      z = *_Image2World_z;
    );
  }

  // ---------------------------------------------------------------------------
  /// Update using single transformation
  void UpdateUsingTransformation(const blocked_range3d<int> &re) const
  {
    _UPDATE0(
      x = i, y = j, z = k;
      _Output->front()->ImageToWorld(x, y, z);
      if (_Invert) _Transformation1->Inverse  (x, y, z, _t, _t0);
      else         _Transformation1->Transform(x, y, z, _t, _t0);
    );
  }

  // ---------------------------------------------------------------------------
  /// Update using lookup table of world coordinates and single transformation
  void UpdateUsingWorldCoordinatesAndTransformation(const blocked_range3d<int> &re) const
  {
    _UPDATE1(Image2WorldMap, _Image2World,
      x = *_Image2World_x;
      y = *_Image2World_y;
      z = *_Image2World_z;
      if (_Invert) _Transformation1->Inverse  (x, y, z, _t, _t0);
      else         _Transformation1->Transform(x, y, z, _t, _t0);
    );
  }

  // ---------------------------------------------------------------------------
  /// Update using single dense displacement field
  void UpdateUsingDisplacement(const blocked_range3d<int> &re) const
  {
    _UPDATE1(DisplacementField, _Displacement1,
      x = i, y = j, z = k;
      _Output->front()->ImageToWorld(x, y, z);
      x += *_Displacement1_x;
      y += *_Displacement1_y;
      z += *_Displacement1_z;
    );
  }

  // ---------------------------------------------------------------------------
  /// Update using lookup table of world coordinates and single dense displacement field
  void UpdateUsingWorldCoordinatesAndDisplacement(const blocked_range3d<int> &re) const
  {
    _UPDATE2(Image2WorldMap, _Image2World, DisplacementField, _Displacement1,
      x = *_Image2World_x + *_Displacement1_x;
      y = *_Image2World_y + *_Displacement1_y;
      z = *_Image2World_z + *_Displacement1_z;
    );
  }

  // ---------------------------------------------------------------------------
  /// Update using fluid composition of transformations
  void UpdateUsingCompositionOfTransformations(const blocked_range3d<int> &re) const
  {
    _UPDATE0(
      x = i, y = j, z = k;
      _Output->front()->ImageToWorld(x, y, z);
      if (_Invert) {
        _Transformation1->Inverse  (x, y, z, _t, _t0);
        _Transformation2->Inverse  (x, y, z, _t, _t0);
      } else {
        _Transformation1->Transform(x, y, z, _t, _t0);
        _Transformation2->Transform(x, y, z, _t, _t0);
      }
    );
  }

  // ---------------------------------------------------------------------------
  /// Update using lookup table of world coordinates and fluid composition of transformations
  void UpdateUsingWorldCoordinatesAndCompositionOfTransformations(const blocked_range3d<int> &re) const
  {
    _UPDATE1(Image2WorldMap, _Image2World,
      x = *_Image2World_x;
      y = *_Image2World_y;
      z = *_Image2World_z;
      if (_Invert) {
        _Transformation1->Inverse  (x, y, z, _t, _t0);
        _Transformation2->Inverse  (x, y, z, _t, _t0);
      } else {
        _Transformation1->Transform(x, y, z, _t, _t0);
        _Transformation2->Transform(x, y, z, _t, _t0);
      }
    );
  }

  // ---------------------------------------------------------------------------
  /// Update using fluid composition of dense displacement field and transformation
  void UpdateUsingCompositionOfDisplacementAndTransformation(const blocked_range3d<int> &re) const
  {
    _UPDATE1(DisplacementField, _Displacement1,
      x = i, y = j, z = k;
      _Output->front()->ImageToWorld(x, y, z);
      x += *_Displacement1_x;
      y += *_Displacement1_y;
      z += *_Displacement1_z;
      if (_Invert) _Transformation2->Inverse  (x, y, z, _t, _t0);
      else         _Transformation2->Transform(x, y, z, _t, _t0);
    );
  }

  // ---------------------------------------------------------------------------
  /// Update using lookup table of world coordinates and fluid composition of
  /// dense displacement field and transformation
  void UpdateUsingWorldCoordinatesAndCompositionOfDisplacementAndTransformation(const blocked_range3d<int> &re) const
  {
    _UPDATE2(Image2WorldMap, _Image2World, DisplacementField, _Displacement1,
      x = *_Image2World_x + *_Displacement1_x;
      y = *_Image2World_y + *_Displacement1_y;
      z = *_Image2World_z + *_Displacement1_z;
      if (_Invert) _Transformation2->Inverse  (x, y, z, _t, _t0);
      else         _Transformation2->Transform(x, y, z, _t, _t0);
    );
  }

  // ---------------------------------------------------------------------------
  /// Update using additive composition of dense displacement fields
  void UpdateUsingSumOfDisplacements(const blocked_range3d<int> &re) const
  {
    _UPDATE2(DisplacementField, _Displacement1, DisplacementField, _Displacement2,
      x = i, y = j, z = k;
      _Output->front()->ImageToWorld(x, y, z);
      x += *_Displacement1_x + *_Displacement2_x;
      y += *_Displacement1_y + *_Displacement2_y;
      z += *_Displacement1_z + *_Displacement2_z;
    );
  }

  // ---------------------------------------------------------------------------
  /// Update using lookup table of world coordinates and additive composition of
  /// dense displacement fields
  void UpdateUsingWorldCoordinatesAndSumOfDisplacements(const blocked_range3d<int> &re) const
  {
    _UPDATE3(Image2WorldMap, _Image2World, DisplacementField, _Displacement1, DisplacementField, _Displacement2,
      x = *_Image2World_x + *_Displacement1_x + *_Displacement2_x;
      y = *_Image2World_y + *_Displacement1_y + *_Displacement2_y;
      z = *_Image2World_z + *_Displacement1_z + *_Displacement2_z;
    );
  }

  // ===========================================================================
  // Entry methods
  // ===========================================================================

public:

  // ---------------------------------------------------------------------------
  /// Update images within given region-of-interest
  /// \note Called for each thread by parallel_for. Do not call directly.
  void operator() (const blocked_range3d<int> &re) const
  {
    // If scalar type is not the same for all output images...
    if (_UniqueScalarType == IRTK_VOXEL_UNKNOWN) {
      // 1. Pre-compute map from output to input voxel indices
      if (_Image2World) {
        // Two displacement maps                         => additive composition
        if      (_Displacement1   && _Displacement2)     TransformUsingWorldCoordinatesAndSumOfDisplacements(re);
        // Two transformations                           => fluid composition
        else if (_Transformation1 && _Transformation2)   TransformUsingWorldCoordinatesAndCompositionOfTransformations(re);
        // Single transformation and/or displacement map => prefer displacement map
        else if (_Displacement1)                         TransformUsingWorldCoordinatesAndDisplacement(re);
        else if (_Transformation1)                       TransformUsingWorldCoordinatesAndTransformation(re);
        // Otherwise                                     => world coordinates map
        else                                             TransformUsingWorldCoordinates(re);
      } else {
        // Two displacement maps                         => additive composition
        if      (_Displacement1   && _Displacement2)     TransformUsingSumOfDisplacements(re);
        // Two transformations                           => fluid composition
        else if (_Transformation1 && _Transformation2)   TransformUsingCompositionOfTransformations(re);
        // Single transformation and/or displacement map => prefer displacement map
        else if (_Displacement1)                         TransformUsingDisplacement(re);
        else if (_Transformation1)                       TransformUsingTransformation(re);
        // Otherwise                                     => nothing to do
      }
      // 2. Update each output separately using the pre-computed indices mapping
      for (int n = _Begin; n < _End; n++) {
        switch ((*_Output)[n]->GetScalarType()) {
          case IRTK_VOXEL_CHAR:           { Update<char          >(n, re); break; }
          case IRTK_VOXEL_UNSIGNED_CHAR:  { Update<unsigned char >(n, re); break; }
          case IRTK_VOXEL_SHORT:          { Update<short         >(n, re); break; }
          case IRTK_VOXEL_UNSIGNED_SHORT: { Update<unsigned short>(n, re); break; }
          case IRTK_VOXEL_INT:            { Update<int           >(n, re); break; }
          case IRTK_VOXEL_UNSIGNED_INT:   { Update<unsigned int  >(n, re); break; }
          case IRTK_VOXEL_FLOAT:          { Update<float         >(n, re); break; }
          case IRTK_VOXEL_DOUBLE:         { Update<double        >(n, re); break; }
          default:                        { DefaultUpdate         (n, re); break; }
        }
      }
      // 3. Update mask using the pre-computed indices mapping
      if (_Mask) UpdateMask(re);
    // Otherwise, perform both steps and update all outputs at the same time to save the
    // additional memory required to store the map from input to output indices
    } else {
      if (_Image2World) {
        // Two displacement maps                         => additive composition
        if      (_Displacement1   && _Displacement2)     UpdateUsingWorldCoordinatesAndSumOfDisplacements(re);
        // Two transformations                           => fluid composition
        else if (_Transformation1 && _Transformation2)   UpdateUsingWorldCoordinatesAndCompositionOfTransformations(re);
        // Single transformation and/or displacement map => prefer displacement map
        else if (_Displacement1)                         UpdateUsingWorldCoordinatesAndDisplacement(re);
        else if (_Transformation1)                       UpdateUsingWorldCoordinatesAndTransformation(re);
        // Otherwise                                     => world coordinates map
        else                                             UpdateUsingWorldCoordinates(re);
      } else {
        // Two displacement maps                         => additive composition
        if      (_Displacement1   && _Displacement2)     UpdateUsingSumOfDisplacements(re);
        // Two transformations                           => fluid composition
        else if (_Transformation1 && _Transformation2)   UpdateUsingCompositionOfTransformations(re);
        // Single transformation and/or displacement map => prefer displacement map
        else if (_Displacement1)                         UpdateUsingDisplacement(re);
        else if (_Transformation1)                       UpdateUsingTransformation(re);
        // Otherwise                                     => nothing to do
      }
    }
  }

  // -------------------------------------------------------------------------
  /// Update images
  void operator() ()
  {
    // Temporal coordinates of target and source images
    _t  = _Input ->front()->ImageToTime(.0);
    _t0 = _Output->front()->ImageToTime(.0);
    // Size of transformed images
    _X = _Output->front()->GetX();
    _Y = _Output->front()->GetY();
    const int Z = _Output->front()->GetZ();
    // Check size of region-of-interest mask
    if (_RegionOfInterest && (_RegionOfInterest->GetX() != _X ||
                              _RegionOfInterest->GetY() != _Y ||
                              _RegionOfInterest->GetZ() !=  Z)) {
      cerr << "irtkMultipleImageTransformation::Run/Update: Region-of-interest mask must have same size as output image(s)" << endl;
      exit(1);
    }
    // Obtain required information about output images
    _NumberOfVoxels        = _X * _Y * Z;
    _MaxNumberOfComponents = 1;
    _UniqueScalarType      = voxel_info<OutputVoxelType1>::type();
    for (int n = _Begin; n < _End; n++) {
      irtkImage *output = (*_Output)[n];
      // Maximum number of voxel components (i.e., image frames)
      _MaxNumberOfComponents = max(_MaxNumberOfComponents, output->GetT());
      // Determine if scalar type differs from OutputVoxelType1 template argument
      if (output->GetScalarType() != _UniqueScalarType) _UniqueScalarType = IRTK_VOXEL_UNKNOWN;
    }
    // Allocate memory for indices map if required
    if (_UniqueScalarType == IRTK_VOXEL_UNKNOWN) {
      if (!_Output2Input) {
        _Output2Input      = new Image2WorldMap();
        _Output2InputOwner = true;
      }
      _Output2Input->Initialize(_X, _Y, Z, 3);
    }
    // If only one transformation/displacement given, always use first slot
    if (!_Transformation1 && !_Displacement1 && (_Transformation2 || _Displacement2)) {
      _Transformation1 = _Transformation2;
      _Displacement1   = _Displacement2;
      _Transformation2 = NULL;
      _Displacement2   = NULL;
    }
    // Update transformed images
    blocked_range3d<int> voxels(0, Z, 0, _Y, 0, _X);
    parallel_for(voxels, *this);
    // Clean up
    if (_Output2InputOwner) {
      delete _Output2Input;
      _Output2Input      = NULL;
      _Output2InputOwner = false;
    }
  }

}; // irtkMultipleImageTransformationUpdateBody

// =============================================================================
// Transformation
// =============================================================================

// -----------------------------------------------------------------------------
template <class InterpolateImageFunction, class VoxelType>
inline void irtkMultipleImageTransformation::RunWithInterpolator(int b, int e, const irtkBinaryImage *roi)
{
  // Instantiate (multi-threaded) method body
  irtkMultipleImageTransformationUpdateBody<InterpolateImageFunction, VoxelType> body;
  body._Input              = &_Input;
  body._Interpolator       = &_Interpolator;
  body._ScaleFactor        = _ScaleFactor;
  body._Offset             = _Offset;
  body._InputPaddingValue  = _InputPaddingValue;
  body._Output             = &_Output;
  body._RegionOfInterest   = roi ? roi : _RegionOfInterest;
  body._Mask               = _Mask;
  body._CheckIfForeground  = _CheckIfForeground;
  body._Begin              = b;
  body._End                = e;
  // Update/Initialize cached data if required
  if (_Transformation) {
    // Fixed and changing transformation given
    if (_FixedTransformation) {
      // Fluid composition of transformations
      if (_Composition == COMPOSITION_FLUID) {
        if (_Invert) {
          delete _FixedDisplacement;
          _FixedDisplacement = NULL;
          ComputeWorldCoordinates(true);
          ComputeDisplacement();
        } else {
          ApplyFixedTransformation();
          ComputeDisplacement();
        }
      // Additive composition of transformations
      } else {
        ComputeWorldCoordinates(true);
        ComputeFixedDisplacement();
        ComputeDisplacement();
      }
      if (_Invert) {
        if (_Composition == COMPOSITION_FLUID && _CacheDisplacements) {
          body._Displacement1   = _Displacement;
        } else {
          body._Transformation1 = _Transformation;
          body._Displacement1   = _Displacement;
          body._Transformation2 = _FixedTransformation;
          body._Displacement2   = _FixedDisplacement;
        }
      } else {
        body._Transformation1 = _FixedTransformation;
        body._Displacement1   = _FixedDisplacement;
        body._Transformation2 = _Transformation;
        body._Displacement2   = _Displacement;
      }
    // Only changing transformation given
    } else {
      ComputeWorldCoordinates(_FixedTransformationApplied);
      ComputeDisplacement();
      body._Transformation1 = _Transformation;
      body._Displacement1   = _Displacement;
    }
  // Only fixed transformation given
  } else if (_FixedTransformation) {
    if (_CacheWorldCoordinates) ApplyFixedTransformation();
    else                        ComputeFixedDisplacement();
    body._Transformation1 = _FixedTransformation;
    body._Displacement1   = _FixedDisplacement;
  // No transformation given and output has same attributes as input
  } else {
    if (_Output[0]->HasSpatialAttributesOf(_Input[0])) {
      for (int n = body._Begin; n < body._End; n++) (*_Output[n]) = (*_Input[n]);
      // FIXME: Only the outputs in the range [_Begin, _End) have been updated!
      _Update = false;
      return; // no transformation needed
    }
    ComputeWorldCoordinates(true);
  }
  // Set world coordinates map - after it was optionally computed above!
  body._Image2World = _Image2World;
  // Update outputs
  body();
  // Further updates only required if changing transformation set
  // FIXME: Only the outputs in the range [_Begin, _End) have been updated!
  _Update = (_Transformation != NULL);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void irtkMultipleImageTransformation::RunWithScalarType(int b, int e, const irtkBinaryImage *roi)
{
  typedef irtkGenericImage<VoxelType> ImageType;

  const char *name = _Interpolator[0]->NameOfClass() + 4;
  const char *p    = NULL;

  // Linear interpolator
  if ((p = strstr(name, "FastLinear")) != NULL) {
    if (p[34/*=strlen("FastLinearInterpolateImageFunction")*/] == '2') {
      RunWithInterpolator<irtkGenericLinearInterpolateImageFunction2D<ImageType>, VoxelType>(b, e, roi);
    } else {
      RunWithInterpolator<irtkGenericLinearInterpolateImageFunction3D<ImageType>, VoxelType>(b, e, roi);
    }
  // Unknown interpolator
  } else {
    RunWithInterpolator<irtkGenericInterpolateImageFunction<ImageType>, VoxelType>(b, e, roi);
  }
}

// -----------------------------------------------------------------------------
void irtkMultipleImageTransformation::Run(int i, int n, const irtkBinaryImage *roi)
{
  if      (i < 0) i = 0, n = NumberOfImages();
  else if (n < 0)        n = NumberOfImages();
  if (i >= n) return;
  // Update transformed images
  switch (_UniqueScalarType) {
    // Determine type of interpolator to enable compiler to generate faster code
    case IRTK_VOXEL_CHAR:           { RunWithScalarType<char          >(i, n, roi); break; }
    case IRTK_VOXEL_UNSIGNED_CHAR:  { RunWithScalarType<unsigned char >(i, n, roi); break; }
    case IRTK_VOXEL_SHORT:          { RunWithScalarType<short         >(i, n, roi); break; }
    case IRTK_VOXEL_UNSIGNED_SHORT: { RunWithScalarType<unsigned short>(i, n, roi); break; }
    case IRTK_VOXEL_INT:            { RunWithScalarType<int           >(i, n, roi); break; }
    case IRTK_VOXEL_UNSIGNED_INT:   { RunWithScalarType<unsigned int  >(i, n, roi); break; }
    case IRTK_VOXEL_FLOAT:          { RunWithScalarType<float         >(i, n, roi); break; }
    case IRTK_VOXEL_DOUBLE:         { RunWithScalarType<double        >(i, n, roi); break; }
    // No unique voxel type
    default: RunWithInterpolator<irtkInterpolateImageFunction, double>(i, n, roi);
  }
  // Update foreground mask of outputs
  for (int j = i; j < n; j++) {
    irtkBinaryImage *mask = _Output[j]->GetMask();
    if (mask && mask != _Mask) *mask = *_Mask;
  }
}

// -----------------------------------------------------------------------------
void irtkMultipleImageTransformation::Update(int i, int n, const irtkBinaryImage *roi)
{
  if (_Update) Run(i, n, roi);
}

// =============================================================================
// Deprecated - for compatibility with irtkImageTransformation
// =============================================================================

// -----------------------------------------------------------------------------
void irtkMultipleImageTransformation::PutTargetPaddingValue(double value)
{
  _OutputPaddingValue = value;
}

// -----------------------------------------------------------------------------
double irtkMultipleImageTransformation::GetTargetPaddingValue() const
{
  return _OutputPaddingValue;
}

// -----------------------------------------------------------------------------
void irtkMultipleImageTransformation::PutSourcePaddingValue(double value)
{
  _DefaultPaddingValue = value;
}

// -----------------------------------------------------------------------------
double irtkMultipleImageTransformation::GetSourcePaddingValue() const
{
  return _DefaultPaddingValue;
}

// -----------------------------------------------------------------------------
void irtkMultipleImageTransformation::PutInterpolator(irtkImageFunction *f)
{
  irtkInterpolateImageFunction *interpolator = dynamic_cast<irtkInterpolateImageFunction *>(f);
  if (interpolator) {
    _Interpolation      = interpolator->InterpolationMode();
    _UnusedInterpolator = interpolator;
  } else {
    cerr << "irtkMultipleImageTransformation::PutInterpolator: Expected irtkInterpolateImageFunction object! Deprecated method, use SetInterpolation instead." << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
irtkImageFunction *irtkMultipleImageTransformation::GetInterpolator() const
{
  return _UnusedInterpolator;
}

// -----------------------------------------------------------------------------
void irtkMultipleImageTransformation::PutScaleFactorAndOffset(double scale, double offset)
{
  SetScaleFactorAndOffset(scale, offset);
}

// -----------------------------------------------------------------------------
void irtkMultipleImageTransformation::PutInputTimeOffset(double)
{
  // TODO: Implement irtkMultipleImageTransformation::PutInputTimeOffset
  cerr << "irtkMultipleImageTransformation::PutInputTimeOffset: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
void irtkMultipleImageTransformation::PutOutputTimeOffset(double)
{
  // TODO: Implement irtkMultipleImageTransformation::PutOutputTimeOffset
  cerr << "irtkMultipleImageTransformation::PutOutputTimeOffset: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
void irtkMultipleImageTransformation::TwoDOn()
{
  // TODO: Implement irtkMultipleImageTransformation::TwoDOn
  cerr << "irtkMultipleImageTransformation::TwoDOn: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
void irtkMultipleImageTransformation::TwoDOff()
{
  // TODO: irtkMultipleImageTransformation::TwoDOff
  cerr << "irtkMultipleImageTransformation::TwoDOff: Not implemented" << endl;
  exit(1);
}
