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

#include <irtkTransformation.h>

// -----------------------------------------------------------------------------
irtkImageTransformation::irtkImageTransformation()
{
  // Set input and output
  _input  = NULL;
  _output = NULL;

  // Set transformation
  _transformation     = NULL;
  _Cache              = NULL;
  _CacheOwner         = false;
  _DisplacementField  = NULL;
  _CacheInterpolation = Interpolation_Linear;
  _CacheExtrapolation = Extrapolation_Const;

  // Set interpolator
  _interpolator = NULL;

  // Set padding value
  _TargetPaddingValue = voxel_limits<double>::min();

  // Set padding value
  _SourcePaddingValue = 0;

  // Scale factor
  _ScaleFactor = 1;
  
  // Offset
  _Offset = 0;

  // Temporal offsets
  _OutputTimeOffset = .0;
  _InputTimeOffset  = .0;

  // Set invert mode
  _Invert = false;

  // Set 2D mode
  _2D = false;

  _NumberOfSingularPoints = 0;
}

// -----------------------------------------------------------------------------
irtkImageTransformation::~irtkImageTransformation()
{
  // Set input and output
  _input  = NULL;
  _output = NULL;

  // Set transformation
  _transformation = NULL;
  if (_CacheOwner) delete _Cache;
  _Cache          = NULL;
  Delete(_DisplacementField);

  // Set interpolator
  _interpolator = NULL;

  // Set padding value
  _TargetPaddingValue = -std::numeric_limits<double>::max();

  // Set padding value
  _SourcePaddingValue = 0;

  // Set invert mode
  _Invert = false;
}

// -----------------------------------------------------------------------------
irtkImageTransformation *irtkImageTransformation::New(irtkTransformation *transformation)
{
  irtkImageTransformation  *imagetransformation = NULL;
  if (strcmp(transformation->NameOfClass(), "irtkHomogeneousTransformation") == 0) {
    imagetransformation = new irtkImageHomogeneousTransformation;
  } else if (strcmp(transformation->NameOfClass(), "irtkRigidTransformation") == 0) {
    imagetransformation = new irtkImageHomogeneousTransformation;
  } else if (strcmp(transformation->NameOfClass(), "irtkAffineTransformation") == 0) {
    imagetransformation = new irtkImageHomogeneousTransformation;
  }
  if (!imagetransformation) imagetransformation = new irtkImageTransformation;
  imagetransformation->SetTransformation(transformation);
  return imagetransformation;
}

// -----------------------------------------------------------------------------
void irtkImageTransformation::SetInput(irtkImage *image)
{
  if (image != NULL) {
    _input = image;
  } else {
    cerr << "irtkImageTransformation::SetInput: Input is NULL\n";
    exit(1);
  }
}

// -----------------------------------------------------------------------------
void irtkImageTransformation::SetInput(irtkImage *image, irtkTransformation *transformation)
{
  if ((image != NULL) && (transformation != NULL)) {
    _input          = image;
    if (_transformation != transformation) {
      _transformation = transformation;
      if (_Cache) _Cache->Modified(true);
    }
  } else {
    cerr << "irtkImageTransformation::SetInput: Input is NULL\n";
    exit(1);
  }
}

// -----------------------------------------------------------------------------
void irtkImageTransformation::SetOutput(irtkImage *image)
{
  if (image != NULL) {
    _output = image;
  } else {
    cerr << "irtkImageTransformation::SetOutput: Output is NULL\n";
    exit(1);
  }
}

// -----------------------------------------------------------------------------
void irtkImageTransformation::SetTransformation(irtkTransformation *transformation)
{
  if (transformation != NULL) {
    if (_transformation != transformation) {
      _transformation = transformation;
      if (_Cache) _Cache->Modified(true);
    }
  } else {
    cerr << "irtkImageTransformation::SetInput: Transformation is NULL\n";
    exit(1);
  }
}

// -----------------------------------------------------------------------------
void irtkImageTransformation::SetCache(irtkImageTransformationCache *cache)
{
  if (_Cache != cache) {
    if (_CacheOwner) delete _Cache;
    Delete(_DisplacementField);
    _Cache      = cache;
    _CacheOwner = false;
  }
}

// -----------------------------------------------------------------------------
void irtkImageTransformation::PutTargetPaddingValue(double PaddingValue)
{
  _TargetPaddingValue = PaddingValue;
}

// -----------------------------------------------------------------------------
double irtkImageTransformation::GetTargetPaddingValue()
{
  return _TargetPaddingValue;
}

// -----------------------------------------------------------------------------
void irtkImageTransformation::PutSourcePaddingValue(double PaddingValue)
{
  _SourcePaddingValue = PaddingValue;
}

// -----------------------------------------------------------------------------
double irtkImageTransformation::GetSourcePaddingValue()
{
  return _SourcePaddingValue;
}

// -----------------------------------------------------------------------------
irtkImageFunction *irtkImageTransformation::GetInterpolator()
{
  return _interpolator;
}

// -----------------------------------------------------------------------------
void irtkImageTransformation::PutInterpolator(irtkImageFunction *interpolator)
{
  _interpolator = interpolator;
}

// -----------------------------------------------------------------------------
void irtkImageTransformation::PutScaleFactorAndOffset(double ScaleFactor, double Offset)
{
	_ScaleFactor = ScaleFactor;
	_Offset      = Offset;
}

// -----------------------------------------------------------------------------
void irtkImageTransformation::PutInputTimeOffset(double Offset)
{
  _InputTimeOffset = Offset;
}

// -----------------------------------------------------------------------------
void irtkImageTransformation::PutOutputTimeOffset(double Offset)
{
  _OutputTimeOffset = Offset;
}

// -----------------------------------------------------------------------------
void irtkImageTransformation::InvertOn()
{
  if (!_Invert) {
    if (_Cache) _Cache->Modified(true);
    _Invert = true;
  }
}

// -----------------------------------------------------------------------------
void irtkImageTransformation::InvertOff()
{
  if (_Invert) {
    if (_Cache) _Cache->Modified(true);
    _Invert = false;
  }
}

// -----------------------------------------------------------------------------
void irtkImageTransformation::TwoDOn()
{
  _2D = true;
}

// -----------------------------------------------------------------------------
void irtkImageTransformation::TwoDOff()
{
  _2D = false;
}

// -----------------------------------------------------------------------------
void irtkImageTransformation::Initialize()
{
  // Clean up previous run
  Delete(_DisplacementField);

  // Check inputs and outputs
  if (_input == NULL) {
    cerr << "irtkImageTransformation::Run: Filter has no input" << endl;
    exit(1);
  }
  if (_output == NULL) {
    cerr << "irtkImageTransformation::Run: Filter has no output" << endl;
    exit(1);
  }
  if (_transformation == NULL) {
    cerr << "irtkImageTransformation::Run: Filter has no transformation" << endl;
    exit(1);
  }
  if (_interpolator == NULL) {
    cerr << "irtkImageTransformation::Run: Filter has no interpolator" << endl;
    exit(1);
  }
  if (_input->IsEmpty() == true) {
    cerr << "irtkImageTransformation::Run: Input is empty" << endl;
    exit(1);
  }
  if (_input == _output) {
    cerr << "irtkImageTransformation::Run: Input equals output" << endl;
    exit(1);
  }

  _NumberOfSingularPoints = 0;

  if (_transformation->RequiresCachingOfDisplacements() && !_Cache) {
    _CacheOwner = true;
    _Cache      = new irtkImageTransformationCache();
    _Cache->Initialize(_output->GetImageAttributes(), 3);
  }

  if (_Cache && !_Cache->IsEmpty()) {

    // Compute and cache displacements
    if (_Cache->Modified()) {
      // TODO: Cache displacements for multiple time intervals.
      const double t0 = _Cache->GetTOrigin();
      const double t  = _input->GetTOrigin();

      _Cache->Initialize();
      if (_Invert) {
        _NumberOfSingularPoints = _transformation->InverseDisplacement(*_Cache, t, t0);
      } else {
        _transformation->Displacement(*_Cache, t, t0);
      }

      _Cache->Modified(false);
    }

    // Initialize cache interpolator
    _DisplacementField = irtkInterpolateImageFunction::New(_CacheInterpolation,
                                                           _CacheExtrapolation,
                                                           _Cache);
    _DisplacementField->SetInput(_Cache);
    _DisplacementField->Initialize();
  }

  // Initialize input interpolator
  _interpolator->SetInput(_input);
  _interpolator->Initialize();
}

// -----------------------------------------------------------------------------
class irtkImageTransformationRun
{
public:

  // ---------------------------------------------------------------------------
  // Data members
  irtkImage                    *_Input;
  irtkImageFunction            *_Interpolator;
  irtkTransformation           *_Transformation;
  irtkInterpolateImageFunction *_DisplacementField;
  irtkImage                    *_Output;
  bool                          _Invert;
  bool                          _2D;
  double                        _InputMin;
  double                        _InputMax;
  double                        _ScaleFactor;
  double                        _Offset;
  double                        _TargetPaddingValue;
  double                        _SourcePaddingValue;
  double                        _InputTimeOffset;
  double                        _OutputTimeOffset;
  mutable int                   _OutputFrame;
  mutable double                _OutputTime;
  mutable int                   _InputFrame;
  mutable double                _InputTime;
  int                           _NumberOfSingularPoints;

  // ---------------------------------------------------------------------------
  /// Default constructor
  irtkImageTransformationRun()
  :
    _Input             (NULL),
    _Interpolator      (NULL),
    _Transformation    (NULL),
    _DisplacementField (NULL),
    _Output            (NULL),
    _Invert            (false),
    _2D                (false),
    _InputMin          (.0),
    _InputMax          (-1),
    _ScaleFactor       (1),
    _Offset            (0),
    _TargetPaddingValue(-1),
    _SourcePaddingValue(-1),
    _InputTimeOffset   (0),
    _OutputTimeOffset  (0),
    _OutputFrame       (0),
    _OutputTime        (-1),
    _InputFrame        (0),
    _InputTime         (0),
    _NumberOfSingularPoints(0)
  {}

  // ---------------------------------------------------------------------------
  /// Copy constructor
  irtkImageTransformationRun(const irtkImageTransformationRun &other)
  :
    _Input             (other._Input),
    _Interpolator      (other._Interpolator),
    _Transformation    (other._Transformation),
    _DisplacementField (other._DisplacementField),
    _Output            (other._Output),
    _Invert            (other._Invert),
    _2D                (other._2D),
    _InputMin          (other._InputMin),
    _InputMax          (other._InputMax),
    _ScaleFactor       (other._ScaleFactor),
    _Offset            (other._Offset),
    _TargetPaddingValue(other._TargetPaddingValue),
    _SourcePaddingValue(other._SourcePaddingValue),
    _InputTimeOffset   (other._InputTimeOffset),
    _OutputTimeOffset  (other._OutputTimeOffset),
    _OutputFrame       (other._OutputFrame),
    _OutputTime        (other._OutputTime),
    _InputFrame        (other._InputFrame),
    _InputTime         (other._InputTime),
    _NumberOfSingularPoints(0)
  {}

  // ---------------------------------------------------------------------------
  /// Split constructor
  irtkImageTransformationRun(const irtkImageTransformationRun &other, split)
  :
    _Input             (other._Input),
    _Interpolator      (other._Interpolator),
    _Transformation    (other._Transformation),
    _DisplacementField (other._DisplacementField),
    _Output            (other._Output),
    _Invert            (other._Invert),
    _2D                (other._2D),
    _InputMin          (other._InputMin),
    _InputMax          (other._InputMax),
    _ScaleFactor       (other._ScaleFactor),
    _Offset            (other._Offset),
    _TargetPaddingValue(other._TargetPaddingValue),
    _SourcePaddingValue(other._SourcePaddingValue),
    _InputTimeOffset   (other._InputTimeOffset),
    _OutputTimeOffset  (other._OutputTimeOffset),
    _OutputFrame       (other._OutputFrame),
    _OutputTime        (other._OutputTime),
    _InputFrame        (other._InputFrame),
    _InputTime         (other._InputTime),
    _NumberOfSingularPoints(0)
  {}

  // ---------------------------------------------------------------------------
  void join(const irtkImageTransformationRun &other)
  {
    _NumberOfSingularPoints += other._NumberOfSingularPoints;
  }

  // ---------------------------------------------------------------------------
  double Clamp(double v) const
  {
    if (_InputMin <= _InputMax) {
      if      (v < _InputMin) return _InputMin;
      else if (v > _InputMax) return _InputMax;
    }
    return v;
  }

  // ---------------------------------------------------------------------------
  /// Applies the given transformation within given output region for frame _T
  void operator() (const blocked_range3d<int> &r)
  {
    double value, x, y, z, u, v, w, disp[3] = {.0, .0, .0};

    for (int k = r.pages().begin(); k != r.pages().end(); ++k)
    for (int j = r.rows ().begin(); j != r.rows ().end(); ++j)
    for (int i = r.cols ().begin(); i != r.cols ().end(); ++i) {
      if (_Output->GetAsDouble(i, j, k, _OutputFrame) > _TargetPaddingValue) {
        // Transform point into world coordinates
        x = i, y = j, z = k;
        _Output->ImageToWorld(x, y, z);
        // Transform point
        if (_DisplacementField) {
          u = x, v = y, w = z;
          _DisplacementField->WorldToImage  (u, v, w);
          _DisplacementField->Evaluate(disp, u, v, w);
          x += disp[0], y += disp[1], z += disp[2];
        } else {
          if (_Invert) {
            if (!_Transformation->Inverse(x, y, z, _InputTime, _OutputTime)) {
              ++_NumberOfSingularPoints;
            }
          } else {
            _Transformation->Transform(x, y, z, _InputTime, _OutputTime);
          }
        }
        // Transform point into image coordinates
        _Input->WorldToImage(x, y, z);
        // Check whether transformed point is in FOV of input
        if (-0.5 < x && x < static_cast<double>(_Input->GetX()) - 0.5 &&
            -0.5 < y && y < static_cast<double>(_Input->GetY()) - 0.5) {
          if (_2D) {
            value = _Interpolator->Evaluate(x, y, k, _InputFrame);
            value = _ScaleFactor * Clamp(value) + _Offset;
          } else if (-0.5 < z && z < static_cast<double>(_Input->GetZ()) - 0.5) {
            value = _Interpolator->Evaluate(x, y, z, _InputFrame);
            value = _ScaleFactor * Clamp(value) + _Offset;
          } else {
            value = _SourcePaddingValue;
          }
        } else {
          value = _SourcePaddingValue;
        }
      } else {
        value = _SourcePaddingValue;
      }
      _Output->PutAsDouble(i, j, k, _OutputFrame, value);
    }
  }

  // ---------------------------------------------------------------------------
  /// Applies the given transformation for the given range of output frames
  void operator() ()
  {
    _NumberOfSingularPoints = 0;
    for (_OutputFrame = 0; _OutputFrame < _Output->GetT(); ++_OutputFrame) {
      _OutputTime = _Output->ImageToTime(_OutputFrame);
      _InputFrame = ((_Input->GetT() > 1) ? static_cast<int>(round(_Input->TimeToImage(_OutputTime))) : 0);
      _InputTime  = _Input->ImageToTime(_InputFrame);

      _OutputTime += _OutputTimeOffset;
      _InputTime  += _InputTimeOffset;

      if (0 <= _InputFrame && _InputFrame < _Input->GetT()) {
        blocked_range3d<int> voxels(0, _Output->GetZ(),
                                    0, _Output->GetY(),
                                    0, _Output->GetX());
        irtkImageTransformationRun body(*this);
        parallel_reduce(voxels, body);
        _NumberOfSingularPoints += body._NumberOfSingularPoints;
      } else {
        for (int k = 0; k < _Output->GetZ(); ++k)
        for (int j = 0; j < _Output->GetY(); ++j)
        for (int i = 0; i < _Output->GetX(); ++i) {
          _Output->PutAsDouble(i, j, k, _OutputFrame, _SourcePaddingValue);
        }
      }
    }
  }

}; // irtkApplyImageTransformation

// -----------------------------------------------------------------------------
void irtkImageTransformation::Run()
{
  IRTK_START_TIMING();

  this->Initialize();

  irtkImageTransformationRun body;
  body._Input              = _input;
  body._Interpolator       = _interpolator;
  body._Transformation     = _transformation;
  body._DisplacementField  = _DisplacementField;
  body._Output             = _output;
  body._Invert             = _Invert;
  body._2D                 = _2D;
  body._ScaleFactor        = _ScaleFactor;
  body._Offset             = _Offset;
  body._TargetPaddingValue = _TargetPaddingValue;
  body._SourcePaddingValue = _SourcePaddingValue;
  body._InputTimeOffset    = _InputTimeOffset;
  body._OutputTimeOffset   = _OutputTimeOffset;

  // Range of output intensities should be the same as the one of the input
  _input->GetMinMaxAsDouble(body._InputMin, body._InputMax);

  body();

  if (_Invert && !_DisplacementField) {
    _NumberOfSingularPoints = body._NumberOfSingularPoints;
  }

  IRTK_DEBUG_TIMING(1, "irtkImageTransformation::Run");
}

