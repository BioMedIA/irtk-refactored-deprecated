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

#include <irtkBalloonForce.h>

#include <vtkMath.h>
#include <vtkImageData.h>
#include <vtkImageStencilData.h>
#include <vtkImageStencilIterator.h>
#include <vtkPolyDataToImageStencil.h>

#include <irtkEdgeTable.h>

using namespace irtk::polydata;


// =============================================================================
// Auxiliary functions
// =============================================================================

namespace irtkBalloonForceUtils {

typedef irtkGenericLinearInterpolateImageFunction<
           irtkGenericImage<irtkRegisteredImage::VoxelType>
         > ImageFunction;

// -----------------------------------------------------------------------------
/// Convert IRTK image to VTK image
///
/// Note that vtkImageData has no implicit orientation. Therefore we just
/// convert the image in voxel coordinates (origin at 0 and voxel size 1x1x1).
vtkSmartPointer<vtkImageData> ConvertImage(const irtkRegisteredImage *image)
{
  vtkSmartPointer<vtkImageData> imagedata = vtkSmartPointer<vtkImageData>::New();
  imagedata->SetOrigin(.0, .0, .0);
  imagedata->SetDimensions(image->X(), image->Y(), image->Z());
  imagedata->SetSpacing(1.0, 1.0, 1.0);
#if VTK_MAJOR_VERSION >= 6
  imagedata->AllocateScalars(VTK_FLOAT, 1);
#else
  imagedata->SetScalarType(VTK_FLOAT);
  imagedata->AllocateScalars();
#endif
  const int nvox = image->NumberOfSpatialVoxels();
  const irtkRegisteredImage::VoxelType *ptr1 = image->Data();
  float *ptr2 = reinterpret_cast<float *>(imagedata->GetScalarPointer());
  for (int i = 0; i < nvox; ++i, ++ptr1, ++ptr2) {
    *ptr2 = static_cast<float>(*ptr1);
  }
  return imagedata;
}

// -----------------------------------------------------------------------------
/// Map surface points to voxel coordinates
vtkSmartPointer<vtkPolyData> WorldToImage(vtkPolyData *surface, const irtkRegisteredImage *image)
{
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(surface->GetNumberOfPoints());
  double p[3];
  for (vtkIdType ptId = 0; ptId < points->GetNumberOfPoints(); ++ptId) {
    surface->GetPoint(ptId, p);
    image->WorldToImage(p[0], p[1], p[2]);
    points->SetPoint(ptId, p);
  }
  vtkSmartPointer<vtkPolyData> output = vtkSmartPointer<vtkPolyData>::New();
  output->ShallowCopy(surface);
  output->SetPoints(points);
  return output;
}

// -----------------------------------------------------------------------------
/// Get inside surface image stencil
vtkSmartPointer<vtkImageStencilData> ImageStencil(vtkPolyData *surface, vtkImageData *image)
{
  vtkSmartPointer<vtkPolyDataToImageStencil> filter;
  filter = vtkSmartPointer<vtkPolyDataToImageStencil>::New();
  SetVTKInput(filter, surface);
  filter->SetOutputOrigin(image->GetOrigin());
  filter->SetOutputSpacing(image->GetSpacing());
  filter->SetOutputWholeExtent(image->GetExtent());
  filter->Update();
  return filter->GetOutput();
}

// -----------------------------------------------------------------------------
/// Convert surface image stencil to binary mask, e.g., for debugging purposes
irtkBinaryImage ImageStencilToMask(const irtkImageAttributes &attr,
                                   vtkImageData              *data,
                                   vtkImageStencilData       *stencil)
{
  irtkBinaryImage mask(attr);
  const float * const start = reinterpret_cast<float *>(data->GetScalarPointer());
  vtkImageStencilIterator<float> it;
  it.Initialize(data, stencil, data->GetExtent());
  while (!it.IsAtEnd()) {
    if (it.IsInStencil()) {
      for (const float *cur = it.BeginSpan(); cur != it.EndSpan(); ++cur) {
        mask(cur - start) = true;
      }
    }
    it.NextSpan();
  }
  return mask;
}

// -----------------------------------------------------------------------------
/// Convert surface image stencil to binary mask, e.g., for debugging purposes
irtkBinaryImage ImageStencilToMask(const irtkImageAttributes &attr,
                                   vtkImageData              *data,
                                   vtkImageStencilData       *stencil,
                                   double                     p[3],
                                   double                     radius)
{
  irtkBinaryImage mask(attr);

  double rx = radius / attr._dx;
  double ry = radius / attr._dy;
  double rz = radius / attr._dz;

  double x = p[0], y = p[1], z = p[2];
  attr.WorldToLattice(x, y, z);

  int extent[6];
  extent[0] = floor(x - rx);
  extent[1] = ceil (x + rx);
  extent[2] = floor(y - ry);
  extent[3] = ceil (y + ry);
  extent[4] = floor(z - rz);
  extent[5] = ceil (z + rz);

  int i, i2, j, k, iter;
  for (k = extent[4]; k <= extent[5]; ++k)
  for (j = extent[2]; j <= extent[3]; ++j) {
    iter = 0;
    while (stencil->GetNextExtent(i, i2, extent[0], extent[1], j, k, iter)) {
      for (; i <= i2; ++i) {
        mask(i, j, k) = true;
      }
    }
  }

  return mask;
}

// -----------------------------------------------------------------------------
/// Compute point intensity thresholds based on local image statistics
struct ComputeLocalIntensityThresholds
{
  vtkPoints           *_Points;
  irtkImage           *_Image;
  vtkImageStencilData *_Stencil;
  vtkDataArray        *_LowerIntensity;
  vtkDataArray        *_UpperIntensity;
  double               _SigmaFactor;
  double               _RadiusX;
  double               _RadiusY;
  double               _RadiusZ;

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    int    n, extent[6], i, i2, j, k, iter;
    double p[3], value, mu, var, delta, sigma;
    vtkSmartPointer<vtkImageStencilData> roi = vtkSmartPointer<vtkImageStencilData>::New();
    vtkImageStencilIterator<float> it;

    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId) {
      n  = 0;
      mu = var = .0;

      _Points->GetPoint(ptId, p);
      _Image->WorldToImage(p[0], p[1], p[2]);

      extent[0] = floor(p[0] - _RadiusX);
      extent[1] = ceil (p[0] + _RadiusX);
      extent[2] = floor(p[1] - _RadiusY);
      extent[3] = ceil (p[1] + _RadiusY);
      extent[4] = floor(p[2] - _RadiusZ);
      extent[5] = ceil (p[2] + _RadiusZ);

      for (k = extent[4]; k <= extent[5]; ++k)
      for (j = extent[2]; j <= extent[3]; ++j) {
        iter = 0;
        while (_Stencil->GetNextExtent(i, i2, extent[0], extent[1], j, k, iter)) {
          for (; i <= i2; ++i) {
            ++n;
            value = _Image->GetAsDouble(i, j, k);
            delta = value - mu;
            mu  += delta / n;
            var += delta * (value - mu);
          }
        }
      }
      if (n > 2) var /= n - 1;
      else       var = .0;

      sigma = sqrt(var);
      _LowerIntensity->SetComponent(ptId, 0, mu - _SigmaFactor * sigma);
      _UpperIntensity->SetComponent(ptId, 0, mu + _SigmaFactor * sigma);
//      _LowerIntensity->SetComponent(ptId, 0, mu - 50);
//      _UpperIntensity->SetComponent(ptId, 0, mu + 200);
    }
  }
};

// -----------------------------------------------------------------------------
/// Compute local statistics of intensities inside/outside surface mesh
struct ComputeLocalIntensityStatistics
{
  vtkPoints             *_Points;
  const irtkImage       *_Image;
  const irtkBinaryImage *_ForegroundMask;
  const irtkBinaryImage *_BackgroundMask;
  vtkDataArray          *_ForegroundStatistics;
  vtkDataArray          *_BackgroundStatistics;
  double                 _RadiusX;
  double                 _RadiusY;
  double                 _RadiusZ;

  enum Component { BG, FG };

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    int       num[2], i, j, k, i1, i2, j1, j2, k1, k2;
    double    p[3], value, delta, mean[2], var[2];
    bool      is_bg, is_fg;
    Component c;

    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId) {
      num[0] = num[1] = 0;
      mean[0] = mean[1] = var[0] = var[1] = .0;

      _Points->GetPoint(ptId, p);
      _Image->WorldToImage(p[0], p[1], p[2]);

      i1 = floor(p[0] - _RadiusX);
      i2 = ceil (p[0] + _RadiusX);
      j1 = floor(p[1] - _RadiusY);
      j2 = ceil (p[1] + _RadiusY);
      k1 = floor(p[2] - _RadiusZ);
      k2 = ceil (p[2] + _RadiusZ);

      for (k = k1; k <= k2; ++k)
      for (j = j1; j <= j2; ++j)
      for (i = i1; i <= i2; ++i) {
        if (_ForegroundMask->IsInside(i, j, k)) {
          is_bg = (_BackgroundMask ? _BackgroundMask->Get(i, j, k) != 0 : true);
          is_fg =                    _ForegroundMask->Get(i, j, k) != 0;
          if (is_fg || is_bg) {
            c = static_cast<Component>(is_fg);
            value = _Image->GetAsDouble(i, j, k);
            delta = value - mean[c];
            num [c] += 1;
            mean[c] += delta / num[c];
            var [c] += delta * (value - mean[c]);
          }
        }
      }
      if (num[0] > 2) var[0] /= num[0] - 1;
      else            var[0] = .0;
      if (num[1] > 2) var[1] /= num[1] - 1;
      else            var[1] = .0;

      _BackgroundStatistics->SetComponent(ptId, 0, mean[BG]);
      _BackgroundStatistics->SetComponent(ptId, 1, sqrt(var[BG]));
      _ForegroundStatistics->SetComponent(ptId, 0, mean[FG]);
      _ForegroundStatistics->SetComponent(ptId, 1, sqrt(var[FG]));
    }
  }
};

// -----------------------------------------------------------------------------
/// Update balloon force magnitude and direction
struct UpdateMagnitude
{
  vtkPoints           *_Points;
  const ImageFunction *_Image;
  double               _LowerIntensity;
  double               _UpperIntensity;
  vtkDataArray        *_LocalLowerIntensity;
  vtkDataArray        *_LocalUpperIntensity;
  vtkDataArray        *_BackgroundStatistics;
  double               _BackgroundSigmaFactor;
  vtkDataArray        *_ForegroundStatistics;
  double               _ForegroundSigmaFactor;
  vtkDataArray        *_Magnitude;
  double               _MagnitudeDamping;
  double               _MagnitudeThreshold;

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    bool   inside;
    double w, p[3], v, mean, sigma, bgPb, fgPb;

    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId) {
      w = _Magnitude->GetComponent(ptId, 0);
      if (w == .0) continue;
      // Get intensity at current node position
      _Points->GetPoint(ptId, p);
      _Image->WorldToImage(p[0], p[1], p[2]);
      v = _Image->Evaluate(p[0], p[1], p[2]);
      // Check global intensity thresholds
      inside = (_LowerIntensity <= v && v <= _UpperIntensity);
      // Check local intensity thresholds
      if (inside) {
        if (_LocalLowerIntensity && v < _LocalLowerIntensity->GetComponent(ptId, 0)) {
          inside = false;
        }
        if (_LocalUpperIntensity && v > _LocalUpperIntensity->GetComponent(ptId, 0)) {
          inside = false;
        }
      }
      // Check background/foreground statistics
      if (inside && _BackgroundStatistics && _ForegroundStatistics) {
        mean = _BackgroundStatistics->GetComponent(ptId, 0);
        if (mean == .0 || v < mean) {
          inside = false;
        } else {
          sigma = _BackgroundStatistics->GetComponent(ptId, 1) * _BackgroundSigmaFactor;
          if (sigma > .0) bgPb = exp(- pow(v - mean, 2) / (2.0 * sigma * sigma)) / sigma;
          else            bgPb = .0;
          mean  = _ForegroundStatistics->GetComponent(ptId, 0);
          sigma = _ForegroundStatistics->GetComponent(ptId, 1) * _ForegroundSigmaFactor;
          if (sigma > .0) {
            fgPb = exp(- pow(v - mean, 2) / (2.0 * sigma * sigma)) / sigma;
            if (bgPb > fgPb) {
              inside = false;
            }
          } else if (v != mean) {
            inside = false;
          }
        }
      }
      // Adjust sign of balloon force and damp magnitude if direction changes
      if (inside) {
        if (w < .0) {
          w = _MagnitudeDamping * (-w);
          if (w < _MagnitudeThreshold) w = .0;
          _Magnitude->SetComponent(ptId, 0, w);
        }
      } else {
        if (w > .0) {
          w = _MagnitudeDamping * (-w);
          if (-w < _MagnitudeThreshold) w = .0;
          _Magnitude->SetComponent(ptId, 0, w);
        }
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Average balloon force magnitude across adjacent nodes
struct SmoothMagnitude
{
  const irtkEdgeTable *_EdgeTable;
  vtkDataArray        *_Normals;
  vtkDataArray        *_Input;
  vtkDataArray        *_Output;
  double               _MinMagnitude;

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    double     v, m, n0[3], n1[3], w, W;
    int        numAdjPts;
    const int *adjPts;

    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId) {
      m = fabs(_Input->GetComponent(ptId, 0));
//      if (m != .0) {
        W = 1.0;
        _Normals->GetTuple(ptId, n0);
        _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPts);
        for (int i = 0; i < numAdjPts; ++i) {
          v = fabs(_Input->GetComponent(adjPts[i], 0));
          _Normals->GetTuple(adjPts[i], n1);
          w = vtkMath::Dot(n0, n1);
          if (w > .5) {
            m += w * v;
            W += w;
          }
        }
        m /= W;
        m = copysign(m, _Input->GetComponent(ptId, 0));
//        m = copysign(max(m, _MinMagnitude), _Input->GetComponent(ptId, 0));
//      }
      _Output->SetComponent(ptId, 0, m);
    }
  }
};

// -----------------------------------------------------------------------------
/// Compute balloon force gradient (i.e., negative force)
struct ComputeGradient
{
  typedef irtkBalloonForce::GradientType Force;

  vtkDataArray *_Normals;
  vtkDataArray *_Magnitude;
  Force        *_Gradient;

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    double w, n[3];
    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId) {
      w = _Magnitude->GetComponent(ptId, 0);
      if (w == .0) continue;
      _Normals->GetTuple(ptId, n);
      vtkMath::MultiplyScalar(n, w);
      _Gradient[ptId] = -Force(n[0], n[1], n[2]);
    }
  }
};


} // namespace irtkBalloonForceUtils
using namespace irtkBalloonForceUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkBalloonForce::irtkBalloonForce(const char *name, double weight)
:
  irtkSurfaceForce(name, weight),
  _LowerIntensity(-numeric_limits<double>::infinity()),
  _UpperIntensity( numeric_limits<double>::infinity()),
  _SigmaFactor(5.0),
  _ForegroundSigmaFactor(1.0),
  _BackgroundSigmaFactor(1.0),
  _Radius(-1.0),
  _DampingFactor(.67),
  _MagnitudeThreshold(.1)
{
}

// -----------------------------------------------------------------------------
void irtkBalloonForce::Copy(const irtkBalloonForce &other)
{
  _LowerIntensity        = other._LowerIntensity;
  _UpperIntensity        = other._UpperIntensity;
  _SigmaFactor           = other._SigmaFactor;
  _ForegroundSigmaFactor = other._ForegroundSigmaFactor;
  _BackgroundSigmaFactor = other._BackgroundSigmaFactor;
  _Radius                = other._Radius;
  _DampingFactor         = other._DampingFactor;
  _MagnitudeThreshold    = other._MagnitudeThreshold;
}

// -----------------------------------------------------------------------------
irtkBalloonForce::irtkBalloonForce(const irtkBalloonForce &other)
:
  irtkSurfaceForce(other)
{
  Copy(other);
}

// -----------------------------------------------------------------------------
irtkBalloonForce &irtkBalloonForce::operator =(const irtkBalloonForce &other)
{
  if (this != &other) {
    irtkSurfaceForce::operator =(other);
    Copy(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
irtkBalloonForce::~irtkBalloonForce()
{
}

// =============================================================================
// Configuration
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkBalloonForce::Set(const char *param, const char *value)
{
  const string name = ParameterNameWithoutPrefix(param);

  if (name == "Minimum intensity" || name == "Lower threshold" || name == "Lower intensity" || name == "Lower intensity threshold") {
    return FromString(value, _LowerIntensity);
  }
  if (name == "Maximum intensity" || name == "Upper threshold" || name == "Upper intensity" || name == "Upper intensity threshold") {
    return FromString(value, _UpperIntensity);
  }
  if (name == "Local intensity sigma" || name == "Local intensity sigma factor") {
    return FromString(value, _SigmaFactor);
  }
  if (name == "Background intensity sigma" || name == "Background intensity sigma factor") {
    return FromString(value, _BackgroundSigmaFactor);
  }
  if (name == "Foreground intensity sigma" || name == "Foreground intensity sigma factor") {
    return FromString(value, _ForegroundSigmaFactor);
  }
  if (name == "Local window radius") {
    return FromString(value, _Radius);
  }
  if (name == "Damping factor" || name == "Magnitude damping factor" || name == "Magnitude damping") {
    return FromString(value, _DampingFactor);
  }
  if (name == "Magnitude threshold") {
    return FromString(value, _MagnitudeThreshold);
  }

  return irtkSurfaceForce::Set(param, value);
}

// -----------------------------------------------------------------------------
irtkParameterList irtkBalloonForce::Parameter() const
{
  irtkParameterList params = irtkSurfaceForce::Parameter();
  InsertWithPrefix(params, "Lower intensity",            _LowerIntensity);
  InsertWithPrefix(params, "Upper intensity",            _UpperIntensity);
  InsertWithPrefix(params, "Local intensity sigma",      _SigmaFactor);
  InsertWithPrefix(params, "Background intensity sigma", _BackgroundSigmaFactor);
  InsertWithPrefix(params, "Foreground intensity sigma", _ForegroundSigmaFactor);
  InsertWithPrefix(params, "Local window radius",        _Radius);
  InsertWithPrefix(params, "Damping factor",             _DampingFactor);
  InsertWithPrefix(params, "Magnitude threshold",        _MagnitudeThreshold);
  return params;
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
void irtkBalloonForce::ComputeLocalIntensityAttributes(bool thresholds, bool bg_fg_stats)
{
  if (!thresholds && !bg_fg_stats) return;

  vtkSmartPointer<vtkImageData>        imagedata = ConvertImage(_Image);
  vtkSmartPointer<vtkPolyData>         surface   = WorldToImage(_PointSet->InputSurface(), _Image);
  vtkSmartPointer<vtkImageStencilData> stencil   = ImageStencil(surface, imagedata);

  if (thresholds) {
    ComputeLocalIntensityThresholds eval;
    eval._Points         = _PointSet->Points();
    eval._Image          = _Image;
    eval._Stencil        = stencil;
    eval._RadiusX        = _Radius / _Image->GetXSize();
    eval._RadiusY        = _Radius / _Image->GetYSize();
    eval._RadiusZ        = _Radius / _Image->GetZSize();
    eval._SigmaFactor    = _SigmaFactor;
    eval._LowerIntensity = GetPointData("Lower intensity");
    eval._UpperIntensity = GetPointData("Upper intensity");
    parallel_for(blocked_range<vtkIdType>(0, _NumberOfPoints), eval);
  }

  if (bg_fg_stats) {
    irtkBinaryImage fg_mask = ImageStencilToMask(_Image->Attributes(), imagedata, stencil);
    irtkBinaryImage bg_mask("CC00081XX08-gm-mask.nii.gz"); // FIXME
    ComputeLocalIntensityStatistics eval;
    eval._Points               = _PointSet->Points();
    eval._Image                = _Image;
    eval._BackgroundMask       = &bg_mask;
    eval._ForegroundMask       = &fg_mask;
    eval._RadiusX              = _Radius / _Image->GetXSize();
    eval._RadiusY              = _Radius / _Image->GetYSize();
    eval._RadiusZ              = _Radius / _Image->GetZSize();
    eval._BackgroundStatistics = GetPointData("Background statistics");
    eval._ForegroundStatistics = GetPointData("Foreground statistics");
    parallel_for(blocked_range<vtkIdType>(0, _NumberOfPoints), eval);
  }
}

// -----------------------------------------------------------------------------
void irtkBalloonForce::Initialize()
{
  // Initialize base class
  irtkSurfaceForce::Initialize();
  if (_NumberOfPoints == 0) return;

  // Initial magnitude and direction of balloon force
  AddPointData("Signed magnitude")->FillComponent(0, 1.0);
  _DampingFactor      = max(.0, min(_DampingFactor,      1.0));
  _MagnitudeThreshold = max(.0, min(_MagnitudeThreshold, 1.0));

  // Radius of box window used for local intensity statistics
  // If zero, only the global intensity thresholds are considered
  if (_Radius < .0) {
    _Radius = 7.0 * max(max(_Image->GetXSize(),
                            _Image->GetYSize()),
                            _Image->GetZSize());
  }

  // Add/remove point data arrays for local intensity thresholds
  // If sigma multiplier is negative, no local thresholds are used
  if (_Radius == .0 || _SigmaFactor < .0) {
    RemovePointData("Lower intensity");
    RemovePointData("Upper intensity");
  } else {
    AddPointData("Lower intensity");
    AddPointData("Upper intensity");
  }

  // Add/remove point data arrays for local inside/outside intensity statistics
  // If sigma multiplier is negative, no local inside/outside statistics are used
  if (_Radius == .0 || _ForegroundSigmaFactor < .0 || _BackgroundSigmaFactor < .0) {
    RemovePointData("Background statistics");
    RemovePointData("Foreground statistics");
  } else {
    AddPointData("Background statistics", 2);
    AddPointData("Foreground statistics", 2);
  }
}

// -----------------------------------------------------------------------------
void irtkBalloonForce::Update(bool gradient)
{
  // Note that _InitialUpdate is set to false by irtkPointSetForce::Update
  const bool initial = _InitialUpdate;

  // Update base class
  irtkSurfaceForce::Update(gradient);

  // Delayed initialization of local intensity thresholds
  // and update of local background/foreground statistics
  if (_Radius > .0) {
    const bool thresholds  = initial && _SigmaFactor > .0;
    const bool bg_fg_stats = _BackgroundSigmaFactor > .0 && _ForegroundSigmaFactor > .0;
    ComputeLocalIntensityAttributes(thresholds, bg_fg_stats);
  }

  // Get (optional) point data (possibly interpolated during remeshing)
  vtkDataArray *magnitude       = GetPointData("Signed magnitude");
  vtkDataArray *lower_intensity = GetPointData("Lower intensity",       true);
  vtkDataArray *upper_intensity = GetPointData("Upper intensity",       true);
  vtkDataArray *bg_statistics   = GetPointData("Background statistics", true);
  vtkDataArray *fg_statistics   = GetPointData("Foreground statistics", true);

  // Initialize continuous intensity function
  ImageFunction image;
  image.Input(_Image);
  image.Initialize();

  // Update force magnitude and direction
  UpdateMagnitude update;
  update._Image                 = &image;
  update._Points                = _PointSet->SurfacePoints();
  update._LowerIntensity        = _LowerIntensity;
  update._UpperIntensity        = _UpperIntensity;
  update._LocalLowerIntensity   = lower_intensity;
  update._LocalUpperIntensity   = upper_intensity;
  update._BackgroundStatistics  = bg_statistics;
  update._BackgroundSigmaFactor = _BackgroundSigmaFactor;
  update._ForegroundStatistics  = fg_statistics;
  update._ForegroundSigmaFactor = _ForegroundSigmaFactor;
  update._Magnitude             = magnitude;
  update._MagnitudeDamping      = _DampingFactor;
  update._MagnitudeThreshold    = _MagnitudeThreshold;
  parallel_for(blocked_range<vtkIdType>(0, _NumberOfPoints), update);

  // Smooth magnitude such that adjacent nodes move coherently
  const int _MagnitudeSmoothing = 0;
  if (_MagnitudeSmoothing > 0) {
    vtkSmartPointer<vtkDataArray> smoothed_magnitude;
    smoothed_magnitude = magnitude->NewInstance();
    smoothed_magnitude->SetNumberOfComponents(magnitude->GetNumberOfComponents());
    smoothed_magnitude->SetNumberOfTuples(magnitude->GetNumberOfTuples());
    irtkEdgeTable edgeTable(_PointSet->Surface());
    SmoothMagnitude smooth;
    smooth._Input        = magnitude;
    smooth._Output       = smoothed_magnitude;
    smooth._EdgeTable    = &edgeTable;
    smooth._Normals      = _PointSet->SurfaceNormals();
    smooth._MinMagnitude = _MagnitudeThreshold;
    blocked_range<vtkIdType> ptIds(0, magnitude->GetNumberOfTuples());
    for (int iter = 0; iter < _MagnitudeSmoothing; ++iter) {
      parallel_for(ptIds, smooth);
      swap(smooth._Input, smooth._Output);
    }
    if (smooth._Output != magnitude) {
      magnitude->CopyComponent(0, smoothed_magnitude, 0);
    }
  }
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void irtkBalloonForce::EvaluateGradient(double *gradient, double step, double weight)
{
  if (_NumberOfPoints == 0) return;

  memset(_Gradient, 0, _NumberOfPoints * sizeof(GradientType));

  ComputeGradient eval;
  eval._Normals   = _PointSet->SurfaceNormals();
  eval._Magnitude = GetPointData("Signed magnitude");
  eval._Gradient  = _Gradient;
  parallel_for(blocked_range<vtkIdType>(0, _NumberOfPoints), eval);

  irtkSurfaceForce::EvaluateGradient(gradient, step, weight / _NumberOfPoints);
}
