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

#include <irtkImageEdgeForce.h>

#include <irtkGaussianBlurring.h>
#include <irtkConvolutionFunction.h>
#include <irtkImageGradientFunction.h>

#include <vtkMath.h>


// =============================================================================
// Auxiliary functions
// =============================================================================

namespace irtkImageEdgeForceUtils {


// Type of edge map gradient interpolation/evaluation function
typedef irtkGenericFastLinearImageGradientFunction<
           irtkGenericImage<irtkRealPixel>
        > EdgeGradient;
//typedef irtkGenericLinearImageGradientFunction<
//           irtkGenericImage<irtkRealPixel>
//        > EdgeGradient;

// -----------------------------------------------------------------------------
/// Detect image edges using the Sobel operator
void DetectEdges(irtkRealImage &im, double sigma = .0)
{
  typedef irtkConvolutionFunction::ConvolveInX<irtkRealPixel> ConvX;
  typedef irtkConvolutionFunction::ConvolveInY<irtkRealPixel> ConvY;
  typedef irtkConvolutionFunction::ConvolveInZ<irtkRealPixel> ConvZ;

  irtkRealPixel h[3] = { .25, .5, .25};
  irtkRealPixel g[3] = {-.5,  .0, .5};

  const irtkImageAttributes &attr = im.Attributes();
  irtkRealImage gx(attr), gy(attr), gz(attr), gm(attr);

  if (sigma > .0) {
    irtkGaussianBlurring<irtkRealPixel> blur(sigma);
    blur.SetInput (&im);
    blur.SetOutput(&gm);
    blur.Run();
    im.CopyFrom(gm);
  }

  ParallelForEachVoxel(ConvX(&im, g, 3, false), attr, im, gx);
  ParallelForEachVoxel(ConvY(&gx, h, 3, false), attr, gx, gm);
  ParallelForEachVoxel(ConvZ(&gm, h, 3, false), attr, gm, gx);

  ParallelForEachVoxel(ConvX(&im, h, 3, false), attr, im, gy);
  ParallelForEachVoxel(ConvY(&gy, g, 3, false), attr, gy, gm);
  ParallelForEachVoxel(ConvZ(&gm, h, 3, false), attr, gm, gy);

  ParallelForEachVoxel(ConvX(&im, h, 3, false), attr, im, gz);
  ParallelForEachVoxel(ConvY(&gz, h, 3, false), attr, gz, gm);
  ParallelForEachVoxel(ConvZ(&gm, g, 3, false), attr, gm, gz);

  for (int i = 0; i < gm.NumberOfVoxels(); ++i) {
    im(i) = sqrt(gx(i) * gx(i) + gy(i) * gy(i) + gz(i) * gz(i));
  }
}

// -----------------------------------------------------------------------------
/// Evaluate edge force
struct EvaluateGradient
{
  typedef irtkImageEdgeForce::GradientType Force;

  vtkPoints    *_Points;
  vtkDataArray *_Normals;
  EdgeGradient *_EdgeGradient;
  Force        *_Gradient; // Note: Opposite direction than force vector!
  double        _MaxNorm;

  EvaluateGradient() : _MaxNorm(.0) {}

  EvaluateGradient(const EvaluateGradient &other, split)
  :
    _Points(other._Points),
    _Normals(other._Normals),
    _EdgeGradient(other._EdgeGradient),
    _Gradient(other._Gradient),
    _MaxNorm(other._MaxNorm)
  {}

  void join(const EvaluateGradient &other)
  {
    if (other._MaxNorm > _MaxNorm) _MaxNorm = other._MaxNorm;
  }

  void operator ()(const blocked_range<vtkIdType> &re)
  {
    double p[3], n[3], g[3], norm;
    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId) {
      _Points->GetPoint(ptId, p);
      _EdgeGradient->WorldToImage(p[0], p[1], p[2]);
      _EdgeGradient->Evaluate(g,  p[0], p[1], p[2]);
      if (_Normals) {
        _Normals->GetTuple(ptId, n);
        norm = -vtkMath::Dot(g, n);
        vtkMath::MultiplyScalar(n, norm);
        _Gradient[ptId] = -Force(n[0], n[1], n[2]);
      } else {
//        vtkMath::Normalize(g);
        norm = vtkMath::Norm(g);
        _Gradient[ptId] = -Force(g[0], g[1], g[2]);
      }
      if (norm > _MaxNorm) _MaxNorm = norm;
    }
  }
};


} // namespace irtkImageEdgeForceUtils
using namespace irtkImageEdgeForceUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkImageEdgeForce::irtkImageEdgeForce(const char *name, double weight)
:
  irtkExternalForce(name, weight),
  _Sigma(-1.0),
  _InNormalDirection(true)
{
  _ParameterPrefix.push_back("Image edge force ");
  _ParameterPrefix.push_back("Intensity edge force ");
  _ParameterPrefix.push_back("Edge force ");
}

// -----------------------------------------------------------------------------
void irtkImageEdgeForce::Copy(const irtkImageEdgeForce &other)
{
  _Sigma             = other._Sigma;
  _InNormalDirection = other._InNormalDirection;
  _EdgeField         = other._EdgeField;
}

// -----------------------------------------------------------------------------
irtkImageEdgeForce::irtkImageEdgeForce(const irtkImageEdgeForce &other)
:
  irtkExternalForce(other)
{
  Copy(other);
}

// -----------------------------------------------------------------------------
irtkImageEdgeForce &irtkImageEdgeForce::operator =(const irtkImageEdgeForce &other)
{
  if (this != &other) {
    irtkExternalForce::operator =(other);
    Copy(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
irtkImageEdgeForce::~irtkImageEdgeForce()
{
}

// =============================================================================
// Configuration
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkImageEdgeForce::Set(const char *param, const char *value)
{
  const string name = ParameterNameWithoutPrefix(param);

  if (name == "Blurring") {
    return FromString(value, _Sigma);
  }
  if (name == "In normal direction") {
    return FromString(value, _InNormalDirection);
  }

  return irtkExternalForce::Set(param, value);
}

// -----------------------------------------------------------------------------
irtkParameterList irtkImageEdgeForce::Parameter() const
{
  irtkParameterList params = irtkExternalForce::Parameter();
  InsertWithPrefix(params, "Blurring",            _Sigma);
  InsertWithPrefix(params, "In normal direction", _InNormalDirection);
  return params;
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
void irtkImageEdgeForce::Initialize()
{
  // Initialize base class
  irtkExternalForce::Initialize();
  if (_NumberOfPoints == 0) return;

  // Default smoothing sigma
  if (_Sigma < .0) {
    _Sigma = .7 * max(_Image->GetXSize(), _Image->GetYSize());
  }

  // Compute edge field
  _EdgeField = *_Image;
  DetectEdges(_EdgeField, _Sigma);

  // FIXME: Move this to irtkImageEdgeForce::WriteDataSets
  if (debug) _EdgeField.Write("irtkImageEdgeForce_edges.nii.gz");
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void irtkImageEdgeForce::EvaluateGradient(double *gradient, double step, double weight)
{
  if (_NumberOfPoints == 0) return;

  memset(_Gradient, 0, _NumberOfPoints * sizeof(GradientType));

  EdgeGradient edge_gradient;
  edge_gradient.WrtWorld(true);
  edge_gradient.Input(&_EdgeField);
  edge_gradient.Initialize();

  irtkImageEdgeForceUtils::EvaluateGradient eval;
  eval._Points              = _PointSet->SurfacePoints();
  eval._Normals             = _InNormalDirection ? _PointSet->SurfaceNormals() : NULL;
  eval._EdgeGradient        = &edge_gradient;
  eval._Gradient            = _Gradient;
  parallel_reduce(blocked_range<vtkIdType>(0, _NumberOfPoints), eval);

  if (eval._MaxNorm > .0) {
    for (int i = 0; i < _NumberOfPoints; ++i) {
      _Gradient[i] /= eval._MaxNorm;
    }
  }

  irtkExternalForce::EvaluateGradient(gradient, step, weight / _NumberOfPoints);
}
