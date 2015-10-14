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

#include <irtkPolyDataSmoothing.h>

#include <irtkCommon.h>
#include <irtkEdgeTable.h>
#include <irtkMatrix3x3.h>
#include <irtkPolyDataUtils.h>

#include <vtkMath.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkPolyDataNormals.h>


namespace irtk { namespace polydata {


// =============================================================================
// Auxiliary functors
// =============================================================================

namespace irtkPolyDataSmoothingUtils {

// -----------------------------------------------------------------------------
/// Uniform node weight
struct UniformWeightKernel
{
  double operator ()(vtkIdType, double [3], double [3]) const
  {
    return 1.0;
  }
};

// -----------------------------------------------------------------------------
/// Inverse Euclidean distance node weight
class InverseDistanceKernel
{
  double _Sigma;

public:

  InverseDistanceKernel(double sigma = .0) : _Sigma(sigma) {}

  double operator ()(vtkIdType, double p0[3], double p1[3]) const
  {
    const double d = vtkMath::Distance2BetweenPoints(p0, p1) + _Sigma;
    return (d == .0 ? .0 : 1.0 / d);
  }
};

// -----------------------------------------------------------------------------
/// Isotropic Gaussian node weight
class GaussianKernel
{
  double _Scale;

public:

  GaussianKernel(double sigma = 1.0) : _Scale(- .5 / (sigma * sigma)) {}

  double operator ()(vtkIdType, double p0[3], double p1[3]) const
  {
    return exp(_Scale * vtkMath::Distance2BetweenPoints(p0, p1));
  }
};

// -----------------------------------------------------------------------------
/// Anisotropic Gaussian node weight
///
/// https://tschumperle.users.greyc.fr/publications/tschumperle_ijcv06.pdf
class AnisotropicGaussianKernel
{
  vtkDataArray *_Tensors;
  vtkDataArray *_Normals;
  vtkDataArray *_MinimumDirection;
  vtkDataArray *_MaximumDirection;
  double        _Scale[3];

  /// Component indices of tensor output array
  ///
  /// This order is compatible with ParaView, which automatically labels the
  /// entries in this order. Also used by irtkPolyDataCurvature.
  enum TensorIndex
  {
    XX = 0, YY = 1, ZZ = 2, XY = 3, YX = XY, YZ = 4, ZY = YZ, XZ = 5, ZX = XZ
  };

public:

  /// Default constructor
  AnisotropicGaussianKernel()
  :
    _Tensors(NULL), _Normals(NULL), _MinimumDirection(NULL), _MaximumDirection(NULL)
  {
    _Scale[0] = _Scale[1] = _Scale[2] = .0;
  }

  /// Construct Gaussian kernel from local geometry tensors
  /// and isotropic standard deviation in locally oriented coordinate system
  AnisotropicGaussianKernel(vtkDataArray *tensors, double sigma = 1.0)
  :
    _Tensors(tensors), _Normals(NULL), _MinimumDirection(NULL), _MaximumDirection(NULL)
  {
    // Tensor contains individual scaling factors along each local axis
    _Scale[0] = _Scale[1] = _Scale[2] = - .5 / (sigma * sigma);
  }

  /// Construct Gaussian kernel from orthonormal geometry tensor with anisotropic
  /// standard deviation along directions of minimum/maximum change
  /// (i.e., second and third local coordinate axes)
  AnisotropicGaussianKernel(vtkDataArray *tensors, double sigma1, double sigma2)
  :
    _Tensors(tensors), _Normals(NULL), _MinimumDirection(NULL), _MaximumDirection(NULL)
  {
    // Tensor should be orthonormal basis
    _Scale[1] = -.5 / (sigma1 * sigma1);
    _Scale[2] = -.5 / (sigma2 * sigma2);
    _Scale[0] = -.5 / pow(min(sigma1, sigma2), 2);
  }

  /// Construct Gaussian kernel from orthonormal local basis vectors with
  /// half the extent along the normal direction and the direction of maximum change
  AnisotropicGaussianKernel(vtkDataArray *n, vtkDataArray *e1, vtkDataArray *e2, double sigma = 1.0)
  :
    _Tensors(NULL), _Normals(n), _MinimumDirection(e1), _MaximumDirection(e2)
  {
    _Scale[1] = -.5   / (sigma * sigma); // sigma_kmin   =   sigma
    _Scale[2] = -.125 / (sigma * sigma); // sigma_kmax   = 2 sigma
    _Scale[0] = _Scale[1];               // sigma_normal = sigma_kmin
  }

  /// Construct Gaussian kernel from orthonormal local basis vectors with
  /// specified standard deviation in directions of minimum and maximum change,
  /// respectively. The standard deviation in the normal direction is equal
  /// the minimum standard deviation in either of the other orthogonal directions.
  AnisotropicGaussianKernel(vtkDataArray *n, vtkDataArray *e1, vtkDataArray *e2,
                            double sigma1, double sigma2)
  :
    _Tensors(NULL), _Normals(n), _MinimumDirection(e1), _MaximumDirection(e2)
  {
    _Scale[1] = -.5 / (sigma1 * sigma1);
    _Scale[2] = -.5 / (sigma2 * sigma2);
    _Scale[0] = -.5 / pow(min(sigma1, sigma2), 2);
  }

  /// Evaluate anisotropic Gaussian kernel centered at \c p0 at \c x=p1
  double operator ()(vtkIdType ptId, double p0[3], double p1[3]) const
  {
    irtkMatrix3x3 T; // local geometry tensor, e.g., curvature tensor
    irtkVector3   x(p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]);

    if (_Tensors) {
      double m[6];
      _Tensors->GetTuple(ptId, m);
      T = irtkMatrix3x3(m[XX], m[XY], m[XZ],
                        m[YX], m[YY], m[YZ],
                        m[ZX], m[ZY], m[ZZ]);
    } else {
      double n[3], e1[3], e2[3];
      if (_MinimumDirection && _MaximumDirection) {
        _MinimumDirection->GetTuple(ptId, e1);
        _MaximumDirection->GetTuple(ptId, e2);
        vtkMath::Normalize(e1);
        vtkMath::Normalize(e2);
        vtkMath::Cross(e1, e2, n);
      } else {
        _Normals->GetTuple(ptId, n);
        if (_MinimumDirection) {
          _MinimumDirection->GetTuple(ptId, e1);
          vtkMath::Normalize(e1);
          vtkMath::Cross(n, e1, e2);
        } else {
          _MaximumDirection->GetTuple(ptId, e2);
          vtkMath::Normalize(e2);
          vtkMath::Cross(e2, n, e1);
        }
      }

      T = irtkMatrix3x3(n[0], e1[0], e2[0],
                        n[1], e1[1], e2[1],
                        n[2], e1[2], e2[2]);
    }

    x = T * x;

    x[0] *= x[0], x[1] *= x[1], x[2] *= x[2];
    return exp(_Scale[0] * x[0] + _Scale[1] * x[1] + _Scale[2] * x[2]);
  }
};

// -----------------------------------------------------------------------------
/// Smooth node position and/or data using the given node weighting kernel function
template <class TKernel>
struct SmoothData
{
  typedef irtkPolyDataSmoothing::DataArrays DataArrays;

  const irtkEdgeTable *_EdgeTable;
  vtkPoints           *_InputPoints;
  vtkPoints           *_OutputPoints;
  const DataArrays    *_InputArrays;
  const DataArrays    *_OutputArrays;
  TKernel              _WeightFunction;
  double               _Lambda;
  bool                 _InclNodeItself;

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    double        p[3] = {.0}, p0[3], p1[3], v[3], w, norm, alpha, beta;
    vtkDataArray *ia, *oa;
    const int    *adjPtIt, *adjPtEnd;
    vtkIdType     adjPtId;
    size_t        i;
    int           j;

    const bool smooth_points = _OutputPoints != NULL;
    const bool smooth_data   = _OutputArrays && !_InputArrays->empty();

    double **data = NULL, *sum;
    if (smooth_data) {
      data = new double *[_InputArrays->size()];
      for (i = 0; i < _InputArrays->size(); ++i) {
        data[i] = new double[(*_InputArrays)[i]->GetNumberOfComponents()];
      }
    }

    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId) {
      _InputPoints->GetPoint(ptId, p0);

      // Initialize sums
      if (_InclNodeItself) {
        norm = (w = _WeightFunction(ptId, p0, p0));
        if (smooth_points) {
          p[0] = w * p0[0];
          p[1] = w * p0[1];
          p[2] = w * p0[2];
        }
        if (smooth_data) {
          for (i = 0; i < _InputArrays->size(); ++i) {
            ia  = (*_InputArrays )[i];
            sum = data[i];
            for (j = 0; j < ia->GetNumberOfComponents(); ++j, ++sum) {
              (*sum) = w * ia->GetComponent(ptId, j);
            }
          }
        }
      } else {
        norm = p[0] = p[1] = p[2] = .0;
        if (smooth_data) {
          for (i = 0; i < _InputArrays->size(); ++i) {
            sum = data[i];
            ia  = (*_InputArrays )[i];
            for (j = 0; j < ia->GetNumberOfComponents(); ++j, ++sum) (*sum) = .0;
          }
        }
      }

      // Weighted sum of input data
      for (_EdgeTable->GetAdjacentPoints(ptId, adjPtIt, adjPtEnd); adjPtIt != adjPtEnd; ++adjPtIt) {
        adjPtId = static_cast<vtkIdType>(*adjPtIt);
        _InputPoints->GetPoint(adjPtId, p1);
        norm += (w = _WeightFunction(ptId, p0, p1));
        if (smooth_points) {
          p[0] += w * p1[0];
          p[1] += w * p1[1];
          p[2] += w * p1[2];
        }
        if (smooth_data) {
          for (i = 0; i < _InputArrays->size(); ++i) {
            ia  = (*_InputArrays )[i];
            sum = data[i];
            if (ia->GetNumberOfComponents() == 3) {
              ia->GetTuple(adjPtId, v);
              if (vtkMath::Dot(sum, v) < .0) {
                vtkMath::MultiplyScalar(v, -1.0);
              }
              vtkMath::MultiplyScalar(v, w);
              vtkMath::Add(sum, v, sum);
            } else {
              for (j = 0; j < ia->GetNumberOfComponents(); ++j, ++sum) {
                (*sum) += w * ia->GetComponent(adjPtId, j);
              }
            }
          }
        }
      }

      // Normalize weights and compute output data
      if (norm > .0) alpha = 1.0 - _Lambda, beta = _Lambda / norm;
      else           alpha = 1.0,           beta = .0;
      if (alpha) {
        if (smooth_points) {
          _OutputPoints->SetPoint(ptId, alpha * p0[0] + beta * p[0],
                                        alpha * p0[1] + beta * p[1],
                                        alpha * p0[2] + beta * p[2]);
        }
        if (smooth_data) {
          for (i = 0; i < _InputArrays->size(); ++i) {
            ia  = (*_InputArrays )[i];
            oa  = (*_OutputArrays)[i];
            sum = data[i];
            for (j = 0; j < oa->GetNumberOfComponents(); ++j, ++sum) {
              oa->SetComponent(ptId, j, alpha * ia->GetComponent(ptId, j) + beta * (*sum));
            }
          }
        }
      } else {
        if (smooth_points) {
          _OutputPoints->SetPoint(ptId, beta * p[0], beta * p[1], beta * p[2]);
        }
        if (smooth_data) {
          for (i = 0; i < _InputArrays->size(); ++i) {
            oa  = (*_OutputArrays)[i];
            sum = data[i];
            for (j = 0; j < oa->GetNumberOfComponents(); ++j, ++sum) {
              oa->SetComponent(ptId, j, beta * (*sum));
            }
          }
        }
      }
    }

    if (smooth_data) {
      for (i = 0; i < _InputArrays->size(); ++i) delete[] data[i];
      delete[] data;
    }
  }

  static void Run(const irtkEdgeTable *edgeTable,
                  vtkPoints           *input_points,
                  vtkPoints           *output_points,
                  const DataArrays    &input_arrays,
                  const DataArrays    &output_arrays,
                  TKernel              kernel,
                  double               lambda,
                  bool                 incl_node)
  {
    SmoothData<TKernel> body;
    body._EdgeTable      = edgeTable;
    body._InputPoints    = input_points;
    body._OutputPoints   = output_points;
    body._InputArrays    = &input_arrays;
    body._OutputArrays   = &output_arrays;
    body._WeightFunction = kernel;
    body._Lambda         = lambda;
    body._InclNodeItself = incl_node;
    blocked_range<vtkIdType> ptIds(0, input_points->GetNumberOfPoints());
    parallel_for(ptIds, body);
  }
};


} // namespace irtkPolyDataSmoothingUtils
using namespace irtkPolyDataSmoothingUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkPolyDataSmoothing::irtkPolyDataSmoothing()
:
  _EdgeTable(NULL),
  _EdgeTableOwner(false),
  _NumberOfIterations(1),
  _Lambda(1.0),
  _Sigma(.0),
  _MaximumDirectionSigma(.0),
  _Weighting(InverseDistance),
  _AdjacentValuesOnly(false),
  _SmoothPoints(false)
{
}

// -----------------------------------------------------------------------------
void irtkPolyDataSmoothing::Copy(const irtkPolyDataSmoothing &other)
{
  if (_EdgeTableOwner) Delete(_EdgeTable);
  _EdgeTable      = NULL;
  _EdgeTableOwner = false;

  if (other._EdgeTable) {
    _EdgeTableOwner = other._EdgeTableOwner;
    _EdgeTable      = (_EdgeTableOwner ? new irtkEdgeTable(*other._EdgeTable) : other._EdgeTable);
  }

  _NumberOfIterations    = other._NumberOfIterations;
  _Lambda                = other._Lambda;
  _Sigma                 = other._Sigma;
  _MaximumDirectionSigma = other._MaximumDirectionSigma;
  _Weighting             = other._Weighting;
  _GeometryTensorName    = other._GeometryTensorName;
  _MinimumDirectionName  = other._MinimumDirectionName;
  _MaximumDirectionName  = other._MaximumDirectionName;
  _AdjacentValuesOnly    = other._AdjacentValuesOnly;
  _SmoothPoints          = other._SmoothPoints;
  _SmoothArrays          = other._SmoothArrays;

  _InputArrays .clear();
  _OutputArrays.clear();
}

// -----------------------------------------------------------------------------
irtkPolyDataSmoothing::irtkPolyDataSmoothing(const irtkPolyDataSmoothing &other)
:
  irtkPolyDataFilter(other)
{
  Copy(other);
}

// -----------------------------------------------------------------------------
irtkPolyDataSmoothing &irtkPolyDataSmoothing::operator =(const irtkPolyDataSmoothing &other)
{
  irtkPolyDataFilter::operator =(other);
  Copy(other);
  return *this;
}

// -----------------------------------------------------------------------------
irtkPolyDataSmoothing::~irtkPolyDataSmoothing()
{
  if (_EdgeTableOwner) Delete(_EdgeTable);
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void irtkPolyDataSmoothing::Initialize()
{
  // Initialize base class -- makes shallow copy of input
  irtkPolyDataFilter::Initialize();

  // Construct edge table if none provided
  if (_EdgeTable == NULL) {
    _EdgeTable      = new irtkEdgeTable(_Input);
    _EdgeTableOwner = true;
  }

  // Check input point data arrays required for anisotropic smoothing
  if (_Weighting == AnisotropicGaussian) {
    vtkDataArray *array;
    if (_GeometryTensorName.empty()) {
      if (_MinimumDirectionName.empty() && _MaximumDirectionName.empty()) {
        _Weighting = Gaussian;
      } else {
        if (!_MinimumDirectionName.empty()) {
          array = _Input->GetPointData()->GetArray(_MinimumDirectionName.c_str());
          if (!array) {
            cerr << "irtkPolyDataSmoothing::Initialize: Missing input point data array named " << _MinimumDirectionName << endl;
            exit(1);
          }
          if (array->GetNumberOfComponents() != 3) {
            cerr << "irtkPolyDataSmoothing::Initialize: Invalid direction array. Must have 3 components." << endl;
            exit(1);
          }
        }
        if (!_MaximumDirectionName.empty()) {
          array = _Input->GetPointData()->GetArray(_MaximumDirectionName.c_str());
          if (!array) {
            cerr << "irtkPolyDataSmoothing::Initialize: Missing input point data array named " << _MaximumDirectionName << endl;
            exit(1);
          }
          if (array->GetNumberOfComponents() != 3) {
            cerr << "irtkPolyDataSmoothing::Initialize: Invalid direction array. Must have 3 components." << endl;
            exit(1);
          }
        }
        if (_MinimumDirectionName.empty() || _MaximumDirectionName.empty()) {
          vtkSmartPointer<vtkPolyDataNormals> filter;
          filter = vtkSmartPointer<vtkPolyDataNormals>::New();
          SetVTKInput(filter, _Input);
          filter->ComputeCellNormalsOff();
          filter->ComputePointNormalsOn();
          filter->SplittingOff();
          filter->AutoOrientNormalsOn();
          filter->ConsistencyOn();
          filter->Update();
          _Output = filter->GetOutput();
        }
      }
    } else {
      array = _Input->GetPointData()->GetArray(_GeometryTensorName.c_str());
      if (!array) {
        cerr << "irtkPolyDataSmoothing::Initialize: Missing input point data array named " << _GeometryTensorName << endl;
        exit(1);
      }
      if (array->GetNumberOfComponents() != 6 &&
          array->GetNumberOfComponents() != 9) {
        cerr << "irtkPolyDataSmoothing::Initialize: Invalid local geometry tensor array. Must have either 6 or 9 components." << endl;
        exit(1);
      }
    }
  }

  // Set default kernel parameter value
  if (_Weighting == Gaussian || _Weighting == AnisotropicGaussian) {
    if (_Sigma <= .0 || _MaximumDirectionSigma < .0) {
      double mean = AverageEdgeLength(_Input->GetPoints(), *_EdgeTable);
      if      (_Sigma == .0) _Sigma = mean;
      else if (_Sigma <  .0) _Sigma = fabs(_Sigma) * mean;
      if (_MaximumDirectionSigma <  .0) {
        _MaximumDirectionSigma = fabs(_MaximumDirectionSigma) * mean;
      }
    }
    if (_MaximumDirectionSigma == .0) {
      if (_GeometryTensorName.empty()) {
        _MaximumDirectionSigma = 2.0 * _Sigma;
      } else {
        _MaximumDirectionSigma = _Sigma;
      }
    }
  }

  // Smooth node positions by default
  if (!_SmoothPoints && _SmoothArrays.empty()) _SmoothPoints = true;

  // Get input point data arrays
  _InputArrays.resize(_SmoothArrays.size());
  for (size_t i = 0; i < _SmoothArrays.size(); ++i) {
    if (_SmoothArrays[i].empty()) {
      cerr << "irtkPolyDataSmoothing::Initialize: Empty input point data array name" << endl;
      exit(1);
    }
    _InputArrays[i] = _Input->GetPointData()->GetArray(_SmoothArrays[i].c_str());
    if (_InputArrays[i] == NULL) {
      cerr << "irtkPolyDataSmoothing::Initialize: Missing input point data array named " << _SmoothArrays[i] << endl;
      exit(1);
    }
  }

  // Allocate output points
  if (_SmoothPoints) {
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(_Input->GetNumberOfPoints());
    _Output->SetPoints(points);
  }

  // Allocate output arrays
  _OutputArrays.resize(_InputArrays.size());
  for (size_t i = 0; i < _InputArrays.size(); ++i) {
    _OutputArrays[i] = _InputArrays[i]->NewInstance();
    _OutputArrays[i]->SetName(_InputArrays[i]->GetName());
    _OutputArrays[i]->SetNumberOfComponents(_InputArrays[i]->GetNumberOfComponents());
    _OutputArrays[i]->SetNumberOfTuples(_InputArrays[i]->GetNumberOfTuples());
  }
}

// -----------------------------------------------------------------------------
void irtkPolyDataSmoothing::Execute()
{
  vtkSmartPointer<vtkPoints> ip = _Input->GetPoints();
  DataArrays                 ia = _InputArrays;
  vtkPoints * const          op = _SmoothPoints ? _Output->GetPoints() : NULL;
  const DataArrays          &oa = _OutputArrays;
  const bool incl_node = !_AdjacentValuesOnly;

  for (int iter = 0; iter < _NumberOfIterations; ++iter) {
    if (iter > 0) {
      if (_SmoothPoints) {
        if (iter == 1) ip = ip->NewInstance();
        ip->DeepCopy(_Output->GetPoints());
      }
      for (size_t i = 0; i < ia.size(); ++i) {
        if (iter == 1) ia[i] = ia[i]->NewInstance();
        ia[i]->DeepCopy(_OutputArrays[i]);
      }
    }
    switch (_Weighting) {
      case Combinatorial: {
        typedef UniformWeightKernel Kernel;
        SmoothData<Kernel>::Run(_EdgeTable, ip, op, ia, oa, Kernel(), _Lambda, incl_node);
      } break;
      case InverseDistance: {
        typedef InverseDistanceKernel Kernel;
        SmoothData<Kernel>::Run(_EdgeTable, ip, op, ia, oa, Kernel(_Sigma), _Lambda, incl_node);
      } break;
      case Gaussian: {
        typedef GaussianKernel Kernel;
        SmoothData<Kernel>::Run(_EdgeTable, ip, op, ia, oa, Kernel(_Sigma), _Lambda, incl_node);
      } break;
      case AnisotropicGaussian: {
        typedef AnisotropicGaussianKernel Kernel;
        vtkPointData * const pd = _Input->GetPointData();
        if (!_GeometryTensorName.empty()) {
          Kernel kernel(pd->GetArray(_GeometryTensorName.c_str()), _Sigma, _MaximumDirectionSigma);
          SmoothData<Kernel>::Run(_EdgeTable, ip, op, ia, oa, kernel, _Lambda, incl_node);
        } else {
          Kernel kernel(pd->GetNormals(),
                        pd->GetArray(_MinimumDirectionName.c_str()),
                        pd->GetArray(_MaximumDirectionName.c_str()),
                        _Sigma, _MaximumDirectionSigma);
          SmoothData<Kernel>::Run(_EdgeTable, ip, op, ia, oa, kernel, _Lambda, incl_node);
        }
      } break;
    }
  }
}

// -----------------------------------------------------------------------------
void irtkPolyDataSmoothing::Finalize()
{
  // Set output point data
  for (size_t i = 0; i < _OutputArrays.size(); ++i) {
    _Output->GetPointData()->RemoveArray(_OutputArrays[i]->GetName());
    _Output->GetPointData()->AddArray(_OutputArrays[i]);
  }
  _InputArrays .clear();
  _OutputArrays.clear();

  // Destroy edge table
  if (_EdgeTableOwner) {
    Delete(_EdgeTable);
    _EdgeTableOwner = false;
  }

  // Finalize base class
  irtkPolyDataFilter::Finalize();
}


} } // namespace irtk::polydata
