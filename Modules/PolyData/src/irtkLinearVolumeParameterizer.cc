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

#include <irtkLinearVolumeParameterizer.h>

#include <irtkParallel.h>
#include <irtkMatrix3x3.h>

#include <vtkIdList.h>
#include <vtkTetra.h>
#include <vtkMath.h>

#include <Eigen/IterativeLinearSolvers>

using namespace std;


namespace irtk { namespace polydata {


// =============================================================================
// Auxiliary functors
// =============================================================================

namespace irtkLinearVolumeParameterizerUtils {


// -----------------------------------------------------------------------------
template <class Scalar>
class LinearSystem
{
  const irtkLinearVolumeParameterizer      *_Filter;
  const irtkLinearVolumeParameterizer      *_Operator;
  std::vector<Eigen::Triplet<Scalar> >      _Coefficients;
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1>  _RightHandSide;

public:

  // ---------------------------------------------------------------------------
  LinearSystem(const irtkLinearVolumeParameterizer *filter,
               const irtkLinearVolumeParameterizer *map, int n)
  :
    _Filter(filter), _Operator(map ? map : filter), _RightHandSide(n)
  {
    _RightHandSide.setZero();
  }

  LinearSystem(const LinearSystem &other, split)
  :
    _Filter(other._Filter), _Operator(other._Operator),
    _RightHandSide(other._RightHandSide.rows())
  {
    _RightHandSide.setZero();
  }

  // ---------------------------------------------------------------------------
  void join(const LinearSystem &other)
  {
    _Coefficients.insert(_Coefficients.end(), other._Coefficients.begin(), other._Coefficients.end());
    _RightHandSide += other._RightHandSide;
  }

  // ---------------------------------------------------------------------------
  void AddWeight(vtkIdType ptId0, bool isBoundary0, vtkIdType ptId1, bool isBoundary1, const irtkMatrix3x3 &weight)
  {
    if (isBoundary0 && isBoundary1) {

      // Unused coefficients

    } else if (isBoundary0) {

      // Point variables base index
      const int c = _Filter->InteriorPointPos()[ptId1];

      // Pre-multiply coefficient by constant boundary coordinates
      for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        const double &w = weight[i][j];
        if (w == .0) continue;
        _RightHandSide(c + j) -= w * _Filter->Coords()->GetComponent(ptId0, i);
        _Coefficients.push_back(Eigen::Triplet<Scalar>(c + j, c + i, -w));
      }

    } else if (isBoundary1) {

      // Point variables base index
      const int r = _Filter->InteriorPointPos()[ptId0];

      // Pre-multiply coefficient by constant boundary coordinates
      for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        const double &w = weight[i][j];
        if (w == .0) continue;
        _RightHandSide(r + i) -= w * _Filter->Coords()->GetComponent(ptId1, j);
        _Coefficients.push_back(Eigen::Triplet<Scalar>(r + i, r + j, -w));
      }

    } else {

      // Point variables base indices
      const int r = _Filter->InteriorPointPos()[ptId0];
      const int c = _Filter->InteriorPointPos()[ptId1];

      // Add symmetric coefficients
      for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        const double &w = weight[i][j];
        if (w == .0) continue;
        _Coefficients.push_back(Eigen::Triplet<Scalar>(r + i, c + j,  w));
        _Coefficients.push_back(Eigen::Triplet<Scalar>(r + i, r + j, -w));
        _Coefficients.push_back(Eigen::Triplet<Scalar>(c + j, r + i,  w));
        _Coefficients.push_back(Eigen::Triplet<Scalar>(c + j, c + i, -w));
      }

    }
  }

  // ---------------------------------------------------------------------------
  void operator ()(const blocked_range<vtkIdType> &cellIds)
  {
    vtkIdType i0, i1, i2, i3;
    bool      b0, b1, b2, b3;
    double    v0[3], v1[3], v2[3], v3[3], volume;

    vtkPointSet * const pointset = _Filter->Volume();
    vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();

    for (vtkIdType cellId = cellIds.begin(); cellId != cellIds.end(); ++cellId) {
      pointset->GetCellPoints(cellId, ptIds);

      i0 = ptIds->GetId(0);
      i1 = ptIds->GetId(1);
      i2 = ptIds->GetId(2);
      i3 = ptIds->GetId(3);

      b0 = _Filter->IsBoundaryPoint()[i0];
      b1 = _Filter->IsBoundaryPoint()[i1];
      b2 = _Filter->IsBoundaryPoint()[i2];
      b3 = _Filter->IsBoundaryPoint()[i3];

      pointset->GetPoint(i0, v0);
      pointset->GetPoint(i1, v1);
      pointset->GetPoint(i2, v2);
      pointset->GetPoint(i3, v3);

      volume = vtkTetra::ComputeVolume(v0, v1, v2, v3);

      AddWeight(i0, b0, i1, b1, _Operator->GetWeight(cellId, v0, v1, v2, v3, volume));
      AddWeight(i0, b0, i2, b2, _Operator->GetWeight(cellId, v0, v2, v3, v1, volume));
      AddWeight(i0, b0, i3, b3, _Operator->GetWeight(cellId, v0, v3, v1, v2, volume));
      AddWeight(i1, b1, i2, b2, _Operator->GetWeight(cellId, v1, v2, v0, v3, volume));
      AddWeight(i1, b1, i3, b3, _Operator->GetWeight(cellId, v1, v3, v2, v0, volume));
      AddWeight(i2, b2, i3, b3, _Operator->GetWeight(cellId, v2, v3, v0, v1, volume));
    }
  }

  // ---------------------------------------------------------------------------
  static void Build(const irtkLinearVolumeParameterizer      *filter,
                    const irtkLinearVolumeParameterizer      *mapop,
                    Eigen::SparseMatrix<Scalar>              &A,
                    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &b,
                    int                                       n)
  {
    LinearSystem problem(filter, mapop, n);
    blocked_range<vtkIdType> cellIds(0, filter->Volume()->GetNumberOfCells());
    parallel_reduce(cellIds, problem);
    A.resize(n, n);
    A.setFromTriplets(problem._Coefficients.begin(), problem._Coefficients.end());
    b = problem._RightHandSide;
  }
};


} // namespace irtkLinearVolumeParameterizerUtils
using namespace irtkLinearVolumeParameterizerUtils;

// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void irtkLinearVolumeParameterizer::Copy(const irtkLinearVolumeParameterizer &other)
{
  _NumberOfIterations = other._NumberOfIterations;
  _Tolerance          = other._Tolerance;
  _RelaxationFactor   = other._RelaxationFactor;
  _InteriorPointId    = other._InteriorPointId;
  _InteriorPointPos   = other._InteriorPointPos;
}

// -----------------------------------------------------------------------------
irtkLinearVolumeParameterizer::irtkLinearVolumeParameterizer()
:
  _NumberOfIterations(200),
  _Tolerance(1e-8),
  _RelaxationFactor(1.0)
{
}

// -----------------------------------------------------------------------------
irtkLinearVolumeParameterizer::irtkLinearVolumeParameterizer(const irtkLinearVolumeParameterizer &other)
:
  irtkVolumeParameterizer(other)
{
  Copy(other);
}

// -----------------------------------------------------------------------------
irtkLinearVolumeParameterizer &irtkLinearVolumeParameterizer::operator =(const irtkLinearVolumeParameterizer &other)
{
  irtkVolumeParameterizer::operator =(other);
  Copy(other);
  return *this;
}

// -----------------------------------------------------------------------------
irtkLinearVolumeParameterizer::~irtkLinearVolumeParameterizer()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void irtkLinearVolumeParameterizer::Initialize()
{
  const int d = 3; // Dimension of output domain

  // Initialize base class
  irtkVolumeParameterizer::Initialize();

  // Pre-compute maps from interior point index to global point ID as well as
  // the position of interior points (given its global ID) in the linear system
  _InteriorPointId .resize(_NumberOfInteriorPoints, -1);
  _InteriorPointPos.resize(_NumberOfPoints,         -1);
  for (int ptId = 0, i = 0; ptId < _NumberOfPoints; ++ptId) {
    if (_IsBoundaryPoint[ptId]) continue;
    _InteriorPointId [i   ] = ptId;
    _InteriorPointPos[ptId] = d * i;
    ++i;
  }
}

// -----------------------------------------------------------------------------
void irtkLinearVolumeParameterizer::Parameterize()
{
  Solve(this);
}

// -----------------------------------------------------------------------------
void irtkLinearVolumeParameterizer::Solve(const irtkLinearVolumeParameterizer *mapop)
{
  typedef double                                   Scalar;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
  typedef Eigen::SparseMatrix<Scalar>              Matrix;
  typedef Eigen::DiagonalPreconditioner<Scalar>    Preconditioner;

  const int d = 3;                           // Dimension of output domain
  const int n = d * _NumberOfInteriorPoints; // Size of linear system

  Matrix A;
  Vector x, b;

  // Use current parameterization of interior points as initial guess
  x.resize(n);
  for (int i = 0, r = 0; i < _NumberOfInteriorPoints; ++i, r += d) {
    for (int j = 0; j < d; ++j) {
      x(r + j) = static_cast<Scalar>(_Coords->GetComponent(_InteriorPointId[i], j));
    }
  }

  // Build linear system
  LinearSystem<Scalar>::Build(this, mapop, A, b, n);

  // Solve linear system
  Eigen::ConjugateGradient<Matrix, Eigen::Upper|Eigen::Lower, Preconditioner> solver;
  solver.setMaxIterations(_NumberOfIterations);
  solver.setTolerance(_Tolerance);
  x = solver.compute(A).solveWithGuess(b, x);

  // Update parameterization of interior points
  for (int i = 0, r = 0; i < _NumberOfInteriorPoints; ++i, r += d) {
    for (int j = 0; j < d; ++j) {
      _Coords->SetComponent(_InteriorPointId[i], j, static_cast<double>(x(r + j)));
    }
  }
}


} } // namespace irtk::polydata
