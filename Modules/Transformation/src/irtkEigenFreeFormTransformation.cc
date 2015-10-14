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

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkEigenFreeFormTransformation::irtkEigenFreeFormTransformation()
:
  _Label(NULL)
{
}

// -----------------------------------------------------------------------------
irtkEigenFreeFormTransformation
::irtkEigenFreeFormTransformation(const irtkBSplineFreeFormTransformation3D &ffd,
                                  const irtkVector &EigenValues,
                                  const irtkMatrix &EigenVectors,
                                  const irtkVector &ShapeVector)
:
  irtkBSplineFreeFormTransformation3D(ffd),
  _EigenValues (EigenValues),
  _EigenVectors(EigenVectors),
  _ShapeVector (ShapeVector)
{
  // Initialize element labels
  const int n = this->NumberOfElements();
  _Label = new int[n];
  memset(_Label, 0, n * sizeof(int));
}

// -----------------------------------------------------------------------------
irtkEigenFreeFormTransformation
::irtkEigenFreeFormTransformation(const irtkEigenFreeFormTransformation &ffd)
:
  irtkBSplineFreeFormTransformation3D(ffd),
  _EigenValues (ffd._EigenValues),
  _EigenVectors(ffd._EigenVectors),
  _ShapeVector (ffd._ShapeVector)
{
  // Copy element labels
  const int n = ffd.NumberOfElements();
  _Label = new int [n];
  memcpy(_Label, ffd._Label, n * sizeof(int));
}

// -----------------------------------------------------------------------------
irtkEigenFreeFormTransformation::~irtkEigenFreeFormTransformation()
{
  delete[] _Label;
}

// =============================================================================
// Transformation parameters
// =============================================================================

// -----------------------------------------------------------------------------
void irtkEigenFreeFormTransformation::Put(int index, double shape)
{
  // Only update control points for nonrigid modes
  if (this->GetStatus(index) == Active) {

    // Compute lambda
    double sqrt_lambda = sqrt(_EigenValues.Get(index));

    // Compute shape difference weight
    double shape_diff = (shape - _ShapeVector.Get(index)) * sqrt_lambda;

    // Update control points
    double x, y, z;
    int    n = 0;
    for (int k = 0; k < _z; ++k) {
      for (int j = 0; j < _y; ++j) {
        for (int i = 0; i < _x; ++i) {
          irtkBSplineFreeFormTransformation3D::Get(i, j, k, x, y, z);
          // Add difference between new and old shape parameter
          irtkBSplineFreeFormTransformation3D::Put(i, j, k, x + shape_diff * _EigenVectors.Get(3*n, index), y + shape_diff * _EigenVectors.Get(3*n+1, index), z + shape_diff * _EigenVectors.Get(3*n+2, index));
          ++n;
        }
      }
    }

    // Change shape vector
    _ShapeVector.Put(index, shape);
  }
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void irtkEigenFreeFormTransformation::Print(irtkIndent indent) const
{
  // Print keyword
  cout << indent << "Eigen FFD:" << endl;
  ++indent;

  // Print FFD
  irtkFreeFormTransformation3D::Print(indent);

  // Print labels
  const int n = this->NumberOfElements();
  cout << indent << "Labels:" << endl;
  for (int i = 0; i < n; ++i) {
    if (i % 40 == 0) cout << (indent + 1);
    cout << _Label[i] << " ";
  }
  cout << endl;

  // Print EigenSystem
  _EigenValues .Print(indent);
  _EigenVectors.Print(indent);
  _ShapeVector .Print(indent);
}

// -----------------------------------------------------------------------------
bool irtkEigenFreeFormTransformation::CanRead(irtkTransformationType format) const
{
  switch (format) {
    case IRTKTRANSFORMATION_EIGEN_FFD_3D_v1:
    case IRTKTRANSFORMATION_EIGEN_FFD_3D_v2:
    case IRTKTRANSFORMATION_EIGEN_FFD_3D_v3:
      return true;
    default:
      return false;
  }
}

// -----------------------------------------------------------------------------
irtkCifstream &irtkEigenFreeFormTransformation::ReadDOFs(irtkCifstream &from, irtkTransformationType format)
{
  if (format < IRTKTRANSFORMATION_EIGEN_FFD_3D_v2) {
    irtkBSplineFreeFormTransformation3D::ReadDOFs(from, IRTKTRANSFORMATION_BSPLINE_FFD_3D_v2);
  } else {
    irtkBSplineFreeFormTransformation3D::ReadDOFs(from, IRTKTRANSFORMATION_BSPLINE_FFD_3D_v3);
  }

  // Get number of elements
  const int n = this->NumberOfElements();

  // Allocate memory for labels
  delete[] _Label;
  _Label = new int[n];

  // Read labels as unsigned char
  unsigned char *label = new unsigned char[n];
  from.ReadAsUChar(label, n);
  for (int i = 0; i < n; i++) {
    _Label[i] = static_cast<int>(label[i]);
  }
  delete[] label;

  // Read EigenSystem
  from >> _EigenValues;
  from >> _EigenVectors;
  from >> _ShapeVector;

  return from;
}

// -----------------------------------------------------------------------------
irtkCofstream &irtkEigenFreeFormTransformation::WriteDOFs(irtkCofstream &to) const
{
  irtkBSplineFreeFormTransformation3D::WriteDOFs(to);

  // Get number of elements
  const int n = this->NumberOfElements();

  // Write labels as unsigned char
  unsigned char *label = new unsigned char[n];
  for (int i = 0; i < n; i++) {
    label[i] = static_cast<unsigned char>(_Label[i]);
  }
  to.WriteAsUChar(label, n);
  delete[] label;

  // Write EigenSystem
  to << _EigenValues;
  to << _EigenVectors;
  to << _ShapeVector;

  return to;
}
