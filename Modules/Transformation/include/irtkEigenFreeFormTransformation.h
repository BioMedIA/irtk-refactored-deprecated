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

#ifndef _IRTKEIGENFREEFORMTRANSFORMATION_H
#define _IRTKEIGENFREEFORMTRANSFORMATION_H


#undef NORMAL


/**
 * Adaptive FFD class definition based on eigen parameterization.
 *
 * This class implements adaptive hierarchical FFD's using
 * the eigenmodes (e.g. derived from biomechanics or statistics).
 *
 * \todo The correct functionality of this transformation has to be verified
 *       yet after a complete refactoring of the transformation package.
 *       -as12312
 */

class irtkEigenFreeFormTransformation : public irtkBSplineFreeFormTransformation3D
{
  irtkTransformationMacro(irtkEigenFreeFormTransformation);

  // ---------------------------------------------------------------------------
  // Data members

  /// Eigenvalues
  irtkPublicAttributeMacro(irtkVector, EigenValues);

  /// Eigenvectors
  irtkPublicAttributeMacro(irtkMatrix, EigenVectors);

  /// Shape vector
  irtkPublicAttributeMacro(irtkVector, ShapeVector);

protected:

  /// Label of each element (0=unlabeled by default)
  int *_Label;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  irtkEigenFreeFormTransformation();

  /// Constructor (FFD, eigen system, shape vector)
  irtkEigenFreeFormTransformation(const irtkBSplineFreeFormTransformation3D &,
                                  const irtkVector &,
                                  const irtkMatrix &,
                                  const irtkVector &);

  /// Copy Constructor
  irtkEigenFreeFormTransformation(const irtkEigenFreeFormTransformation &);

  /// Destructor
  virtual ~irtkEigenFreeFormTransformation();

  // ---------------------------------------------------------------------------
  // Transformation parameters

  /// Returns the number of DOFs (rows of _ShapeVector)
  virtual int NumberOfDOFs() const;

  /// Puts a shape parameter
  virtual void Put(int, double);

  /// Gets a shape parameter
  virtual double Get(int) const;

  /// Gets a DOF status (overloaded, always returns Active)
  virtual DOFStatus GetStatus(int) const;

  /// Returns the label of the indexed control point
  virtual int GetLabel(int) const;

  /// Sets the label of the indexed control point
  void PutLabel(int, int);

  /// Returns the number of elements of the transformation
  virtual int NumberOfElements() const;

  /// Gets an element by filling a pre-allocated array with node indices
  virtual void GetElement(int, int*) const;

  /// Gets an element by allocating array and performing overloaded call
  int *GetElement(int) const;

  // ---------------------------------------------------------------------------
  // I/O

  /// Prints the parameters of the transformation
  virtual void Print(irtkIndent = 0) const;

  /// Whether this transformation can read a file of specified type (i.e. format)
  virtual bool CanRead(irtkTransformationType) const;

protected:

  /// Reads a transformation parameters from a file stream
  virtual irtkCifstream &ReadDOFs(irtkCifstream &, irtkTransformationType);

  /// Writes a transformation parameters to a file stream
  virtual irtkCofstream &WriteDOFs(irtkCofstream &) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Transformation parameters
// =============================================================================

// -----------------------------------------------------------------------------
inline int irtkEigenFreeFormTransformation::NumberOfDOFs() const
{
  return _ShapeVector.Rows();
}

// -----------------------------------------------------------------------------
inline double irtkEigenFreeFormTransformation::Get(int index) const
{
  return _ShapeVector.Get(index);
}

// -----------------------------------------------------------------------------
inline DOFStatus irtkEigenFreeFormTransformation::GetStatus(int) const
{
  return Active;
}

// -----------------------------------------------------------------------------
inline void irtkEigenFreeFormTransformation::PutLabel(int index, int label)
{
  _Label[index] = label;
}

// -----------------------------------------------------------------------------
inline int irtkEigenFreeFormTransformation::GetLabel(int index) const
{
  return _Label[index];
}

// -----------------------------------------------------------------------------
inline int irtkEigenFreeFormTransformation::NumberOfElements() const
{
  if (_z > 1) {
    return (_x - 1) * (_y - 1) * (_z - 1);
  } else {
    return (_x - 1) * (_y - 1);
  }
}

// -----------------------------------------------------------------------------
inline void irtkEigenFreeFormTransformation::GetElement(int index, int *element) const
{
  int i, j, k;

  // Find top left node of element (=N1)
  // See vtkStructuredGrid::GetCell()
  if (_z > 1) {
    // 3D (case VTK_XYZ_GRID) - still needs checking!
    i =  index % ( _x - 1);
    j = (index / ( _x - 1)) % (_y - 1);
    k =  index / ((_x - 1)  * (_y - 1));
  } else {
    // 2D (case VTK_XY_PLANE)
    i = index % ( _x - 1);
    j = index / ( _x - 1);
    k = 0;
  }

  // In-plane element node indices
  element[0] = this->LatticeToIndex(i+1, j,   k);
  element[1] = this->LatticeToIndex(i,   j,   k);
  element[2] = this->LatticeToIndex(i,   j+1, k);
  element[3] = this->LatticeToIndex(i+1, j+1, k);

  // Through-plane element node indices, if applicable
  if (_z > 1) {
    element[4] = this->LatticeToIndex(i+1, j,   k+1);
    element[5] = this->LatticeToIndex(i,   j,   k+1);
    element[6] = this->LatticeToIndex(i,   j+1, k+1);
    element[7] = this->LatticeToIndex(i+1, j+1, k+1);
  }
}

// -----------------------------------------------------------------------------
inline int *irtkEigenFreeFormTransformation::GetElement(int index) const
{
  int *element = NULL;

  // Allocate 2D or 3D element
  if (_z > 1) element = new int[8];
  else        element = new int[4];

  // Overloaded call
  this->GetElement(index, element);

  // Return element
  return element;
}


#endif
