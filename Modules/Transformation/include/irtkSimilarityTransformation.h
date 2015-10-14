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

#ifndef _IRTKSIMILARITYTRANSFORMATION_H
#define _IRTKSIMILARITYTRANSFORMATION_H


/**
 * Class for similarity transformations.
 *
 * This class defines and implements similarity transformations. In addition to
 * the rigid body transformation parameters, similarity transformations are
 * parameterized by a global scaling parameter. The scaling parameter defines
 * the scaling along all axis of the coordinate transformations.
 *
 */

class irtkSimilarityTransformation : public irtkRigidTransformation
{
  irtkTransformationMacro(irtkSimilarityTransformation);

protected:

  /// Update transformation matrix after change of parameter
  virtual void UpdateMatrix();

  /// Update transformation parameters after change of matrix
  virtual void UpdateDOFs();

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor with given number of parameters
  irtkSimilarityTransformation(int);

  /// Copy constructor with given number of parameters
  irtkSimilarityTransformation(const irtkRigidTransformation &, int);

  /// Copy constructor with given number of parameters
  irtkSimilarityTransformation(const irtkSimilarityTransformation &, int);

public:

  /// Default constructor
  irtkSimilarityTransformation();

  /// Copy constructor
  irtkSimilarityTransformation(const irtkRigidTransformation &);

  /// Copy constructor
  irtkSimilarityTransformation(const irtkSimilarityTransformation &);

  /// Destructor
  virtual ~irtkSimilarityTransformation();

  // ---------------------------------------------------------------------------
  // Approximation

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds !new! parameters such that the resulting
  /// transformation approximates the displacements as good as possible.
  virtual void ApproximateDOFs(const double *, const double *, const double *, const double *,
                               const double *, const double *, const double *, int);

  // ---------------------------------------------------------------------------
  // Transformation parameters

  /// Copy active transformation parameters (DoFs) from given
  /// transformation if possible and return \c false, otherwise
  virtual bool CopyFrom(const irtkTransformation *);

  /// Puts scaling factor
  virtual void PutScale(double);

  /// Gets scaling factor
  virtual double GetScale() const;

  /// Construct a matrix based on parameters passed in the array
  static irtkMatrix DOFs2Matrix(const double *);

  // ---------------------------------------------------------------------------
  // Derivatives

  // Do not overwrite other base class overloads
  using irtkTransformation::JacobianDOFs;

  /// Calculates the Jacobian of the transformation w.r.t the parameters
  virtual void JacobianDOFs(double [3], int, double, double, double, double = 0, double = -1) const;

  /// Calculates the derivative of the Jacobian of the transformation (w.r.t. world coordinates) w.r.t. a transformation parameter
  virtual void DeriveJacobianWrtDOF(irtkMatrix &, int, double, double, double, double = 0, double = -1) const;

  // ---------------------------------------------------------------------------
  // I/O
  using irtkRigidTransformation::Write;

  /// Prints the parameters of the transformation
  virtual void Print(irtkIndent = 0) const;

  /// Whether this transformation can read a file of specified type (i.e. format)
  virtual bool CanRead(irtkTransformationType) const;

  /// Writes transformation to a file stream
  virtual irtkCofstream &Write(irtkCofstream &) const;

protected:

  /// Reads transformation from a file stream
  virtual irtkCifstream &ReadDOFs(irtkCifstream &, irtkTransformationType);

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Transformation parameters
// =============================================================================

// -----------------------------------------------------------------------------
inline void irtkSimilarityTransformation::PutScale(double s)
{
  Put(SG, s);
}

// -----------------------------------------------------------------------------
inline double irtkSimilarityTransformation::GetScale() const
{
  return Get(SG);
}


#endif
