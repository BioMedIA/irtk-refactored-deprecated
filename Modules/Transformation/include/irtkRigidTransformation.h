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

#ifndef _IRTKRIGIDTRANSFORMATION_H
#define _IRTKRIGIDTRANSFORMATION_H


/**
 * Class for rigid transformations.
 *
 * This class defines and implements rigid body transformations. The rigid
 * body transformations are parameterized by three rotations around the axes
 * of the coordinate system followed by three translations along the axes of
 * the coordinate system. Note that the order of rotations is defined as a
 * rotation around the z-axis, the y-axis and finally around the x-axis. In
 * total, the transformation is parameterized by six degrees of freedom.
 *
 */

class irtkRigidTransformation : public irtkHomogeneousTransformation
{
  irtkTransformationMacro(irtkRigidTransformation);

  // ---------------------------------------------------------------------------
  // Data members

protected:

  /// Sine of rotation angle rx
  double _sinrx;

  /// Sine of rotation angle ry
  double _sinry;

  /// Sine of rotation angle rz
  double _sinrz;

  /// Cosine of rotation angle rx
  double _cosrx;
  
  /// Cosine of rotation angle ry
  double _cosry;
  
  /// Cosine of rotation angle rz
  double _cosrz;

  /// Update cached sine and cosine of rotation angles
  void UpdateRotationSineCosine();

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor with given number of parameters
  irtkRigidTransformation(int);

  /// Copy constructor with given number of parameters
  irtkRigidTransformation(const irtkRigidTransformation &, int);

public:

  /// Default constructor
  irtkRigidTransformation();

  /// Copy constructor
  irtkRigidTransformation(const irtkRigidTransformation &);

  /// Destructor
  virtual ~irtkRigidTransformation();

  // ---------------------------------------------------------------------------
  // Approximation

  /// Approximate displacements: This function takes a set of points and a set
  /// of displacements and finds !new! parameters such that the resulting
  /// transformation approximates the displacements as good as possible.
  virtual void ApproximateDOFs(const double *, const double *, const double *, const double *,
                               const double *, const double *, const double *, int);

  // ---------------------------------------------------------------------------
  // Transformation parameters

  /// Construct a matrix based on parameters passed in the array
  static irtkMatrix DOFs2Matrix(const double *);

  /// Return an array with parameters corresponding to a given matrix
  static void Matrix2DOFs(const irtkMatrix &, double *);

  /// Updates transformation matrix after change of parameter
  virtual void UpdateMatrix();

  /// Updates transformation parameters after change of matrix
  virtual void UpdateDOFs();

  /// Puts translation along the x-axis (transformation matrix is updated)
  void PutTranslationX(double);

  /// Gets translation along the x-axis
  double GetTranslationX() const;

  /// Puts translation along the y-axis (transformation matrix is updated)
  void PutTranslationY(double);

  /// Gets translation along the y-axis
  double GetTranslationY() const;

  /// Puts translation along the z-axis (transformation matrix is updated)
  void PutTranslationZ(double);

  /// Gets translation along the z-axis
  double GetTranslationZ() const;

  /// Puts rotation angle around the x-axis (transformation matrix is updated)
  void PutRotationX(double);

  /// Gets rotation angle around the x-axis
  double GetRotationX() const;

  /// Puts rotation angle around the y-axis (transformation matrix is updated)
  void PutRotationY(double);

  /// Gets rotation angle around the y-axis
  double GetRotationY() const;

  /// Puts rotation angle around the z-axis (transformation matrix is updated)
  void PutRotationZ(double);

  /// Gets rotation angle around the z-axis
  double GetRotationZ() const;

public:

  // ---------------------------------------------------------------------------
  // Point transformation

  /// Transforms a single point by the translation part of the rigid transformation
  virtual void Translate(double& x, double& y, double& z) const;

  /// Transforms a single point by the rotation part of the rigid transformation
  virtual void Rotate(double& x, double& y, double& z) const;

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

  /// Prints the parameters of the transformation
  virtual void Print(irtkIndent = 0) const;

public:

  // ---------------------------------------------------------------------------
  // Deprecated

  /// Set transformation parameters (DoFs)
  /// \deprecated Use Put(params) instead.
  void SetParameters(double *params);

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Transformation parameters
// =============================================================================

// -----------------------------------------------------------------------------
inline void irtkRigidTransformation::PutRotationX(double rx)
{
  Put(RX, rx);
}

// -----------------------------------------------------------------------------
inline double irtkRigidTransformation::GetRotationX() const
{
  return Get(RX);
}

// -----------------------------------------------------------------------------
inline void irtkRigidTransformation::PutRotationY(double ry)
{
  Put(RY, ry);
}

// -----------------------------------------------------------------------------
inline double irtkRigidTransformation::GetRotationY() const
{
  return Get(RY);
}

// -----------------------------------------------------------------------------
inline void irtkRigidTransformation::PutRotationZ(double rz)
{
  Put(RZ, rz);
}

// -----------------------------------------------------------------------------
inline double irtkRigidTransformation::GetRotationZ() const
{
  return Get(RZ);
}

// -----------------------------------------------------------------------------
inline void irtkRigidTransformation::PutTranslationX(double tx)
{
  Put(TX, tx);
}

// -----------------------------------------------------------------------------
inline double irtkRigidTransformation::GetTranslationX() const
{
  return Get(TX);
}

// -----------------------------------------------------------------------------
inline void irtkRigidTransformation::PutTranslationY(double ty)
{
  Put(TY, ty);
}

// -----------------------------------------------------------------------------
inline double irtkRigidTransformation::GetTranslationY() const
{
  return Get(TY);
}

// -----------------------------------------------------------------------------
inline void irtkRigidTransformation::PutTranslationZ(double tz)
{
  Put(TZ, tz);
}

// -----------------------------------------------------------------------------
inline double irtkRigidTransformation::GetTranslationZ() const
{
  return Get(TZ);
}

// =============================================================================
// Point transformation
// =============================================================================

// -----------------------------------------------------------------------------
inline void irtkRigidTransformation::Translate(double& x, double& y, double& z) const
{
  x += _matrix(0, 3);
  y += _matrix(1, 3);
  z += _matrix(2, 3);
}

// -----------------------------------------------------------------------------
inline void irtkRigidTransformation::Rotate(double& x, double& y, double& z) const
{
  double a = _matrix(0, 0) * x + _matrix(0, 1) * y + _matrix(0, 2) * z;
  double b = _matrix(1, 0) * x + _matrix(1, 1) * y + _matrix(1, 2) * z;
  double c = _matrix(2, 0) * x + _matrix(2, 1) * y + _matrix(2, 2) * z;

  x = a;
  y = b;
  z = c;
}

// =============================================================================
// Deprecated
// =============================================================================

// -----------------------------------------------------------------------------
inline void irtkRigidTransformation::SetParameters(double *params)
{
  Put(params);
}


#endif
