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

#ifndef _IRTKFREEFORMTRANSFORMATION3D_H
#define _IRTKFREEFORMTRANSFORMATION3D_H


/**
 * Base class for 3D free-form transformations.
 */

class irtkFreeFormTransformation3D : public irtkFreeFormTransformation
{
  irtkAbstractTransformationMacro(irtkFreeFormTransformation3D);

  friend class Evaluate2ndOrderBSplineFFDDerivatives;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  // Import other overloads from base class
  using irtkFreeFormTransformation::DefaultAttributes;

  /// Default attributes of free-form transformation lattice
  static irtkImageAttributes DefaultAttributes(double, double, double,
                                               double, double, double,
                                               double, double, double,
                                               const double *,
                                               const double *,
                                               const double *);

protected:

  /// Default constructor
  irtkFreeFormTransformation3D(CPInterpolator &, CPInterpolator * = NULL);

  /// Copy Constructor
  irtkFreeFormTransformation3D(const irtkFreeFormTransformation3D &,
                               CPInterpolator &, CPInterpolator * = NULL);

public:

  /// Destructor
  virtual ~irtkFreeFormTransformation3D();

  // Import other Initialize overloades from base class
  using irtkFreeFormTransformation::Initialize;

  /// Initialize free-form transformation
  virtual void Initialize(const irtkImageAttributes &);

  // ---------------------------------------------------------------------------
  // Derivatives
  using irtkFreeFormTransformation::JacobianDOFs;
  using irtkFreeFormTransformation::ParametricGradient;

  /// Calculates the Jacobian of the transformation w.r.t a control point
  virtual void JacobianDOFs(irtkMatrix &, int, int, int, double, double, double) const;

  /// Calculates the Jacobian of the transformation w.r.t a control point
  virtual void JacobianDOFs(irtkMatrix &, int, double, double, double, double = 0, double = -1) const;

  /// Calculates the Jacobian of the transformation w.r.t a control point
  virtual void JacobianDOFs(double [3], int, int, int, double, double, double) const;

  /// Calculates the Jacobian of the transformation w.r.t a transformation parameter
  virtual void JacobianDOFs(double [3], int, double, double, double, double = 0, double = -1) const;

  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation
  virtual void ParametricGradient(const irtkGenericImage<double> *, double *,
                                  const irtkWorldCoordsImage *,
                                  const irtkWorldCoordsImage *,
                                  double = -1, double = 1) const;

  /// Applies the chain rule to convert point-wise non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  virtual void ParametricGradient(const irtkPointSet &, const irtkVector3D<double> *,
                                  double *, double = 0, double = -1, double = 1) const;

  // ---------------------------------------------------------------------------
  // I/O

protected:

  /// Reads transformation parameters from a file stream
  virtual irtkCifstream &ReadDOFs(irtkCifstream &, irtkTransformationType);

  /// Writes transformation parameters to a file stream
  virtual irtkCofstream &WriteDOFs(irtkCofstream &) const;

public:

  // ---------------------------------------------------------------------------
  // Deprecated

  /// \deprecated Use PutStatus instead.
  void PutStatusCP(int, int, int, DOFStatus, DOFStatus, DOFStatus);
  
  /// \deprecated Use GetStatus instead.
  void GetStatusCP(int, int, int, DOFStatus &, DOFStatus &, DOFStatus &) const;

  /// \deprecated Use overloaded BoundingBox method instead.
  ///             Note, however, that the new function swaps coordinates to
  ///             enforce p1[i] <= p2[i]. This is not done by this function.
  void BoundingBoxCP(int cp, irtkPoint &, irtkPoint &, double = 1) const;

  /// \deprecated Use overloaded BoundingBox member function instead.
  void BoundingBoxImage(const irtkImage *, int, int &, int &, int &,
                                                int &, int &, int &, double = 1) const;

  /// \deprecated Use overloaded BoundingBox member function instead.
  void MultiBoundingBoxImage(const irtkImage *,
                             int, int &, int &, int &,
                                  int &, int &, int &, double = 1) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation3D::JacobianDOFs(double jac[3], int, int, int, double, double, double) const
{
  cerr << this->NameOfClass() << "::JacobianDOFs: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation3D::JacobianDOFs(double jac[3], int dof, double x, double y, double z, double, double) const
{
  int ci, cj, ck;
  this->IndexToLattice(dof / 3, ci, cj, ck);
  this->JacobianDOFs(jac, ci, cj, ck, x, y, z);
  const int c = dof % 3;
  for (int i = 0; i < 3; ++i) {
    if (i != c) jac[i] = .0;
  }
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation3D::JacobianDOFs(irtkMatrix &jac, int ci, int cj, int ck, double x, double y, double z) const
{
  double tmp[3];
  this->JacobianDOFs(tmp, ci, cj, ck, x, y, z);
  jac.Initialize(3, 3);
  jac(0, 0) = tmp[0];
  jac(1, 1) = tmp[1];
  jac(2, 2) = tmp[2];
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation3D::JacobianDOFs(irtkMatrix &jac, int cp, double x, double y, double z, double, double) const
{
  int ci, cj, ck;
  this->IndexToLattice(cp, ci, cj, ck);
  this->JacobianDOFs(jac, ci, cj, ck, x, y, z);
}

// =============================================================================
// Deprecated
// =============================================================================

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation3D::PutStatusCP(int i, int j, int k, DOFStatus sx, DOFStatus sy, DOFStatus sz)
{
  PutStatus(i, j, k, sx, sy, sz);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation3D::GetStatusCP(int i, int j, int k, DOFStatus &sx, DOFStatus &sy, DOFStatus &sz) const
{
  GetStatus(i, j, k, sx, sy, sz);
}

// -----------------------------------------------------------------------------
inline void irtkFreeFormTransformation3D::BoundingBoxCP(int cp, irtkPoint &p1, irtkPoint &p2, double fraction) const
{
  int i, j, k;
  IndexToLattice(cp, i, j, k);
  p1 = irtkPoint(i - 2 * fraction, j - 2 * fraction, k - 2 * fraction);
  p2 = irtkPoint(i + 2 * fraction, j + 2 * fraction, k + 2 * fraction);
  this->LatticeToWorld(p1);
  this->LatticeToWorld(p2);
}

// -----------------------------------------------------------------------------
inline void
irtkFreeFormTransformation3D
::BoundingBoxImage(const irtkImage *image, int cp, int &i1, int &j1, int &k1,
                                                   int &i2, int &j2, int &k2, double fraction) const
{
  // Calculate bounding box in world coordinates
  irtkPoint p1, p2;
  this->BoundingBoxCP(cp, p1, p2, fraction);

  // Transform world coordinates to image coordinates
  image->WorldToImage(p1._x, p1._y, p1._z);
  image->WorldToImage(p2._x, p2._y, p2._z);

  // Calculate bounding box in image coordinates
  i1 = (p1._x < 0) ? 0 : int(p1._x)+1;
  j1 = (p1._y < 0) ? 0 : int(p1._y)+1;
  i2 = (int(p2._x) >= image->GetX()) ? image->GetX()-1 : int(p2._x);
  j2 = (int(p2._y) >= image->GetY()) ? image->GetY()-1 : int(p2._y);

  if (image->GetZ() == 1) {
    k1 = k2 = 0;
  } else {
    k1 = (p1._z < 0) ? 0 : int(p1._z)+1;
    k2 = (int(p2._z) >= image->GetZ()) ? image->GetZ()-1 : int(p2._z);
  }
}

// -----------------------------------------------------------------------------
inline void
irtkFreeFormTransformation3D
::MultiBoundingBoxImage(const irtkImage *image, int cp, int &i1, int &j1, int &k1,
                                                        int &i2, int &j2, int &k2, double fraction) const
{
  BoundingBoxImage(image, cp, i1, j1, k1, i2, j2, k2, 2.0 * fraction);
}


#endif
