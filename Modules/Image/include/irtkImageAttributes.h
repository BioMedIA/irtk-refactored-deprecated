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

#ifndef IRTKIMAGEATTRIBUTES_H
#define IRTKIMAGEATTRIBUTES_H

#include <irtkVector.h>
#include <irtkMatrix.h>


/**
 * Class which defines the attributes of the imaging geometry
 */

class irtkImageAttributes : public irtkObject
{
  irtkObjectMacro(irtkImageAttributes);

public:

  /// Image x-dimension (in voxels)
  int _x;
  
  /// Image y-dimension (in voxels)
  int _y;
  
  /// Image z-dimension (in voxels)
  int _z;
  
  /// Image t-dimension (in voxels)
  int _t;

  /// Voxel x-dimensions (in mm)
  double _dx;
  
  /// Voxel y-dimensions (in mm)
  double _dy;
  
  /// Voxel z-dimensions (in mm)
  double _dz;
  
  /// Voxel t-dimensions (in ms)
  double _dt;

  /// Image x-origin (in mm)
  double _xorigin;

  /// Image y-origin (in mm)
  double _yorigin;
  
  /// Image z-origin (in mm)
  double _zorigin;

  /// Image t-origin (in ms)
  double _torigin;

  /// Direction of x-axis
  double _xaxis[3];

  /// Direction of y-axis
  double _yaxis[3];

  /// Direction of z-axis
  double _zaxis[3];

  /// Affine transformation matrix
  irtkMatrix _smat;

  /// Constructor
  irtkImageAttributes();

  /// Constructor
  irtkImageAttributes(int, int, double = 1.0, double = 1.0);

  /// Constructor
  irtkImageAttributes(int, int, int, double = 1.0, double = 1.0, double = 1.0);

  /// Constructor
  irtkImageAttributes(int, int, int, int, double = 1.0, double = 1.0, double = 1.0, double = 1.0);

  /// Copy constructor
  irtkImageAttributes(const irtkImageAttributes &);

  /// Copy operator
  irtkImageAttributes& operator= (const irtkImageAttributes &);

  /// Check if attributes are valid
  operator bool() const;

  /// Whether given lattice is fully contained by this image lattice
  bool ContainsInSpace(const irtkImageAttributes &attr) const;

  /// Whether spatial attributes are equal
  bool EqualInSpace(const irtkImageAttributes &attr) const;

  /// Whether temporal attributes are equal
  bool EqualInTime(const irtkImageAttributes &attr) const;

  /// Equality operator
  bool operator==(const irtkImageAttributes &attr) const;

  /// Inequality operator
  bool operator!=(const irtkImageAttributes &attr) const;

  /// Number of lattice points in image domain
  int NumberOfLatticePoints() const;

  /// Number of spatial lattice points
  int NumberOfSpatialPoints() const;

  /// Number of lattice points
  /// (i.e., NumberOfLatticePoints if _dt != 0, otherwise NumberOfSpatialPoints)
  int NumberOfPoints() const;

  /// Get Index from Lattice
  int LatticeToIndex(int, int, int = 0, int = 0) const;

  /// Get Index from Lattice
  void IndexToLattice(int, int *, int *, int * = NULL, int * = NULL) const;

  /// Get Index from Lattice
  void IndexToLattice(int, int &, int &) const;

  /// Get Index from Lattice
  void IndexToLattice(int, int &, int &, int &) const;

  /// Get Index from Lattice
  void IndexToLattice(int, int &, int &, int &, int &) const;

  /// Get world coordinates (in mm) of lattice point
  void IndexToWorld(int, double &, double &) const;

  /// Get world coordinates (in mm) of lattice point
  void IndexToWorld(int, double &, double &, double &) const;

  /// Get world coordinates (in mm) of lattice point
  void IndexToWorld(int, irtkPoint &) const;

  /// Get world coordinates (in mm) of lattice point
  irtkPoint IndexToWorld(int) const;

  /// Convert lattice to world coordinate
  void LatticeToWorld(double &, double &, double &) const;

  /// Compute spatial world coordinates of lattice points (_x * _y * _z)
  void LatticeToWorld(double *, double *, double *) const;

  /// Compute world coordinates of all lattice points (_x * _y * _z * _t)
  void LatticeToWorld(double *, double *, double *, double *) const;

  /// Convert world to lattice coordinate
  void WorldToLattice(double &, double &, double &) const;

  /// Convert lattice to time coordinate
  double LatticeToTime(double) const;

  /// Convert time to lattice coordinate
  double TimeToLattice(double) const;

  /// Put affine world coordinate transformation which is applied
  /// after the image to world coordinate transformation derived from the
  /// imaging geometry when mapping voxel indices to world coordinates.
  /// This transformation can be the inverse of the affine transformation
  /// obtained by an affine registration with this image as source.
  ///
  /// \param[in] m     Homogeneous transformation matrix.
  /// \param[in] apply Whether to apply the translation, rotation, and
  ///                  scaling directly to the image attributes and store
  ///                  only the shearing (if any) as additional transformation.
  void PutAffineMatrix(const irtkMatrix &m, bool apply = false);

  /// Return transformation matrix for lattice to world coordinates
  irtkMatrix GetLatticeToWorldMatrix() const;

  /// Return transformation matrix for world to lattice coordinates
  irtkMatrix GetWorldToLatticeMatrix() const;

  /// Return orientation part of lattice to world coordinate transformation
  irtkMatrix GetLatticeToWorldOrientation() const;

  /// Return orientation part of lattice to world coordinate transformation
  irtkMatrix GetWorldToLatticeOrientation() const;

  /// Alias for GetLatticeToWorldMatrix
  inline irtkMatrix GetImageToWorldMatrix()      const { return GetLatticeToWorldMatrix(); }
  /// Alias for GetWorldToLatticeMatrix
  inline irtkMatrix GetWorldToImageMatrix()      const { return GetWorldToLatticeMatrix(); }
  /// Alias for GetLatticeToWorldOrientation
  inline irtkMatrix GetImageToWorldOrientation() const { return GetLatticeToWorldOrientation(); }
  /// Alias for GetWorldToLatticeOrientation
  inline irtkMatrix GetWorldToImageOrientation() const { return GetWorldToLatticeOrientation(); }

  /// Print attributes
  void Print(irtkIndent = 0) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline irtkImageAttributes::operator bool() const
{
  // Note: _dz may be zero for 2D images, _dt may even be negative!
  return _x > 0 && _y > 0 && _z > 0 && _t > 0 && _dx > .0 && _dy > .0 && _dz >= .0;
}

// -----------------------------------------------------------------------------
inline bool irtkImageAttributes::operator==(const irtkImageAttributes &attr) const
{
  return EqualInSpace(attr) && EqualInTime(attr);
}

// -----------------------------------------------------------------------------
inline bool irtkImageAttributes::operator!=(const irtkImageAttributes &attr) const
{
  return !(*this == attr);
}

// -----------------------------------------------------------------------------
inline int irtkImageAttributes::LatticeToIndex(int i, int j, int k, int l) const
{
  return i + _x * (j + _y * (k + _z * l));
}

// -----------------------------------------------------------------------------
inline void irtkImageAttributes::IndexToLattice(int index, int *i, int *j, int *k, int *l) const
{
  int n = _x * _y * _z;
  if (l) *l = index / n;
  index = index % n;
  n = _x * _y;
  if (k) *k = index / n;
  *j = index % n / _x;
  *i = index % n % _x;
}

// -----------------------------------------------------------------------------
inline int irtkImageAttributes::NumberOfLatticePoints() const
{
  return _x * _y * _z * _t;
}

// -----------------------------------------------------------------------------
inline int irtkImageAttributes::NumberOfSpatialPoints() const
{
  return _x * _y * _z;
}

// -----------------------------------------------------------------------------
inline int irtkImageAttributes::NumberOfPoints() const
{
  return _dt ? NumberOfLatticePoints() : NumberOfSpatialPoints();
}

// -----------------------------------------------------------------------------
inline void irtkImageAttributes::IndexToLattice(int index, int &i, int &j) const
{
  IndexToLattice(index, &i, &j);
}

// -----------------------------------------------------------------------------
inline void irtkImageAttributes::IndexToLattice(int index, int &i, int &j, int &k) const
{
  IndexToLattice(index, &i, &j, &k);
}

// -----------------------------------------------------------------------------
inline void irtkImageAttributes::IndexToLattice(int index, int &i, int &j, int &k, int &l) const
{
  IndexToLattice(index, &i, &j, &k, &l);
}

// -----------------------------------------------------------------------------
inline void irtkImageAttributes::IndexToWorld(int idx, double &x, double &y) const
{
  int i, j, k;
  IndexToLattice(idx, i, j, k);
  x = i, y = j;
  double z = k;
  LatticeToWorld(x, y, z);
}

// -----------------------------------------------------------------------------
inline void irtkImageAttributes::IndexToWorld(int idx, double &x, double &y, double &z) const
{
  int i, j, k;
  IndexToLattice(idx, i, j, k);
  x = i, y = j, z = k;
  LatticeToWorld(x, y, z);
}

// -----------------------------------------------------------------------------
inline void irtkImageAttributes::IndexToWorld(int idx, irtkPoint &p) const
{
  IndexToWorld(idx, p._x, p._y, p._z);
}

// -----------------------------------------------------------------------------
inline irtkPoint irtkImageAttributes::IndexToWorld(int idx) const
{
  irtkPoint p;
  IndexToWorld(idx, p);
  return p;
}

// -----------------------------------------------------------------------------
inline void irtkImageAttributes::LatticeToWorld(double &x, double &y, double &z) const
{
  irtkMatrix m = GetImageToWorldMatrix();
  double a = m(0, 0) * x + m(0, 1) * y + m(0, 2) * z + m(0, 3);
  double b = m(1, 0) * x + m(1, 1) * y + m(1, 2) * z + m(1, 3);
  double c = m(2, 0) * x + m(2, 1) * y + m(2, 2) * z + m(2, 3);
  x = a, y = b, z = c;
}

// -----------------------------------------------------------------------------
inline double irtkImageAttributes::LatticeToTime(double t) const
{
  return _torigin + t * _dt;
}

// -----------------------------------------------------------------------------
inline void irtkImageAttributes::WorldToLattice(double &x, double &y, double &z) const
{
  irtkMatrix m = GetWorldToImageMatrix();
  double a = m(0, 0) * x + m(0, 1) * y + m(0, 2) * z + m(0, 3);
  double b = m(1, 0) * x + m(1, 1) * y + m(1, 2) * z + m(1, 3);
  double c = m(2, 0) * x + m(2, 1) * y + m(2, 2) * z + m(2, 3);
  x = a, y = b, z = c;
}

// -----------------------------------------------------------------------------
inline double irtkImageAttributes::TimeToLattice(double t) const
{
  return (_dt ? ((t - _torigin) / _dt) : .0);
}

// -----------------------------------------------------------------------------
inline irtkMatrix irtkImageAttributes::GetLatticeToWorldOrientation() const
{
  irtkMatrix R(3, 3);
  R(0, 0) = _xaxis[0];
  R(1, 0) = _xaxis[1];
  R(2, 0) = _xaxis[2];
  R(0, 1) = _yaxis[0];
  R(1, 1) = _yaxis[1];
  R(2, 1) = _yaxis[2];
  R(0, 2) = _zaxis[0];
  R(1, 2) = _zaxis[1];
  R(2, 2) = _zaxis[2];
  return R;
}

// -----------------------------------------------------------------------------
inline irtkMatrix irtkImageAttributes::GetWorldToLatticeOrientation() const
{
  irtkMatrix R(3, 3);
  R(0, 0) = _xaxis[0];
  R(0, 1) = _xaxis[1];
  R(0, 2) = _xaxis[2];
  R(1, 0) = _yaxis[0];
  R(1, 1) = _yaxis[1];
  R(1, 2) = _yaxis[2];
  R(2, 0) = _zaxis[0];
  R(2, 1) = _zaxis[1];
  R(2, 2) = _zaxis[2];
  return R;
}

////////////////////////////////////////////////////////////////////////////////
// Image domain helpers
////////////////////////////////////////////////////////////////////////////////

/// Ensure that image domain has orthogonal basis vectors, i.e., no additional
/// affine transformation which contains any shearing
///
/// \returns Image grid which fully contains the input image domain but without
///          additional affine transformation (_smat is identity matrix).
///
/// \remarks The returned image domain need not be an axis-aligned bounding box!
///          When the Geometric Tools Engine (GTEngine) library is available,
///          the minimum-volume bounding box is returned.
irtkImageAttributes OrthogonalFieldOfView(const irtkImageAttributes &);

/// This method implements a method to find a common image grid which fully
/// contains all given image grids in world space
///
/// The voxel size of the resulting grid corresponds to the average voxel size
/// of the input images. The orientation and origin are chosen such that the
/// resulting image domain overlaps all input images in world space with a near
/// to minimal image volume. The current implementation only uses an heuristic
/// approach to find such minimal oriented bounding box of the input domains.
/// When the Geometric Tools Engine (GTEngine) library is available, the
/// minimum-volume bounding box is returned.
///
/// Additional information regarding algorithms to find the optimal
/// minimum-volume oriented bounding box (OBB) can be found at:
/// - http://stackoverflow.com/questions/7282805/where-can-i-find-a-c-c-implementation-of-the-minimum-bound-box-algorithm
/// - http://www.geometrictools.com/LibMathematics/Containment/Containment.html
/// - http://link.springer.com/article/10.1007%2FBF00991005?LI=true
///
/// \note The final overall image grid computed by this function must be
///       independent of the order in which the input attributes are given
///       in the input list. Otherwise an inverse consistent registration
///       might still depend on the order of the input images.
irtkImageAttributes OverallFieldOfView(const vector<irtkImageAttributes> &);


#endif
