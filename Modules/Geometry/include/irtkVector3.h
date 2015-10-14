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

#ifndef IRTKVECTOR3_H

#define IRTKVECTOR3_H

#include <irtkArith.h>

class irtkVector3
{
public:
  // construction
  irtkVector3 ();
  irtkVector3 (double fX, double fY, double fZ);
  irtkVector3 (double afCoordinate[3]);
  irtkVector3 (const irtkVector3& rkVector);

  // member access (allows V.x or V[0], V.y or V[1], V.z or V[2])
  double x, y, z;
  double& operator[] (int i) const;
  operator double* ();

  // assignment and comparison
  irtkVector3& operator= (const irtkVector3& rkVector);
  bool operator== (const irtkVector3& rkVector) const;
  bool operator!= (const irtkVector3& rkVector) const;

  // arithmetic operations
  irtkVector3 operator+ (const irtkVector3& rkVector) const;
  irtkVector3 operator- (const irtkVector3& rkVector) const;
  irtkVector3 operator* (double fScalar) const;
  irtkVector3 operator/ (double fScalar) const;
  irtkVector3 operator- () const;
  friend irtkVector3 operator* (double fScalar, const irtkVector3& rkVector);

  // arithmetic updates
  irtkVector3& operator+= (const irtkVector3& rkVector);
  irtkVector3& operator-= (const irtkVector3& rkVector);
  irtkVector3& operator*= (double fScalar);
  irtkVector3& operator/= (double fScalar);

  // vector operations
  double Length () const;
  double SquaredLength () const;
  double Dot (const irtkVector3& rkVector) const;
  double Unitize (double fTolerance = 1e-06);
  irtkVector3 Cross (const irtkVector3& rkVector) const;
  irtkVector3 UnitCross (const irtkVector3& rkVector) const;

  // Gram-Schmidt orthonormalization.
  static void Orthonormalize (irtkVector3 akVector[3]);

  // Input W must be initialize to a nonzero vector, output is {U,V,W}
  // an orthonormal basis.  A hint is provided about whether or not W
  // is already unit length.
  static void GenerateOrthonormalBasis (irtkVector3& rkU, irtkVector3& rkV,
                                        irtkVector3& rkW, bool bUnitLengthW = true);

  // special points
  static const irtkVector3 ZERO;
  static const irtkVector3 UNIT_X;
  static const irtkVector3 UNIT_Y;
  static const irtkVector3 UNIT_Z;
};

#endif
