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

#ifndef _IRTKPARTIALMULTILEVELSTATIONARYVELOCITYTRANSFORMATION_H
#define _IRTKPARTIALMULTILEVELSTATIONARYVELOCITYTRANSFORMATION_H


class irtkMultiLevelStationaryVelocityTransformation;


/**
 * Decorator for SV MFFD transformation.
 *
 * This decorator wraps a SV MFFD transformation but defines its own attributes
 * regarding integration and whether or not to invert the transformation by
 * default. As the SV MFFD is parameterized by a sum of stationary velocities,
 * inverting the transformation is simply achieved by a negative upper
 * integration limit.
 */

class irtkPartialMultiLevelStationaryVelocityTransformation : public irtkMultiLevelTransformation
{
  irtkTransformationMacro(irtkPartialMultiLevelStationaryVelocityTransformation);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Pointer to decorated transformation
  irtkPublicAggregateMacro(irtkMultiLevelStationaryVelocityTransformation, Transformation);

  /// Fraction of decorated transformation, negative value corresponds to inverse
  irtkPublicAttributeMacro(double, Fraction);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  irtkPartialMultiLevelStationaryVelocityTransformation(irtkMultiLevelStationaryVelocityTransformation * = NULL, double = 1.0);

  /// Copy constructor
  irtkPartialMultiLevelStationaryVelocityTransformation(const irtkPartialMultiLevelStationaryVelocityTransformation &);

  /// Destructor
  virtual ~irtkPartialMultiLevelStationaryVelocityTransformation();

  // ---------------------------------------------------------------------------
  // Levels

  /// Returns the number of levels
  int NumberOfLevels() const;

  /// Gets global transformation
  irtkAffineTransformation *GetGlobalTransformation();

  // Get global transformation
  const irtkAffineTransformation *GetGlobalTransformation() const;

  /// Gets local transformation
  irtkFreeFormTransformation *GetLocalTransformation(int);

  /// Gets local transformation
  const irtkFreeFormTransformation *GetLocalTransformation(int) const;

  /// Put local transformation and return pointer to previous one (needs to be deleted if not used)
  irtkFreeFormTransformation *PutLocalTransformation(irtkFreeFormTransformation *, int);

  /// Push local transformation on stack (append transformation)
  void PushLocalTransformation(irtkFreeFormTransformation *);

  /// Insert local transformation
  void InsertLocalTransformation(irtkFreeFormTransformation *, int = 0);

  /// Pop local transformation from stack (remove last transformation)
  irtkFreeFormTransformation *PopLocalTransformation();

  /// Remove local transformation and return the pointer (need to be deleted if not used)
  irtkFreeFormTransformation *RemoveLocalTransformation(int = 0);

  /// Combine local transformations on stack
  virtual void CombineLocalTransformation();

  /// Convert the global transformation from a matrix representation to a
  /// FFD and incorporate it with any existing local transformation
  virtual void MergeGlobalIntoLocalDisplacement();

  // ---------------------------------------------------------------------------
  // Transformation parameters (DoFs)

  /// Copy active transformation parameters (DoFs) from given
  /// transformation if possible and return \c false, otherwise
  virtual bool CopyFrom(const irtkTransformation *);

  /// Get number of transformation parameters
  virtual int NumberOfDOFs() const;

  /// Put value of transformation parameter
  virtual void Put(int, double);

  /// Put values of transformation parameters
  virtual void Put(const DOFValue *);

  /// Add change to transformation parameters
  virtual void Add(const DOFValue *);

  /// Update transformation parameters given parametric gradient
  virtual double Update(const DOFValue *);

  /// Get value of transformation parameter
  virtual double Get(int) const;

  /// Get values of transformation parameters
  virtual void Get(DOFValue *) const;

  /// Put status of transformation parameter
  virtual void PutStatus(int, DOFStatus);

  /// Get status of transformation parameter
  virtual DOFStatus GetStatus(int) const;

  /// Checks whether transformation depends on the same vector of parameters
  virtual bool HasSameDOFsAs(const irtkTransformation *) const;

  /// Checks whether the transformation is an identity mapping
  virtual bool IsIdentity() const;

  // ---------------------------------------------------------------------------
  // Parameters (non-DoFs)

  // Import other overloads
  using irtkTransformation::Parameter;

  /// Set named (non-DoF) parameter from value as string
  virtual bool Set(const char *, const char *);

  /// Get (non-DoF) parameters as key/value as string map
  virtual irtkParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Point transformation
  using irtkMultiLevelTransformation::LocalTransform;
  using irtkMultiLevelTransformation::Transform;
  using irtkMultiLevelTransformation::LocalInverse;
  using irtkMultiLevelTransformation::Inverse;
  using irtkMultiLevelTransformation::Displacement;
  using irtkMultiLevelTransformation::InverseDisplacement;

  /// Whether the caching of the transformation displacements is required
  /// (or preferred) by this transformation. For some transformations such as
  /// those parameterized by velocities, caching of the displacements for
  /// each target voxel results in better performance or is needed for example
  /// for the scaling and squaring method.
  virtual bool RequiresCachingOfDisplacements() const;

  /// Transforms a single point using the global transformation component only
  virtual void GlobalTransform(double &, double &, double &, double = 0, double = 1) const;

  /// Transforms a single point using the local transformation component only
  virtual void LocalTransform(int, int, double &, double &, double &, double = 0, double = 1) const;

  /// Transforms a single point
  virtual void Transform(int, int, double &, double &, double &, double = 0, double = 1) const;

  /// Transforms a single point using the inverse of the global transformation only
  virtual void GlobalInverse(double &, double &, double &, double = 0, double = 1) const;

  /// Transforms a single point using the inverse of the local transformation only
  virtual bool LocalInverse(int, int, double &, double &, double &, double = 0, double = 1) const;

  /// Transforms a single point using the inverse of the transformation
  virtual bool Inverse(int, int, double &, double &, double &, double = 0, double = 1) const;

  /// Calculates the displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  virtual void Displacement(int, int, irtkGenericImage<double> &, double, double = 1, const irtkWorldCoordsImage * = NULL) const;

  /// Calculates the displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  virtual void Displacement(int, int, irtkGenericImage<float> &, double, double = 1, const irtkWorldCoordsImage * = NULL) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  ///
  /// \returns Always zero.
  virtual int InverseDisplacement(int, int, irtkGenericImage<double> &, double, double = 1, const irtkWorldCoordsImage * = NULL) const;

  /// Calculates the inverse displacement vectors for a whole image domain
  ///
  /// \attention The displacements are computed at the positions after applying the
  ///            current displacements at each voxel. These displacements are then
  ///            added to the current displacements. Therefore, set the input
  ///            displacements to zero if only interested in the displacements of
  ///            this transformation at the voxel positions.
  ///
  /// \returns Always zero.
  virtual int InverseDisplacement(int, int, irtkGenericImage<float> &, double, double = 1, const irtkWorldCoordsImage * = NULL) const;

  // ---------------------------------------------------------------------------
  // Derivatives
  using irtkMultiLevelTransformation::ParametricGradient;

  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  virtual void ParametricGradient(const irtkGenericImage<double> *, double *,
                                  const irtkWorldCoordsImage *,
                                  const irtkWorldCoordsImage *,
                                  double = 1, double = 1) const;


  /// Applies the chain rule to convert spatial non-parametric gradient
  /// to a gradient w.r.t the parameters of this transformation.
  virtual void ParametricGradient(const irtkGenericImage<double> **, int, double *,
                                  const irtkWorldCoordsImage *,
                                  const irtkWorldCoordsImage *,
                                  const double * = NULL, double = 1) const;

  // ---------------------------------------------------------------------------
  // I/O
  using irtkMultiLevelTransformation::Read;
  using irtkMultiLevelTransformation::Write;

  /// Prints information about the transformation
  virtual void Print(irtkIndent = 0) const;

  /// Reads a transformation from a file stream
  virtual irtkCifstream &Read(irtkCifstream &);

  /// Writes a transformation to a file stream
  virtual irtkCofstream &Write(irtkCofstream &) const;

};


#endif
