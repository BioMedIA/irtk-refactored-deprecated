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

#ifndef _IRTKNONSELFINTERSECTINGSURFACECONSTRAINT_H
#define _IRTKNONSELFINTERSECTINGSURFACECONSTRAINT_H

#include <irtkSurfaceConstraint.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkDataArray.h>


/**
 * Non-self intersecting surface forces
 *
 * This internal force makes nearby non-adjacent faces to repel each other
 * in order to avoid self-intersection of the surface.
 *
 * Park et al., A non-self-intersecting adaptive deformable surface for
 * complex boundary extraction from volumetric images, 25, 421–440 (2001).
 */
class irtkNonSelfIntersectionConstraint : public irtkSurfaceConstraint
{
  irtkObjectMacro(irtkNonSelfIntersectionConstraint);

  // ---------------------------------------------------------------------------
  // Types

public:

  /// Pair of candidate triangles which are too close to each other
  struct CellPair
  {
    vtkIdType _CellId1;
    vtkIdType _CellId2;
    double    _Point1[3];
    double    _Point2[3];
    double    _Distance;
  };

  // ---------------------------------------------------------------------------
  // Attributes

protected:

  /// Minimum distance
  /// Set to non-positive value to use initial average edge length.
  irtkPublicAttributeMacro(double, MinDistance);

  /// Candidates for self-intersection forces
  irtkAttributeMacro(vector<CellPair>, Candidates);

  /// Computed vertex normals
  irtkAttributeMacro(vtkSmartPointer<vtkDataArray>, Normals);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  irtkNonSelfIntersectionConstraint(const char * = "", double = 1.0);

  /// Copy constructor
  irtkNonSelfIntersectionConstraint(const irtkNonSelfIntersectionConstraint &);

  /// Assignment operator
  irtkNonSelfIntersectionConstraint &operator =(const irtkNonSelfIntersectionConstraint &);

  /// Destructor
  virtual ~irtkNonSelfIntersectionConstraint();

  // ---------------------------------------------------------------------------
  // Configuration

  // Import other overloads
  using irtkSurfaceConstraint::Parameter;

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  /// Get parameter name/value pairs
  virtual irtkParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Initialize internal force term once input and parameters have been set
  virtual void Initialize();

  /// Reinitialize internal force term after change of input topology
  virtual void Reinitialize();

  /// Update internal state upon change of input
  virtual void Update(bool = true);

protected:

  /// Common (re-)initialization steps of this class (non-virtual function!)
  void Init();

  /// Compute penalty for current transformation estimate
  virtual double Evaluate();

  /// Compute internal force w.r.t. transformation parameters
  virtual void EvaluateGradient(double *, double, double);

};


#endif
