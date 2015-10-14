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

#ifndef _IRTKEDGELENGTHCONSTRAINT_H
#define _IRTKEDGELENGTHCONSTRAINT_H

#include <irtkPointSetConstraint.h>

#include <irtkEdgeTable.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkDataArray.h>


/**
 * Edge length / metric distortion constraint
 *
 * This spring force keeps the edge lengths of a simplicial complex either
 * similar to a specified average edge length or the initial length in the
 * undeformed complex. It can also be used to minimize metric distortion
 * while inflating a surface mesh.
 *
 * Lachaud and Montanvert, Deformable meshes with automated topology changes
 * for coarse-to-fine three-dimensional surface extraction. Medical Image Analysis,
 * 3(2), 187–207, doi:10.1016/S1361-8415(99)80012-7 (1999).
 *
 * Fischl et al.. Cortical Surface-Based Analysis II: Inflation, Flattening,
 * and a Surface-Based Coordinate System. NeuroImage, 9(2), 195–207 (1999).
 *
 * Park et al., A non-self-intersecting adaptive deformable surface for
 * complex boundary extraction from volumetric images, 25, 421–440 (2001).
 */
class irtkEdgeLengthConstraint : public irtkPointSetConstraint
{
  irtkObjectMacro(irtkEdgeLengthConstraint);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Length of edges at rest
  ///
  /// -1: Use average edge length of undeformed mesh
  /// -2: Use initial lengths of individual edges
  irtkPublicAttributeMacro(double, RestLength);

  /// Initial point positions
  irtkAttributeMacro(vtkSmartPointer<vtkPoints>, InitialPoints);

  /// Extracted edges
  irtkAttributeMacro(irtk::polydata::irtkEdgeTable, EdgeTable);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  irtkEdgeLengthConstraint(const char * = "", double = 1.0);

  /// Copy constructor
  irtkEdgeLengthConstraint(const irtkEdgeLengthConstraint &);

  /// Assignment operator
  irtkEdgeLengthConstraint &operator =(const irtkEdgeLengthConstraint &);

  /// Destructor
  virtual ~irtkEdgeLengthConstraint();

  // ---------------------------------------------------------------------------
  // Configuration

  // Import other overloads
  using irtkPointSetConstraint::Parameter;

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

protected:

  /// Common (re-)initialization steps of this internal force term
  /// \note Must be a non-virtual function!
  void Init();

  /// Compute penalty for current transformation estimate
  virtual double Evaluate();

  /// Compute internal force w.r.t. transformation parameters
  virtual void EvaluateGradient(double *, double, double);

};


#endif
