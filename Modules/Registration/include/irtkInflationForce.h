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

#ifndef _IRTKINFLATIONFORCE_H
#define _IRTKINFLATIONFORCE_H

#include <irtkSurfaceConstraint.h>

#include <vtkSmartPointer.h>
#include <vtkDataArray.h>


/**
 * Surface inflation force
 *
 * The inflation is driven by the average convexity or concavity of a region.
 * That is, points which lie in concave regions move outwards over time,
 * while points in convex regions move inwards.
 *
 * Fischl et al., Cortical Surface-Based Analysis II: Inflation, Flattening,
 * and a Surface-Based Coordinate System. NeuroImage, 9(2), 195–207 (1999).
 */
class irtkInflationForce : public irtkSurfaceConstraint
{
  irtkObjectMacro(irtkInflationForce);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Weight of inflation force in convex regions
  irtkPublicAttributeMacro(double, ConvexityWeight);

  /// Weight of inflation force in concave regions
  irtkPublicAttributeMacro(double, ConcavityWeight);

  /// Sum of concavity/convexity weights for each node
  irtkAttributeMacro(vtkSmartPointer<vtkDataArray>, Norm);

  /// Precomputed energy value (side effect of computing weight norms
  irtkAttributeMacro(double, EnergyValue);

  /// Copy attributes of this class from another instance
  void Copy(const irtkInflationForce &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  irtkInflationForce(const char * = "", double = 1.0);

  /// Copy constructor
  irtkInflationForce(const irtkInflationForce &);

  /// Assignment operator
  irtkInflationForce &operator =(const irtkInflationForce &);

  /// Destructor
  virtual ~irtkInflationForce();

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Initialize internal force term
  virtual void Initialize();

  /// Reinitialize internal force term after change of input topology
  virtual void Reinitialize();

  /// Update internal force data structures
  virtual void Update(bool);

protected:

  /// Common (re-)initialization code of this class only (non-virtual function!)
  void Init();

  /// Compute penalty for current transformation estimate
  virtual double Evaluate();

  /// Compute internal force w.r.t. transformation parameters
  virtual void EvaluateGradient(double *, double, double);

};


#endif
