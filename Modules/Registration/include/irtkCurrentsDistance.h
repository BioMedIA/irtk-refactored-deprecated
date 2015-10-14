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

#ifndef _IRTKCURRENTSDISTANCE_H
#define _IRTKCURRENTSDISTANCE_H

#include <vtkSmartPointer.h>
#include <vtkPointSet.h>
#include <vtkPolyData.h>

#include <irtkPointSetDistance.h>


/**
 * Currents distance measure of
 * points (0-currents), curves (1-currents), or surfaces (2-currents)
 *
 */
class irtkCurrentsDistance : public irtkPointSetDistance
{
  irtkObjectMacro(irtkCurrentsDistance);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Current representation of target data set
  irtkAttributeMacro(vtkSmartPointer<vtkPolyData>, TargetCurrent);

  /// Current representation of source data set
  irtkAttributeMacro(vtkSmartPointer<vtkPolyData>, SourceCurrent);

  /// Sigma value of currents kernel
  irtkPublicAttributeMacro(double, Sigma);

  /// Whether to ensure symmetry of currents dot product
  irtkPublicAttributeMacro(bool, Symmetric);

  /// Sum of squared norm of fixed (i.e., untransformed) data set(s)
  irtkAttributeMacro(double, TargetNormSquared);

  // ---------------------------------------------------------------------------
  // Currents representation
protected:

  /// Convert data set to current
  static vtkSmartPointer<vtkPolyData> ToCurrent(vtkPointSet *);

  /// Convert surface mesh to current
  static vtkSmartPointer<vtkPolyData> SurfaceToCurrent(vtkPolyData *);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  irtkCurrentsDistance(const char * = "", double = 1.0);

  /// Copy constructor
  irtkCurrentsDistance(const irtkCurrentsDistance &);

  /// Assignment operator
  irtkCurrentsDistance &operator =(const irtkCurrentsDistance &);

  /// Destructor
  virtual ~irtkCurrentsDistance();

  // ---------------------------------------------------------------------------
  // Initialization

protected:

  /// Common (re-)initialization code of this class (must be non-virtual function!)
  void Init();

public:

  /// Initialize distance measure after input and parameters were set
  virtual void Initialize();

  /// Reinitialize distance measure after input topology changed
  virtual void Reinitialize();

  // ---------------------------------------------------------------------------
  // Parameters

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  // Import other overloads
  using irtkPointSetDistance::Parameter;

  /// Get parameter key/value as string map
  virtual irtkParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Update moving input points and internal state of distance measure
  virtual void Update(bool);

protected:

  /// Evaluate unweighted energy term
  virtual double Evaluate();

  /// Compute non-parametric gradient w.r.t points of given data set
  ///
  /// \param[in]  target   Transformed data set.
  /// \param[out] gradient Non-parametric gradient of polydata distance measure.
  virtual void NonParametricGradient(const irtkRegisteredPointSet *target,
                                     GradientType                 *gradient);

  // ---------------------------------------------------------------------------
  // Debugging
public:

  /// Write input of data fidelity term
  virtual void WriteDataSets(const char *, const char *, bool = true) const;

};


#endif
