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

#ifndef _IRTKTRANSFORMATIONAPPROXIMATIONERROR_H
#define _IRTKTRANSFORMATIONAPPROXIMATIONERROR_H

#include <irtkObjectiveFunction.h>

#include <irtkImage.h>
#include <irtkTransformation.h>


/**
 * Mean-squared-error of transformation approximation
 *
 * This objective function is minimized by irtkTransformation::ApproximateDOFs
 * to find the transformation parameters which minimize the mean squared error
 * of the approximation.
 */
class irtkTransformationApproximationError : public irtkObjectiveFunction
{
  irtkObjectMacro(irtkTransformationApproximationError);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Transformation
  irtkReadOnlyAggregateMacro(irtkTransformation, Transformation);

  /// Number of distinct time points
  irtkReadOnlyAttributeMacro(int, NumberOfTimePoints);

  /// Number of points
  irtkReadOnlyAttributeMacro(int, NumberOfPoints);

  /// Centroid of target points
  irtkReadOnlyAttributeMacro(irtkPoint, TargetCenter);

  /// Centroid of source points
  irtkReadOnlyAttributeMacro(irtkPoint, SourceCenter);

protected:

  vector<irtkPointSet>  _Target;
  vector<double>        _TargetTime;
  vector<irtkPointSet>  _Current;
  vector<irtkPointSet>  _Source;
  irtkVector3D<double> *_Gradient;

public:

  // ---------------------------------------------------------------------------
  // Construction/destruction

  /// Constructor
  irtkTransformationApproximationError(irtkTransformation *,
      const double *, const double *, const double *, const double *,
      const double *, const double *, const double *, int);

  /// Destructor
  virtual ~irtkTransformationApproximationError();

  /// Subtract centroid from each point set
  void CenterPoints();

  // ---------------------------------------------------------------------------
  // Function parameters (DoFs)

  /// Get number of DoFs
  ///
  /// \returns Number of free function parameters.
  virtual int NumberOfDOFs() const;

  /// Set function parameter values
  ///
  /// This is function can be used to set the parameters of the objective function
  /// to particular values. In particular, it can be used to restore the function
  /// parameters after a failed incremental update which did not result in the
  /// desired improvement.
  ///
  /// \param[in] x Function parameter (DoF) values.
  virtual void Put(const double *x);

  /// Get function parameter value
  ///
  /// \param[in] i Function parameter (DoF) index.
  ///
  /// \returns Value of specified function parameter (DoF).
  virtual double Get(int i) const;

  /// Get function parameter values
  ///
  /// This function can be used to store a backup of the current funtion parameter
  /// values before an update such that these can be restored using the Put
  /// member function if the update did not result in the desired change of the
  /// overall objective function value.
  ///
  /// \param[in] x Function parameter (DoF) values.
  virtual void Get(double *x) const;

  /// Add change (i.e., scaled gradient) to each parameter value
  ///
  /// This function updates each DoF of the objective function given a vector
  /// of corresponding changes, i.e., the computed gradient of the objective
  /// function w.r.t. these parameters or a desired change computed otherwise.
  ///
  /// \param[in] dx Change of each function parameter (DoF) as computed by the
  ///               Gradient member function and scaled by a chosen step length.
  ///
  /// \returns Maximum change of function parameter.
  virtual double Step(const double *dx);

  /// Update internal state after change of parameters
  ///
  /// \param[in] gradient Update also internal state required for evaluation of
  ///                     gradient of objective function.
  virtual void Update(bool gradient = true);

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Evaluate objective function value
  virtual double Value();

  /// Evaluate gradient of objective function w.r.t its DoFs
  ///
  /// \param[int] step    Step length for finite differences.
  /// \param[out] dx      Gradient of objective function.
  /// \param[out] sgn_chg Whether function parameter value is allowed to
  ///                     change sign when stepping along the computed gradient.
  virtual void Gradient(double *dx, double step = .0, bool *sgn_chg = NULL);

  /// Compute norm of gradient of objective function
  ///
  /// This norm is used to define a unit for the step length used by gradient
  /// descent methods. It is, for example, the maximum absolute value norm for
  /// linear transformations and the maximum control point displacement for FFDs.
  /// The computation of the norm may be done after conjugating the gradient
  /// vector obtained using the Gradient member function.
  ///
  /// \param[in] dx Gradient of objective function.
  virtual double GradientNorm(const double *dx) const;

};


#endif
