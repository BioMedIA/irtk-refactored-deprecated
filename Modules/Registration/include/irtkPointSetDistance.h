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

#ifndef _IRTKPOINTSETDISTANCE_H
#define _IRTKPOINTSETDISTANCE_H

#include <irtkDataFidelity.h>
#include <irtkRegisteredPointSet.h>

/**
 * Base class for point set distance measures
 */
class irtkPointSetDistance : public irtkDataFidelity
{
  irtkAbstractMacro(irtkPointSetDistance);

  // ---------------------------------------------------------------------------
  // Types
public:

  /// Type of gradient w.r.t a single transformed data point
  typedef irtkVector3D<double> GradientType;

  // ---------------------------------------------------------------------------
  // Attributes
protected:

  /// First point set
  irtkPublicAggregateMacro(irtkRegisteredPointSet, Target);

  /// Second point set
  irtkPublicAggregateMacro(irtkRegisteredPointSet, Source);

  /// Memory for (non-parametric) gradient w.r.t points of target
  irtkComponentMacro(GradientType, GradientWrtTarget);

  /// Memory for (non-parametric) gradient w.r.t points of source
  irtkComponentMacro(GradientType, GradientWrtSource);

  /// Whether Update has not been called since initialization
  irtkAttributeMacro(bool, InitialUpdate);

  /// Allocate memory for (non-parametric) gradient w.r.t points of target
  void AllocateGradientWrtTarget(int);

  /// Allocate memory for (non-parametric) gradient w.r.t points of source
  void AllocateGradientWrtSource(int);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
protected:

  /// Constructor
  irtkPointSetDistance(const char * = "", double = 1.0);

  /// Copy constructor
  irtkPointSetDistance(const irtkPointSetDistance &, int = -1, int = -1);

  /// Assignment operator
  irtkPointSetDistance &operator =(const irtkPointSetDistance &);

  /// Copy attributes from other point set distance measure
  void Copy(const irtkPointSetDistance &, int = -1, int = -1);

public:

  /// Instantiate specified similarity measure
  static irtkPointSetDistance *New(irtkPointSetDistanceMeasure);

  /// Destructor
  virtual ~irtkPointSetDistance();

  // ---------------------------------------------------------------------------
  // Initialization

protected:

  /// Initialize distance measure once input and parameters have been set
  void Initialize(int, int);

  /// Reinitialize distance measure after change of input topology
  ///
  /// This function is called in particular when an input surface has been
  /// reparameterized, e.g., by a local remeshing filter.
  void Reinitialize(int, int);

public:

  /// Initialize distance measure once input and parameters have been set
  virtual void Initialize();

  /// Reinitialize distance measure after change of input topology
  ///
  /// This function is called in particular when an input surface has been
  /// reparameterized, e.g., by a local remeshing filter.
  virtual void Reinitialize();

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Update moving input points and internal state of distance measure
  virtual void Update(bool = true);

protected:

  /// Compute non-parametric gradient w.r.t points of given data set
  ///
  /// \param[in]  target   Transformed point set.
  /// \param[out] gradient Non-parametric gradient of polydata distance measure.
  virtual void NonParametricGradient(const irtkRegisteredPointSet *target,
                                     GradientType                 *gradient) = 0;

  /// Convert non-parametric gradient of polydata distance measure into
  /// gradient w.r.t transformation parameters
  ///
  /// This function calls irtkTransformation::ParametricGradient of the
  /// transformation to apply the chain rule in order to obtain the gradient
  /// of the distance measure w.r.t the transformation parameters.
  /// It adds the weighted gradient to the final registration energy gradient.
  ///
  /// \param[in]     target      Transformed point set.
  /// \param[in]     np_gradient Point-wise non-parametric gradient.
  /// \param[in,out] gradient    Gradient to which the computed parametric gradient
  ///                            is added, after multiplication by the given \p weight.
  /// \param[in]     weight      Weight of polydata distance measure.
  virtual void ParametricGradient(const irtkRegisteredPointSet *target,
                                  const GradientType           *np_gradient,
                                  double                       *gradient,
                                  double                        weight);

  /// Evaluate gradient of polydata distance measure
  ///
  /// This function calls the virtual NonParametricGradient function to be
  /// implemented by subclasses for each transformed input data set to obtain
  /// the gradient of the polydata distance measure. It then converts this gradient
  /// into a gradient w.r.t the transformation parameters using the ParametricGradient.
  ///
  /// If both target and source data sets are transformed by different transformations,
  /// the resulting gradient vector contains first the derivative values w.r.t
  /// the parameters of the target transformation followed by those computed
  /// w.r.t the parameters of the source transformation. If both data sets are
  /// transformed by the same transformation, the sum of the derivative values
  /// is added to the resulting gradient vector. This is in particular the case
  /// for a velocity based transformation model which is applied to deform both
  /// data sets "mid-way". Otherwise, only one input data set is transformed
  /// (usually the target) and the derivative values of only the respective
  /// transformation parameters added to the gradient vector.
  ///
  /// \sa NonParametricGradient, ParametricGradient
  ///
  /// \param[in,out] gradient Gradient to which the computed gradient of the
  ///                         polydata similarity is added after multiplying by
  ///                         the given similarity \p weight.
  /// \param[in]     step     Step length for finite differences (unused).
  /// \param[in]     weight   Weight of polydata similarity.
  virtual void EvaluateGradient(double *gradient, double step, double weight);

  // ---------------------------------------------------------------------------
  // Debugging
public:

  /// Write input of data fidelity term
  virtual void WriteDataSets(const char *, const char *, bool = true) const;

  /// Write gradient of data fidelity term w.r.t each transformed input
  virtual void WriteGradient(const char *, const char *) const;

protected:

  /// Write gradient of data fidelity term w.r.t each transformed input
  virtual void WriteGradient(const char *,
                             const irtkRegisteredPointSet *,
                             const GradientType *,
                             const vector<int> * = NULL) const;

};


#endif
