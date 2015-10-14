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

#ifndef _IRTKPOINTCORRESPONDENCEDISTANCE_H
#define _IRTKPOINTCORRESPONDENCEDISTANCE_H

#include <irtkPointSetDistance.h>
#include <irtkPointCorrespondence.h>
#include <irtkRadialErrorFunction.h>

#include <irtkTestProd.h>
class irtkEventDelegate;


/**
 * Distance error of established/known point correspondences
 *
 * This distance term evaluates the residual registration error of each
 * transformed target point relative to its corresponding source point.
 * It is a generic point set distance measure which can be used for many
 * types of point sets, including in particular point clouds, curves, and
 * surface meshes.
 *
 * A corresponding point locator is utilized to find the corresponding point
 * in the source data set. This can be, for example, the closest point in the
 * source data set, the closest point on the source surface, or the matched
 * source data set point which is closest in terms of some other feature than
 * just Euclidean distance of the points. A special case of point correspondence
 * is one which remains fixed and maps the i-th point in the target data set to
 * the j-th point in the source data set. This is the case for manually labeled
 * fiducial markers or pre-computed correspondence maps using an external
 * point/surface matching algorithm.
 *
 * The Euclidean distance of corresponding points is further weighted by an
 * error function which may non-linearly penalize established correspondences
 * and reduce influence of outliers or incorrect matches in order to improve
 * the robustness of the point set registration.
 */
class irtkPointCorrespondenceDistance : public irtkPointSetDistance
{
  irtkObjectMacro(irtkPointCorrespondenceDistance);

  FRIEND_TEST(irtkPointCorrespondenceDistance, FuzzyMatch);

  // ---------------------------------------------------------------------------
  // Attributes
protected:

  /// Forwards correspondence map event messages to observers of energy term
  irtkEventDelegate _EventDelegate;

  /// Approximate distance between sampled target points
  irtkPublicAttributeMacro(double, TargetSampleDistance);

  /// Approximate distance between sampled source points
  irtkPublicAttributeMacro(double, SourceSampleDistance);

  /// Maximum number of target samples used
  irtkPublicAttributeMacro(int, NumberOfTargetSamples);

  /// Maximum number of target samples used
  irtkPublicAttributeMacro(int, NumberOfSourceSamples);

  /// Indices of sampled target points
  irtkAttributeMacro(vector<int>, TargetSample);

  /// Indices of sampled source points
  irtkAttributeMacro(vector<int>, SourceSample);

  /// Number of Update calls between reevaluation of correspondences
  irtkPublicAttributeMacro(int, UpdatePeriod);

  /// Number of invocations of Update modulo _UpdatePeriod
  irtkAttributeMacro(int, NumberOfUpdates);

  /// Point correspondence map
  irtkComponentMacro(irtkPointCorrespondence, Correspondence);

  /// Polydata registration error weight function
  irtkComponentMacro(irtkRadialErrorFunction, ErrorFunction);

  /// Whether to evaluate error of target points and corresponding source points
  irtkPublicAttributeMacro(bool, EvaluateTargetError);

  /// Whether to evaluate error of source points and corresponding target points
  irtkPublicAttributeMacro(bool, EvaluateSourceError);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Construct distance with given default point correspondence type and error function
  irtkPointCorrespondenceDistance(const char *, double,
                                  irtkPointCorrespondence *,
                                  irtkRadialErrorFunction * = NULL);

public:

  /// Constructor
  irtkPointCorrespondenceDistance(const char * = "", double = 1.0);

  /// Copy constructor
  irtkPointCorrespondenceDistance(const irtkPointCorrespondenceDistance &);

  /// Assignment operator
  irtkPointCorrespondenceDistance &operator =(const irtkPointCorrespondenceDistance &);

  /// Destructor
  virtual ~irtkPointCorrespondenceDistance();

  // ---------------------------------------------------------------------------
  // Parameters

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  // Do not hide other overload
  using irtkPointSetDistance::Parameter;

  /// Get parameter key/value as string map
  virtual irtkParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Initialization/update

protected:

  /// Sample point sets
  void SamplePoints();

public:

  /// Initialize error measure once input and parameters have been set
  virtual void Initialize();

  /// Reinitialize error measure after change of input topology
  virtual void Reinitialize();

  /// Update moving input points and internal state of distance measure
  virtual void Update(bool = true);

  /// Update energy term after convergence
  virtual bool Upgrade();

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Whether to evaluate target to source error
  bool DoEvaluateTargetError() const;

  /// Whether to evaluate source to target error
  bool DoEvaluateSourceError() const;

protected:

  /// Evaluate unweighted energy term
  virtual double Evaluate();

  /// Compute non-parametric gradient w.r.t the given point set
  ///
  /// \param[in]  source   Set of transformed fiducial points.
  /// \param[out] gradient Non-parametric fiducial registration error gradient.
  virtual void NonParametricGradient(const irtkRegisteredPointSet *source,
                                     GradientType                 *gradient);

  /// Convert non-parametric gradient of polydata distance measure into
  /// gradient w.r.t transformation parameters
  ///
  /// This function calls irtkTransformation::ParametricGradient of the
  /// transformation to apply the chain rule in order to obtain the gradient
  /// of the distance measure w.r.t the transformation parameters.
  /// It adds the weighted gradient to the final registration energy gradient.
  ///
  /// \param[in]    source      Transformed point set.
  /// \param[in]    np_gradient Point-wise non-parametric gradient.
  /// \param[inout] gradient    Gradient to which the computed parametric gradient
  ///                           is added, after multiplication by the given \p weight.
  /// \param[in]    weight      Weight of polydata distance measure.
  virtual void ParametricGradient(const irtkRegisteredPointSet *source,
                                  const GradientType           *np_gradient,
                                  double                       *gradient,
                                  double                        weight);

  /// Evaluate gradient of point distance measure
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
  void EvaluateGradient(double *gradient, double step, double weight);

  /// Forward point correspondence map event
  void ForwardEvent(irtkObservable *, irtkEvent, const void *);

  // ---------------------------------------------------------------------------
  // Debugging
public:

  /// Write input of data fidelity term
  virtual void WriteDataSets(const char *, const char *, bool = true) const;

  /// Write gradient of data fidelity term w.r.t each transformed input
  virtual void WriteGradient(const char *, const char *) const;

protected:

  /// Write given input data set to specified file
  virtual void WriteDataSet(const char *,
                            const irtkRegisteredPointSet *,
                            const vector<int> &,
                            const irtkPointCorrespondence *) const;

};


#endif
