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

#ifndef _IRTKROBUSTPOINTMATCH_H
#define _IRTKROBUSTPOINTMATCH_H

#include <irtkFuzzyCorrespondence.h>


/**
 * Robust point match (RPM) correspondences
 *
 * The point correspondences returned by an instance of this class correspond
 * to the cluster centers of the Robust Point Matching (RPM) algorithm of
 *
 *   Chui and Rangarajan, "A new point matching algorithm for non-rigid registration",
 *   Computer Vision and Image Understanding, 89(2-3), pp. 114–141, 2003.
 *
 * where the smoothness regularization term is separate. Therefore, the weight
 * of the smoothness term is not directly coupled to the annealing temperature
 * used by this class, unlike the original RPM proposal.
 *
 * \sa irtkFiducialRegistrationError
 */

class irtkRobustPointMatch : public irtkFuzzyCorrespondence
{
  irtkObjectMacro(irtkRobustPointMatch);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Initial temperature of the deterministic annealing process
  irtkPublicAttributeMacro(double, InitialTemperature);

  /// Annealing rate by which the temperature is multiplied or number of NN
  ///
  /// If set to a negative value, the average mean distance to the
  /// corresponding number of nearest neighbors is used as next temperature.
  irtkPublicAttributeMacro(double, AnnealingRate);

  /// Final temperature which stops the deterministic annealing process
  irtkPublicAttributeMacro(double, FinalTemperature);

  /// Current temperature of the deterministic annealing process
  irtkPublicAttributeMacro(double, Temperature);

  /// Variance of extra features (resp. of their differences)
  irtkPublicAttributeMacro(double, VarianceOfFeatures);

  /// Cluster for outliers in target (i.e., centroid of source points!)
  irtkReadOnlyAttributeMacro(irtkPoint, TargetOutlierCluster);

  /// Cluster for outliers in source (i.e., centroid of target points!)
  irtkReadOnlyAttributeMacro(irtkPoint, SourceOutlierCluster);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Default constructor
  irtkRobustPointMatch();

  /// Construct correspondence map and initialize it
  irtkRobustPointMatch(const irtkRegisteredPointSet *,
                       const irtkRegisteredPointSet *);

  /// Copy constructor
  irtkRobustPointMatch(const irtkRobustPointMatch &);

  /// Copy construct a new instance
  virtual irtkPointCorrespondence *NewInstance() const;

  /// Destructor
  virtual ~irtkRobustPointMatch();

  /// Type enumeration value
  virtual TypeId Type() const;

  // ---------------------------------------------------------------------------
  // Parameters

  // Import other overloads
  using irtkFuzzyCorrespondence::Parameter;

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  /// Get parameter key/value as string map
  virtual irtkParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Correspondences

  /// Initialize correspondence map
  virtual void Initialize();

  /// Update correspondence map after convergence
  virtual bool Upgrade();

protected:

  /// Initialize annealing process
  virtual void InitializeAnnealing();

  /// (Re-)calculate weights of correspondence links
  virtual void CalculateWeights();

};


#endif
