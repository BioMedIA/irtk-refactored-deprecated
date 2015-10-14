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

#ifndef _IRTKROBUSTCLOSESTPOINT_H
#define _IRTKROBUSTCLOSESTPOINT_H

#include <irtkFuzzyCorrespondence.h>


/**
 * Improved (iterative) closest point (ICP) correspondences
 *
 * The point correspondences returned by an instance of this class implement
 * an improvied iterative closest point matching with outlier rejection, the
 * reference method used by
 *
 *   Chui and Rangarajan, "A new point matching algorithm for non-rigid registration",
 *   Computer Vision and Image Understanding, 89(2-3), pp. 114–141, 2003.
 *
 * for comparison to the proposed robust point matching (RPM) algorithm.
 * This algorithm is realized by irtkRobustPointMatch.
 *
 * \sa irtkFiducialRegistrationError
 */

class irtkRobustClosestPoint : public irtkFuzzyCorrespondence
{
  irtkObjectMacro(irtkRobustClosestPoint);

  // ---------------------------------------------------------------------------
  // Attributes
protected:

  /// Factor by which standard deviation of point distances is multiplied
  /// for outlier rejection during (iterative) closest point matching
  irtkPublicAttributeMacro(double, Sigma);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Default constructor
  irtkRobustClosestPoint();

  /// Construct correspondence map and initialize it
  irtkRobustClosestPoint(const irtkRegisteredPointSet *,
                         const irtkRegisteredPointSet *);

  /// Copy constructor
  irtkRobustClosestPoint(const irtkRobustClosestPoint &);

  /// Copy construct a new instance
  virtual irtkPointCorrespondence *NewInstance() const;

  /// Destructor
  virtual ~irtkRobustClosestPoint();

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
protected:

  /// (Re-)calculate weights of correspondence links
  virtual void CalculateWeights();

};


#endif
