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

#ifndef _IRTKSPECTRALMATCH_H
#define _IRTKSPECTRALMATCH_H

#include <irtkPointCorrespondence.h>


/**
 * Correspondence map based on diffeomorphic spectral matching
 */

class irtkSpectralMatch : public irtkPointCorrespondence
{
  irtkObjectMacro(irtkSpectralMatch);

  // ---------------------------------------------------------------------------
  // Attributes
protected:

  /// Number of spectral nearest neighbors to consider
  irtkPublicAttributeMacro(int, NumberOfNeighbors);

  /// Standard deviation of Gaussian interpolation kernel
  irtkPublicAttributeMacro(double, Sigma);

  /// Target points corresponding to source points
  irtkAttributeMacro(irtkPointSet, TargetPoints);

  /// Source points corresponding to target points
  irtkAttributeMacro(irtkPointSet, SourcePoints);

  irtkPointLocator *_TargetLocator;
  irtkPointLocator *_SourceLocator;

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  irtkSpectralMatch();

  /// Copy constructor
  irtkSpectralMatch(const irtkSpectralMatch &);

  /// Copy construct a new instance
  virtual irtkPointCorrespondence *NewInstance() const;

  /// Destructor
  virtual ~irtkSpectralMatch();

  /// Type enumeration value
  virtual TypeId Type() const;

  // ---------------------------------------------------------------------------
  // Parameters

  // Import other overloads
  using irtkPointCorrespondence::Parameter;

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  /// Get parameter key/value as string map
  virtual irtkParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Correspondences

  /// Initialize correspondence map
  virtual void Initialize();

  /// Update correspondence map
  virtual void Update();

  /// Update correspondence map after convergence
  virtual bool Upgrade();

  /// Get untransformed target point corresponding to i-th source (sample) point
  virtual bool GetInputTargetPoint(int, irtkPoint &) const;

  /// Get untransformed source point corresponding to i-th target (sample) point
  virtual bool GetInputSourcePoint(int, irtkPoint &) const;

  /// Get (transformed) target point corresponding to i-th source (sample) point
  virtual bool GetTargetPoint(int, irtkPoint &) const;

  /// Get (transformed) source point corresponding to i-th target (sample) point
  virtual bool GetSourcePoint(int, irtkPoint &) const;

};


#endif
