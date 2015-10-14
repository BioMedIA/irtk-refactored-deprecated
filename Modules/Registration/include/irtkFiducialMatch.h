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

#ifndef _IRTKFIDUCIALMATCH_H
#define _IRTKFIDUCIALMATCH_H

#include <irtkPointCorrespondence.h>


/**
 * Point correspondence based on the fiducial indices
 */

class irtkFiducialMatch : public irtkPointCorrespondence
{
  irtkObjectMacro(irtkFiducialMatch);

  /// Optional map of corresponding target point indices
  irtkPublicAttributeMacro(vector<int>, TargetIndex);

  /// Optional map of corresponding source point indices
  irtkPublicAttributeMacro(vector<int>, SourceIndex);

  /// Name of text file listing corresponding point indices
  irtkPublicAttributeMacro(string, CorrespondenceMap);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  irtkFiducialMatch();

  /// Copy constructor
  irtkFiducialMatch(const irtkFiducialMatch &);

  /// Copy construct a new instance
  virtual irtkPointCorrespondence *NewInstance() const;

  /// Destructor
  virtual ~irtkFiducialMatch();

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
protected:

  /// Check if correspondence map has a valid entry for every target point
  void ValidateCorrespondenceMap(const irtkRegisteredPointSet *,
                                 const irtkRegisteredPointSet *,
                                 const vector<int> &,
                                 const char * = NULL) const;

  /// Compute inverse correspondence map
  void InvertCorrespondenceMap(const irtkRegisteredPointSet *,
                               const irtkRegisteredPointSet *,
                               const vector<int> &, vector<int> &) const;

public:

  /// Initialize correspondence map
  virtual void Initialize();

  /// Get untransformed target point corresponding to i-th source (sample) point
  virtual bool GetInputTargetPoint(int, irtkPoint &) const;

  /// Get untransformed source point corresponding to i-th target (sample) point
  virtual bool GetInputSourcePoint(int, irtkPoint &) const;

  /// Get (transformed) target point corresponding to i-th source (sample) point
  virtual bool GetTargetPoint(int, irtkPoint &) const;

  /// Get (transformed) source point corresponding to i-th target (sample) point
  virtual bool GetSourcePoint(int, irtkPoint &) const;

  /// Get index of target point corresponding to i-th source (sample) point
  ///
  /// \returns Index of corresponding target point and -1 if point is an outlier
  ///          or its corresponding point is not a target vertex position.
  virtual int GetTargetIndex(int) const;

  /// Get index of source point corresponding to i-th target (sample) point
  ///
  /// \returns Index of corresponding source point and -1 if point is an outlier
  ///          or its corresponding point is not a source vertex position.
  virtual int GetSourceIndex(int) const;

};


#endif
