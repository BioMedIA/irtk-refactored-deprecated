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

#ifndef _IRTKCLOSESTCELL_H
#define _IRTKCLOSESTCELL_H

#include <irtkPointCorrespondence.h>
#include <irtkPointSet.h>


/**
 * Closest point on line/surface correspondence map
 */

class irtkClosestCell : public irtkPointCorrespondence
{
  irtkObjectMacro(irtkClosestCell);

  // ---------------------------------------------------------------------------
  // Types
public:

  /// Enumeration value of supported cell locators
  enum LocatorType { Default, CellTree, BSPTree, OBBTree };

  // ---------------------------------------------------------------------------
  // Attributes
protected:

  /// Type of point locator
  irtkPublicAttributeMacro(enum LocatorType, LocatorType);

  /// Number of cells per node
  irtkPublicAttributeMacro(int, NumberOfCellsPerNode);

  /// Mark correspondences with distance greater than the mean distance plus
  /// _Sigma times the standard deviation of the current point distances as outliers
  irtkPublicAttributeMacro(double, Sigma);

  /// Mark correspondences with distance greater than this as outliers
  /// If _Sigma is positive, it is used to set this attribute automatically
  irtkPublicAttributeMacro(double, MaxDistance);

  /// Target points corresponding to source points
  irtkAttributeMacro(irtkPointSet, TargetPoints);

  /// Distance of target points from source points
  irtkAttributeMacro(vector<double>, TargetDistance);

  /// Source points corresponding to target points
  irtkAttributeMacro(irtkPointSet, SourcePoints);

  /// Distance of source points from target points
  irtkAttributeMacro(vector<double>, SourceDistance);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  irtkClosestCell();

  /// Copy constructor
  irtkClosestCell(const irtkClosestCell &);

  /// Copy construct a new instance
  virtual irtkPointCorrespondence *NewInstance() const;

  /// Destructor
  virtual ~irtkClosestCell();

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
