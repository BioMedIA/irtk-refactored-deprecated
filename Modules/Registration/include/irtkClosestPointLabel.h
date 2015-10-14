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

#ifndef _IRTKCLOSESTPOINTLABEL_H
#define _IRTKCLOSESTPOINTLABEL_H

#include <irtkClosestPoint.h>


/**
 * Closest point correspondence map
 *
 * This class establishes point correspondence based on a known parcellation
 * of the surfaces, i.e., the discrete "Labels" point data array
 * (cf. IRTK polydataassignlabels tool). A source point is said to correspond
 * to a target point if it is the closest point on the source surface which
 * has the same label as the target point. If no such point is found,
 * the target point is marked as outlier.
 */

class irtkClosestPointLabel : public irtkClosestPoint
{
  irtkObjectMacro(irtkClosestPointLabel);

  typedef map<int, int> LabelToComponentIndexMap;

  // ---------------------------------------------------------------------------
  // Attributes

  /// IDs of points on the target surface closest to each respective surface label
  irtkAttributeMacro(vtkSmartPointer<vtkIdTypeArray>, TargetIds);

  /// Map of target surface label to target ID array component index
  irtkAttributeMacro(LabelToComponentIndexMap, TargetComponent);

  /// IDs of points on the source surface closest to each respective surface label
  irtkAttributeMacro(vtkSmartPointer<vtkIdTypeArray>, SourceIds);

  /// Map of target surface label to target ID array component index
  irtkAttributeMacro(LabelToComponentIndexMap, SourceComponent);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  irtkClosestPointLabel();

  /// Copy constructor
  irtkClosestPointLabel(const irtkClosestPointLabel &);

  /// Copy construct a new instance
  virtual irtkPointCorrespondence *NewInstance() const;

  /// Destructor
  virtual ~irtkClosestPointLabel();

  /// Type enumeration value
  virtual TypeId Type() const;

  // ---------------------------------------------------------------------------
  // Correspondences

protected:

  /// Common (re-)initialization steps of this class
  /// \note Must be a non-virtual function!
  void Init();

public:

  /// Initialize correspondence map after input and parameters are set
  virtual void Initialize();

  /// Reinitialize correspondence map after change of input topology
  virtual void Reinitialize();

  /// Update correspondence map
  virtual void Update();

};


#endif
