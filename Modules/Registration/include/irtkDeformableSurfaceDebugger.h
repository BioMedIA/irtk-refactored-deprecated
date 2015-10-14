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

#ifndef _IRTKDEFORMABLESURFACEDEBUGGER_H
#define _IRTKDEFORMABLESURFACEDEBUGGER_H

#include <irtkObserver.h>

class irtkImageSimilarity;
class irtkDeformableSurfaceModel;


/**
 * Writes intermediate surfaces to the current working directory
 *
 * This object observes the optimizer of the deformable surface model and must
 * therefore be attached to the respective irtkLocalOptimizer instance
 * (typically irtkEulerMethod or a subclass of it).
 */
class irtkDeformableSurfaceDebugger : public irtkObserver
{
  irtkObjectMacro(irtkDeformableSurfaceDebugger);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Prefix for output file names
  irtkPublicAttributeMacro(string, Prefix);

  /// Reference to the deformable surface model
  irtkPublicAggregateMacro(const irtkDeformableSurfaceModel, Model);

  /// Write intermediate results only every n gradient steps
  irtkPublicAttributeMacro(int, Interval);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
private:

  /// Copy construction
  /// \note Intentionally not implemented.
  irtkDeformableSurfaceDebugger(const irtkDeformableSurfaceDebugger &);

  /// Assignment operator
  /// \note Intentionally not implemented.
  irtkDeformableSurfaceDebugger &operator =(const irtkDeformableSurfaceDebugger &);

public:

  /// Constructor
  irtkDeformableSurfaceDebugger(const irtkDeformableSurfaceModel * = NULL, const char * = "");

  /// Destructor
  ~irtkDeformableSurfaceDebugger();

  /// Handle event and print message to output stream
  void HandleEvent(irtkObservable *, irtkEvent, const void *);

};


#endif
