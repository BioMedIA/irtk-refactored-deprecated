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

#ifndef _IRTKGENERICREGISTRATIONDEBUGGER_H
#define _IRTKGENERICREGISTRATIONDEBUGGER_H

#include <irtkObserver.h>

class irtkImageSimilarity;
class irtkGenericRegistrationFilter;


/**
 * Writes intermediate registration data to the current working directory
 *
 * Usage:
 * \code
 * irtkImageRegistrationFilter   registration;
 * irtkImageRegistrationDebugger debugger;
 * registration.AddObserver(debugger);
 * \endcode
 */
class irtkGenericRegistrationDebugger : public irtkObserver
{
  irtkObjectMacro(irtkGenericRegistrationDebugger);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Prefix for output file names
  irtkPublicAttributeMacro(string, Prefix);

  /// Whether to use level specific file name prefix
  irtkPublicAttributeMacro(bool, LevelPrefix);

  int _Level;                                 ///< Current level
  int _Iteration;                             ///< Current iteration
  int _LineIteration;                         ///< Current line iteration
  vector<irtkImageSimilarity *> _Similarity;  ///< Similarity terms

  /// Reference to the registration filter object
  irtkGenericRegistrationFilter *_Registration;

  // ---------------------------------------------------------------------------
  // Construction/Destruction
private:

  /// Copy construction
  /// \note Intentionally not implemented.
  irtkGenericRegistrationDebugger(const irtkGenericRegistrationDebugger &);

  /// Assignment operator
  /// \note Intentionally not implemented.
  irtkGenericRegistrationDebugger &operator =(const irtkGenericRegistrationDebugger &);

public:

  /// Constructor
  irtkGenericRegistrationDebugger(const char * = "");

  /// Destructor
  ~irtkGenericRegistrationDebugger();

  /// Handle event and print message to output stream
  void HandleEvent(irtkObservable *, irtkEvent, const void *);

};


#endif
