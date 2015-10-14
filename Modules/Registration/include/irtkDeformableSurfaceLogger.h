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

#ifndef _IRTKDEFORMABLESURFACELOGGER_H
#define _IRTKDEFORMABLESURFACELOGGER_H

#include <irtkObserver.h>


/**
 * Prints progress of deformable surface to output stream
 */
class irtkDeformableSurfaceLogger : public irtkObserver
{
  irtkObjectMacro(irtkDeformableSurfaceLogger);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Verbosity level
  irtkPublicAttributeMacro(int, Verbosity);

  /// Output stream for progress report
  irtkPublicAggregateMacro(ostream, Stream);

  /// Whether to use SGR commands for colored terminal output
  irtkPublicAttributeMacro(bool, Color);

  /// Whether to flush stream buffer after each printed message
  irtkPublicAttributeMacro(bool, FlushBuffer);

  bool _Converged;             ///< Track whether line search converged
  int  _NumberOfIterations;    ///< Number of actual line search iterations
  int  _NumberOfSteps;         ///< Number of iterative line search steps
  int  _NumberOfGradientSteps; ///< Number of gradient descent steps

  // ---------------------------------------------------------------------------
  // Construction/Destruction
private:

  /// Copy construction
  /// \note Intentionally not implemented.
  irtkDeformableSurfaceLogger(const irtkDeformableSurfaceLogger &);

  /// Assignment operator
  /// \note Intentionally not implemented.
  irtkDeformableSurfaceLogger &operator =(const irtkDeformableSurfaceLogger &);

public:

  /// Constructor
  irtkDeformableSurfaceLogger(ostream * = &cout);

  /// Destructor
  virtual ~irtkDeformableSurfaceLogger();

  /// Handle event and print message to output stream
  void HandleEvent(irtkObservable *, irtkEvent, const void *);

};


#endif
