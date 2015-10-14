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

#ifndef _IRTKIMAGEEDGEFORCE_H
#define _IRTKIMAGEEDGEFORCE_H

#include <irtkExternalForce.h>


/**
 *
 */
class irtkImageEdgeForce : public irtkExternalForce
{
  irtkObjectMacro(irtkImageEdgeForce);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Standard deviation of Gaussian smoothing kernel
  irtkPublicAttributeMacro(double, Sigma);

  /// Whether to project edge field gradient onto face normal
  irtkPublicAttributeMacro(bool, InNormalDirection);

  /// Edge field
  irtkAttributeMacro(irtkRealImage, EdgeField);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Copy attributes of this class from another instance
  void Copy(const irtkImageEdgeForce &);

public:

  /// Constructor
  irtkImageEdgeForce(const char * = "", double = 1.0);

  /// Copy constructor
  irtkImageEdgeForce(const irtkImageEdgeForce &);

  /// Assignment operator
  irtkImageEdgeForce &operator =(const irtkImageEdgeForce &);

  /// Destructor
  virtual ~irtkImageEdgeForce();

  // ---------------------------------------------------------------------------
  // Configuration

  // Import other overloads
  using irtkExternalForce::Parameter;

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  /// Get parameter name/value pairs
  virtual irtkParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Initialization

  /// Initialize external force once input and parameters have been set
  virtual void Initialize();

  // ---------------------------------------------------------------------------
  // Evaluation

protected:

  /// Evaluate external force
  virtual void EvaluateGradient(double *, double, double);

};


#endif
