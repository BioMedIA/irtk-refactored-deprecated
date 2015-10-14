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

#ifndef _IRTKIMAGEFUNCTION_H
#define _IRTKIMAGEFUNCTION_H

#include <irtkImage.h>


/**
 * Abstract base class for any general image interpolation function filter.
 *
 * This is the abstract base class which defines a common interface for all
 * filters which take an image as input and sample that image at arbitrary
 * location. Each derived class has to implement all abstract member functions.
 */

class irtkImageFunction : public irtkObject
{
  irtkAbstractMacro(irtkImageFunction);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Default value to return
  irtkPublicAttributeMacro(double, DefaultValue);

  /// Debugging flag
  /// \todo Use global \c debug flag instead.
  irtkPublicAttributeMacro(bool, DebugFlag);

protected:

  /// Input image for filter
  irtkImage *_input;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  irtkImageFunction();

  /// Copy constructor
  irtkImageFunction(const irtkImageFunction &);

public:

  /// Destructor
  virtual ~irtkImageFunction();

  // ---------------------------------------------------------------------------
  // Initialization

  /// Set input image for filter
  virtual void SetInput(irtkImage *);

  /// Get input image for filter
  irtkImage *GetInput();

  /// Get input image for filter
  const irtkImage *GetInput() const;

  /** Initialize the filter. This function must be called by any derived
   *  filter class to perform some initialize tasks. */
  virtual void Initialize();

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Evaluate the filter at an arbitrary image location (in pixels)
  virtual double Evaluate(double, double, double, double = 0) /*const*/ = 0;

  // ---------------------------------------------------------------------------
  // Debugging

protected:

  /// Print debugging messages if debugging is enabled
  virtual void Debug(const char *);

  // ---------------------------------------------------------------------------
  // Deprecated

public:

  /// Set default value
  /// \deprecated Use DefaultValue(double) instead.
  irtkSetMacro(DefaultValue, double);

  /// Get default value
  /// \deprecated Use DefaultValue() instead.
  irtkGetMacro(DefaultValue, double);
};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline irtkImage *irtkImageFunction::GetInput()
{
  return _input;
}

// -----------------------------------------------------------------------------
inline const irtkImage *irtkImageFunction::GetInput() const
{
  return _input;
}

////////////////////////////////////////////////////////////////////////////////
// Image functions
////////////////////////////////////////////////////////////////////////////////

#include <irtkInterpolateImageFunction.h>
#include <irtkExtrapolateImageFunction.h>


#endif
