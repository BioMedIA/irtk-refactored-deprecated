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

#ifndef _IRTKHOMOGENEOUSIMAGETRANSFORMATION_H

#define _IRTKHOMOGENEOUSIMAGETRANSFORMATION_H

#include <irtkTransformation.h>

class irtkImageHomogeneousTransformation  : public irtkImageTransformation
{
  irtkObjectMacro(irtkImageHomogeneousTransformation);

public:

  /** Constructor. This constructs an transformation filter with a given
   *  interpolation mode and padding value. By default the interpolation
   *  mode is set to trilinear.
   */
  irtkImageHomogeneousTransformation();

  /// Destructor
  virtual ~irtkImageHomogeneousTransformation();

  /// Sets transformation
  virtual void SetTransformation(irtkTransformation *);

  /// Runs the filter
  virtual void Run();

};

#endif
