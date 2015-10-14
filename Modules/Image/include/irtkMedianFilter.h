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

#ifndef _IRTKMEDIAN_H

#define _IRTKMEDIAN_H

#include <irtkImageToImage.h>

/**
 * Class for median filtering an image
 *
 */

template <class VoxelType>
class irtkMedianFilter : public irtkImageToImage<VoxelType>
{
  irtkImageFilterMacro(irtkMedianFilter);

protected:

  /// Initialize the filter
  virtual void Initialize();

  /// What connectivity to assume when running the filter.
  int _kernelRadius;

  irtkRealImage* _mask;

public:

  /// Constructor
  irtkMedianFilter();

  /// Destructor
  ~irtkMedianFilter();

  /// Run dilation
  virtual void Run();

  /// Set input image for filter
  void SetMask (irtkRealImage*);

  SetMacro(kernelRadius, int);

  GetMacro(kernelRadius, int);
};

#endif
