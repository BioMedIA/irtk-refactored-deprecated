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

#ifndef _IRTKCONVOLUTION_H

#define _IRTKCONVOLUTION_H

#include <irtkImageToImage.h>

template <class VoxelType>
class irtkConvolution : public irtkImageToImage<VoxelType>
{
  irtkImageFilterMacro(irtkConvolution);

protected:

  /// Flag whether to normalize convolution
  bool _Normalization;

public:

  /// Constructor
  irtkConvolution(bool = false);

  /// Set normalization on/off
  SetMacro(Normalization, bool);

  /// Set normalization on/off
  GetMacro(Normalization, bool);

};

// Convolution filters without padding
#include <irtkConvolution_1D.h>
#include <irtkConvolution_2D.h>
#include <irtkConvolution_3D.h>

// Convolution filters with padding
#include <irtkConvolutionWithPadding_1D.h>
#include <irtkConvolutionWithPadding_2D.h>
#include <irtkConvolutionWithPadding_3D.h>

#endif
