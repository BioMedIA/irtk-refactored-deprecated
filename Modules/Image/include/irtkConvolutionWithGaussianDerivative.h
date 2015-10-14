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

#ifndef _IRTKCONVOLUTIONWITHGAUSSIANDERIVATIVE_H

#define _IRTKCONVOLUTIONWITHGAUSSIANDERIVATIVE_H

#include <irtkImageToImage.h>

/** 
 * Class for convolution with 1st order Gaussian derivative 
 * 
 * This class defines and implements the 1st order gaussian derivative filtering of images. 
 */

template <class VoxelType>
class irtkConvolutionWithGaussianDerivative : public irtkImageToImage<VoxelType>
{
  irtkImageFilterMacro(irtkConvolutionWithGaussianDerivative);

protected:

  /// Sigma (standard deviation of Gaussian kernel)
  double _Sigma;

public:

  /// Constructor
  irtkConvolutionWithGaussianDerivative(double);

  /// Destructor
  ~irtkConvolutionWithGaussianDerivative();

  /// Compute derivatives
  void Ix();

  /// Compute derivatives
  void Iy();

  /// Compute derivatives
  void Iz();
  
  /// Set sigma
  SetMacro(Sigma, double);
 
  /// Get sigma
  GetMacro(Sigma, double);

};


#endif
