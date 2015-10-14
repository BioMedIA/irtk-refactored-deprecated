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

#ifndef _IRTKCONVOLUTIONWITHGAUSSIANDERIVATIVE2_H

#define _IRTKCONVOLUTIONWITHGAUSSIANDERIVATIVE2_H

#include <irtkImageToImage.h>

/** 
 * Class for convolution with 2nd order Gaussian derivative 
 * 
 * This class defines and implements the 2nd order gaussian derivative filtering of images. 
 */

template <class VoxelType>
class irtkConvolutionWithGaussianDerivative2 : public irtkImageToImage<VoxelType>
{
  irtkImageFilterMacro(irtkConvolutionWithGaussianDerivative2);

protected:

  /// Sigma (standard deviation of Gaussian kernel)
  double _Sigma;

public:

  /// Constructor
  irtkConvolutionWithGaussianDerivative2(double);

  /// Destructor
  ~irtkConvolutionWithGaussianDerivative2();

  /// Compute derivatives
  void Ixx();

  /// Compute derivatives
  void Iyy();

  /// Compute derivatives
  void Izz();

  /// Compute derivatives
  void Ixy();

  /// Compute derivatives
  void Ixz();

  /// Compute derivatives
  void Iyz();
  
  /// Set sigma
  SetMacro(Sigma, double);
 
  /// Get sigma
  GetMacro(Sigma, double);

};


#endif
