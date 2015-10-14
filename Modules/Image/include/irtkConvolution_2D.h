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

#ifndef _IRTKCONVOLUTION_2D_H

#define _IRTKCONVOLUTION_2D_H

/**
 * Class for two-dimensional convolution.
 *
 * This class defines and implements two-dimensional convolutions of an image
 * with a filter kernel. The convolution is computed along the x- and y-axis.
 * This class assumes that the filter kernel is two-dimensional and its size
 * along the z-axis must be 1.
 */

template <class VoxelType>
class irtkConvolution_2D : public irtkConvolution<VoxelType>
{
  irtkImageFilterMacro(irtkConvolution_2D);

protected:

  /// Second input, i.e. the filter kernel
  irtkGenericImage<irtkRealPixel> *_input2;

  /** Run the convolution filter. This method is protected and should only
   *  be called from within public member function Run().
   */
  virtual double Run(int, int, int, int);

public:

  /// Constructor
  irtkConvolution_2D(bool = false);

  /// Set second input, i.e. the filter kernel
  virtual void SetInput2(irtkGenericImage<irtkRealPixel> *);

  /// Initialize the convolution filter
  virtual void Initialize();
};

#endif
