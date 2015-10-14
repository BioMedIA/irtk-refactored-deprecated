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

#ifndef _IRTKDIFFERENCEOFCOMPOSITIONLIEBRACKETIMAGEFILTER3D_H
#define _IRTKDIFFERENCEOFCOMPOSITIONLIEBRACKETIMAGEFILTER3D_H


/**
 * Image filter for computation of Lie bracket of two 3D vector fields.
 *
 * This filter implements the definition of the Lie bracket as the difference
 * of the composition of the first vector field with the second and vice versa,
 * i.e., [X,Y] = X(Y) - Y(X).
 */

template <class VoxelType>
class irtkDifferenceOfCompositionLieBracketImageFilter3D
: public irtkLieBracketImageFilter<VoxelType>
{
  irtkImageFilterMacro(irtkDifferenceOfCompositionLieBracketImageFilter3D);

  /// Vector field interpolation mode
  irtkPublicAttributeMacro(irtkInterpolationMode, Interpolation);

  /// Vector field extrapolation mode
  irtkPublicAttributeMacro(irtkExtrapolationMode, Extrapolation);

  /// Whether to compute interpolation coefficients from the given input
  /// or if the input images contain the coefficients already
  irtkPublicAttributeMacro(bool, ComputeInterpolationCoefficients);

protected:

  typedef irtkGenericImage<VoxelType>                      InputType;
  typedef irtkGenericInterpolateImageFunction<InputType>   InterpolatorType;

  InterpolatorType *_Interpolator[2]; /// Input vector field interpolators
  double            _Scaling     [2]; /// Scaling of input vector fields

  /// Initialize filter
  virtual void Initialize();

public:

  /// Constructor
  irtkDifferenceOfCompositionLieBracketImageFilter3D();

  /// Destructor
  virtual ~irtkDifferenceOfCompositionLieBracketImageFilter3D();

  /// Set scaling of n-th input vector field
  void SetScaling(int, double);

  /// Run filter on every voxel
  virtual void Run();

  /// Run filter on single voxel
  virtual void Run(double [3], int, int, int);

  /// Run filter on single voxel and component
  virtual double Run(int, int, int, int);

};


#endif
