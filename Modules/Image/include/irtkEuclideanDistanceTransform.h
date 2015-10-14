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

#ifndef _IRTKEUCLIDEANDISTANCETRANSFORM_H

#define _IRTKEUCLIDEANDISTANCETRANSFORM_H

#define EDT_MAX_IMAGE_DIMENSION 26754
#define EDT_MAX_DISTANCE_SQUARED 2147329548
#define EDT_MAX_DISTANCE_SQUARED_ANISOTROPIC 2147329548

#include <irtkImageToImage.h>

template <class VoxelType>
class irtkEuclideanDistanceTransform : public irtkImageToImage<VoxelType>
{
  irtkInPlaceImageFilterMacro(irtkEuclideanDistanceTransform);

public:

  /// 2D or 3D distance transform
  enum irtkDistanceTransformMode { irtkDistanceTransform2D, irtkDistanceTransform3D };

protected:

  /// 2D or 3D distance transform
  irtkDistanceTransformMode _distanceTransformMode;

  /// Calculate the Vornoi diagram
  int edtVornoiEDT(long *, long);

  /// Calculate 2D distance transform
  void edtComputeEDT_2D(char *, long *, long, long);

  /// Calculate 3D distance transform
  void edtComputeEDT_3D(char *, long *, long, long, long);

  /// Calculate the Vornoi diagram for anisotripic voxel sizes
  int edtVornoiEDT_anisotropic(VoxelType *, long, double);

  /// Calculate 2D distance transform for anisotripic voxel sizes
  void edtComputeEDT_2D_anisotropic(VoxelType *, VoxelType *, long, long, double, double);

  /// Calculate 3D distance transform for anisotripic voxel sizes
  void edtComputeEDT_3D_anisotropic(VoxelType *, VoxelType *, long, long, long, double, double, double);

public:

  /// Default constructor
  irtkEuclideanDistanceTransform(irtkDistanceTransformMode = irtkDistanceTransform3D);

  /// Destructor (empty).
  ~irtkEuclideanDistanceTransform() {};

  // Run distance transform
  virtual void Run();

  // Get Radial
  virtual void Radial();

  // Get Radial+Thickness
  virtual void TRadial();
};


#endif
