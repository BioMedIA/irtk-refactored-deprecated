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

#ifndef _IRTKREGISTRATIONUTILS_H
#define _IRTKREGISTRATIONUTILS_H

#include <irtkCxxLib.h>

class irtkBaseImage;
class irtkImageAttributes;
class irtkPoint;


/// Get average interval length given order set of (abscissa) values
double AverageInterval(const set<double> &);

/// Determine center of mass of image foreground
///
/// \param[in]  image   Intensity image.
/// \param[in]  padding Background value.
/// \param[out] center  Centroid of image foreground.
///
/// \returns Number of foreground intensities encountered or zero if the image does
///          not contain any foreground, in which case the \p center is invalid.
int CenterOfForeground(const irtkBaseImage *image, double padding, irtkPoint &center);

/// Determine minimal axes-aligned foreground bounding region
///
/// \param[in] image      Intensity image.
/// \param[in] padding    Background value.
/// \param[in] sigma      Standard deviation of Gaussian blurring filter that will
///                       be applied to the intensity \p image. Used to increase
///                       the extent of the foreground region such that also background
///                       voxels which are averaged with foreground voxels are included.
/// \param[in] orthogonal Whether to orthogonalize image axes. Useful only in case
///                       of input images where a previous affine (12 DoFs) alignment
///                       has been applied to the attributes.
///
/// \returns Attributes of foreground region.
irtkImageAttributes ImageDomain(const irtkBaseImage *image, double padding, double sigma = .0, bool orthogonal = true);


#ifdef HAS_VTK
class vtkPointSet;
class vtkPolyData;
template <typename> struct irtkVector3D;


/// Determine bounding box of point set
///
/// \param[in] data Data set.
/// \param[in] dx   Desired lattice spacing along x axis.
///                 If non-positive, the average interval between the x
///                 coordinates of the input points is used.
/// \param[in] dy   Desired lattice spacing along y axis.
///                 If non-positive, the average interval between the y
///                 coordinates of the input points is used.
/// \param[in] dz   Desired lattice spacing along z axis.
///                 If non-positive, the average interval between the z
///                 coordinates of the input points is used.
///
/// \returns Attributes of oriented minimum-volume bounding box.
irtkImageAttributes PointSetDomain(vtkPointSet *data, double dx = -1, double dy = -1, double dz = -1);

/// Determine bounding box of point set
///
/// \param[in] data Data set.
/// \param[in] ds   Desired lattice spacing.
///                 If an entry is non-positive, the average interval between
///                 the coordinates of the input points along this axis is used.
///
/// \returns Attributes of oriented minimum-volume bounding box.
irtkImageAttributes PointSetDomain(vtkPointSet *data, const irtkVector3D<double> &ds);

/// Determine bounding box of polydata points
///
/// \param[in] data Data set.
/// \param[in] dx   Desired lattice spacing along x axis.
///                 If non-positive, the average interval between the x
///                 coordinates of the input points is used.
/// \param[in] dy   Desired lattice spacing along y axis.
///                 If non-positive, the average interval between the y
///                 coordinates of the input points is used.
/// \param[in] dz   Desired lattice spacing along z axis.
///                 If non-positive, the average interval between the z
///                 coordinates of the input points is used.
///
/// \returns Attributes of oriented minimum-volume bounding box.
irtkImageAttributes PolyDataDomain(vtkPolyData *data, double dx = -1, double dy = -1, double dz = -1);

/// Determine bounding box of polydata points
///
/// \param[in] data Data set.
/// \param[in] ds   Desired lattice spacing.
///                 If an entry is non-positive, the average interval between
///                 the coordinates of the input points along this axis is used.
///
/// \returns Attributes of oriented minimum-volume bounding box.
irtkImageAttributes PolyDataDomain(vtkPolyData *data, const irtkVector3D<double> &ds);


#endif
#endif
