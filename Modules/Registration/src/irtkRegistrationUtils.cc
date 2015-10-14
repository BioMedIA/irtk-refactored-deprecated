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

#include <irtkRegistrationUtils.h>

#include <irtkPoint.h>
#include <irtkImageAttributes.h>
#include <irtkBaseImage.h>
#include <irtkGaussianBlurring.h>
#include <irtkVector3D.h>

#ifdef HAS_VTK
#include <vtkPolyData.h>
#include <irtkVector3.h>
#include <irtkMatrix3x3.h>
#endif


// -----------------------------------------------------------------------------
double AverageInterval(const set<double> &values)
{
  double avg = .0;
  if (values.size() > 1) {
    set<double>::const_iterator j = values.begin();
    set<double>::const_iterator i = j++;
    for (; j != values.end(); ++i, ++j) avg += (*j) - (*i);
    avg /= (values.size() - 1);
  }
  return avg;
}

// -----------------------------------------------------------------------------
int CenterOfForeground(const irtkBaseImage *image, double padding, irtkPoint &center)
{
  int n = 0;
  center._x = .0, center._y = .0, center._z = .0;
  for (int k = 0; k < image->GetZ(); ++k) {
    for (int j = 0; j < image->GetY(); ++j) {
      for (int i = 0; i < image->GetX(); ++i) {
        if (image->GetAsDouble(i, j, k) > padding) {
          center._x += i;
          center._y += j;
          center._z += k;
          ++n;
        }
      }
    }
  }
  if (n > 0) center /= n;
  image->ImageToWorld(center);
  return n;
}

// -----------------------------------------------------------------------------
irtkImageAttributes ImageDomain(const irtkBaseImage *image, double bg, double sigma, bool orthogonal)
{
  irtkImageAttributes attr = image->Attributes();
  if (attr._z == 1) attr._dz = .0;
  if (attr._t == 1) attr._dt = .0;
  // Determine lower bound along x axis: i1
  int i1 = attr._x;
  for (int k = 0; k < attr._z; ++k) {
  for (int j = 0; j < attr._y; ++j) {
  for (int i = 0; i < attr._x; ++i) {
    if (image->GetAsDouble(i, j, k) > bg) {
      if (i < i1) i1 = i;
      break;
    }
  }
  }
  }
  // Determine upper bound along x axis: i2
  int i2 = -1;
  for (int k = 0; k < attr._z; ++k) {
  for (int j = 0; j < attr._y; ++j) {
  for (int i = attr._x - 1; i >= i1; --i) {
    if (image->GetAsDouble(i, j, k) > bg) {
      if (i > i2) i2 = i;
      break;
    }
  }
  }
  }
  // Determine lower bound along y axis: j1
  int j1 = attr._y;
  for (int k = 0; k < attr._z; ++k) {
  for (int i = i1; i <= i2; ++i) {
  for (int j = 0; j < attr._y; ++j) {
    if (image->GetAsDouble(i, j, k) > bg) {
      if (j < j1) j1 = j;
      break;
    }
  }
  }
  }
  // Determine upper bound along y axis: j2
  int j2 = -1;
  for (int k = 0; k < attr._z; ++k) {
  for (int i = i1; i <= i2; ++i) {
  for (int j = attr._y - 1; j >= j1; --j) {
    if (image->GetAsDouble(i, j, k) > bg) {
      if (j > j2) j2 = j;
      break;
    }
  }
  }
  }
  // Determine lower bound along z axis: k1
  int k1 = attr._z;
  for (int j = j1; j <= j2; ++j) {
  for (int i = i1; i <= i2; ++i) {
  for (int k = 0; k < attr._z; ++k) {
    if (image->GetAsDouble(i, j, k) > bg) {
      if (k < k1) k1 = k;
      break;
    }
  }
  }
  }
  // Determine upper bound along z axis: k2
  int k2 = -1;
  for (int j = j1; j <= j2; ++j) {
  for (int i = i1; i <= i2; ++i) {
  for (int k = attr._z - 1; k >= k1; --k) {
    if (image->GetAsDouble(i, j, k) > bg) {
      if (k > k2) k2 = k;
      break;
    }
  }
  }
  }
  // "Crop" image region
  if (i1 <= i2 && j1 <= j2 && k1 <= k2) {
    // The following adjustment accounts for the blurring of the input image
    // before the downsampling which smears the foreground into background
    if (sigma > .0) {
      const int di = (attr._dx ? irtkGaussianBlurring<double>::KernelSize(sigma / attr._dx) : 0);
      const int dj = (attr._dy ? irtkGaussianBlurring<double>::KernelSize(sigma / attr._dy) : 0);
      const int dk = (attr._dz ? irtkGaussianBlurring<double>::KernelSize(sigma / attr._dz) : 0);
      i1 -= di, i2 += di;
      j1 -= dj, j2 += dj;
      k1 -= dk, k2 += dk;
      if (i1 <        0) i1 = 0;
      if (j1 <        0) j1 = 0;
      if (k1 <        0) k1 = 0;
      if (i2 >= attr._x) i2 = attr._x - 1;
      if (j2 >= attr._y) j2 = attr._y - 1;
      if (k2 >= attr._z) k2 = attr._z - 1;
    }
    // Convert upper index bounds to margin widths
    i2 = (attr._x - 1) - i2;
    j2 = (attr._y - 1) - j2;
    k2 = (attr._z - 1) - k2;
    // Adjust image lattice -- Attention: Order matters!
    attr._xorigin = 0.5 * ((attr._x - 1) + (i1 - i2));
    attr._yorigin = 0.5 * ((attr._y - 1) + (j1 - j2));
    attr._zorigin = 0.5 * ((attr._z - 1) + (k1 - k2));
    attr._x -= i1 + i2;
    attr._y -= j1 + j2;
    attr._z -= k1 + k2;
    image->ImageToWorld(attr._xorigin, attr._yorigin, attr._zorigin);
  }
  // Orthogonalize coordinate system; required in case of input images where
  // a previous affine (12 DoFs) alignment has been applied to the attributes
  return (orthogonal ? OrthogonalFieldOfView(attr) : attr);
}

#ifdef HAS_VTK


// -----------------------------------------------------------------------------
irtkImageAttributes PointSetDomain(vtkPointSet *data, double dx, double dy, double dz)
{
  irtkImageAttributes attr;
  attr._dx = dx;
  attr._dy = dy >= .0 ? dy : dx;
  attr._dz = dz >= .0 ? dz : dx;
  attr._dt = .0;
  // Compute eigenvectors of covariance matrix
  double c[3], p[3];
  data->GetCenter(c);
  irtkMatrix3x3 covar(.0);
  for (vtkIdType i = 0; i < data->GetNumberOfPoints(); ++i) {
    data->GetPoint(i, p);
    for (int d = 0; d < 3; ++d) p[d] -= c[d];
    for (int r = 0; r < 3; ++r) {
      for (int c = 0; c < 3; ++c) {
        covar[r][c] += p[r] * p[c];
      }
    }
  }
  double      eigen[3];
  irtkVector3 axis [3];
  covar.EigenSolveSymmetric(eigen, axis);
  irtkVector3::GenerateOrthonormalBasis(axis[0], axis[1], axis[2]);
  // Set output origin and orientation
  attr._xorigin = c[0];
  attr._yorigin = c[1];
  attr._zorigin = c[2];
  for (int d = 0; d < 3; ++d) {
    attr._xaxis[d] = axis[0][d];
    attr._yaxis[d] = axis[1][d];
    attr._zaxis[d] = axis[2][d];
  }
  // Determine bounds of reoriented data set
  double x, y, z;
  double xmin = numeric_limits<double>::max(), ymin = xmin, zmin = xmin;
  double xmax = -xmin, ymax = -ymin, zmax = -zmin;
  set<double> xs, ys, zs;
  for (vtkIdType i = 0; i < data->GetNumberOfPoints(); ++i) {
    data->GetPoint(i, p);
    x = attr._xaxis[0] * p[0] + attr._xaxis[1] * p[1] + attr._xaxis[2] * p[2];
    y = attr._yaxis[0] * p[0] + attr._yaxis[1] * p[1] + attr._yaxis[2] * p[2];
    z = attr._zaxis[0] * p[0] + attr._zaxis[1] * p[1] + attr._zaxis[2] * p[2];
    if (x < xmin) xmin = x;
    if (x > xmax) xmax = x;
    if (y < ymin) ymin = y;
    if (y > ymax) ymax = y;
    if (z < zmin) zmin = z;
    if (z > zmax) zmax = z;
    if (attr._dx <= .0) xs.insert(x);
    if (attr._dy <= .0) ys.insert(y);
    if (attr._dz <= .0) zs.insert(z);
  }
  // Set output resolution and size
  const double extent[3]  = {xmax - xmin, ymax - ymin, zmax - zmin};
  const double avg_extent = (extent[0] + extent[1] + extent[2]) / 3.0;
  if (attr._dx <= .0) attr._dx = ((extent[0] / avg_extent > 1e-3) ? AverageInterval(xs) : .0);
  if (attr._dy <= .0) attr._dy = ((extent[1] / avg_extent > 1e-3) ? AverageInterval(ys) : .0);
  if (attr._dz <= .0) attr._dz = ((extent[2] / avg_extent > 1e-3) ? AverageInterval(zs) : .0);
  attr._x  = (attr._dx > .0 ? round(extent[0] / attr._dx) : 0) + 1;
  attr._y  = (attr._dy > .0 ? round(extent[1] / attr._dy) : 0) + 1;
  attr._z  = (attr._dz > .0 ? round(extent[2] / attr._dz) : 0) + 1;
  attr._dx = (attr._x  >  1 ? extent[0] / (attr._x - 1) : .0);
  attr._dy = (attr._y  >  1 ? extent[1] / (attr._y - 1) : .0);
  attr._dz = (attr._z  >  1 ? extent[2] / (attr._z - 1) : .0);
  return attr;
}

// -----------------------------------------------------------------------------
irtkImageAttributes PointSetDomain(vtkPointSet *data, const irtkVector3D<double> &ds)
{
  return PointSetDomain(data, ds._x, ds._y, ds._z);
}

// -----------------------------------------------------------------------------
irtkImageAttributes PolyDataDomain(vtkPolyData *data, double dx, double dy, double dz)
{
  return PointSetDomain(data, dx, dy, dz);
}

// -----------------------------------------------------------------------------
irtkImageAttributes PolyDataDomain(vtkPolyData *data, const irtkVector3D<double> &ds)
{
  return PointSetDomain(data, ds);
}


#endif // defined(HAS_VTK)
