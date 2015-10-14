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

#include <irtkGeometry.h>
#include <irtkScalarGaussianDxDz.h>


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkScalarGaussianDxDz
::irtkScalarGaussianDxDz()
:
  irtkScalarGaussian()
{
}

// -----------------------------------------------------------------------------
irtkScalarGaussianDxDz
::irtkScalarGaussianDxDz(double sigma)
:
  irtkScalarGaussian(sigma)
{
}

// -----------------------------------------------------------------------------
irtkScalarGaussianDxDz
::irtkScalarGaussianDxDz(double sigma, double x_0, double y_0, double z_0, double t_0)
:
  irtkScalarGaussian(sigma, x_0, y_0, z_0, t_0)
{
}

// -----------------------------------------------------------------------------
irtkScalarGaussianDxDz
::irtkScalarGaussianDxDz(double sigma_x, double sigma_y, double sigma_z,
                         double x_0,     double y_0,     double z_0)
:
  irtkScalarGaussian(sigma_x, sigma_y, sigma_z, x_0, y_0, z_0)
{
}

// -----------------------------------------------------------------------------
irtkScalarGaussianDxDz
::irtkScalarGaussianDxDz(double sigma_x, double sigma_y, double sigma_z, double sigma_t,
                         double x_0,     double y_0,     double z_0,     double t_0)
:
  irtkScalarGaussian(sigma_x, sigma_y, sigma_z, sigma_t, x_0, y_0, z_0, t_0)
{
}

// -----------------------------------------------------------------------------
irtkScalarGaussianDxDz::~irtkScalarGaussianDxDz()
{
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
double irtkScalarGaussianDxDz::Evaluate(double)
{
  return .0;
}

// -----------------------------------------------------------------------------
double irtkScalarGaussianDxDz::Evaluate(double, double)
{
  return .0;
}

// -----------------------------------------------------------------------------
double irtkScalarGaussianDxDz::Evaluate(double x, double y, double z)
{
  x -= _Center._x, y -= _Center._y, z -= _Center._z;
  return _Norm[2] * exp(- 0.5 * (  (x * x) / _Variance._x
                                 + (y * y) / _Variance._y
                                 + (z * z) / _Variance._z))
         * (-x / _Variance._x) * (-z / _Variance._z);
}

// -----------------------------------------------------------------------------
double irtkScalarGaussianDxDz::Evaluate(double x, double y, double z, double t)
{
  x -= _Center._x, y -= _Center._y, z -= _Center._z, t -= _Center._t;
  return _Norm[3] * exp(- 0.5 * (  (x * x) / _Variance._x
                                 + (y * y) / _Variance._y
                                 + (z * z) / _Variance._z
                                 + (t * t) / _Variance._t))
         * (-x / _Variance._x) * (-z / _Variance._z);
}
