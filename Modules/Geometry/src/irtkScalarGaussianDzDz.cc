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
#include <irtkScalarGaussianDzDz.h>


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkScalarGaussianDzDz
::irtkScalarGaussianDzDz()
:
  irtkScalarGaussian()
{
}

// -----------------------------------------------------------------------------
irtkScalarGaussianDzDz
::irtkScalarGaussianDzDz(double sigma)
:
  irtkScalarGaussian(sigma)
{
}

// -----------------------------------------------------------------------------
irtkScalarGaussianDzDz
::irtkScalarGaussianDzDz(double sigma, double x_0, double y_0, double z_0, double t_0)
:
  irtkScalarGaussian(sigma, x_0, y_0, z_0, t_0)
{
}

// -----------------------------------------------------------------------------
irtkScalarGaussianDzDz
::irtkScalarGaussianDzDz(double sigma_x, double sigma_y, double sigma_z,
                         double x_0,     double y_0,     double z_0)
:
  irtkScalarGaussian(sigma_x, sigma_y, sigma_z, x_0, y_0, z_0)
{
}

// -----------------------------------------------------------------------------
irtkScalarGaussianDzDz
::irtkScalarGaussianDzDz(double sigma_x, double sigma_y, double sigma_z, double sigma_t,
                         double x_0,     double y_0,     double z_0,     double t_0)
:
  irtkScalarGaussian(sigma_x, sigma_y, sigma_z, sigma_t, x_0, y_0, z_0, t_0)
{
}

// -----------------------------------------------------------------------------
irtkScalarGaussianDzDz::~irtkScalarGaussianDzDz()
{
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
double irtkScalarGaussianDzDz::Evaluate(double)
{
  return .0;
}

// -----------------------------------------------------------------------------
double irtkScalarGaussianDzDz::Evaluate(double, double)
{
  return .0;
}

// -----------------------------------------------------------------------------
double irtkScalarGaussianDzDz::Evaluate(double x, double y, double z)
{
  x -= _Center._x, y -= _Center._y, z -= _Center._z;
  return _Norm[2] * exp(- 0.5 * (  (x * x) / _Variance._x
                                 + (y * y) / _Variance._y
                                 + (z * z) / _Variance._z))
         * ((z * z) / (_Variance._z * _Variance._z) - 1.0 / _Variance._z);
}

// -----------------------------------------------------------------------------
double irtkScalarGaussianDzDz::Evaluate(double x, double y, double z, double t)
{
  x -= _Center._x, y -= _Center._y, z -= _Center._z, t -= _Center._t;
  return _Norm[3] * exp(- 0.5 * (  (x * x) / _Variance._x
                                 + (y * y) / _Variance._y
                                 + (z * z) / _Variance._z
                                 + (t * t) / _Variance._t))
         * ((z * z) / (_Variance._z * _Variance._z) - 1.0 / _Variance._z);
}
