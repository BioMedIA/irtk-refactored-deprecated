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

// MUST be included first due to nasty IRTK "round" macro
#include <boost/random.hpp>

#include <irtkPointSamples.h>

typedef boost::mt19937 RandomNumberGeneratorType;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkPointSamples::irtkPointSamples(int n, int seed)
:
  irtkPointSet(n), _RandomNumberGenerator(new RandomNumberGeneratorType())
{
  RandomNumberGeneratorType &rng = *(RandomNumberGeneratorType*)_RandomNumberGenerator;
  if (seed < 0) rng.seed(static_cast<unsigned int>(std::time(0)));
  else          rng.seed(seed);
}

// -----------------------------------------------------------------------------
irtkPointSamples::~irtkPointSamples()
{
  delete (RandomNumberGeneratorType*)_RandomNumberGenerator;
}

// =============================================================================
// Uniform grid sampling
// =============================================================================

// -----------------------------------------------------------------------------
void irtkPointSamples::SampleGrid(const irtkPoint &p1, const irtkPoint &p2,
                                  int nx, int ny, int nz)
{
  SampleGrid(p1._x, p1._y, p1._z, p2._x, p2._y, p2._z, nx, ny, nz);
}

// -----------------------------------------------------------------------------
void irtkPointSamples::SampleGrid(const irtkPoint &p1, const irtkPoint &p2,
                                  double dx, double dy, double dz)
{
  SampleGrid(p1._x, p1._y, p1._z, p2._x, p2._y, p2._z, dx, dy, dz);
}

// -----------------------------------------------------------------------------
void irtkPointSamples::SampleGrid(double x1, double y1, double z1,
                                  double x2, double y2, double z2,
                                  int    nx, int    ny, int    nz)
{
  if (x2 < x1) swap(x1, x2);
  if (y2 < y1) swap(y1, y2);
  if (z2 < z1) swap(z1, z2);
  if (nx <= 0) nx = 1;
  if (ny <= 0) ny = 1;
  if (nz <= 0) nz = 1;
  const double sx = x2 - x1;
  const double sy = x2 - x1;
  const double sz = x2 - x1;
  const double dx = sx / nx;
  const double dy = sy / ny;
  const double dz = sz / nz;
  Size(nx * ny * nz);
  int n = 0;
  for (int k = 0; k < nz; ++k)
  for (int j = 0; j < ny; ++j)
  for (int i = 0; i < nx; ++i, ++n) {
    _data[n] = irtkPoint(x1 + i * dx, y1 + j * dy, z1 + k * dz);
  }
}

// -----------------------------------------------------------------------------
void irtkPointSamples::SampleGrid(double x1, double y1, double z1,
                                  double x2, double y2, double z2,
                                  double dx, double dy, double dz)
{
  if (x2 < x1) swap(x1, x2);
  if (y2 < y1) swap(y1, y2);
  if (z2 < z1) swap(z1, z2);
  if (dx < 0) dx = .0;
  if (dy < 0) dy = .0;
  if (dz < 0) dz = .0;
  const double sx = x2 - x1;
  const double sy = x2 - x1;
  const double sz = x2 - x1;
  const int    nx = (dx > .0 ? round(sx / dx) : 1);
  const int    ny = (dy > .0 ? round(sy / dy) : 1);
  const int    nz = (dz > .0 ? round(sz / dz) : 1);
  Size(nx * ny * nz);
  int n = 0;
  for (int k = 0; k < nz; ++k)
  for (int j = 0; j < ny; ++j)
  for (int i = 0; i < nx; ++i, ++n) {
    _data[n] = irtkPoint(x1 + i * dx, y1 + j * dy, z1 + k * dz);
  }
}

// =============================================================================
// Uniform spherical distribution
// =============================================================================

// -----------------------------------------------------------------------------
void irtkPointSamples::SampleSphere(double r)
{
  SampleSphere(.0, r);
}

// -----------------------------------------------------------------------------
void irtkPointSamples::SampleSphere(double c, double r)
{
  SampleSphere(c, c, c, r, r, r);
}

// -----------------------------------------------------------------------------
void irtkPointSamples::SampleSphere(const irtkPoint & c, double r)
{
  SampleSphere(c._x, c._y, c._z, r, r, r);
}

// -----------------------------------------------------------------------------
void irtkPointSamples::SampleSphere(const irtkPoint & c, double rx, double ry, double rz)
{
  SampleSphere(c._x, c._y, c._z, rx, ry, rz);
}

// -----------------------------------------------------------------------------
void irtkPointSamples::SampleSphere(double cx, double cy, double cz, double r)
{
  SampleSphere(cx, cy, cz, r, r, r);
}

// -----------------------------------------------------------------------------
void irtkPointSamples::SampleSphere(double cx, double cy, double cz,
                                    double rx, double ry, double rz)
{
  typedef boost::uniform_on_sphere<double>                                   Distribution;
  typedef boost::variate_generator<RandomNumberGeneratorType&, Distribution> Generator;

  RandomNumberGeneratorType &rng = *(RandomNumberGeneratorType*)_RandomNumberGenerator;

  Distribution dist(3);
  Generator next(rng, dist);

  vector<double> p(3);
  for (int i = 0; i < _n; ++i) {
    p = next();
    _data[i]._x = cx + rx * p[0];
    _data[i]._y = cy + ry * p[1];
    _data[i]._z = cz + rz * p[2];
  }
}

// =============================================================================
// Normal distribution
// =============================================================================

// -----------------------------------------------------------------------------
void irtkPointSamples::SampleGaussian(double s)
{
  SampleGaussian(.0, s);
}

// -----------------------------------------------------------------------------
void irtkPointSamples::SampleGaussian(double m, double s)
{
  SampleGaussian(m, m, m, s, s, s);
}

// -----------------------------------------------------------------------------
void irtkPointSamples::SampleGaussian(const irtkPoint & m, double s)
{
  SampleGaussian(m._x, m._y, m._z, s, s, s);
}

// -----------------------------------------------------------------------------
void irtkPointSamples::SampleGaussian(const irtkPoint & m, double sx, double sy, double sz)
{
  SampleGaussian(m._x, m._y, m._z, sx, sy, sz);
}

// -----------------------------------------------------------------------------
void irtkPointSamples::SampleGaussian(double mx, double my, double mz, double s)
{
  SampleGaussian(mx, my, mz, s, s, s);
}

// -----------------------------------------------------------------------------
void irtkPointSamples::SampleGaussian(double mx, double my, double mz,
                                      double sx, double sy, double sz)
{
  typedef boost::normal_distribution<double>                                 Distribution;
  typedef boost::variate_generator<RandomNumberGeneratorType&, Distribution> Generator;

  RandomNumberGeneratorType &rng = *(RandomNumberGeneratorType*)_RandomNumberGenerator;

  Distribution distx(mx, sx);
  Distribution disty(my, sy);
  Distribution distz(mz, sz);

  Generator x(rng, distx);
  Generator y(rng, disty);
  Generator z(rng, distz);

  for (int i = 0; i < _n; ++i) {
    _data[i]._x = x();
    _data[i]._y = y();
    _data[i]._z = z();
  }
}
