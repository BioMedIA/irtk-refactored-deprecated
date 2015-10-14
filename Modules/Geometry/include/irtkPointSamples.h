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

#ifndef IRTKPOINTSAMPLES_H_
#define IRTKPOINTSAMPLES_H_

#include <irtkPointSet.h>


/**
 * Auxiliary class for generation of point samples
 */
class irtkPointSamples : public irtkPointSet
{
  irtkObjectMacro(irtkPointSamples);

  /// Random number generator
  void *_RandomNumberGenerator;

public:

  // ---------------------------------------------------------------------------
  // Construction/destruction

  /// Constructor
  ///
  /// \param[in] n    Number of point samples.
  /// \param[in] seed Seed of random number generator.
  irtkPointSamples(int n = 0, int seed = 0);

  /// Destructor
  virtual ~irtkPointSamples();

  // ---------------------------------------------------------------------------
  // Uniform grid sampling

  /// Sample axes-aligned uniform grid
  void SampleGrid(const irtkPoint &p1, const irtkPoint &p2,
                  int nx, int ny, int nz);

  /// Sample axes-aligned uniform grid
  void SampleGrid(const irtkPoint &p1, const irtkPoint &p2,
                  double dx, double dy, double dz);

  /// Sample axes-aligned uniform grid
  void SampleGrid(double x1, double y1, double z1,
                  double x2, double y2, double z2,
                  int    nx, int    ny, int    nz);

  /// Sample axes-aligned uniform grid
  void SampleGrid(double x1, double y1, double z1,
                  double x2, double y2, double z2,
                  double dx, double dy, double dz);

  // ---------------------------------------------------------------------------
  // Uniform spherical distribution

  /// Add uniform spherical point samples
  void SampleSphere(double r = 1.0);

  /// Add uniform spherical point samples
  ///
  /// \param[in] c Center in each dimension.
  /// \param[in] r Radius in each dimension.
  void SampleSphere(double c, double r);

  /// Add uniform spherical point samples
  ///
  /// \param[in] c Center point.
  /// \param[in] r Radius in each dimension.
  void SampleSphere(const irtkPoint & c, double r = 1.0);

  /// Add uniform spherical point samples
  ///
  /// \param[in] c  Center point.
  /// \param[in] rx Radius in x direction.
  /// \param[in] ry Radius in x direction.
  /// \param[in] rz Radius in x direction.
  void SampleSphere(const irtkPoint & c, double rx, double ry, double rz);

  /// Add uniform spherical point samples
  ///
  /// \param[in] cx Center in x direction.
  /// \param[in] cy Center in y direction.
  /// \param[in] cz Center in z direction.
  /// \param[in] r  Radius in each dimension.
  void SampleSphere(double cx, double cy, double cz, double r);

  /// Add uniform spherical point samples
  ///
  /// \param[in] cx Center in x direction.
  /// \param[in] cy Center in y direction.
  /// \param[in] cz Center in z direction.
  /// \param[in] rx Radius in x direction.
  /// \param[in] ry Radius in x direction.
  /// \param[in] rz Radius in x direction.
  void SampleSphere(double cx, double cy, double cz,
                    double rx, double ry, double rz);

  // ---------------------------------------------------------------------------
  // Normal distribution

  /// Add normally distributed point samples
  ///
  /// \param[in] s Standard deviation in each dimension.
  void SampleGaussian(double s = 1.0);

  /// Add normally distributed point samples
  ///
  /// \param[in] m Mean in each dimension.
  /// \param[in] s Standard deviation in each dimension.
  void SampleGaussian(double m, double s);

  /// Add normally distributed point samples
  ///
  /// \param[in] mx Mean in x direction.
  /// \param[in] my Mean in y direction.
  /// \param[in] mz Mean in z direction.
  /// \param[in] s  Standard deviation in each dimension.
  void SampleGaussian(double mx, double my, double mz, double s);

  /// Add normally distributed point samples
  ///
  /// \param[in] m Mean of normal distribution.
  /// \param[in] s Standard deviation in each dimension.
  void SampleGaussian(const irtkPoint & m, double s);

  /// Add normally distributed point samples
  ///
  /// \param[in] m  Mean of normal distribution.
  /// \param[in] sx Standard deviation in x direction.
  /// \param[in] sy Standard deviation in y direction.
  /// \param[in] sz Standard deviation in z direction.
  void SampleGaussian(const irtkPoint & m, double sx, double sy, double sz);

  /// Add normally distributed point samples
  ///
  /// \param[in] mx Mean in x direction.
  /// \param[in] my Mean in y direction.
  /// \param[in] mz Mean in z direction.
  /// \param[in] sx Standard deviation in x direction.
  /// \param[in] sy Standard deviation in y direction.
  /// \param[in] sz Standard deviation in z direction.
  void SampleGaussian(double mx, double my, double mz,
                      double sx, double sy, double sz);

};


#endif
