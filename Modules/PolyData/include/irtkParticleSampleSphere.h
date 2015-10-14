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

//References
// [1] Thomson problem: http://en.wikipedia.org/wiki/Thomson_problem
// [2] Saff, E.B., Kuijlaars, A.B.J. (1997) ‘Distributing many points on a sphere’
// [3] Shao M.-Z., Badler, N. (1996),'Spherical sampling by Archimedes' theorem'

// Create an approximately uniform triangular tessellation of the unit
// sphere by minimizing generalized electrostatic potential energy
// (aka Reisz s-energy) of the system of charged particles. Effectively,
// this function produces a locally optimal solution to the problem that
// involves finding a minimum Reisz s-energy configuration of N equal
// charges confined to the surface of the unit sphere (s=1 corresponds to
// the problem originally posed by J. J. Thomson).

#ifndef IRTK_PARTICLESAMPLESPHERE_H_
#define IRTK_PARTICLESAMPLESPHERE_H_

#include <vtkMath.h>
#include <vtkTriangle.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkPlatonicSolidSource.h>
#include <vtkPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>

class irtkParticleSampleSphere
{
private:

  double _stepSize,_minStepSize,_maxStepSize,_energyTolerance;
  int _maxIterations,_numPoints,_subdivisions;
  bool _sphericalCoordsOn;
  vtkSmartPointer<vtkPolyData> _output;
  void Initialize();
  vtkSmartPointer<vtkPoints> StratifiedSample();
  void OptimizeParticlePositions(vtkSmartPointer<vtkPoints>);
  void AddSphericalCoordScalars(vtkSmartPointer<vtkPolyData>);

public:

  irtkParticleSampleSphere();
  ~irtkParticleSampleSphere();

  void SetNumberOfPoints(int);
  void SetNumberOfSubdivisions(int);

  void SphericalCoordsOn();
  void SphericalCoordsOff();

  vtkSmartPointer<vtkPolyData> GetSphere();
};


#endif
