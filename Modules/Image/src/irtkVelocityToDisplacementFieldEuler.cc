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

#include <irtkVelocityToDisplacementFieldEuler.h>

#include <irtkVoxelFunction.h>

using irtkNaryVoxelFunction::ExpVelocityFieldEuler2D;
using irtkNaryVoxelFunction::ExpVelocityFieldEuler3D;


// ---------------------------------------------------------------------------
// Forward declaration of CUDA accelerated implementation
#ifdef USE_CUDA
#include <irtkCUImage.h>
namespace irtkCUVelocityToDisplacementFieldEuler {
  template <class VoxelType>
  void Run(irtkCUGenericImage<VoxelType>       *,
           const irtkCUGenericImage<VoxelType> *,
           const irtkCUGenericImage<VoxelType> *, int, double);
}
#endif


// ===========================================================================
// Construction/Destruction
// ===========================================================================

// ---------------------------------------------------------------------------
template <class VoxelType>
irtkVelocityToDisplacementFieldEuler<VoxelType>::irtkVelocityToDisplacementFieldEuler()
:
  irtkVelocityToDisplacementField<VoxelType>(),
  _VelocityInterpolator(NULL)
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
irtkVelocityToDisplacementFieldEuler<VoxelType>::~irtkVelocityToDisplacementFieldEuler()
{
  delete _VelocityInterpolator;
}

// ===========================================================================
// Filter implementation
// ===========================================================================

// ---------------------------------------------------------------------------
template <class VoxelType>
void irtkVelocityToDisplacementFieldEuler<VoxelType>::Initialize()
{
  // Initialize base class
  irtkVelocityToDisplacementField<VoxelType>::Initialize();

  // Initialize interpolator
  if (_VelocityInterpolator) delete _VelocityInterpolator;
  _VelocityInterpolator = irtkInterpolateImageFunction::New(this->Interpolation(),
                                                            this->Extrapolation(),
                                                            this->GetInput());
  _VelocityInterpolator->SetInput(this->GetInput());
  _VelocityInterpolator->Initialize();
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void irtkVelocityToDisplacementFieldEuler<VoxelType>::Run()
{
  // Do the initial set up
  this->Initialize();

  // Perform forward Euler integration
#ifdef USE_CUDA
  irtkCUGenericImage<VoxelType> *v    = dynamic_cast<irtkCUGenericImage<VoxelType> *>(this->GetInput());
  irtkCUGenericImage<VoxelType> *din  = dynamic_cast<irtkCUGenericImage<VoxelType> *>(this->GetInput(1));
  irtkCUGenericImage<VoxelType> *dout = dynamic_cast<irtkCUGenericImage<VoxelType> *>(this->GetOutput());

  if (use_gpu && v && dout && din == this->GetInput(1)) {
    const double dt = this->T() / static_cast<double>(this->NumberOfSteps());
    irtkCUVelocityToDisplacementFieldEuler::Run(dout, din, v, this->NumberOfSteps(), dt);
  } else
#endif
  if (this->GetOutput()->GetZ() > 1) {
    const irtkImageAttributes &grid = this->GetOutput()->GetImageAttributes();
    ExpVelocityFieldEuler3D<> exp(_VelocityInterpolator, this->NumberOfSteps(), this->T());
    if (this->GetInput(1)) ParallelForEachVoxel(grid, this->GetInput(1), this->GetOutput(), exp);
    else                   ParallelForEachVoxel(grid,                    this->GetOutput(), exp);
  } else {
    const irtkImageAttributes &grid = this->GetOutput()->GetImageAttributes();
    ExpVelocityFieldEuler2D<> exp(_VelocityInterpolator, this->NumberOfSteps(), this->T());
    if (this->GetInput(1)) ParallelForEachVoxel(grid, this->GetInput(1), this->GetOutput(), exp);
    else                   ParallelForEachVoxel(grid,                    this->GetOutput(), exp);
  }

  // Do the final cleaning up
  this->Finalize();
}

// ===========================================================================
// Explicit template instantiations
// ===========================================================================

template class irtkVelocityToDisplacementFieldEuler<float>;
template class irtkVelocityToDisplacementFieldEuler<double>;
