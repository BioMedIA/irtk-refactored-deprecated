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

#include <irtkVelocityToDisplacementField.h>


// ===========================================================================
// Construction/Destruction
// ===========================================================================

// ---------------------------------------------------------------------------
template <class VoxelType>
irtkVelocityToDisplacementField<VoxelType>::irtkVelocityToDisplacementField()
:
  _Interpolation                   (Interpolation_Linear),
  _Extrapolation                   (Extrapolation_NN),
  _ComputeInterpolationCoefficients(true),
  _ComputeInverse                  (false),
  _NumberOfSteps                   (64), // i.e., 6 squaring steps in case of SS method
  _T                               (1.0),
  _InputDisplacementField          (NULL)
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
irtkVelocityToDisplacementField<VoxelType>::~irtkVelocityToDisplacementField()
{
}

// ===========================================================================
// Filter implementation
// ===========================================================================

// --------------------------------------------------------------------------
template <class VoxelType>
void irtkVelocityToDisplacementField<VoxelType>::SetInput(int i, irtkGenericImage<VoxelType> *image)
{
  if      (i == 0) irtkImageToImage<VoxelType>::SetInput(image);
  else if (i == 1) _InputDisplacementField = image;
  else {
    cerr << this->NameOfClass() << "::SetInput: Input index out of range: " << i << endl;
    exit(1);
  }
}

// --------------------------------------------------------------------------
template <class VoxelType>
irtkGenericImage<VoxelType> *irtkVelocityToDisplacementField<VoxelType>::GetInput(int i)
{
  if      (i == 0) return this->_input;
  else if (i == 1) return _InputDisplacementField;
  else {
    cerr << this->NameOfClass() << "::GetInput: Input index out of range: " << i << endl;
    exit(1);
  }
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void irtkVelocityToDisplacementField<VoxelType>::Initialize()
{
  // Initialize base class
  irtkImageToImage<VoxelType>::Initialize();

  // Check input
  if (this->GetInput()->GetZ() < 1) {
    // not allowed as otherwise a for loop over z from 0 to _z-1 would never
    // be executed, not even once as it should be in case of 2 dimensions
    cerr << this->NameOfClass() << "::Initialize: Size of z dimension must be 1 in case of 2D vector field" << endl;
    exit(1);
  }

  if ((this->GetInput()->GetZ() <= 1 && this->GetInput()->GetT() != 2 && this->GetInput()->GetT() != 3) ||
      (this->GetInput()->GetZ() >  1 && this->GetInput()->GetT() != 3)) {
    cerr << this->NameOfClass() << "::Initialize: Input must be a 2D or 3D vector field" << endl;
    exit(1);
  }

  // Check parameters
  if (_T == 0) {
    cerr << this->NameOfClass() << "::Initialize: Upper integration limit (T) must not be zero" << endl;
    exit(1);
  }
  if (_NumberOfSteps < 1) {
    cerr << this->NameOfClass() << "::Initialize: Number of integration steps must be positive" << endl;
    exit(1);
  }
  if (_ComputeInverse) _T = -_T; // Inverse <=> backward integration

  // Output is vector field
  this->GetOutput()->PutTSize(.0);
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void irtkVelocityToDisplacementField<VoxelType>::Finalize()
{
  // Reset T attribute
  if (_ComputeInverse) _T = -_T;
  // Finalize base class
  irtkImageToImage<VoxelType>::Finalize();
}

// ===========================================================================
// Explicit template instantiations
// ===========================================================================

template class irtkVelocityToDisplacementField<float>;
template class irtkVelocityToDisplacementField<double>;
