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

#include <irtkVoxelFunction.h>
#include <irtkImageFunction.h>
#include <irtkGaussianBlurring2D.h>
#include <irtkLieBracketImageFilter.h>
#include <irtkVelocityToDisplacementFieldEuler.h>

#ifdef USE_CUDA
#  include <irtkCUImage.h>
namespace irtkCUDisplacementToVelocityFieldBCH {
  template <class VoxelType>
  void Update(irtkCUGenericImage<VoxelType> *,
              const irtkCUGenericImage<VoxelType> *,
              const irtkCUGenericImage<VoxelType> *,
              const irtkCUGenericImage<VoxelType> *,
              const irtkCUGenericImage<VoxelType> *,
              const irtkCUGenericImage<VoxelType> *, int);
}
#endif

#include <irtkDisplacementToVelocityField.h>


// ===========================================================================
// Construction/Destruction
// ===========================================================================

// ---------------------------------------------------------------------------
template <class VoxelType>
irtkDisplacementToVelocityFieldBCH<VoxelType>::irtkDisplacementToVelocityFieldBCH()
:
  _dv(NULL), _l1(NULL), _l2(NULL), _l3(NULL), _l4(NULL),
  _ExponentialFilter(new irtkVelocityToDisplacementFieldSS<VoxelType>),
  _CustomExponentialFilter(false),
  _NumberOfIterations(8),
  _NumberOfTerms(3),
  _UseJacobian(false),
  _SmoothVelocities(false)
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
irtkDisplacementToVelocityFieldBCH<VoxelType>::~irtkDisplacementToVelocityFieldBCH()
{
  if (!_CustomExponentialFilter) delete _ExponentialFilter;
}

// ===========================================================================
// Filter implementation
// ===========================================================================

// ----------------------------------------------------------------------------
template <class VoxelType>
void irtkDisplacementToVelocityFieldBCH<VoxelType>::SetExponentialFilter(irtkVelocityToDisplacementField<VoxelType> *filter)
{
  if (!_CustomExponentialFilter) delete _ExponentialFilter;
  _ExponentialFilter       = filter;
  _CustomExponentialFilter = (_ExponentialFilter != NULL);
}

// ----------------------------------------------------------------------------
template <class VoxelType>
irtkVelocityToDisplacementField<VoxelType> *irtkDisplacementToVelocityFieldBCH<VoxelType>::GetExponentialFilter()
{
  return _ExponentialFilter;
}

// ----------------------------------------------------------------------------
template <class VoxelType>
const irtkVelocityToDisplacementField<VoxelType> *irtkDisplacementToVelocityFieldBCH<VoxelType>::GetExponentialFilter() const
{
  return _ExponentialFilter;
}

// ----------------------------------------------------------------------------
template <class VoxelType>
void irtkDisplacementToVelocityFieldBCH<VoxelType>::SetT(double t)
{
  _ExponentialFilter->T(t);
}

// ----------------------------------------------------------------------------
template <class VoxelType>
double irtkDisplacementToVelocityFieldBCH<VoxelType>::GetT()
{
  return _ExponentialFilter->T();
}

// ----------------------------------------------------------------------------
template <class VoxelType>
void irtkDisplacementToVelocityFieldBCH<VoxelType>::SetNumberOfSteps(int n)
{
  _ExponentialFilter->NumberOfSteps(n);
}

// ----------------------------------------------------------------------------
template <class VoxelType>
int irtkDisplacementToVelocityFieldBCH<VoxelType>::GetNumberOfSteps()
{
  return _ExponentialFilter->NumberOfSteps();
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void irtkDisplacementToVelocityFieldBCH<VoxelType>::Initialize()
{
  // Initialize base class
  irtkDisplacementToVelocityField<VoxelType>::Initialize();

  // Get attributes of input vector field
  const irtkImageAttributes &grid = this->GetInput()->GetImageAttributes();

  // Check input
  if ((grid._z <= 1 && grid._t != 2 && grid._t != 3) || (grid._z > 1 && grid._t != 3)) {
    cerr << this->NameOfClass() << "::Initialize: Input must be a valid vector field" << endl;
    exit(1);
  }

  // Check parameters
  if (_NumberOfTerms < 2 || _NumberOfTerms > 6) {
    cerr << this->NameOfClass() << "::Initialize: Number of BCH terms must be 2-6" << endl;
    exit(1);
  }

  // Ensure that exponential filter is set
  if (_ExponentialFilter == NULL) {
    // Developer may by mistake passed a NULL  pointer to SetExponentialFilter,
    // so let them know about this mistake rather than instantiating a default
    // exponential filter here. The default filter is created in the constructor.
    cerr << this->NameOfClass() << "::Initialize: No filter for exponentiation of velocity field set" << endl;
    exit(1);
  }

  // Allocate intermediate images
  switch (_NumberOfTerms) {
    // Attention: Must be in descending order and without break statements!
    case 6: _l4 = new ImageType(grid);
    case 5: _l3 = new ImageType(grid);
    case 4: _l2 = new ImageType(grid);
    case 3: _l1 = new ImageType(grid);
    case 2: _dv = new ImageType(grid);
  };

  // Initialize output velocity field
  // v_0 = 0
  // v_1 = exp(-v_0) ° phi = phi
  this->GetOutput()->CopyFrom(this->GetInput()->GetPointerToVoxels());

  // Initialize exponential filter
  _ExponentialFilter->SetInput(0, this->GetOutput());
  _ExponentialFilter->SetInput(1, this->GetInput ());
  _ExponentialFilter->SetOutput(_dv);
  _ExponentialFilter->ComputeInverse(true);
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void irtkDisplacementToVelocityFieldBCH<VoxelType>::Finalize()
{
  // Deallocate intermediate images
  Delete(_dv);
  Delete(_l1);
  Delete(_l2);
  Delete(_l3);
  Delete(_l4);
  // Finalize base class
  irtkDisplacementToVelocityField<VoxelType>::Finalize();
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void irtkDisplacementToVelocityFieldBCH<VoxelType>::Run()
{
  IRTK_START_TIMING();

  // Do the initial set up
  this->Initialize();

  // Get pointer to output (AFTER initialization of base class!)
  irtkGenericImage<VoxelType> *v = this->GetOutput();

  #ifdef USE_CUDA
  irtkCUGenericImage<VoxelType> *d_dv = dynamic_cast<irtkCUGenericImage<VoxelType> *>(_dv);
  irtkCUGenericImage<VoxelType> *d_l1 = dynamic_cast<irtkCUGenericImage<VoxelType> *>(_l1);
  irtkCUGenericImage<VoxelType> *d_l2 = dynamic_cast<irtkCUGenericImage<VoxelType> *>(_l2);
  irtkCUGenericImage<VoxelType> *d_l3 = dynamic_cast<irtkCUGenericImage<VoxelType> *>(_l3);
  irtkCUGenericImage<VoxelType> *d_l4 = dynamic_cast<irtkCUGenericImage<VoxelType> *>(_l4);
  irtkCUGenericImage<VoxelType> *d_v  = dynamic_cast<irtkCUGenericImage<VoxelType> *>( v);
  const bool gpu = (use_gpu && d_v && d_dv) && (_NumberOfTerms < 3 || d_l1)
                                            && (_NumberOfTerms < 4 || d_l2)
                                            && (_NumberOfTerms < 5 || d_l3)
                                            && (_NumberOfTerms < 6 || d_l4)
  #endif

  // Iteratively update velocity field using Baker-Campbell-Hausdorff formula
  for (int n = 0; n < _NumberOfIterations; n++) {
    // Compute dv = exp(-v) ° phi
    _ExponentialFilter->Run();
    // Smooth to stabilize computation
    // (see Insight Journal article of Vercauteren et al. at http://hdl.handle.net/10380/3060)
    if (_SmoothVelocities) {
      irtkGaussianBlurring<VoxelType> blur(2.0 * v->GetXSize(), 2.0 * v->GetYSize(),
                                           v->GetZ() > 1 ? 2.0 * v->GetZSize() : .0);
      //blur.SetInput (v);
      //blur.SetOutput(v);
      //blur.Run();
      blur.SetInput (_dv);
      blur.SetOutput(_dv);
      blur.Run();
    }
    // Calculate required Lie brackets
    if (_NumberOfTerms > 2) liebracket(_l1,  v,  _dv, _UseJacobian); //          [v, dv]
    if (_NumberOfTerms > 3) liebracket(_l2,  v,  _l1, _UseJacobian); // [v,      [v, dv]]
    if (_NumberOfTerms > 4) liebracket(_l3, _dv, _l1, _UseJacobian); // [dv,     [v, dv]]
    if (_NumberOfTerms > 5) liebracket(_l4, _dv, _l2, _UseJacobian); // [dv, [v, [v, dv]]]
    // Compute update using the BCH formula
    #ifdef USE_CUDA
    if (gpu) {
      irtkCUDisplacementToVelocityFieldBCH::Update(d_v, d_dv, d_l1, d_l2, d_l3, d_l4, _NumberOfTerms);
    } else
    #endif
    {
      irtkNaryVoxelFunction::EvaluateBCHFormula bch;
      if      (_NumberOfTerms == 2) ParallelForEachScalar(v, _dv,                     v, bch);
      else if (_NumberOfTerms == 3) ParallelForEachScalar(v, _dv, _l1,                v, bch);
      else if (_NumberOfTerms == 4) ParallelForEachScalar(v, _dv, _l1, _l2,           v, bch);
      else if (_NumberOfTerms == 5) ParallelForEachScalar(v, _dv, _l1, _l2, _l3,      v, bch);
      else if (_NumberOfTerms == 6) ParallelForEachScalar(v, _dv, _l1, _l2, _l3, _l4, v, bch);
    }
  }

  // Do the final cleaning up
  this->Finalize();

  IRTK_DEBUG_TIMING(1, "irtkDisplacementToVelocityFieldBCH");
}

// ===========================================================================
// Explicit template instantiations
// ===========================================================================

template class irtkDisplacementToVelocityFieldBCH<float>;
template class irtkDisplacementToVelocityFieldBCH<double>;
