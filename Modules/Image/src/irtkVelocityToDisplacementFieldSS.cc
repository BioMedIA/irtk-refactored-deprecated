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

#include <irtkImage.h>
#include <irtkImageFunction.h>
#include <irtkVoxelFunction.h>
#include <irtkInterpolateImageFunction.h>
#include <irtkGaussianPyramidFilter.h>
#include <irtkResampling.h>

#include <irtkVelocityToDisplacementField.h>

using irtkBinaryVoxelFunction::ComposeDisplacementFields2D;
using irtkBinaryVoxelFunction::ComposeDisplacementFields3D;


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkVelocityToDisplacementFieldSS<VoxelType>::irtkVelocityToDisplacementFieldSS()
:
  _NumberOfSquaringSteps   (0),
  _Upsample                (false),
  _SmoothBeforeDownsampling(false),
  _MaxScaledVelocity       (.0),
  _ExternalCache           (NULL),
  _Displacement            (NULL),
  _Interpolator            (NULL)
{
}

// -----------------------------------------------------------------------------
template <class VoxelType>
irtkVelocityToDisplacementFieldSS<VoxelType>::~irtkVelocityToDisplacementFieldSS()
{
  if (_Displacement != _ExternalCache) Delete(_Displacement);
  Delete(_Interpolator);
}

// =============================================================================
// Filter implementation
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
bool irtkVelocityToDisplacementFieldSS<VoxelType>::RequiresBuffering()
{
  // During initialization, the input velocity field is scaled and the resulting
  // scaled velocities stored in a temporary image which in turn is used during
  // the squaring step where the result of each squaring is written to the output
  // image. At this point, the input image is already no longer needed.
  return false;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkVelocityToDisplacementFieldSS<VoxelType>::Initialize()
{
  // Cleanup previous initialization if Finalize not called before
  if (_Displacement != _ExternalCache) Delete(_Displacement);
  Delete(_Interpolator);

  // Check parameters
  if (this->NumberOfSteps() < 1 && _NumberOfSquaringSteps && _MaxScaledVelocity <= .0) {
    cerr << this->NameOfClass() << "::Run: Either number of integration or squaring steps or maximum scaled velocity must be positive" << endl;
    exit(1);
  }

  // Initialize base class
  //
  // Circumvent failure during initialization of base class in case of
  // a non-positive number of steps which was accepted by this filter
  const int n = this->NumberOfSteps();
  this->NumberOfSteps(1);
  irtkVelocityToDisplacementField<VoxelType>::Initialize();
  this->NumberOfSteps(n);

  // Get number of voxels (incl. vector field components)
  const int nvox = this->GetInput()->GetNumberOfVoxels();

  // Sample input vector field
  //
  // According to
  //   Bossa, M., Zacur, E., & Olmos, S. (2008). Algorithms for
  //   computing the group exponential of diffeomorphisms: Performance evaluation.
  //   In 2008 IEEE Computer Society Conference on Computer Vision and Pattern
  //   Recognition Workshops (pp. 1–8). IEEE. doi:10.1109/CVPRW.2008.4563005
  // an upsampling of the input velocity field by a factor of r * 2^d, where
  // r is the interpolation radius (i.e., 1 in case of linear interpolation)
  // and d is the dimension of the vector field (i.e., 3 for 3D) would reduce
  // the error made by using the approximation that the exponential of the
  // scaled velocity field is equal to the scaled velocity itself and especially
  // the larger sampling error at later stages of the squaring steps.
  irtkImageAttributes attr = this->GetInput()->GetImageAttributes();
  if (_Upsample) {
    attr._x *= 2; attr._dx /= 2;
    attr._y *= 2; attr._dy /= 2;
    if (attr._z > 1) { attr._z *= 2; attr._dz /= 2; }
  }
  if (_ExternalCache) _Displacement = _ExternalCache;
  else                _Displacement = new irtkGenericImage<VoxelType>();
  _Displacement->Initialize(this->GetInput()->GetImageAttributes());
  _Interpolator = irtkInterpolateImageFunction::New(this->Interpolation(),
                                                    this->Extrapolation(),
                                                    this->GetInput());
  _Interpolator->SetInput  (this->GetInput());
  _Interpolator->Initialize(this->_ComputeInterpolationCoefficients == false);
  _Interpolator->Evaluate  (*_Displacement);

  // During squaring, use interpolator to sample intermediate vector field.
  // Note that the interpolator must be re-initialized before each squaring
  // step in case it has to convert the vector field to interpolation coefficients.
  // Therefore, this is done during the squaring iterations.
  _Interpolator->SetInput(_Displacement);

  // Scale sampled input velocities by minimum number of steps
  _NumberOfSquaringSteps = this->NumberOfSquaringSteps();
  if (_NumberOfSquaringSteps <= 0) _NumberOfSquaringSteps = static_cast<int>(ceil(log(static_cast<double>(n)) / log(2.0)));
  if (_NumberOfSquaringSteps <  0) _NumberOfSquaringSteps = 0;

  double max   = .0;
  double scale = this->T() / pow(2.0, _NumberOfSquaringSteps);

  VoxelType *s = _Displacement->GetPointerToVoxels();
  for (int idx = 0; idx < nvox; ++idx, ++s) {
    (*s) *= scale;
    if (fabs(*s) > max) max = fabs(*s);
  }

  // Continue halfing input velocities as long as maximum absolute velocity
  // exceeds the specified maximum; skip if fixed number of steps
  if (_MaxScaledVelocity > .0) {
    scale = 1.0;
    while ((max * scale) > _MaxScaledVelocity) {
      scale *= 0.5;
      _NumberOfSquaringSteps++;
    }
    if (scale != 1.0) {
      VoxelType *s = _Displacement->GetPointerToVoxels();
      for (int idx = 0; idx < nvox; ++idx, ++s) {
        (*s) *= scale;
      }
    }
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkVelocityToDisplacementFieldSS<VoxelType>::Finalize()
{
  // Finalize base class
  irtkVelocityToDisplacementField<VoxelType>::Finalize();

  // Destroy temporary vector field
  if (_Displacement != _ExternalCache) Delete(_Displacement);

  // Destroy interpolator
  Delete(_Interpolator);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void irtkVelocityToDisplacementFieldSS<VoxelType>::Run()
{
  // Do the initial set up and scaling step
  // Note: MUST be done before getting pointer to output image as it may
  //       be replaced by a temporary buffer if input == output.
  this->Initialize();

  const irtkGenericImage<VoxelType> *din  = this->GetInput(1);
  irtkGenericImage<VoxelType>       *dout = this->GetOutput();

  // Copy input displacement field if it is used as output
  if (dout == din) din = new irtkGenericImage<VoxelType>(*din);

  // Get attributes of vector fields
  irtkImageAttributes attr = dout->GetImageAttributes();
  const int nbytes = dout->GetNumberOfVoxels() * sizeof(VoxelType);

  // 2D
  if (attr._t == 2) {

    // Squaring steps
    int n = _NumberOfSquaringSteps;
    while (n--) {
      _Interpolator->Initialize(false);
      ComposeDisplacementFields2D<VoxelType> square(_Interpolator, _Displacement);
      ParallelForEachVoxel(attr, _Displacement, dout, square);
      memcpy(_Displacement->GetPointerToVoxels(), dout->GetPointerToVoxels(), nbytes);
    }
    // Either compose resulting displacement field with input displacement field
    if (din) {
      _Interpolator->Initialize(false);
      dout->Initialize(din->GetImageAttributes(), 2);
      attr = dout->GetImageAttributes();
      ComposeDisplacementFields2D<VoxelType> compose(_Interpolator, din);
      ParallelForEachVoxel(attr, din, dout, compose);
    // or downsample it if scaled velocity field was upsampled before
    } else if (_Upsample) {
      if (_SmoothBeforeDownsampling) {
        irtkGaussianPyramidFilter<VoxelType> downsampler(1);
        downsampler.SetInput (_Displacement);
        downsampler.SetOutput(dout);
        downsampler.Run();
      } else {
        _Interpolator->Initialize();
        irtkResampling<VoxelType> downsampler(this->GetInput()->GetX(),
                                              this->GetInput()->GetY(),
                                              this->GetInput()->GetZ());
        downsampler.SetInterpolator(_Interpolator);
        downsampler.SetInput       (_Displacement);
        downsampler.SetOutput      (dout);
        downsampler.Run();
      }
    }

  // 3D
  } else {

    // Squaring steps
    int n = _NumberOfSquaringSteps;
    while (n--) {
      _Interpolator->Initialize(false);
      ComposeDisplacementFields3D<VoxelType> square(_Interpolator, _Displacement);
      ParallelForEachVoxel(attr, _Displacement, dout, square);
      memcpy(_Displacement->GetPointerToVoxels(), dout->GetPointerToVoxels(), nbytes);
    }
    // Either compose resulting displacement field with input displacement field
    if (din) {
      _Interpolator->Initialize(false);
      dout->Initialize(din->GetImageAttributes(), 3);
      attr = dout->GetImageAttributes();
      ComposeDisplacementFields3D<VoxelType> compose(_Interpolator, din);
      ParallelForEachVoxel(attr, din, dout, compose);
    // or downsample it if scaled velocity field was upsampled before
    } else if (_Upsample) {
      if (_SmoothBeforeDownsampling) {
        irtkGaussianPyramidFilter<VoxelType> downsampler(1);
        downsampler.SetInput       (_Displacement);
        downsampler.SetOutput      (dout);
        downsampler.Run();
      } else {
        _Interpolator->Initialize();
        irtkResampling<VoxelType> downsampler(this->GetInput()->GetX(),
                                              this->GetInput()->GetY(),
                                              this->GetInput()->GetZ());
        downsampler.SetInterpolator(_Interpolator);
        downsampler.SetInput       (_Displacement);
        downsampler.SetOutput      (dout);
        downsampler.Run();
      }
    }

  }

  // Do the final cleaning up
  if (din != this->GetInput(1)) delete din;
  this->Finalize();
}

// =============================================================================
// Explicit template instantiations
// =============================================================================

template class irtkVelocityToDisplacementFieldSS<float>;
template class irtkVelocityToDisplacementFieldSS<double>;
