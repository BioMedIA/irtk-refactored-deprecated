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
#include <irtkResampling.h>

// ---------------------------------------------------------------------------
template <class VoxelType>
irtkResampling<VoxelType>::irtkResampling(double dx, double dy, double dz)
{
  _X     = 0;
  _Y     = 0;
  _Z     = 0;
  _XSize = dx;
  _YSize = dy;
  _ZSize = dz;
  _Interpolator = NULL;
}

// ---------------------------------------------------------------------------
template <class VoxelType>
irtkResampling<VoxelType>::irtkResampling(int x, int y, int z)
{
  _X     = x;
  _Y     = y;
  _Z     = z;
  _XSize = 0;
  _YSize = 0;
  _ZSize = 0;
  _Interpolator = NULL;
}

// ---------------------------------------------------------------------------
template <class VoxelType>
irtkResampling<VoxelType>::irtkResampling(int x, int y, int z, double dx, double dy, double dz)
{
  _X     = x;
  _Y     = y;
  _Z     = z;
  _XSize = dx;
  _YSize = dy;
  _ZSize = dz;
  _Interpolator = NULL;
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void irtkResampling<VoxelType>::Initialize()
{
  // Initialize base class
  irtkImageToImage<VoxelType>::Initialize();

  // Set up interpolator
  if (_Interpolator == NULL) {
    cerr << "irtkResampling::Initialize: No interpolator found!" << endl;
    exit(1);
  }
  _Interpolator->SetInput(this->_input);
  _Interpolator->Initialize();

  // Initialize output image
  this->InitializeOutput();
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void irtkResampling<VoxelType>::InitializeOutput()
{
  irtkImageAttributes attr = this->_input->GetImageAttributes();

  if (_X > 0 && _XSize > 0) {
    attr._x  = _X;
    attr._dx = _XSize;
  } else if (_X > 0) {
    attr._x  = _X;
    attr._dx = this->_input->GetX() * this->_input->GetXSize() / static_cast<double>(_X);
  } else if (_XSize > 0) {
    attr._x  = round(this->_input->GetX() * this->_input->GetXSize() / this->_XSize);
    if (attr._x < 1) attr._x  = 1;
    else             attr._dx = _XSize;
  }

  if (_Y > 0 && _YSize > 0) {
    attr._y  = _Y;
    attr._dy = _YSize;
  } else if (_Y > 0) {
    attr._y  = _Y;
    attr._dy = this->_input->GetY() * this->_input->GetYSize() / static_cast<double>(_Y);
  } else if (_YSize > 0) {
    attr._y  = round(this->_input->GetY() * this->_input->GetYSize() / this->_YSize);
    if (attr._y < 1) attr._y  = 1;
    else             attr._dy = _YSize;
  }
  
  if (_Z > 0 && _ZSize > 0) {
    attr._z  = _Z;
    attr._dz = _ZSize;
  } else if (_Z > 0) {
    attr._z  = _Z;
    attr._dz = this->_input->GetZ() * this->_input->GetZSize() / static_cast<double>(_Z);
  } else if (_ZSize > 0) {
    attr._z  = round(this->_input->GetZ() * this->_input->GetZSize() / this->_ZSize);
    if (attr._z < 1) attr._z  = 1;
    else             attr._dz = _ZSize;
  }

  this->_output->Initialize(attr);
}

// ---------------------------------------------------------------------------
// Serial/Multi-threaded body of irtkResampling::Run method
class irtkResamplingRun
{
public:
  // -------------------------------------------------------------------------
  // Data members
  irtkImage         *_Input;
  irtkImage         *_Output;
  irtkImageFunction *_Interpolator;
  mutable int        _T;

  // -------------------------------------------------------------------------
  /// Default constructor
  irtkResamplingRun()
  :
    _Input(NULL),
    _Output(NULL),
    _Interpolator(NULL),
    _T(0)
  {
  }

  // -------------------------------------------------------------------------
  /// Copy constructor
  irtkResamplingRun(const irtkResamplingRun &other)
  :
    _Input(other._Input),
    _Output(other._Output),
    _Interpolator(other._Interpolator),
    _T(other._T)
  {
  }

  // -------------------------------------------------------------------------
  /// Resamples the image within given output region for frame _T
  void operator() (const blocked_range3d<int> &r) const
  {
    double x, y, z;
    for (int k = r.pages().begin(); k != r.pages().end(); k++) {
      for (int j = r.rows().begin(); j != r.rows().end(); j++) {
        for (int i = r.cols().begin(); i != r.cols().end(); i++) {
          x = i;
          y = j;
          z = k;
          _Output->ImageToWorld(x, y, z);
          _Input ->WorldToImage(x, y, z);
          _Output->PutAsDouble(i, j, k, _T, _Interpolator->Evaluate(x, y, z, _T));
        }
      }
    }
  }

  // -------------------------------------------------------------------------
  /// Resamples the image for the given range of output frames
  void operator() (const blocked_range<int> &r) const
  {
    blocked_range3d<int> voxels(0, _Output->GetZ(),
                                0, _Output->GetY(),
                                0, _Output->GetX());

    for (_T = r.begin(); _T != r.end(); _T++) {
      parallel_for(voxels, irtkResamplingRun(*this));
    }
  }

};

// ---------------------------------------------------------------------------
template <class VoxelType>
void irtkResampling<VoxelType>::Run()
{
  IRTK_START_TIMING();

  // Do the initial set up
  this->Initialize();

  // Resample image either serial or in parallel if enabled
  irtkResamplingRun body;
  body._Input        = this->_input;
  body._Output       = this->_output;
  body._Interpolator = this->_Interpolator;

  blocked_range<int> frames(0, this->_output->GetT());
  body(frames);

  // Do the final cleaning up
  this->Finalize();

  IRTK_DEBUG_TIMING(5, "irtkResampling::Run");
}

// ---------------------------------------------------------------------------
// Explicit template instantiations
template class irtkResampling<char>;
template class irtkResampling<unsigned char>;
template class irtkResampling<short>;
template class irtkResampling<unsigned short>;
template class irtkResampling<int>;
template class irtkResampling<float>;
template class irtkResampling<double>;
