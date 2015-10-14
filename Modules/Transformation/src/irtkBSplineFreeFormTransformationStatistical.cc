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

#include <irtkTransformation.h>

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkBSplineFreeFormTransformationStatistical
::irtkBSplineFreeFormTransformationStatistical()
{
}

// -----------------------------------------------------------------------------
irtkBSplineFreeFormTransformationStatistical
::irtkBSplineFreeFormTransformationStatistical(const irtkImageAttributes    &attr,
                                               irtkVector3D<DOFStatus>   ****status,
                                               const irtkMatrix             &mat,
                                               const irtkVector             &vec)
:
  _BasisVectors(mat),
  _MeanVector  (vec)
{
  // Initialize DOFs
  InitializeDOFs(mat.Cols());

  // Initialize control points with dofs flag set to false, i.e.,
  // specifying that in case of this transformation, the CPs are *not* the DOFs
  InitializeCPs(attr, false);

  // Copy control point status
  for (int k = 0; k < attr._z; k++) {
    for (int j = 0; j < attr._y; j++) {
      for (int i = 0; i < attr._x; i++) {
        _CPStatus[0][k][j][i] = status[0][k][j][i];
      }
    }
  }
}

// -----------------------------------------------------------------------------
irtkBSplineFreeFormTransformationStatistical
::irtkBSplineFreeFormTransformationStatistical(const irtkBSplineFreeFormTransformationStatistical &ffd)
:
  _BasisVectors(ffd._BasisVectors),
  _MeanVector  (ffd._MeanVector)
{
  // Initialize DOFs
  InitializeDOFs(ffd._BasisVectors.Cols());
  const irtkImageAttributes &attr = ffd.Attributes();
  if (attr._x > 0 && attr._y > 0 && attr._z > 0) {
    // Initialize control points with dofs flag set to false, i.e.,
    // specifying that in case of this transformation, the CPs are *not* the DOFs
    InitializeCPs(attr, false);
  }
}

// -----------------------------------------------------------------------------
irtkBSplineFreeFormTransformationStatistical::~irtkBSplineFreeFormTransformationStatistical()
{
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationStatistical::Initialize(const irtkImageAttributes &)
{
  // For the moment, initialize the trasformation as one with no DOFs
  // until the statistical deformation model is read from disk
  InitializeDOFs(0);
}

// =============================================================================
// Approximation/Interpolation
// =============================================================================

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationStatistical
::ApproximateDOFs(const double *x,  const double *y,  const double *z,  const double *t,
                  const double *dx, const double *dy, const double *dz, int no)
{
  // Approximate control point displacements
  irtkBSplineFreeFormTransformation3D::ApproximateDOFs(x, y, z, t, dx, dy, dz, no);
  // Update transformation parameters
  UpdateDOFs();
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationStatistical
::ApproximateDOFsGradient(const double *, const double *, const double *, const double *,
                          const double *, const double *, const double *,
                          int, double *, double) const
{
  cerr << this->NameOfClass() << "::ApproximateDOFsGradient: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationStatistical
::Interpolate(const double *dx, const double *dy, const double *dz)
{
  // Interpolate control point displacements
  irtkBSplineFreeFormTransformation3D::Interpolate(dx, dy, dz);
  // Update transformation parameters
  UpdateDOFs();
}

// =============================================================================
// Lattice
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkBSplineFreeFormTransformationStatistical::CropPadPassiveCPs(int, int, int, int, bool)
{
  // Do nothing as the FFD lattice attributes depend on the statistical
  // deformation model and may not be modified once such model was loaded
  return true;
}

// =============================================================================
// Updating
// =============================================================================

// -----------------------------------------------------------------------------
struct irtkBSplineFreeFormTransformationStatisticalCPUpdate
{
  irtkBSplineFreeFormTransformationStatistical *_FFD;
  irtkVector3D<double> *_Data;
  const irtkVector     *_Mean;
  const irtkMatrix     *_Bases;
  double               *_Input;
  int                   _NumberOfDOFs;

  void operator()(const blocked_range<int> &re) const
  {
    int    cp_index_x, cp_index_y, cp_index_z;
    double val_x, val_y, val_z;

    for (int cp = re.begin(); cp != re.end(); ++cp) {
      _FFD->IndexToDOFs(cp, cp_index_x, cp_index_y, cp_index_z);
      val_x = _Mean->Get(cp_index_x);
      val_y = _Mean->Get(cp_index_y);
      val_z = _Mean->Get(cp_index_z);
      for(int q = 0; q < _NumberOfDOFs; q++) {
        val_x += _Bases->Get(cp_index_x, q) * _Input[q];
        val_y += _Bases->Get(cp_index_y, q) * _Input[q];
        val_z += _Bases->Get(cp_index_z, q) * _Input[q];
      }
      _Data[cp]._x = val_x;
      _Data[cp]._y = val_y;
      _Data[cp]._z = val_z;
    }
  }
};

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationStatistical::UpdateCPs()
{
  irtkBSplineFreeFormTransformationStatisticalCPUpdate updater;
  updater._FFD          = this;
  updater._Data         = _CPImage.Data();
  updater._Mean         = &_MeanVector;
  updater._Bases        = &_BasisVectors;
  updater._NumberOfDOFs = _NumberOfDOFs;
  updater._Input        = _Param;

  blocked_range<int> range(0, this->NumberOfCPs());
  parallel_for(range, updater);
}

// -----------------------------------------------------------------------------
struct irtkBSplineFreeFormTransformationStatisticalDOFUpdate
{
  irtkBSplineFreeFormTransformationStatistical *_FFD;
  irtkVector3D<double> *_Data;
  const irtkVector     *_Mean;
  const irtkMatrix     *_Bases;
  double               *_Output;
  int                   _NumberOfCPs;

  void operator()(const blocked_range<int> &re) const
  {
    int cp_index_x, cp_index_y, cp_index_z;

    for (int q = re.begin(); q != re.end(); ++q) {
      for(int cp = 0; cp < _NumberOfCPs; cp++) {
        _FFD->IndexToDOFs(cp, cp_index_x, cp_index_y, cp_index_z);
        _Output[q] += _Bases->Get(cp_index_x, q) * (_Data[cp]._x - _Mean->Get(cp_index_x)) +
                      _Bases->Get(cp_index_y, q) * (_Data[cp]._y - _Mean->Get(cp_index_y)) +
                      _Bases->Get(cp_index_z, q) * (_Data[cp]._z - _Mean->Get(cp_index_z));
      }
    }
  }
};

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationStatistical::UpdateDOFs()
{
  // Clear DOFs first
  memset(_Param, 0, _NumberOfDOFs * sizeof(DOFValue));

  irtkBSplineFreeFormTransformationStatisticalDOFUpdate updater;
  updater._FFD          = this;
  updater._Data         = _CPImage.Data();
  updater._Mean         = &_MeanVector;
  updater._Bases        = &_BasisVectors;
  updater._NumberOfCPs  = this->NumberOfCPs();
  updater._Output       = _Param;

  blocked_range<int> range(0, this->NumberOfDOFs());
  parallel_for(range, updater);
}

// =============================================================================
// Parameters (non-DoFs)
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkBSplineFreeFormTransformationStatistical::Set(const char *name, const char *value)
{
  if (strcmp(name, "Statistical deformation model file") == 0) {
    ReadSDM(value);
    return true;
  }
  return irtkTransformation::Set(name, value);
}

// -----------------------------------------------------------------------------
irtkParameterList irtkBSplineFreeFormTransformationStatistical::Parameter() const
{
  irtkParameterList params = irtkBSplineFreeFormTransformation3D::Parameter();
  if (!_ModelFile.empty()) {
    Insert(params, "Statistical deformation model file", _ModelFile);
  }
  return params;
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
struct irtkBSplineFreeFormTransformationStatisticalCPGradientProjection
{
  const irtkMatrix *_Bases;
  int               _NumberOfCPs;
  double           *_Input;
  double           *_Output;

  void operator()(const blocked_range<int> &re) const
  {
    for (int q = re.begin(); q != re.end(); ++q) {
      for (int p = 0; p != 3 * _NumberOfCPs; ++p) {
        _Output[q] += _Input[p] * _Bases->Get(p, q);
      }
    }
  }
};

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationStatistical
::ParametricGradient(const irtkGenericImage<double> *in,  double *out,
                     const irtkWorldCoordsImage *i2w, const irtkWorldCoordsImage *wc,
                     double t0, double w) const
{
  // Compute parametric gradient w.r.t control point coefficients
  double *tmp_gradient = CAllocate<double>(3 * this->NumberOfCPs());
  memset(tmp_gradient, 0, 3 * this->NumberOfCPs() * sizeof(double));
  irtkBSplineFreeFormTransformation3D::ParametricGradient(in, tmp_gradient, i2w, wc, t0, w);

  // Apply chain rule to obtain gradient w.r.t statistical parameters
  irtkBSplineFreeFormTransformationStatisticalCPGradientProjection proj;
  proj._Bases       = &_BasisVectors;
  proj._NumberOfCPs = this->NumberOfCPs();
  proj._Input       = tmp_gradient;
  proj._Output      = out;

  blocked_range<int> range(0, this->NumberOfDOFs());
  parallel_for(range, proj);

  Deallocate(tmp_gradient);
}

// =============================================================================
// Properties
// =============================================================================

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationStatistical
::BendingEnergyGradient(double *gradient, double w, bool incl_passive, bool wrt_world) const
{
  // Compute bending gradient w.r.t control point coefficients
  double *tmp_gradient = CAllocate<double>(3 * this->NumberOfCPs());
  memset(tmp_gradient, 0, 3 * this->NumberOfCPs() * sizeof(double));
  irtkBSplineFreeFormTransformation3D::BendingEnergyGradient(tmp_gradient, w, incl_passive, wrt_world);

  // Apply chain rule to obtain gradient w.r.t statistical parameters
  irtkBSplineFreeFormTransformationStatisticalCPGradientProjection proj;
  proj._Bases       = &_BasisVectors;
  proj._NumberOfCPs = this->NumberOfCPs();
  proj._Input       = tmp_gradient;
  proj._Output      = gradient;

  blocked_range<int> range(0, this->NumberOfDOFs());
  parallel_for(range, proj);

  Deallocate(tmp_gradient);
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationStatistical::Print(irtkIndent indent) const
{
  cout << indent << "3D Statistical B-spline FFD:" << endl;
  irtkFreeFormTransformation3D::Print(indent + 1);
}

// -----------------------------------------------------------------------------
irtkCofstream &irtkBSplineFreeFormTransformationStatistical::Write(irtkCofstream &to) const
{
  // Write magic no. for transformations
  unsigned int magic_no = IRTKTRANSFORMATION_MAGIC;
  to.WriteAsUInt(&magic_no, 1);

  // Write transformation type (write it as a 3D Bspline FFD)
  unsigned int trans_type = IRTKTRANSFORMATION_BSPLINE_FFD_3D_v3;
  to.WriteAsUInt(&trans_type, 1);

  return this->WriteDOFs(to);
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationStatistical::ReadSDM(const char *file)
{
  char type[8], version[5];

  irtkCifstream from;
  from.Open(file);

  // Read header with file format version information
  from.ReadAsChar(type,    8);
  from.ReadAsChar(version, 5);

  if (strncmp(type, "irtkSDM", 8) != 0) {
    cerr << "The file '" << file << "' is not a valid statistical deformation model!" << endl;
    exit(1);
  }

  // Read lattice attributes
  irtkImageAttributes attr;

  from.ReadAsInt(&attr._x, 1);
  from.ReadAsInt(&attr._y, 1);
  from.ReadAsInt(&attr._z, 1);
  attr._t = 1;

  from.ReadAsDouble(attr._xaxis, 3);
  from.ReadAsDouble(attr._yaxis, 3);
  from.ReadAsDouble(attr._zaxis, 3);

  from.ReadAsDouble(&attr._dx, 1);
  from.ReadAsDouble(&attr._dy, 1);
  from.ReadAsDouble(&attr._dz, 1);
  attr._dt = 0;

  from.ReadAsDouble(&attr._xorigin, 1);
  from.ReadAsDouble(&attr._yorigin, 1);
  from.ReadAsDouble(&attr._zorigin, 1);

  // Initialize control points whith dofs flag in false, i.e.,
  // specifying that in our case CPs are *not* DOFs
  InitializeCPs(attr, false);

  // Read status
  from.ReadAsInt(reinterpret_cast<int *>(_CPStatus[0][0][0]),
                 3 * attr._x * attr._y * attr._z);

  // Read matrix and vector
  from >> _BasisVectors >> _MeanVector;

  if (_BasisVectors.Rows() != _MeanVector.Rows()) {
    cerr << "Invalid statistical deformation model! "
            "Basis vectors and mean vector must have the same number of elements." << endl;
    exit(1);
  }

  if (_MeanVector.Rows() != 3 * attr._x * attr._y * attr._z) {
    cerr << "Incompatible statistical deformation model! "
            "The number of vector elements must be equal to 3 times the number of control points." << endl;
    exit(1);
  }

  // Reinitialize DOFs now that we know how many of them we need
  InitializeDOFs(_BasisVectors.Cols());

  //Initialize the control points in order to be consistent with the newly read SDM
  UpdateCPs();

  // Initialize interpolator
  InitializeInterpolator();

  // Remember model file name
  ModelFile(file);
}

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationStatistical::WriteSDM(const char *file)
{
  irtkCofstream to;
  to.Open(file);

  // Write header with file format version information
  to.WriteAsChar("irtkSDM", 8);
  to.WriteAsChar("v1.0",    5);

  // Write lattice attributes
  to.WriteAsInt(&_attr._x, 1);
  to.WriteAsInt(&_attr._y, 1);
  to.WriteAsInt(&_attr._z, 1);

  to.WriteAsDouble(_attr._xaxis, 3);
  to.WriteAsDouble(_attr._yaxis, 3);
  to.WriteAsDouble(_attr._zaxis, 3);

  to.WriteAsDouble(&_attr._dx, 1);
  to.WriteAsDouble(&_attr._dy, 1);
  to.WriteAsDouble(&_attr._dz, 1);

  to.WriteAsDouble(&_attr._xorigin, 1);
  to.WriteAsDouble(&_attr._yorigin, 1);
  to.WriteAsDouble(&_attr._zorigin, 1);

  to.WriteAsInt(reinterpret_cast<const int *>(_CPStatus[0][0][0]),
                3 * _attr._x * _attr._y * _attr._z);

  // Write statistical deformation model data
  to << _BasisVectors << _MeanVector;
}

// -----------------------------------------------------------------------------
irtkCifstream &irtkBSplineFreeFormTransformationStatistical
::ReadDOFs(irtkCifstream &from, irtkTransformationType format)
{
  from = irtkBSplineFreeFormTransformation3D::ReadDOFs(from, format);
  UpdateDOFs();
  return from;
}

// =============================================================================
// Others
// =============================================================================

// -----------------------------------------------------------------------------
void irtkBSplineFreeFormTransformationStatistical::Verify()
{
  if (_NumberOfDOFs == 0) {
    cerr << "No statistical deformation model was provided! "
            "Make sure 'Statistical deformation model file' is set on the parameter file" << endl;
    exit(1);
  }

  if (_BasisVectors.Norm() == 0) {
    cerr << "Invalid statistical deformation model: "
            "All basis verctors have zero norm!" << endl;
    exit(1);
  }
}
