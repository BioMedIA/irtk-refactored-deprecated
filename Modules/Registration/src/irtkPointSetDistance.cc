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

#include <vtkSmartPointer.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkVertex.h>

#include <irtkPointSetDistance.h>
#include <irtkPolyDataUtils.h>
using namespace irtk::polydata;

// =============================================================================
// Factory
// =============================================================================

#include <irtkFiducialRegistrationError.h>
#include <irtkPointCorrespondenceDistance.h>
#include <irtkCurrentsDistance.h>

// -----------------------------------------------------------------------------
irtkPointSetDistance *irtkPointSetDistance::New(irtkPointSetDistanceMeasure metric)
{
  switch (metric) {
    case PDM_FRE:                    return new irtkFiducialRegistrationError();
    case PDM_CorrespondenceDistance: return new irtkPointCorrespondenceDistance();
    case PDM_CurrentsDistance:       return new irtkCurrentsDistance();
    default:
      cerr << "irtkPointSetDistance::New: Unknown point set distance measure: " << metric << endl;
      exit(1);
  }
  return NULL;
}

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
void irtkPointSetDistance::AllocateGradientWrtTarget(int m)
{
  Deallocate(_GradientWrtTarget);
  if (m > 0) _GradientWrtTarget = Allocate<GradientType>(m);
}

// -----------------------------------------------------------------------------
void irtkPointSetDistance::AllocateGradientWrtSource(int n)
{
  Deallocate(_GradientWrtTarget);
  if (n > 0) _GradientWrtTarget = Allocate<GradientType>(n);
}

// -----------------------------------------------------------------------------
irtkPointSetDistance::irtkPointSetDistance(const char *name, double weight)
:
  irtkDataFidelity(name, weight),
  _Target(NULL),
  _Source(NULL),
  _GradientWrtTarget(NULL),
  _GradientWrtSource(NULL),
  _InitialUpdate    (false)
{
  _ParameterPrefix.push_back("Point set distance ");
  _ParameterPrefix.push_back("Polydata distance "); // legacy
}

// -----------------------------------------------------------------------------
irtkPointSetDistance::irtkPointSetDistance(const irtkPointSetDistance &other, int m, int n)
:
  irtkDataFidelity(other),
  _Target(other._Target),
  _Source(other._Source),
  _GradientWrtTarget(NULL),
  _GradientWrtSource(NULL),
  _InitialUpdate(other._InitialUpdate)
{
  AllocateGradientWrtTarget(m < 0 && _Target && other._GradientWrtTarget ? other._Target->NumberOfPoints() : 0);
  AllocateGradientWrtSource(n < 0 && _Source && other._GradientWrtSource ? other._Source->NumberOfPoints() : 0);
}

// -----------------------------------------------------------------------------
void irtkPointSetDistance::Copy(const irtkPointSetDistance &other, int m, int n)
{
  irtkDataFidelity::operator =(other);
  _Target        = other._Target;
  _Source        = other._Source;
  _InitialUpdate = other._InitialUpdate;
  AllocateGradientWrtTarget(m < 0 && _Target && other._GradientWrtTarget ? _Target->NumberOfPoints() : m);
  AllocateGradientWrtSource(n < 0 && _Source && other._GradientWrtSource ? _Source->NumberOfPoints() : n);
}

// -----------------------------------------------------------------------------
irtkPointSetDistance &irtkPointSetDistance::operator =(const irtkPointSetDistance &other)
{
  Copy(other);
  return *this;
}

// -----------------------------------------------------------------------------
irtkPointSetDistance::~irtkPointSetDistance()
{
  Deallocate(_GradientWrtTarget);
  Deallocate(_GradientWrtSource);
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
void irtkPointSetDistance::Initialize(int m, int n)
{
  // Initialize base class
  irtkDataFidelity::Initialize();

  // Check inputs
  if (_Target == NULL) {
    cerr << "irtkPointSetDistance::Initialize: Target dataset is NULL" << endl;
    exit(1);
  }
  if (_Source == NULL) {
    cerr << "irtkPointSetDistance::Initialize: Source dataset is NULL" << endl;
    exit(1);
  }

  // Force next update if data sets have their SelfUpdate flag set
  _InitialUpdate = true;

  // Allocate memory for non-parametric gradient
  Deallocate(_GradientWrtTarget);
  Deallocate(_GradientWrtSource);
  if (_Target->Transformation()) _GradientWrtTarget = Allocate<GradientType>(m);
  if (_Source->Transformation()) _GradientWrtSource = Allocate<GradientType>(n);
}

// -----------------------------------------------------------------------------
void irtkPointSetDistance::Initialize()
{
  // Check inputs
  if (_Target == NULL) {
    cerr << "irtkPointSetDistance::Initialize: Target dataset is NULL" << endl;
    exit(1);
  }
  if (_Source == NULL) {
    cerr << "irtkPointSetDistance::Initialize: Source dataset is NULL" << endl;
    exit(1);
  }

  // Initialize with allocation of memory for all points
  Initialize(_Target->NumberOfPoints(), _Source->NumberOfPoints());
}

// -----------------------------------------------------------------------------
void irtkPointSetDistance::Reinitialize(int m, int n)
{
  // Allocate memory for non-parametric gradient
  Deallocate(_GradientWrtTarget);
  Deallocate(_GradientWrtSource);
  if (_Target->Transformation()) _GradientWrtTarget = Allocate<GradientType>(m);
  if (_Source->Transformation()) _GradientWrtSource = Allocate<GradientType>(n);
}

// -----------------------------------------------------------------------------
void irtkPointSetDistance::Reinitialize()
{
  // Reinitialize with allocation of memory for all points
  Reinitialize(_Target->NumberOfPoints(), _Source->NumberOfPoints());
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void irtkPointSetDistance::Update(bool)
{
  if (_InitialUpdate || _Target->Transformation()) {
    _Target->Update(_InitialUpdate && _Target->SelfUpdate());
  }
  if (_InitialUpdate || _Source->Transformation()) {
    _Source->Update(_InitialUpdate && _Source->SelfUpdate());
  }
  _InitialUpdate = false;
}

// -----------------------------------------------------------------------------
void irtkPointSetDistance
::ParametricGradient(const irtkRegisteredPointSet *wrt_pset,
                     const irtkVector3D<double>   *np_gradient,
                     double                       *gradient,
                     double                        weight)
{
  const irtkTransformation * const T = wrt_pset->Transformation();
  irtkAssert(T != NULL, "point set is being transformed");
  const double t0 = wrt_pset->InputTime();
  const double t  = wrt_pset->Time();
  T->ParametricGradient(wrt_pset->InputPoints(), np_gradient, gradient, t, t0, weight);
}

// -----------------------------------------------------------------------------
void irtkPointSetDistance::EvaluateGradient(double *gradient, double, double weight)
{
  // Get transformations of input data sets
  const irtkTransformation * const T1 = _Target->Transformation();
  const irtkTransformation * const T2 = _Source->Transformation();
  // Compute parametric gradient w.r.t target transformation
  if (T1) {
    this->NonParametricGradient(_Target, _GradientWrtTarget);
    this->ParametricGradient   (_Target, _GradientWrtTarget, gradient, weight);
  }
  // If target and source are transformed by different transformations,
  // the gradient vector contains first the derivative values w.r.t the
  // parameters of the target transformation followed by those computed
  // w.r.t the parameters of the source transformation. Otherwise, if
  // both point sets are transformed by the same transformation, i.e., a
  // velocity based transformation integrated half way in both directions,
  // the derivative values are summed up instead.
  if (T1 && T2 && !HaveSameDOFs(T1, T2)) gradient += T2->NumberOfDOFs();
  // Compute parametric gradient w.r.t source transformation
  if (T2) {
    this->NonParametricGradient(_Source, _GradientWrtSource);
    this->ParametricGradient   (_Source, _GradientWrtSource, gradient, weight);
  }
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
inline vtkSmartPointer<vtkFloatArray> ToFloatArray(const irtkPointSetDistance::GradientType *v, int n)
{
  vtkSmartPointer<vtkFloatArray> array = vtkSmartPointer<vtkFloatArray>::New();
  array->SetNumberOfComponents(3);
  array->SetNumberOfTuples(n);
  for (int i = 0; i < n; ++i) array->SetTuple3(i, v[i]._x, v[i]._y, v[i]._z);
  return array;
}

// -----------------------------------------------------------------------------
void irtkPointSetDistance::WriteDataSets(const char *p, const char *suffix, bool all) const
{
  const int   sz = 1024;
  char        fname[sz];
  string      _prefix = Prefix(p);
  const char *prefix  = _prefix.c_str();

  if (_Target->Transformation() || all) {
    snprintf(fname, sz, "%starget%s%s", prefix, suffix, _Target->DefaultExtension());
    _Target->Write(fname);
  }
  if (_Source->Transformation() || all) {
    snprintf(fname, sz, "%ssource%s%s", prefix, suffix, _Source->DefaultExtension());
    _Source->Write(fname);
  }
}

// -----------------------------------------------------------------------------
void irtkPointSetDistance::WriteGradient(const char *p, const char *suffix) const
{
  const int   sz = 1024;
  char        fname[sz];
  string      _prefix = Prefix(p);
  const char *prefix  = _prefix.c_str();

  if (_GradientWrtTarget) {
    snprintf(fname, sz, "%sgradient_wrt_target%s.vtp", prefix, suffix);
    this->WriteGradient(fname, _Target, _GradientWrtTarget);
  }
  if (_GradientWrtSource) {
    snprintf(fname, sz, "%sgradient_wrt_source%s.vtp", prefix, suffix);
    this->WriteGradient(fname, _Source, _GradientWrtSource);
  }
}

// -----------------------------------------------------------------------------
void irtkPointSetDistance::WriteGradient(const char                   *fname,
                                         const irtkRegisteredPointSet *data,
                                         const GradientType           *g,
                                         const vector<int>            *sample) const
{
  bool samples_only = sample && !sample->empty();
  const vtkIdType n = (samples_only ? sample->size() : data->NumberOfPoints());

  vtkSmartPointer<vtkPoints>    points;
  vtkSmartPointer<vtkCellArray> vertices;
  vtkSmartPointer<vtkPolyData>  output;

  points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(n);

  vertices = vtkSmartPointer<vtkCellArray>::New();
  vertices->Allocate(n);

  double p[3];
  for (vtkIdType i = 0; i < n; ++i) {
    data->GetInputPoint((samples_only ? (*sample)[i] : i), p);
    points  ->SetPoint(i, p);
    vertices->InsertNextCell(1, &i);
  }

  vtkSmartPointer<vtkDataArray> gradient = ToFloatArray(g, n);
  gradient->SetName("gradient");

  output = vtkSmartPointer<vtkPolyData>::New();
  output->SetPoints(points);
  output->SetVerts(vertices);
  output->GetPointData()->AddArray(gradient);

  WritePolyData(fname, output);
}
