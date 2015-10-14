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
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkVertex.h>

#include <irtkPointSetForce.h>
#include <irtkPolyDataUtils.h>
using namespace irtk::polydata;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
void irtkPointSetForce::AllocateGradient(int n)
{
  if (_GradientSize < n || n <= 0) {
    Deallocate(_Gradient);
    if (n > 0) {
      _Gradient     = Allocate<GradientType>(n);
      _GradientSize = n;
    } else {
      _GradientSize = 0;
    }
  }
}

// -----------------------------------------------------------------------------
void irtkPointSetForce::AllocateCount(int n)
{
  if (_CountSize < n || n <= 0) {
    Deallocate(_Count);
    if (n > 0) {
      _Count     = Allocate<int>(n);
      _CountSize = n;
    } else {
      _CountSize = 0;
    }
  }
}

// -----------------------------------------------------------------------------
void irtkPointSetForce::Copy(const irtkPointSetForce &other)
{
  _PointSet      = other._PointSet;
  _SurfaceForce  = other._SurfaceForce;
  _InitialUpdate = other._InitialUpdate;
  AllocateGradient(other._GradientSize);
  AllocateCount(other._CountSize);
}

// -----------------------------------------------------------------------------
irtkPointSetForce::irtkPointSetForce(const char *name, double weight)
:
  irtkEnergyTerm(name, weight),
  _PointSet(NULL),
  _SurfaceForce(false),
  _Gradient(NULL),
  _GradientSize(0),
  _Count(NULL),
  _CountSize(0),
  _InitialUpdate(false)
{
}

// -----------------------------------------------------------------------------
irtkPointSetForce::irtkPointSetForce(const irtkPointSetForce &other)
:
  irtkEnergyTerm(other),
  _SurfaceForce(false),
  _Gradient(NULL),
  _GradientSize(0),
  _Count(NULL),
  _CountSize(0)
{
  Copy(other);
}

// -----------------------------------------------------------------------------
irtkPointSetForce &irtkPointSetForce::operator =(const irtkPointSetForce &other)
{
  irtkEnergyTerm::operator =(other);
  Copy(other);
  return *this;
}

// -----------------------------------------------------------------------------
irtkPointSetForce::~irtkPointSetForce()
{
  Deallocate(_Gradient);
  Deallocate(_Count);
}

// =============================================================================
// Point set attributes
// =============================================================================

// -----------------------------------------------------------------------------
void irtkPointSetForce::AddPointData(const char *name, vtkSmartPointer<vtkDataArray> &data)
{
  // Remove previously added array if any
  RemovePointData(name);

  // Remove array from point set attributes to prevent duplicate additions
  vtkPointData *pd;
  if (_SurfaceForce) pd = _PointSet->Surface ()->GetPointData();
  else               pd = _PointSet->PointSet()->GetPointData();
  for (int i = 0; i < pd->GetNumberOfArrays(); ++i) {
    if (pd->GetArray(i) == data) pd->RemoveArray(i--);
  }

  // Set unique array name
  string prefix  = ParameterNameWithPrefix(name);
  string unique  = prefix;
  bool is_unique = true;
  for (int i = 0; i < pd->GetNumberOfArrays(); ++i) {
    if (unique == pd->GetArrayName(i)) {
      is_unique = false;
      break;
    }
  }
  if (!is_unique) {
    for (int j = 1; j <= 99; ++j) {
      unique = prefix + ToString(j);
      is_unique = true;
      for (int i = 0; i < pd->GetNumberOfArrays(); ++i) {
        if (unique == pd->GetArrayName(i)) {
          is_unique = false;
          break;
        }
      }
      if (is_unique) break;
    }
    if (!is_unique) unique = prefix + "X";
  }
  data->SetName(unique.c_str());

  // Add point data array
  pd->AddArray(data);
  _PointDataName[name] = data->GetName();
}

// -----------------------------------------------------------------------------
vtkDataArray *irtkPointSetForce::AddPointData(const char *name, int c, int type)
{
  vtkSmartPointer<vtkDataArray> data = GetPointData(name, true);
  if (!data || data->GetArrayType() != type || data->GetNumberOfComponents() != c) {
    switch (type) {
      case VTK_CHAR:           data = vtkSmartPointer<vtkCharArray>::New();
      case VTK_UNSIGNED_CHAR:  data = vtkSmartPointer<vtkUnsignedCharArray>::New();
      case VTK_SHORT:          data = vtkSmartPointer<vtkShortArray>::New();
      case VTK_UNSIGNED_SHORT: data = vtkSmartPointer<vtkUnsignedShortArray>::New();
      case VTK_INT:            data = vtkSmartPointer<vtkIntArray>::New();
      case VTK_UNSIGNED_INT:   data = vtkSmartPointer<vtkUnsignedIntArray>::New();
      case VTK_DOUBLE:         data = vtkSmartPointer<vtkDoubleArray>::New();
      default:                 data = vtkSmartPointer<vtkFloatArray>::New();
    }
    data->SetNumberOfComponents(c);
  }
  if (_NumberOfPoints > 0) data->SetNumberOfTuples(_NumberOfPoints);
  AddPointData(name, data);
  return data;
}

// -----------------------------------------------------------------------------
void irtkPointSetForce::RemovePointData(const char *name)
{
  NameMapIterator it = _PointDataName.find(name);
  if (it == _PointDataName.end()) return;
  vtkPointData *pd;
  if (_SurfaceForce) pd = _PointSet->Surface ()->GetPointData();
  else               pd = _PointSet->PointSet()->GetPointData();
  pd->RemoveArray(it->second.c_str());
  _PointDataName.erase(it);
}

// -----------------------------------------------------------------------------
vtkDataArray *irtkPointSetForce::GetPointData(const char *name, bool optional) const
{
  vtkDataArray *data = NULL;
  NameMapConstIterator it = _PointDataName.find(name);
  if (it != _PointDataName.end()) {
    vtkPointData *pd;
    if (_SurfaceForce) pd = _PointSet->Surface ()->GetPointData();
    else               pd = _PointSet->PointSet()->GetPointData();
    data = pd->GetArray(it->second.c_str());
    if (data) return data;
  }
  if (!optional) {
    cerr << "irtkPointSetForce::GetPointData: Point data array has invalid size!" << endl;
    cerr << "  This indicates that the point data array was not correctly adjusted" << endl;
    cerr << "  during the remeshing of the deformed point set/surface. Please report" << endl;
    cerr << "  this bug or debug the program execution in order to fix this issue." << endl;
    exit(1);
  }
  return NULL;
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPoints> irtkPointSetForce::GetInitialPoints() const
{
  vtkSmartPointer<vtkPoints> points;
  double p[3];

  const irtkMultiLevelTransformation *mffd;
  mffd = dynamic_cast<const irtkMultiLevelTransformation *>(_PointSet->Transformation());

  if (_SurfaceForce) {

    if (mffd) {
      points = vtkSmartPointer<vtkPoints>::New();
      points->SetNumberOfPoints(_PointSet->NumberOfSurfacePoints());
      for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i) {
        _PointSet->GetInputSurfacePoint(i, p);
        mffd->GlobalTransform(p[0], p[1], p[2]);
        points->SetPoint(i, p);
      }
    } else {
      points = vtkSmartPointer<vtkPoints>::New();
      points->DeepCopy(_PointSet->InputSurface()->GetPoints());
    }

  } else {

    if (mffd) {
      points = vtkSmartPointer<vtkPoints>::New();
      points->SetNumberOfPoints(_PointSet->NumberOfPoints());
      for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i) {
        _PointSet->GetInputPoint(i, p);
        mffd->GlobalTransform(p[0], p[1], p[2]);
        points->SetPoint(i, p);
      }
    } else {
      points = vtkSmartPointer<vtkPoints>::New();
      points->DeepCopy(_PointSet->InputPointSet()->GetPoints());
    }

  }

  return points;
}

// -----------------------------------------------------------------------------
void irtkPointSetForce::Init()
{
  // Get number of mesh vertices
  if (_SurfaceForce) {
    _NumberOfPoints = _PointSet->NumberOfSurfacePoints();
  } else {
    _NumberOfPoints = _PointSet->NumberOfPoints();
  }

  // Allocate gradient vector
  // Note: _Count must be first allocated by subclass if required
  AllocateGradient(_NumberOfPoints);
  if (_CountSize > 0) AllocateCount(_NumberOfPoints);
}

// -----------------------------------------------------------------------------
void irtkPointSetForce::Initialize()
{
  // Initialize base class
  irtkEnergyTerm::Initialize();

  // Indicates also that force is inactive
  _NumberOfPoints = 0;

  // Free previously allocated memory
  AllocateGradient(0);
  AllocateCount(0);

  // Check input
  if (!_PointSet) {
    cerr << "irtkPointSetForce::Initialize: Input point set not set" << endl;
    exit(1);
  }
  if (!_PointSet->InputPointSet()) {
    cerr << "irtkPointSetForce::Initialize: Undeformed point set not set" << endl;
    exit(1);
  }

  // Initialize this class
  irtkPointSetForce::Init();

  // Delayed initialization upon next Update call
  _InitialUpdate = true;
}

// -----------------------------------------------------------------------------
void irtkPointSetForce::Reinitialize()
{
  irtkPointSetForce::Init();
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void irtkPointSetForce::Update(bool)
{
  if (_InitialUpdate || _PointSet->Transformation()) {
    _PointSet->Update(_InitialUpdate && _PointSet->SelfUpdate());
  }
  _InitialUpdate = false;
}

// -----------------------------------------------------------------------------
void irtkPointSetForce::EvaluateGradient(double *gradient, double, double weight)
{
  if (_NumberOfPoints > 0) {
    if (_PointSet->Transformation()) {
      const double t0 = _PointSet->InputTime();
      const double t  = _PointSet->Time();
      const irtkPointSet &pos = _SurfaceForce ? _PointSet->InputSurfacePoints()
                                              : _PointSet->InputPoints();
      _PointSet->Transformation()->ParametricGradient(pos, _Gradient, gradient, t, t0, weight);
    } else {
      double         *g;
      vtkIdTypeArray *origPtIds = _PointSet->OriginalSurfacePointIds();
      if (_SurfaceForce && origPtIds) {
        for (int i = 0, j; i < _NumberOfPoints; ++i) {
          j = static_cast<int>(origPtIds->GetComponent(i, 0));
          g = gradient + 3 * j;
          g[0] += weight * _Gradient[i]._x;
          g[1] += weight * _Gradient[i]._y;
          g[2] += weight * _Gradient[i]._z;
        }
      } else {
        double *g = gradient;
        for (int i = 0; i < _NumberOfPoints; ++i, g += 3) {
          g[0] += weight * _Gradient[i]._x;
          g[1] += weight * _Gradient[i]._y;
          g[2] += weight * _Gradient[i]._z;
        }
      }
    }
  }
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
void irtkPointSetForce::WriteDataSets(const char *p, const char *suffix, bool all) const
{
  if (_NumberOfPoints == 0 && !all) return;

  const int   sz = 1024;
  char        fname[sz];
  string      _prefix = Prefix(p);
  const char *prefix  = _prefix.c_str();

  if (_SurfaceForce) {
    snprintf(fname, sz, "%ssurface%s.vtp", prefix, suffix);
    WritePolyData(fname, _PointSet->Surface());
  } else {
    snprintf(fname, sz, "%spointset%s%s", prefix, suffix, _PointSet->DefaultExtension());
    _PointSet->Write(fname);
  }
}

// -----------------------------------------------------------------------------
inline vtkSmartPointer<vtkFloatArray> ToFloatArray(const irtkPointSetForce::GradientType *v, int n)
{
  vtkSmartPointer<vtkFloatArray> array = vtkSmartPointer<vtkFloatArray>::New();
  array->SetNumberOfComponents(3);
  array->SetNumberOfTuples(n);
  for (int i = 0; i < n; ++i) array->SetTuple3(i, v[i]._x, v[i]._y, v[i]._z);
  return array;
}

// -----------------------------------------------------------------------------
void irtkPointSetForce::WriteGradient(const char *p, const char *suffix) const
{
  if (_NumberOfPoints == 0) return;

  const int   sz = 1024;
  char        fname[sz];
  string      _prefix = Prefix(p);
  const char *prefix  = _prefix.c_str();

  vtkSmartPointer<vtkDataArray> gradient = ToFloatArray(_Gradient, _NumberOfPoints);
  gradient->SetName("gradient");

  if (_SurfaceForce) {

    vtkSmartPointer<vtkPoints> points;
    points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(_NumberOfPoints);

    vtkSmartPointer<vtkCellArray> vertices;
    vertices = vtkSmartPointer<vtkCellArray>::New();
    vertices->Allocate(_NumberOfPoints);

    double pt[3];
    for (vtkIdType i = 0; i < _NumberOfPoints; ++i) {
      _PointSet->GetInputSurfacePoint(i, pt);
      points->SetPoint(i, pt);
      vertices->InsertNextCell(1, &i);
    }

    vtkSmartPointer<vtkPolyData> output;
    output = vtkSmartPointer<vtkPolyData>::New();
    output->SetPoints(points);
    output->SetVerts(vertices);
    output->GetPointData()->AddArray(gradient);

    snprintf(fname, sz, "%sgradient%s.vtp", prefix, suffix);
    WritePolyData(fname, output);

  } else {

    snprintf(fname, sz, "%sgradient%s%s", prefix, suffix, _PointSet->DefaultExtension());
    _PointSet->Write(fname, gradient);

  }
}
