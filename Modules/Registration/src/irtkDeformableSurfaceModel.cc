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

#include <irtkDeformableSurfaceModel.h>

#include <irtkPolyDataRemeshing.h>
#include <irtkSurfaceCollisions.h>

#include <vtkCellData.h>

#include <irtkPolyDataUtils.h>
using namespace irtk::polydata;


// =============================================================================
// Auxiliary functor classes for parallel execution
// =============================================================================

namespace irtkDeformableSurfaceModelUtils {


// -----------------------------------------------------------------------------
/// Determine maximum vertex displacement
class MaxVertexDisplacement
{
private:

  const double *_Gradient;
  double        _MaxNorm;

public:

  /// Constructor
  MaxVertexDisplacement(const double *gradient)
  :
    _Gradient(gradient), _MaxNorm(.0)
  {}

  /// Copy constructor
  MaxVertexDisplacement(const MaxVertexDisplacement &other)
  :
    _Gradient(other._Gradient),
    _MaxNorm (other._MaxNorm)
  {}

  /// Split constructor
  MaxVertexDisplacement(const MaxVertexDisplacement &other, split)
  :
    _Gradient(other._Gradient),
    _MaxNorm (other._MaxNorm)
  {}

  /// Join results
  void join(const MaxVertexDisplacement &other)
  {
    if (other._MaxNorm > _MaxNorm) _MaxNorm = other._MaxNorm;
  }

  /// Maximum norm
  double Norm() const { return sqrt(_MaxNorm); }

  /// Determine maximum norm of specified vertex displacements
  void operator()(const blocked_range<int> &re)
  {
    double norm;
    const double *g = _Gradient + 3 * re.begin();
    for (int i = re.begin(); i != re.end(); ++i, g += 3) {
      norm = pow(g[0], 2) + pow(g[1], 2) + pow(g[2], 2);
      if (norm > _MaxNorm) _MaxNorm = norm;
    }
  }
};

// -----------------------------------------------------------------------------
/// Move points of deformable surface model
struct MovePoints
{
  vtkPoints    *_InputPoints;
  vtkPoints    *_OutputPoints;
  const double *_Displacement;

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    double x[3];
    const double *dx = _Displacement + 3 * re.begin();
    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId, dx += 3) {
      _InputPoints->GetPoint(ptId, x);
      x[0] += dx[0], x[1] += dx[1], x[2] += dx[2];
      _OutputPoints->SetPoint(ptId, x);
    }
  }

  static void Run(vtkPoints *input, const double *dx, vtkPoints *output = NULL)
  {
    MovePoints move;
    move._InputPoints  = input;
    move._OutputPoints = output ? output : input;
    move._Displacement = dx;
    move._OutputPoints->SetNumberOfPoints(input->GetNumberOfPoints());
    parallel_for(blocked_range<vtkIdType>(0, input->GetNumberOfPoints()), move);
  }
};

// -----------------------------------------------------------------------------
/// Perform one iteration of gradient smoothing
struct SmoothGradient
{
  const irtkEdgeTable *_EdgeTable;
  double              *_Input;
  double              *_Output;

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    int n;
    const int *adjPts;
    double *in, *out = _Output + 3 * re.begin();

    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId, out += 3) {
      in = _Input + 3 * ptId;
      out[0] = in[0], out[1] = in[1], out[2] = in[2];
      _EdgeTable->GetAdjacentPoints(ptId, n, adjPts);
      for (int i = 0; i < n; ++i) {
        in = _Input + 3 * adjPts[i];
        out[0] += in[0], out[1] += in[1], out[2] += in[2];
      }
      n += 1;
      out[0] /= n, out[1] /= n, out[2] /= n;
    }
  }
};


} // namespace irtkDeformableSurfaceModelUtils
using namespace irtkDeformableSurfaceModelUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkDeformableSurfaceModel::irtkDeformableSurfaceModel()
:
  _Image(NULL),
  _Transformation(NULL),
  _NumberOfTerms(0),
  _GradientSmoothing(0),
  _MinEdgeLength(-1.0),
  _MaxEdgeLength(-1.0),
  _MinFeatureAngle(20.0),
  _MaxFeatureAngle(60.0),
  _RemeshInterval(0),
  _RemeshCounter(0)
{
}

// -----------------------------------------------------------------------------
irtkDeformableSurfaceModel::~irtkDeformableSurfaceModel()
{
  Clear();
}

// =============================================================================
// Energy terms
// =============================================================================

// -----------------------------------------------------------------------------
void irtkDeformableSurfaceModel::Initialize()
{
  // Check input
  if (_Input == NULL || _Input->GetNumberOfPoints() == 0) {
    cerr << "irtkDeformableSurfaceModel::Initialize: Missing initial surface mesh" << endl;
    exit(1);
  }
  if (_NumberOfTerms == 0) {
    cerr << "irtkDeformableSurfaceModel::Initialize: No internal and/or external forces added" << endl;
    exit(1);
  }
  if (!_Transformation) {
    for (size_t i = 0; i < _Constraint.size(); ++i) _Constraint[i]->Weight(.0);
  }

  // Local adaptive remeshing settings
  if (IsTriangularMesh(_Input) && _RemeshInterval > 0) {
    if (_MinEdgeLength < .0 || _MaxEdgeLength <= .0) {
      double mean = AverageEdgeLength(_Input);
      if (_MinEdgeLength <  .0) _MinEdgeLength =  .5 * mean;
      if (_MaxEdgeLength <= .0) _MaxEdgeLength = 2.0 * mean;
    }
  } else {
    _RemeshInterval = 0;
  }
  _RemeshCounter = 0;

  // Initialize output surface mesh
  const bool deep_copy_points = (_Transformation == NULL);
  const bool init_edge_tables = (_GradientSmoothing > 0);
  _Surface.InputPointSet(_Input);
  _Surface.Transformation(_Transformation);
  _Surface.SelfUpdate(false);
  _Surface.Initialize(deep_copy_points, init_edge_tables);
  this->Changed(true);

  // Initialize energy terms
  for (size_t i = 0; i < _ExternalForce.size(); ++i) {
    _ExternalForce[i]->PointSet(&_Surface);
    _ExternalForce[i]->Image(_Image);
  }
  for (size_t i = 0; i < _InternalForce.size(); ++i) {
    _InternalForce[i]->PointSet(&_Surface);
  }
  for (int i = 0; i < _NumberOfTerms; ++i) {
    irtkEnergyTerm *term = Term(i);
    if (term->Weight() != .0) {
      term->Transformation(_Transformation);
      term->Initialize();
    }
  }

  // Initialize cache of computed energy term values
  _Value.resize(_NumberOfTerms);
  for (int i = 0; i < _NumberOfTerms; ++i) {
    _Value[i] = numeric_limits<double>::quiet_NaN();
  }
}

// -----------------------------------------------------------------------------
void irtkDeformableSurfaceModel::Clear()
{
  for (size_t i = 0; i < _ExternalForce.size(); ++i) delete _ExternalForce[i];
  for (size_t i = 0; i < _InternalForce.size(); ++i) delete _InternalForce[i];
  for (size_t i = 0; i < _Constraint   .size(); ++i) delete _Constraint   [i];
  _ExternalForce.clear();
  _InternalForce.clear();
  _Constraint   .clear();
  _NumberOfTerms = 0;
}

// -----------------------------------------------------------------------------
bool irtkDeformableSurfaceModel::Empty() const
{
  return _NumberOfTerms == 0;
}

// -----------------------------------------------------------------------------
int irtkDeformableSurfaceModel::NumberOfForces() const
{
  return NumberOfInternalForces() + NumberOfExternalForces();
}

// -----------------------------------------------------------------------------
int irtkDeformableSurfaceModel::NumberOfInternalForces() const
{
  return static_cast<int>(_InternalForce.size());
}

// -----------------------------------------------------------------------------
int irtkDeformableSurfaceModel::NumberOfExternalForces() const
{
  return static_cast<int>(_ExternalForce.size());
}

// -----------------------------------------------------------------------------
void irtkDeformableSurfaceModel::Add(irtkExternalForce *term)
{
  _ExternalForce.push_back(term);
  ++_NumberOfTerms;
}

// -----------------------------------------------------------------------------
void irtkDeformableSurfaceModel::Sub(irtkExternalForce *term)
{
  vector<irtkExternalForce *>::iterator it = _ExternalForce.begin();
  while (it != _ExternalForce.end()) {
    if (*it == term) {
      _ExternalForce.erase(it);
      --_NumberOfTerms;
      break;
    }
    ++it;
  }
}

// -----------------------------------------------------------------------------
void irtkDeformableSurfaceModel::Add(irtkPointSetConstraint *term)
{
  _InternalForce.push_back(term);
  ++_NumberOfTerms;
}

// -----------------------------------------------------------------------------
void irtkDeformableSurfaceModel::Sub(irtkPointSetConstraint *term)
{
  vector<irtkPointSetConstraint *>::iterator it = _InternalForce.begin();
  while (it != _InternalForce.end()) {
    if (*it == term) {
      _InternalForce.erase(it);
      --_NumberOfTerms;
      break;
    }
    ++it;
  }
}

// -----------------------------------------------------------------------------
void irtkDeformableSurfaceModel::Add(irtkTransformationConstraint *term)
{
  _Constraint.push_back(term);
  ++_NumberOfTerms;
}

// -----------------------------------------------------------------------------
void irtkDeformableSurfaceModel::Sub(irtkTransformationConstraint *term)
{
  vector<irtkTransformationConstraint *>::iterator it = _Constraint.begin();
  while (it != _Constraint.end()) {
    if (*it == term) {
      _Constraint.erase(it);
      --_NumberOfTerms;
      break;
    }
    ++it;
  }
}

// -----------------------------------------------------------------------------
irtkEnergyTerm *irtkDeformableSurfaceModel::Term(int i)
{
  if (i < NumberOfExternalForces()) return _ExternalForce[i];
  i -= NumberOfExternalForces();
  if (i < NumberOfInternalForces()) return _InternalForce[i];
  i -= NumberOfInternalForces();
  return _Constraint[i];
}

// -----------------------------------------------------------------------------
const irtkEnergyTerm *irtkDeformableSurfaceModel::Term(int i) const
{
  if (i < NumberOfExternalForces()) return _ExternalForce[i];
  i -= NumberOfExternalForces();
  if (i < NumberOfInternalForces()) return _InternalForce[i];
  i -= NumberOfInternalForces();
  return _Constraint[i];
}

// -----------------------------------------------------------------------------
irtkExternalForce *irtkDeformableSurfaceModel::ExternalForce(int i)
{
  return _ExternalForce[i];
}

// -----------------------------------------------------------------------------
const irtkExternalForce *irtkDeformableSurfaceModel::ExternalForce(int i) const
{
  return _ExternalForce[i];
}

// -----------------------------------------------------------------------------
irtkPointSetConstraint *irtkDeformableSurfaceModel::InternalForce(int i)
{
  return _InternalForce[i];
}

// -----------------------------------------------------------------------------
const irtkPointSetConstraint *irtkDeformableSurfaceModel::InternalForce(int i) const
{
  return _InternalForce[i];
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkDeformableSurfaceModel::Set(const char *name, const char *value)
{
  if (strcmp(name, "No. of gradient smoothing iterations") == 0) {
    return FromString(value, _GradientSmoothing);
  }
  if (strcmp(name, "Minimum edge length") == 0) {
    return FromString(value, _MinEdgeLength);
  }
  if (strcmp(name, "Maximum edge length") == 0) {
    return FromString(value, _MaxEdgeLength);
  }
  if (strcmp(name, "Remesh interval") == 0) {
    return FromString(value, _RemeshInterval);
  }

  bool known = false;
  for (int i = 0; i < _NumberOfTerms; ++i) {
    known = Term(i)->Set(name, value) || known;
  }
  return known;
}

// -----------------------------------------------------------------------------
irtkParameterList irtkDeformableSurfaceModel::Parameter() const
{
  irtkParameterList params;
  for (int i = 0; i < _NumberOfTerms; ++i) {
    Insert(params, Term(i)->Parameter());
  }
  Insert(params, "No. of gradient smoothing iterations", _GradientSmoothing);
  Insert(params, "Minimum edge length", _MinEdgeLength);
  Insert(params, "Maximum edge length", _MaxEdgeLength);
  Insert(params, "Remesh interval", _RemeshInterval);
  return params;
}

// =============================================================================
// Degrees of freedom
// =============================================================================

// -----------------------------------------------------------------------------
int irtkDeformableSurfaceModel::NumberOfDOFs() const
{
  return _Transformation ? _Transformation->NumberOfDOFs() : 3 * _Surface.NumberOfPoints();
}

// -----------------------------------------------------------------------------
void irtkDeformableSurfaceModel::Put(const double *x)
{
  if (_Transformation) {
    _Transformation->Put(x);
  } else {
    // Set vertex positions of deformable surface mesh
    vtkPoints *points = _Surface.Points();
    for (int i = 0; i < _Surface.NumberOfPoints(); ++i, x += 3) {
      points->SetPoint(i, x);
    }
    _Surface.PointsChanged();
  }
  // Reset cached energy values
  for (int i = 0; i < _NumberOfTerms; ++i) {
    _Value[i] = numeric_limits<double>::quiet_NaN();
  }
  // Mark deformable surface model as changed
  this->Changed(true);
}

// -----------------------------------------------------------------------------
void irtkDeformableSurfaceModel::Get(double *x) const
{
  if (_Transformation) {
    _Transformation->Get(x);
  } else {
    for (int i = 0; i < _Surface.NumberOfPoints(); ++i, x += 3) {
      _Surface.GetPoint(i, x);
    }
  }
}

// -----------------------------------------------------------------------------
double irtkDeformableSurfaceModel::Get(int dof) const
{
  if (_Transformation) {
    return _Transformation->Get(dof);
  } else {
    double p[3];
    _Surface.GetPoint(dof / 3, p);
    return p[dof % 3];
  }
}

// -----------------------------------------------------------------------------
double irtkDeformableSurfaceModel::Step(const double *dx)
{
  double delta;
  if (_Transformation) {
    delta = _Transformation->Update(dx);
  } else {
    // Determine maximum vertex displacement (before change of dx pointer)
    MaxVertexDisplacement max(dx);
    parallel_reduce(blocked_range<int>(0, _Surface.NumberOfPoints()), max);
    delta = max.Norm();
    // Update points of output surface
    MovePoints::Run(_Surface.Points(), dx);
    _Surface.PointsChanged();
  }
  // Reset cached energy values
  for (int i = 0; i < _NumberOfTerms; ++i) {
    _Value[i] = numeric_limits<double>::quiet_NaN();
  }
  // Mark deformable surface model as changed
  if (delta != .0) this->Changed(true);
  // Return maximum vertex displacement
  return delta;
}

// -----------------------------------------------------------------------------
void irtkDeformableSurfaceModel::Update(bool gradient)
{
  if (this->Changed() || gradient) {
    IRTK_START_TIMING();
    // Update edge table needed for smoothing operations
    if (_Transformation) _Surface.Update(true);
    // Update energy terms
    for (int i = 0; i < _NumberOfTerms; ++i) {
      irtkEnergyTerm *term = Term(i);
      if (term->Weight() != .0) term->Update(gradient);
    }
    // Mark deformable surface model as up-to-date
    this->Changed(false);
    IRTK_DEBUG_TIMING(3, "update of energy function");
  }
}

// -----------------------------------------------------------------------------
bool irtkDeformableSurfaceModel::Upgrade()
{
  return false;
}

// -----------------------------------------------------------------------------
bool irtkDeformableSurfaceModel::Remesh()
{
  // Currently only remeshing of a triangulated surface mesh is supported
  if (_RemeshInterval == 0) return false;

  // Increment counter counting the number of iterations since the last remeshing
  ++_RemeshCounter;
  if (_RemeshCounter < _RemeshInterval) return false;
  _RemeshCounter = 0;

  IRTK_START_TIMING();

  // Perform local adaptive remeshing
  irtkPolyDataRemeshing remesher;
  remesher.MeltingOrder(irtkPolyDataRemeshing::SHORTEST_EDGE);
  remesher.MeltNodesOff();
  remesher.MeltTrianglesOn();
  remesher.MinEdgeLength(_MinEdgeLength);
  remesher.MaxEdgeLength(_MaxEdgeLength);
  remesher.MinFeatureAngle(_MinFeatureAngle);
  remesher.MaxFeatureAngle(_MaxFeatureAngle);
  remesher.SkipTriangulationOn();

  if (_Transformation) {
    remesher.Input(_Surface.InputSurface());
    remesher.Transformation(_Transformation);
  } else {
    vtkSmartPointer<vtkPolyData> surface = _Surface.Surface();
    if (!surface->GetPointData()->HasArray("InitialPoints")) {
      vtkSmartPointer<vtkDataArray> initial_points;
      initial_points = vtkSmartPointer<vtkFloatArray>::New();
      initial_points->SetName("InitialPoints");
      initial_points->SetNumberOfComponents(3);
      initial_points->SetNumberOfTuples(surface->GetNumberOfPoints());
      vtkPoints *points = _Surface.InputSurface()->GetPoints();
      for (vtkIdType ptId = 0; ptId < surface->GetNumberOfPoints(); ++ptId) {
        initial_points->SetTuple(ptId, points->GetPoint(ptId));
      }
      surface->GetPointData()->AddArray(initial_points);
    }
    remesher.Input(_Surface.Surface());
  }

  remesher.Run();
  if (remesher.Output() != remesher.Input()) {
    vtkSmartPointer<vtkPolyData> surface = remesher.Output();

    // Update deformable surface mesh
    const bool init_edge_tables = (_GradientSmoothing > 0);

    _Surface.InputPointSet(surface);
    if (_Transformation) {
      _Surface.Initialize(false, init_edge_tables);
      _Surface.Update(true);
    } else {
      // Initialize surface, making deep copy of new deformed surface points
      _Surface.Initialize(true, init_edge_tables);
      // Reset points of input surface mesh such that
      // irtkPointSetForce::GetInitialPoints returns the correct positions
      vtkPoints    *points         = surface->GetPoints();
      vtkDataArray *initial_points = surface->GetPointData()->GetArray("InitialPoints");
      for (vtkIdType ptId = 0; ptId < points->GetNumberOfPoints(); ++ptId) {
        points->SetPoint(ptId, initial_points->GetTuple(ptId));
      }
    }

    // Reinitialize internal and external force terms
    for (size_t i = 0; i < _ExternalForce.size(); ++i) {
      if (_ExternalForce[i]->Weight() != .0) {
        _ExternalForce[i]->Reinitialize();
      }
    }
    for (size_t i = 0; i < _InternalForce.size(); ++i) {
      if (_InternalForce[i]->Weight() != .0) {
        _InternalForce[i]->Reinitialize();
      }
    }

    // Mark deformable surface model as up-to-date
    this->Changed(false);
  }

  IRTK_DEBUG_TIMING(3, "local adaptive remeshing");
  return remesher.Output() != remesher.Input();
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
double irtkDeformableSurfaceModel::RawValue(int i)
{
  irtkEnergyTerm *term = Term(i);
  if (term->Weight() == .0) return .0;
  if (IsNaN(_Value[i])) _Value[i] = term->Value();
  return term->RawValue(_Value[i]);
}

// -----------------------------------------------------------------------------
double irtkDeformableSurfaceModel::InitialValue()
{
  IRTK_START_TIMING();

  double sum = .0;
  for (int i = 0; i < _NumberOfTerms; ++i) {
    irtkEnergyTerm *term = Term(i);
    if (term->Weight() != .0) {
      sum += term->InitialValue();
    } else {
      _Value[i] = .0;
    }
  }

  IRTK_DEBUG_TIMING(3, "initial evaluation of energy function");
  return sum;
}

// -----------------------------------------------------------------------------
double irtkDeformableSurfaceModel::InitialValue(int i)
{
  return Term(i)->InitialValue();
}

// -----------------------------------------------------------------------------
double irtkDeformableSurfaceModel::Value()
{
  IRTK_START_TIMING();

  double sum = .0;
  for (int i = 0; i < _NumberOfTerms; ++i) {
    irtkEnergyTerm *term = Term(i);
    if (term->Weight() != .0) {
      if (IsNaN(_Value[i])) _Value[i] = term->Value();
      sum += _Value[i];
    } else {
      _Value[i] = .0;
    }
  }

  IRTK_DEBUG_TIMING(3, "evaluation of energy function");
  return sum;
}

// -----------------------------------------------------------------------------
double irtkDeformableSurfaceModel::Value(int i)
{
  return _Value[i];
}

// -----------------------------------------------------------------------------
void irtkDeformableSurfaceModel::Gradient(double *gradient, double step, bool *sgn_chg)
{
  IRTK_START_TIMING();

  // Use default step length if none specified
  if (step <= .0) step = _StepLength;

  // Initialize output variables
  const int ndofs = this->NumberOfDOFs();
  memset(gradient, 0, ndofs * sizeof(double));
  if (sgn_chg) {
    for (int dof = 0; dof < ndofs; ++dof) {
      sgn_chg[dof] = true;
    }
  }

  // Sum (weighted) internal and external forces
  for (int i = 0; i < _NumberOfTerms; ++i) {
    irtkEnergyTerm *term = Term(i);
    if (term->Weight() != .0) {
      term->Gradient(gradient, step);
    }
  }

  IRTK_DEBUG_TIMING(3, "evaluation of energy gradient");
}

// -----------------------------------------------------------------------------
double irtkDeformableSurfaceModel::GradientNorm(const double *dx) const
{
  if (_Transformation) return _Transformation->DOFGradientNorm(dx);
  MaxVertexDisplacement max(dx);
  parallel_reduce(blocked_range<int>(0, _Surface.NumberOfPoints()), max);
  return max.Norm();
}

// -----------------------------------------------------------------------------
void irtkDeformableSurfaceModel::GradientStep(const double *dx, double &min, double &max) const
{
  for (int i = 0; i < _NumberOfTerms; ++i) {
    const irtkEnergyTerm *term = Term(i);
    if (term->Weight() != .0) {
      term->GradientStep(dx, min, max);
    }
  }
}

// -----------------------------------------------------------------------------
double irtkDeformableSurfaceModel::Evaluate(double *dx, double step, bool *sgn_chg)
{
  // Update energy function
  if (this->Changed()) this->Update(dx != NULL);

  // Evaluate gradient
  if (dx) this->Gradient(dx, step, sgn_chg);

  // Evaluate energy
  return this->Value();
}

// -----------------------------------------------------------------------------
void irtkDeformableSurfaceModel::EnforceHardConstraints(double *dx) const
{
  // Hard constraints only apply to non-parametric deformable surface models
  if (_Transformation) return;

  // Non-self-intersection constraint
  if (IsSurfaceMesh(_Surface.PointSet())) {
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    MovePoints::Run(_Surface.Points(), dx, points);

    vtkSmartPointer<vtkPolyData> surface = vtkSmartPointer<vtkPolyData>::New();
    surface->ShallowCopy(_Surface.Surface());
    surface->SetPoints(points);

    irtkSurfaceCollisions nsi;
    nsi.Input(surface);
    nsi.AdjacentIntersectionTest(true);
    nsi.NonAdjacentIntersectionTest(true);
    nsi.FrontfaceCollisionTest(false);
    nsi.BackfaceCollisionTest(false);
    nsi.MinDistance(.0);
    nsi.MaxAngle(45.0);
    nsi.Run();

    vtkDataArray *collisions = nsi.GetCollisionTypeArray();
    _Surface.Surface()->GetCellData()->RemoveArray(collisions->GetName());
    _Surface.Surface()->GetCellData()->AddArray(collisions);

    vtkIdType npts, *pts;
    for (vtkIdType cellId = 0; cellId < surface->GetNumberOfCells(); ++cellId) {
      if (collisions->GetComponent(cellId, 0) != .0) {
        surface->GetCellPoints(cellId, npts, pts);
        for (vtkIdType i = 0; i < npts; ++i) {
          double * const d = dx + 3 * pts[i];
          d[0] = d[1] = d[2] = .0;
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
void irtkDeformableSurfaceModel::SmoothGradient(double *dx) const
{
  // Only applies when DoFs are the node positions themselves
  if (_Transformation || _GradientSmoothing <= 0) return;

  // Array of detected self-collisions
  vtkPolyData  *surface    = _Surface.Surface();
  vtkDataArray *collisions = surface->GetCellData()->GetArray("CollisionType");

  // Smooth vertex displacements such that adjacent nodes move coherently
  IRTK_START_TIMING();
  irtkDeformableSurfaceModelUtils::SmoothGradient smooth;
  smooth._Input     = dx;
  smooth._Output    = Allocate<double>(this->NumberOfDOFs());
  smooth._EdgeTable = _Surface.SurfaceEdges();
  blocked_range<vtkIdType> ptIds(0, _Surface.NumberOfPoints());
  for (int iter = 0; iter < _GradientSmoothing; ++iter) {
    parallel_for(ptIds, smooth);
    if (collisions) {
      vtkIdType npts, *pts;
      for (vtkIdType cellId = 0; cellId < surface->GetNumberOfCells(); ++cellId) {
        if (collisions->GetComponent(cellId, 0) != .0) {
          surface->GetCellPoints(cellId, npts, pts);
          for (vtkIdType i = 0; i < npts; ++i) {
            double * const d = smooth._Output + 3 * pts[i];
            d[0] = d[1] = d[2] = .0;
          }
        }
      }
    }
    swap(smooth._Input, smooth._Output);
  }
  if (smooth._Output == dx) {
    delete[] smooth._Input;
  } else {
    memcpy(dx, smooth._Output, this->NumberOfDOFs() * sizeof(double));
    delete[] smooth._Output;
  }
  IRTK_DEBUG_TIMING(3, "smoothing of energy gradient (#iter=" << _GradientSmoothing << ")");
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
void irtkDeformableSurfaceModel::WriteDataSets(const char *prefix, const char *suffix, bool all) const
{
  const int sz = 1024;
  char      fname[sz];
  snprintf(fname, sz, "%soutput%s%s", prefix, suffix, _Surface.DefaultExtension());
  _Surface.Write(fname);

  if (all) {
    for (int i = 0; i < _NumberOfTerms; ++i) {
      const irtkEnergyTerm *term = Term(i);
      if (term->Weight() != .0) {
        term->WriteDataSets(prefix, suffix, all);
      }
    }
  }
}

// -----------------------------------------------------------------------------
void irtkDeformableSurfaceModel::WriteGradient(const char *prefix, const char *suffix) const
{
  for (int i = 0; i < _NumberOfTerms; ++i) {
    const irtkEnergyTerm *term = Term(i);
    if (term->Weight() != .0) {
      term->WriteGradient(prefix, suffix);
    }
  }
}
