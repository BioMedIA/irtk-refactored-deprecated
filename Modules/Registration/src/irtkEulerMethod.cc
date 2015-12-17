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

#include <irtkEulerMethod.h>

#include <irtkDeformableSurfaceModel.h>

#include <vtkPointData.h>


// =============================================================================
// Auxiliary functors
// =============================================================================

namespace irtkEulerMethodUtils {


// -----------------------------------------------------------------------------
/// Negate deformable surface model gradient to obtain sum of forces
class NegateValues
{
  const double *_Input;
  double       *_Output;

public:

  NegateValues(double *b, const double *a = NULL)
  :
    _Input(a ? a : b), _Output(b)
  {}

  void operator ()(const blocked_range<int> &re) const
  {
    for (int i = re.begin(); i != re.end(); ++i) {
      _Output[i] = -_Input[i];
    }
  }
};

// -----------------------------------------------------------------------------
/// Compute/update node velocities
class ComputeVelocities
{
  irtkEulerMethod *_Optimizer;
  const double    *_Force;
  vtkDataArray    *_Velocity;

public:

  ComputeVelocities(irtkEulerMethod *obj, const double *f, vtkDataArray *v)
  :
    _Optimizer(obj), _Force(f), _Velocity(v)
  {}

  void operator ()(const blocked_range<int> &re) const
  {
    double v[3];
    const double *f = _Force + 3 * re.begin();
    for (int ptId = re.begin(); ptId != re.end(); ++ptId, f += 3) {
      _Velocity->GetTuple(ptId, v);
      _Optimizer->ComputeVelocity(v, f);
      _Velocity->SetTuple(ptId, v);
    }
  }
};

// -----------------------------------------------------------------------------
/// Compute displacements from velocities and temporal integration step length
class ComputeDisplacements
{
  vtkDataArray *_Velocity;
  double       *_Displacement;
  double        _StepLength;

public:

  ComputeDisplacements(vtkDataArray *v, double dt, double *d)
  :
    _Velocity(v), _Displacement(d), _StepLength(dt)
  {}

  void operator ()(const blocked_range<int> &re) const
  {
    double v[3], *d = _Displacement + 3 * re.begin();
    for (int ptId = re.begin(); ptId != re.end(); ++ptId, d += 3) {
      _Velocity->GetTuple(ptId, v);
      d[0] = v[0] * _StepLength;
      d[1] = v[1] * _StepLength;
      d[2] = v[2] * _StepLength;
    }
  }
};


} // namespace irtkEulerMethodUtils
using namespace irtkEulerMethodUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
void irtkEulerMethod::Copy(const irtkEulerMethod &other)
{
  Deallocate(_Force);
  Deallocate(_Displacement);
  _NumberOfSteps = other._NumberOfSteps;
  _DampingFactor = other._DampingFactor;
  _StepLength    = other._StepLength;
  _NumberOfDOFs  = other._NumberOfDOFs;
  if (_NumberOfDOFs > 0) {
    _Force = Allocate<double>(_NumberOfDOFs);
    memcpy(_Force, other._Force, _NumberOfDOFs * sizeof(double));
    _Displacement = Allocate<double>(_NumberOfDOFs);
    memcpy(_Displacement, other._Displacement, _NumberOfDOFs * sizeof(double));
  }
}

// -----------------------------------------------------------------------------
irtkEulerMethod::irtkEulerMethod(irtkObjectiveFunction *f)
:
  irtkLocalOptimizer(f),
  _NumberOfSteps(100),
  _DampingFactor(1),
  _StepLength(1.0),
  _Force(NULL),
  _Displacement(NULL),
  _NumberOfDOFs(0)
{
}

// -----------------------------------------------------------------------------
irtkEulerMethod::irtkEulerMethod(const irtkEulerMethod &other)
:
  irtkLocalOptimizer(other),
  _Force(NULL),
  _Displacement(NULL)
{
  Copy(other);
}

// -----------------------------------------------------------------------------
irtkEulerMethod &irtkEulerMethod::operator =(const irtkEulerMethod &other)
{
  irtkLocalOptimizer::operator =(other);
  Copy(other);
  return *this;
}

// -----------------------------------------------------------------------------
irtkEulerMethod::~irtkEulerMethod()
{
  Deallocate(_Force);
  Deallocate(_Displacement);
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkEulerMethod::Set(const char *name, const char *value)
{
  if (strcmp(name, "Maximum no. of iterations")        == 0 ||
      strcmp(name, "Maximum number of iterations")     == 0 ||
      strcmp(name, "Maximum no. of gradient steps")    == 0 ||
      strcmp(name, "Maximum number of gradient steps") == 0 ||
      strcmp(name,    "No. of iterations")             == 0 ||
      strcmp(name, "Number of iterations")             == 0 ||
      strcmp(name,    "No. of gradient steps")         == 0 ||
      strcmp(name, "Number of gradient steps")         == 0) {
    return FromString(value, _NumberOfSteps);
  }
  if (strcmp(name, "Deformable surface damping") == 0) {
    return FromString(value, _DampingFactor);
  }
  if (strcmp(name, "Deformable surface step length") == 0 ||
      strcmp(name, "Length of steps")                == 0 ||
      strcmp(name, "Maximum length of steps")        == 0) {
    return FromString(value, _StepLength);
  }
  return irtkLocalOptimizer::Set(name, value);
}

// -----------------------------------------------------------------------------
irtkParameterList irtkEulerMethod::Parameter() const
{
  irtkParameterList params = irtkLocalOptimizer::Parameter();
  Insert(params, "Maximum no. of iterations",  _NumberOfSteps);
  Insert(params, "Deformable surface damping", _DampingFactor);
  Insert(params, "Length of steps",            _StepLength);
  return params;
}

// =============================================================================
// Optimization
// =============================================================================

// -----------------------------------------------------------------------------
void irtkEulerMethod::Initialize()
{
  // Initialize base class
  irtkLocalOptimizer::Initialize();

  // Cast objective function to deformable surface model
  _Model = dynamic_cast<irtkDeformableSurfaceModel *>(_Function);
  if (_Model == NULL) {
    cerr << "irtkEulerMethod::Initialize: Objective function must be a deformable surface model" << endl;
    exit(1);
  }
  if (_Model->Transformation()) {
    cerr << "irtkEulerMethod::Initialize: Optimizer can only be used for non-parametric deformable surface models" << endl;
    exit(1);
  }

  // Allocate memory for node vectors if not done before
  if (_Model->NumberOfDOFs() > _NumberOfDOFs) {
    Deallocate(_Force);
    Deallocate(_Displacement);
    _NumberOfDOFs = _Model->NumberOfDOFs();
    Allocate(_Force,        _NumberOfDOFs);
    Allocate(_Displacement, _NumberOfDOFs);
  }

  // Add point data array with initial node velocities such that these
  // are interpolated at new node positions during the remeshing
  vtkSmartPointer<vtkDataArray> velocity;
  velocity = _Model->Output()->GetPointData()->GetArray("Velocity");
  if (!velocity) {
    velocity = vtkSmartPointer<vtkFloatArray>::New();
    velocity->SetName("Velocity");
    velocity->SetNumberOfComponents(3);
    _Model->Output()->GetPointData()->AddArray(velocity);
  }
  velocity->SetNumberOfTuples(_Model->NumberOfPoints());
  velocity->FillComponent(0, .0);
  velocity->FillComponent(1, .0);
  velocity->FillComponent(2, .0);
}

// -----------------------------------------------------------------------------
double irtkEulerMethod::Run()
{
  // Initialize
  this->Initialize();

  // Initial update of deformable surface model before start event because
  // the update can trigger some lazy initialization which in turn may
  // broadcast some log events for verbose command output
  _Model->Update(true);

  // Notify observers about start of optimization
  Broadcast(StartEvent);

  // Get initial energy value
  double value = _Model->Value();
  if (IsNaN(value)) {
    cerr << "irtkEulerMethod::Run:" << __LINE__ << ": NaN objective function value!" << endl;
    exit(1);
  }

  // Perform explicit integration steps
  irtkIteration step(0, _NumberOfSteps);
  while (step.Next()) {

    // Notify observers about start of iteration
    Broadcast(IterationStartEvent, &step);

    // Update node velocities
    this->UpdateVelocity();

    // Check equilibrium condition
    //
    // TODO: Determine ratio of nodes with (close to) zero velocity and
    //       break if it falls below a user-defined threshold.
    vtkDataArray *velocity = _Model->Output()->GetPointData()->GetArray("Velocity");
    const double max_delta = velocity->GetMaxNorm();
    if (fequal(max_delta, .0)) break;

    // Move points of surface model
    ComputeDisplacements mul(velocity, _StepLength / max_delta, _Displacement);
    parallel_for(blocked_range<int>(0, _Model->NumberOfPoints()), mul);
    _Model->EnforceHardConstraints(_Displacement);
    _Model->SmoothGradient(_Displacement);
    if (_Model->GradientNorm(_Displacement) <= _Delta) break;
    _Model->Step(_Displacement);

    // Perform local adaptive remeshing
    if (_Model->Remesh()) {
      if (_Model->NumberOfDOFs() > _NumberOfDOFs) {
        Deallocate(_Force);
        Deallocate(_Displacement);
        _NumberOfDOFs = _Model->NumberOfDOFs();
        Allocate(_Force,        _NumberOfDOFs);
        Allocate(_Displacement, _NumberOfDOFs);
      }
    }

    // Update deformable surface model
    _Model->Update(true);

    // Check convergence towards local minimum of energy function
    //
    // For deformable surface models which are not defined as a local minimum
    // of a well-defined energy function but only via the equilibrium of
    // internal and external forces, the energy values corresponding to the
    // external forces are infinite and hence the total energy value
    if (!IsInf(value)) {
      const double prev = value;
      value = _Model->Value();
      if (IsNaN(value)) {
        cerr << "irtkEulerMethod::Run:" << __LINE__ << ": NaN objective function value!" << endl;
        exit(1);
      }
      if (fabs(value - prev) <= _Epsilon) break;
    }

    // Notify observers about end of iteration
    Broadcast(IterationEndEvent, &step);
  }

  // Notify observers about end of optimization
  Broadcast(EndEvent, &value);

  // Finalize
  this->Finalize();

  return value;
}

// -----------------------------------------------------------------------------
void irtkEulerMethod::UpdateVelocity()
{
  // Compute forces
  double * const gradient = _Force;
  _Model->Gradient(gradient);
  NegateValues negate(_Force, gradient);
  parallel_for(blocked_range<int>(0, _Model->NumberOfDOFs()), negate);

  // Update velocities
  vtkDataArray *velocity = _Model->Output()->GetPointData()->GetArray("Velocity");
  ComputeVelocities eval(this, _Force, velocity);
  parallel_for(blocked_range<int>(0, _Model->NumberOfPoints()), eval);
}

// -----------------------------------------------------------------------------
void irtkEulerMethod::ComputeVelocity(double *v, const double *f) const
{
  v[0] = f[0] / _DampingFactor;
  v[1] = f[1] / _DampingFactor;
  v[2] = f[2] / _DampingFactor;
}

// -----------------------------------------------------------------------------
void irtkEulerMethod::Finalize()
{
  Deallocate(_Force);
  Deallocate(_Displacement);
  _Model->Output()->GetPointData()->RemoveArray("Velocity");
}
