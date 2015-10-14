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

#include <irtkLimitedMemoryBFGSDescent.h>

#include <lbfgs.h>


// =============================================================================
// liblbfgs error messages
// =============================================================================

// -----------------------------------------------------------------------------
const char *lbfgs_error_message(int err)
{
  switch (err) {
    case LBFGS_ALREADY_MINIMIZED:            return "The initial variables already minimize the objective function.";
    case LBFGSERR_LOGICERROR:                return "Logic error.";
    case LBFGSERR_OUTOFMEMORY:               return "Insufficient memory.";
    case LBFGSERR_CANCELED:                  return "The minimization process has been canceled.";
    case LBFGSERR_INVALID_N:                 return "Invalid number of variables specified.";
    case LBFGSERR_INVALID_N_SSE:             return "Invalid number of variables (for SSE) specified.";
    case LBFGSERR_INVALID_X_SSE:             return "The array x must be aligned to 16 (for SSE).";
    case LBFGSERR_INVALID_EPSILON:           return "Invalid parameter lbfgs_parameter_t::epsilon specified.";
    case LBFGSERR_INVALID_TESTPERIOD:        return "Invalid parameter lbfgs_parameter_t::past specified.";
    case LBFGSERR_INVALID_DELTA:             return "Invalid parameter lbfgs_parameter_t::delta specified.";
    case LBFGSERR_INVALID_LINESEARCH:        return "Invalid parameter lbfgs_parameter_t::linesearch specified.";
    case LBFGSERR_INVALID_MINSTEP:           return "Invalid parameter lbfgs_parameter_t::max_step specified.";
    case LBFGSERR_INVALID_MAXSTEP:           return "Invalid parameter lbfgs_parameter_t::max_step specified.";
    case LBFGSERR_INVALID_FTOL:              return "Invalid parameter lbfgs_parameter_t::ftol specified.";
    case LBFGSERR_INVALID_WOLFE:             return "Invalid parameter lbfgs_parameter_t::wolfe specified.";
    case LBFGSERR_INVALID_GTOL:              return "Invalid parameter lbfgs_parameter_t::gtol specified.";
    case LBFGSERR_INVALID_XTOL:              return "Invalid parameter lbfgs_parameter_t::xtol specified.";
    case LBFGSERR_INVALID_MAXLINESEARCH:     return "Invalid parameter lbfgs_parameter_t::max_linesearch specified.";
    case LBFGSERR_INVALID_ORTHANTWISE:       return "Invalid parameter lbfgs_parameter_t::orthantwise_c specified.";
    case LBFGSERR_INVALID_ORTHANTWISE_START: return "Invalid parameter lbfgs_parameter_t::orthantwise_start specified.";
    case LBFGSERR_INVALID_ORTHANTWISE_END:   return "Invalid parameter lbfgs_parameter_t::orthantwise_end specified.";
    case LBFGSERR_OUTOFINTERVAL:             return "The line-search step went out of the interval of uncertainty.";
    case LBFGSERR_INCORRECT_TMINMAX:         return "A logic error occurred; alternatively, the interval of uncertainty became too small.";
    case LBFGSERR_ROUNDING_ERROR:            return "A rounding error occurred; alternatively, no line-search step satisfies the sufficient decrease and curvature conditions.";
    case LBFGSERR_MINIMUMSTEP:               return "The line-search step became smaller than lbfgs_parameter_t::min_step.";
    case LBFGSERR_MAXIMUMSTEP:               return "The line-search step became larger than lbfgs_parameter_t::max_step.";
    case LBFGSERR_MAXIMUMLINESEARCH:         return "The line-search routine reaches the maximum number of evaluations.";
    case LBFGSERR_MAXIMUMITERATION:          return "The algorithm routine reaches the maximum number of iterations.";
    case LBFGSERR_WIDTHTOOSMALL:             return "Relative width of the interval of uncertainty is at most lbfgs_parameter_t::xtol.";
    case LBFGSERR_INVALIDPARAMETERS:         return "A logic error (negative line-search step) occurred.";
    case LBFGSERR_INCREASEGRADIENT:          return "The current search direction increases the objective function value.";
    default:                                 return "Unknown error.";
  }
}

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkLimitedMemoryBFGSDescent::irtkLimitedMemoryBFGSDescent(irtkObjectiveFunction *f)
:
  irtkLocalOptimizer(f),
  _NumberOfIterations(100),
  _MinStepLength     (0.01),
  _MaxStepLength     (1.0)
{
}

// -----------------------------------------------------------------------------
irtkLimitedMemoryBFGSDescent::irtkLimitedMemoryBFGSDescent(const irtkLimitedMemoryBFGSDescent &other)
:
  irtkLocalOptimizer(other),
  _NumberOfIterations(other._NumberOfIterations),
  _MinStepLength     (other._MinStepLength),
  _MaxStepLength     (other._MaxStepLength)
{
}

// -----------------------------------------------------------------------------
irtkLimitedMemoryBFGSDescent &irtkLimitedMemoryBFGSDescent::operator =(const irtkLimitedMemoryBFGSDescent &other)
{
  irtkLocalOptimizer::operator =(other);
  _NumberOfIterations = other._NumberOfIterations;
  _MinStepLength      = other._MinStepLength;
  _MaxStepLength      = other._MaxStepLength;
  return *this;
}

// -----------------------------------------------------------------------------
irtkLimitedMemoryBFGSDescent::~irtkLimitedMemoryBFGSDescent()
{
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool irtkLimitedMemoryBFGSDescent::Set(const char *name, const char *value)
{
  if (strcmp(name, "No. of iterations") == 0) {
    return FromString(value, _NumberOfIterations) && _NumberOfIterations > 0;
  } else if (strcmp(name, "Minimum length of steps") == 0) {
    return FromString(value, _MinStepLength) && _MinStepLength > .0;
  } else if (strcmp(name, "Maximum length of steps") == 0) {
    return FromString(value, _MaxStepLength) && _MaxStepLength > .0;
  } else if (strcmp(name, "Length of steps") == 0) {
    if (!FromString(value, _MaxStepLength) || _MaxStepLength <= .0) return false;
    _MinStepLength = _MaxStepLength;
  }
  return irtkLocalOptimizer::Set(name, value);
}

// -----------------------------------------------------------------------------
irtkParameterList irtkLimitedMemoryBFGSDescent::Parameter() const
{
  irtkParameterList params = irtkLocalOptimizer::Parameter();
  Insert(params, "No. of iterations", ToString(_NumberOfIterations));
  if (_MinStepLength == _MaxStepLength) {
    Insert(params, "Length of steps", ToString(_MinStepLength));
  } else {
    Insert(params, "Minimum length of steps", ToString(_MinStepLength));
    Insert(params, "Maximum length of steps", ToString(_MaxStepLength));
  }
  return params;
}

// =============================================================================
// Optimization
// =============================================================================

// L-BFGS helper
struct LBFGSCallback
{
  static lbfgsfloatval_t Evaluate(
      void                  *o,    // Instance of irtkLimitedMemoryBFGSDescent
      const lbfgsfloatval_t *x,    // Current parameters
      lbfgsfloatval_t       *g,    // Gradient of objective function
      const int              n,    // Number of parameters
      const lbfgsfloatval_t  step) // Line search step length
  {
    irtkLimitedMemoryBFGSDescent *_this = reinterpret_cast<irtkLimitedMemoryBFGSDescent *>(o);
    _this->Function()->Put(x);
    lbfgsfloatval_t m = _this->Function()->Evaluate(g, step);
    _this->_CurrentStep._Value = m;
    return m;
  }

  static int Progress(
      void                  *o,        // Instance of irtkLimitedMemoryBFGSDescent
      const lbfgsfloatval_t *x,        // Current parameters
      const lbfgsfloatval_t */*g*/,    // Current gradient
      const lbfgsfloatval_t fx,        // Current value
      const lbfgsfloatval_t xnorm,     // Norm of parameter vector
      const lbfgsfloatval_t /*gnorm*/, // Norm of gradient vector
      const lbfgsfloatval_t step,      // Current step length
      int n,                           // Number of parameters
      int k,                           // Current iteration
      int ls                           // Current line search iteration
      )
  {
    irtkLimitedMemoryBFGSDescent *_this = reinterpret_cast<irtkLimitedMemoryBFGSDescent *>(o);
    string msg = "Current best metric value is ";
    msg += ToString(fx);
    msg += ", step ";
    msg += ToString(step);
    msg += ", after ";
    msg += ToString(ls);
    msg += " line search steps\n";
    _this->Broadcast(LogEvent, msg.c_str());
    /*
    irtkIteration iter  (0, _this->NumberOfIterations());
    irtkIteration lsiter(0, 1);
    iter  ._Iter = k;
    lsiter._Iter = ls;
    if (ls == 0) {
      _this->_CurrentStep._Current   = fx;
      _this->_CurrentStep._MinLength = _this->MinStepLength();
      _this->_CurrentStep._MaxLength = _this->MaxStepLength();
      if (k > 0) _this->Broadcast(LineSearchEndEvent, &_this->_CurrentStep);
      _this->Broadcast(IterationEvent, &iter);
      if (ls == 0) _this->Broadcast(LineSearchStartEvent, &_this->_CurrentStep);
    } else {
      _this->Broadcast(LineSearchIterationStartEvent, &lsiter);
      _this->_CurrentStep._Length = step;
      _this->_CurrentStep._Delta  = xnorm;
      if (fx < _this->_CurrentStep._Current - _this->Epsilon()) {
        _this->Broadcast(AcceptedStepEvent, &_this->_CurrentStep);
        _this->_CurrentStep._Current = fx;
      } else {
        _this->Broadcast(RejectedStepEvent, &_this->_CurrentStep);
      }
      _this->Broadcast(LineSearchIterationEndEvent, &lsiter);
    }
    */
    return 0;
  }
};

// -----------------------------------------------------------------------------
double irtkLimitedMemoryBFGSDescent::Run()
{
  lbfgs_parameter_t param;
  lbfgs_parameter_init(&param);

  param.past           = 1;
  param.delta          = _Epsilon;
  param.epsilon        = _Delta;
  param.max_iterations = _NumberOfIterations;
  param.min_step       = _MinStepLength;
  param.max_step       = _MaxStepLength;

  const int        n = Function()->NumberOfDOFs();
  lbfgsfloatval_t *x = lbfgs_malloc(n);
  Function()->Get(x);

  int ret = lbfgs(n, x, NULL, LBFGSCallback::Evaluate, LBFGSCallback::Progress, this, &param);
  if (ret < 0 && ret != LBFGSERR_MAXIMUMITERATION && ret != LBFGSERR_MAXIMUMLINESEARCH) {
    cerr << "L-BFGS optimization failed: " << lbfgs_error_message(ret) << " (error code: " << ret << ")" << endl;
    exit(1);
  }

  lbfgs_free(x);
  return _CurrentStep._Value;
}
