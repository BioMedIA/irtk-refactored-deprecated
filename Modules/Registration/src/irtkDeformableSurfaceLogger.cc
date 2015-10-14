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

#include <irtkDeformableSurfaceLogger.h>
#include <irtkDeformableSurfaceModel.h>
#include <irtkLocalOptimizer.h>


// -----------------------------------------------------------------------------
ostream &PrintNumber(ostream &os, double value)
{
  if (value != .0 && (fabs(value) < 1.e-5 || fabs(value) >= 1.e5)) {
    os << scientific << setprecision(5) << setw(16) << value; // e-0x
  } else os << fixed << setprecision(5) << setw(12) << value << "    ";
  return os;
}

// -----------------------------------------------------------------------------
ostream &PrintNormalizedNumber(ostream &os, double value)
{
  os << fixed << setprecision(5) << setw(8) << value;
  return os;
}

// -----------------------------------------------------------------------------
ostream &PrintWeight(ostream &os, double weight, int nterms)
{
  const int w = round(fabs(weight) * 100.0);
  if      (w ==   0)               os << "< 1";
  else if (w == 100 && nterms > 1) os << ">99";
  else                             os << fixed << setw(3) << w;
  return os;
}

// -----------------------------------------------------------------------------
irtkDeformableSurfaceLogger::irtkDeformableSurfaceLogger(ostream *stream)
:
  _Verbosity  (0),
  _Stream     (stream),
  _Color      (stdout_color),
  _FlushBuffer(true)
{
}

// -----------------------------------------------------------------------------
irtkDeformableSurfaceLogger::~irtkDeformableSurfaceLogger()
{
}

// -----------------------------------------------------------------------------
void irtkDeformableSurfaceLogger::HandleEvent(irtkObservable *obj, irtkEvent event, const void *data)
{
  if (!_Stream) return;
  ostream &os = *_Stream;

  // Change/Remember stream output format
  const streamsize w = os.width(0);
  const streamsize p = os.precision(6);

  // Deformable surface model (i.e., energy function)
  irtkLocalOptimizer *optimizer = dynamic_cast<irtkLocalOptimizer *>(obj);
  if (optimizer == NULL) return;
  irtkDeformableSurfaceModel *model = dynamic_cast<irtkDeformableSurfaceModel *>(optimizer->Function());
  if (model == NULL) return;

  // Typed pointers to event data (which of these is valid depends on the type of event)
  const irtkIteration      *iter = reinterpret_cast<const irtkIteration      *>(data);
  const irtkLineSearchStep *step = reinterpret_cast<const irtkLineSearchStep *>(data);
  const char               *msg  = reinterpret_cast<const char               *>(data);

  // Note: std::endl forces flush of stream buffer! Use \n instead and only flush buffer
  //       at end of this function if _FlushBuffer flag is set.

  switch (event) {

    // Start of optimization
    case StartEvent: {
      if (_Verbosity == 0) {
        if (_Color) os << xboldblack;
        os << "          Energy     ";
        string name;
        int    i;
        size_t w;
        for (i = 0; i < model->NumberOfTerms(); ++i) {
          const irtkEnergyTerm *term = model->Term(i);
          if (term->Weight() != .0) {
            w = 16;
            if (model->NumberOfTerms() < 5 &&
                term->DivideByInitialValue()) w += 10;
            name                   = term->Name();
            if (name.empty()) name = term->NameOfClass();
            if (name.length() > w-3) {
              name = name.substr(0, w-3), name += "...";
            }
            os << string(3 + (w - name.length())/2, ' ');
            os << name;
            os << string((w - name.length() + 1)/2, ' ');
          }
        }
        if (_Color) os << xreset;
      }
      _NumberOfGradientSteps = 0;
    } break;

    // Next gradient step
    case IterationEvent:
    case IterationStartEvent: {
      ++_NumberOfGradientSteps;
      os << "\n";
      if (debug_time || _Verbosity > 0) {
        os << "Iteration " << setw(3) << left << _NumberOfGradientSteps << right << " ";
      } else {
        os << " " << setw(3) << _NumberOfGradientSteps << " ";
      }
      if (debug_time) os << "\n";
      // If this event is followed by an EndEvent without an IterationEndEvent
      // or LineSearchIterationStartEvent+LineSearchIterationEndEvent before,
      // one of the convergence criteria must have been fullfilled
      _NumberOfIterations = 0;
      _NumberOfSteps      = 0;
      _Converged          = true;
    } break;

    // Start of line search given current gradient direction
    case LineSearchStartEvent: {

      model->Value(); // update cached value of individual terms if necessary

      if (_Verbosity > 0) {
        if (debug_time) os << "\n              ";
        os << "Line search with ";
        if (step->_Info) os << step->_Info << " ";
        os << "step length in [";
        os.unsetf(ios_base::floatfield);
        os << setprecision(2) << step->_MinLength << ", " << step->_MaxLength << "]";
      }
      string name;
      int i, nterms = 0;
      bool divide_by_init = false;
      for (i = 0; i < model->NumberOfTerms(); ++i) {
        const irtkEnergyTerm *term = model->Term(i);
        if (term->Weight() != .0) {
          ++nterms;
          if (term->DivideByInitialValue()) divide_by_init = true;
        }
      }
      if (_Verbosity > 0) {
        os << "\n\n  Energy = ";
        if (_Color) os << xboldblue;
        PrintNumber(os, step->_Current);
        if (_Color) os << xreset;
        for (i = 0; i < model->NumberOfTerms(); ++i) {
          const irtkEnergyTerm *term = model->Term(i);
          if (term->Weight() != .0) {
            if (nterms > 1) {
              os << " = ";
              if (term->DivideByInitialValue()) {
                PrintNormalizedNumber(os, model->Value(i));
                os << " (";
                PrintNumber(os, model->RawValue(i));
                os << ")";
              } else {
                if (divide_by_init) os << "          ";
                PrintNumber(os, model->Value(i));
                if (divide_by_init) os << " ";
              }
            }
            os << "    ";
            PrintWeight(os, term->Weight(), nterms);
            name                   = term->Name();
            if (name.empty()) name = term->NameOfClass();
            os << "% " << name << "\n";
            break;
          }
        }
        for (i = i + 1; i < model->NumberOfTerms(); ++i) {
          const irtkEnergyTerm *term = model->Term(i);
          if (term->Weight() != .0) {
            os << "                            + ";
            if (term->DivideByInitialValue()) {
              PrintNormalizedNumber(os, model->Value(i));
              os << " (";
              PrintNumber(os, model->RawValue(i));
              os << ")";
            } else {
              if (divide_by_init) os << "          ";
              PrintNumber(os, model->Value(i));
              if (divide_by_init) os << " ";
            }
            os << "    ";
            PrintWeight(os, term->Weight(), nterms);
            name                   = term->Name();
            if (name.empty()) name = term->DefaultName();
            os << "% " << name << "\n";
          }
        }
        if (_Color) os << xboldblack;
        os << "\n                 Energy        Step Length        Max. Delta\n\n";
        if (_Color) os << xreset;
      } else {
        if (_Color) os << xboldblue;
        PrintNumber(os, step->_Current);
        if (_Color) os << xreset;
        nterms = 0;
        for (i = 0; i < model->NumberOfTerms(); ++i) {
          const irtkEnergyTerm *term = model->Term(i);
          if (term->Weight() != .0) {
            if (model->NumberOfTerms() < 5) {
              if (nterms++ == 0) os << " = ";
              else               os << " + ";
              if (term->DivideByInitialValue()) {
                PrintNormalizedNumber(os, model->Value(i));
                os << " (";
                PrintNumber(os, model->RawValue(i));
                os << ")";
              } else {
                os << " ";
                PrintNumber(os, model->Value(i));
              }
            } else {
              os << "    ";
              PrintNumber(os, model->RawValue(i));
            }
          }
        }
      }
      _NumberOfIterations = 0;
      _NumberOfSteps      = 0;
      _Converged          = false;
      break;
    }

    case LineSearchIterationStartEvent:
      if (_Verbosity > 0) {
        if (debug_time) os << "\nStep " << left << setw(3) << iter->Count() << "\n";
        else            os <<   "     " << setw(3) << iter->Count() << "   ";
      }
      // If this event is followed by LineSearchEndEvent without a
      // LineSearchIterationEndEvent before, one of the convergence
      // criteria must have been fullfilled
      _Converged = true;
      break;

    case AcceptedStepEvent:
    case RejectedStepEvent: {
      if (_Verbosity > 0) {
        if (debug_time) os << "\n           ";
        if (_Color) os << (event == AcceptedStepEvent ? xgreen : xbrightred);
        PrintNumber(os, step->_Value ) << "  ";
        PrintNumber(os, step->_Length) << "  ";
        PrintNumber(os, step->_Delta );
        if (_Color) os << xreset;
        else os << "    " << ((event == AcceptedStepEvent) ? "Accepted" : "Rejected");
        os << "\n";
      }
      if (event == AcceptedStepEvent) ++_NumberOfSteps;
      // Increment counter of actual line search iterations performed
      // Note that the Delta convergence criterium on the minimum maximum
      // DoF change can cause the line search to stop immediately
      ++_NumberOfIterations;
      break;
    }

    case LineSearchIterationEndEvent: {
      if (_Verbosity > 0 && iter->Count() == iter->Total()) {
        if (_Color) os << xboldblack;
        os << "\n              Maximum number of iterations exceeded\n";
        if (_Color) os << xreset;
      }
      _Converged = false;
      break;
    }

    // End of line search
    case LineSearchEndEvent: {
      if (_Verbosity > 0) {
        // The minimum maximum DoF change convergence criterium kicked in immediately...
        if (_NumberOfIterations == 0) {
          os << "                              ";
          if (_Color) os << xboldred;
          PrintNumber(os, step->_Delta) << "\n";
        }
        // Report if line search was successful or no improvement
        if (_Color) os << xboldblack;
        if (_NumberOfSteps > 0) {
          if (_Verbosity > 0) {
            os << "\n               Step length = ";
            PrintNumber(os, step->_TotalLength) << " / ";
            PrintNumber(os, step->_Unit)        << "\n";
          }
        } else {
          os << "\n         No further improvement within search range\n";
        }
        if (_Color) os << xreset;
        os << "\n";
      }
      break;
    }

    // End of iteration without line search
    case IterationEndEvent: {
      if (_Color) os << xboldblue;
      PrintNumber(os, model->Value());
      if (_Color) os << xreset;
      int nterms = 0;
      for (int i = 0; i < model->NumberOfTerms(); ++i) {
        const irtkEnergyTerm *term = model->Term(i);
        if (term->Weight() != .0) {
          if (model->NumberOfTerms() < 5) {
            if (nterms++ == 0) os << " = ";
            else               os << " + ";
            if (term->DivideByInitialValue()) {
              PrintNormalizedNumber(os, model->Value(i));
              os << " (";
              PrintNumber(os, model->RawValue(i));
              os << ")";
            } else {
              os << " ";
              PrintNumber(os, model->Value(i));
            }
          } else {
            os << "    ";
            PrintNumber(os, model->RawValue(i));
          }
        }
      }
      ++_NumberOfSteps;
      _Converged = false;
    } break;

    // Energy function modified after convergence and optimization restarted
    case RestartEvent: {
      os << "\n";
      if (_Verbosity > 0) {
        os << "\nContinue with modified energy function\n";
      }
    } break;

    // End of optimization
    case EndEvent: {
      // IterationStartEvent followed directly by EndEvent (i.e. no line search used)
      if (_Converged && _NumberOfSteps == 0 && _NumberOfIterations == 0) {
        if (_Color) os << xboldblue;
        PrintNumber(os, model->Value());
        if (_Color) os << xboldblack;
        os << "No further improvement\n";
        if (_Color) os << xreset;
      }
      os << "\n";
      break;
    }

    // Status message broadcasted by registration filter
    case StatusEvent: {
      if (_Color) os << xboldblack;
      os << msg;
      if (_Color) os << xreset;
      break;
    }

    // Log message broadcasted by registration filter
    case LogEvent: {
      if (_Verbosity > 0) os << msg;
      break;
    }

    default: break;
  }

  // Flush output buffer
  if (_FlushBuffer) os.flush();

  // Reset stream output format
  os.width(w);
  os.precision(p);
}
