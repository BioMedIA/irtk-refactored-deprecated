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
#include <irtkDeformableSurfaceDebugger.h>


// -----------------------------------------------------------------------------
irtkDeformableSurfaceDebugger
::irtkDeformableSurfaceDebugger(const irtkDeformableSurfaceModel *model, const char *prefix)
:
  _Prefix(prefix),
  _Model(model),
  _Interval(1)
{
}

// -----------------------------------------------------------------------------
irtkDeformableSurfaceDebugger::~irtkDeformableSurfaceDebugger()
{
}

// -----------------------------------------------------------------------------
void irtkDeformableSurfaceDebugger::HandleEvent(irtkObservable *obj, irtkEvent event, const void *data)
{
  if (_Model == NULL) return;

  const int sz = 8;
  char suffix[sz];
  switch (event) {

    // -------------------------------------------------------------------------
    // Write intermediate results after each gradient step
    case IterationEvent:
    case IterationStartEvent: {
      const int iter = static_cast<const irtkIteration *>(data)->Iter();
      if (_Interval == 1 || (iter % _Interval) == 1) {
        snprintf(suffix, sz, "_%03d", iter);
        _Model->WriteDataSets(_Prefix.c_str(), suffix, debug >= 3);
      }
    } break;

    case IterationEndEvent: {
      if (debug >= 2) {
        const int iter = static_cast<const irtkIteration *>(data)->Iter();
        if (_Interval == 1 || (iter % _Interval == 1)) {
          snprintf(suffix, sz, "_%03d", iter);
          _Model->WriteGradient(_Prefix.c_str(), suffix);
        }
      }
    } break;

    // -------------------------------------------------------------------------
    // Unhandled event
    default: break;
  }
}
