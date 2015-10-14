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

#include <irtkImageFunction.h>

// -----------------------------------------------------------------------------
irtkExtrapolateImageFunction::irtkExtrapolateImageFunction()
{
}

// -----------------------------------------------------------------------------
irtkExtrapolateImageFunction::~irtkExtrapolateImageFunction()
{
}

// -----------------------------------------------------------------------------
irtkExtrapolateImageFunction *
irtkExtrapolateImageFunction::New(irtkExtrapolationMode mode, const irtkBaseImage *image)
{
  irtkExtrapolateImageFunction *p = NULL;
  switch (mode) {
    case Extrapolation_None:   { p = NULL;                                              break; }
    case Extrapolation_Const:  { p = new irtkConstExtrapolateImageFunction();           break; }
    case Extrapolation_NN:     { p = new irtkNearestNeighborExtrapolateImageFunction(); break; }
    case Extrapolation_Repeat: { p = new irtkRepeatExtrapolateImageFunction();          break; }
    case Extrapolation_Mirror: { p = new irtkMirrorExtrapolateImageFunction();          break; }
    case Extrapolation_ConstWithPeriodicTime:
      p = new irtkConstExtrapolateImageFunctionWithPeriodicTime();
      break;
    default:
   	  cerr << "irtkExtrapolateImageFunction::New: Unknwon extrapolation mode: " << mode  << endl;
      exit(1);
  }
  if (p) p->Input(image);
  return p;
}
