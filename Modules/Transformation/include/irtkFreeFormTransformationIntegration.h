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

#ifndef _IRTKFREEFORMTRANSFORMATIONINTEGRATION_H
#define _IRTKFREEFORMTRANSFORMATIONINTEGRATION_H

#include <irtkEnums.h> // FFDIntegrationMethods


// ============================================================================
// Auxiliary macros
// ============================================================================

// ---------------------------------------------------------------------------
/// Define integration method for particular FFD
#define IRTK_FFDIM2(NAME, FFD) \
  typedef irtkFreeFormTransformationIntegration##NAME<FFD> NAME;

// ---------------------------------------------------------------------------
/// Define integration method for particular FFD
#define IRTK_FFDIM3(NAME, FFD, METHOD) \
  typedef irtkFreeFormTransformationIntegration##METHOD<FFD> NAME;

#ifdef USE_CUDA

// ---------------------------------------------------------------------------
/// Define integration method for particular FFD
#define IRTKCU_FFDIM2(NAME, FFD) \
  typedef irtkCUFreeFormTransformationIntegration##NAME<FFD> NAME;

// ---------------------------------------------------------------------------
/// Define integration method for particular FFD
#define IRTKCU_FFDIM3(NAME, FFD, METHOD) \
  typedef irtkCUFreeFormTransformationIntegration##METHOD<FFD> NAME;

#endif

// ============================================================================
// Integration methods
// ============================================================================

// Runge-Kutta integration methods
#include <irtkFreeFormTransformationRungeKutta.h>


#endif
