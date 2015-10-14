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

#ifndef _IRTKDEFINES_H
#define _IRTKDEFINES_H


// ===========================================================================
// General
// ===========================================================================

/// Whether to build for execution on Microsoft Windows
#define WINDOWS (defined (_WIN32) || defined (WIN32) || defined (_WINDOWS))

/// Precision of floating point types to use by default
/// 0: single-precision 1: double-precision
#define USE_FLOAT_BY_DEFAULT 0

// ===========================================================================
// CUDA
// ===========================================================================

// ---------------------------------------------------------------------------
#ifndef IRTKCU_API
#  if __CUDACC__
#    define IRTKCU_API __device__ __host__
#  else
#    define IRTKCU_API
#  endif
#endif

// ---------------------------------------------------------------------------
#ifndef IRTKCU_HOST_API
#  if __CUDACC__
#    define IRTKCU_HOST_API __host__
#  else
#    define IRTKCU_HOST_API
#  endif
#endif

// ---------------------------------------------------------------------------
#ifndef IRTKCU_DEVICE_API
#  if __CUDACC__
#    define IRTKCU_DEVICE_API __device__
#  else
#    define IRTKCU_DEVICE_API
#  endif
#endif

// =============================================================================
// irtkAssert
// =============================================================================

// -----------------------------------------------------------------------------
#ifndef NDEBUG
#  define irtkAssert(condition, message)                                       \
     do {                                                                      \
       if (!(condition)) {                                                     \
         std::cerr << "Assertion `" #condition "` failed in " << __FILE__      \
                   << " line " << __LINE__ << ": " << message << std::endl;    \
         std::exit(EXIT_FAILURE);                                              \
       }                                                                       \
     } while (false)
#else
#  define irtkAssert(condition, message) do { } while (false)
#endif

// =============================================================================
// VTK 5/6 transition
// =============================================================================
#ifdef HAS_VTK
#include <vtkConfigure.h>

// Auxiliary macro to set VTK filter input
#if VTK_MAJOR_VERSION >= 6
#  define SetVTKInput(filter, dataset) (filter)->SetInputData(dataset);
#else
#  define SetVTKInput(filter, dataset) (filter)->SetInput(dataset);
#endif
// Auxiliary macro to set VTK filter input connection
#if VTK_MAJOR_VERSION >= 6
#  define SetVTKConnection(filter2, filter1) (filter2)->SetInputConnection((filter1)->GetOutputPort());
#else
#  define SetVTKConnection(filter2, filter1) (filter2)->SetInput((filter1)->GetOutput());
#endif

#endif // HAS_VTK


#endif
