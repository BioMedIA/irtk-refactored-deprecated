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

#ifndef _IRTKSINC_H
#define _IRTKSINC_H

#include <cmath>
#include <irtkCuda.h>


/**
 * Sinc function
 */

template <class TReal>
class irtkSinc
{
public:

  /// Radius of Sinc kernel
  /// (for use as Kernel::Radius, where Kernel is a typedef of irtkSinc<T>)
  IRTKCU_API static const int Radius = 6;

  /// Size of Sinc kernel
  IRTKCU_API static const int KernelSize = 2 * Radius + 1;

  /// Size of Sinc function lookup table
  /// Note that the actual array size is: LookupTableSize * (Radius + 0.5)
  IRTKCU_API static const int LookupTableSize = 1000000;

  /// Lookup table of Sinc function values
  IRTKCU_API static TReal *LookupTable;

  /// Initialize lookup table of Sinc function values
  IRTKCU_API static void Initialize();

  /// Lookup Sinc function value
  IRTKCU_API static TReal Lookup(TReal);

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
template <class TReal>
inline TReal irtkSinc<TReal>::Lookup(TReal x)
{
  return LookupTable[static_cast<int>(round(fabs(x) * LookupTableSize))];
}


#endif
