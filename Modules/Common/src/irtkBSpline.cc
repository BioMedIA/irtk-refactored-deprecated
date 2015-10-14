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

#include <irtkBSpline.h>

// Lookup table of B-spline function values
template <class TReal> TReal irtkBSpline<TReal>::WeightLookupTable[irtkBSpline<TReal>::LookupTableSize];

// Lookup table of B-spline basis function values
template <class TReal> TReal irtkBSpline<TReal>::LookupTable   [irtkBSpline<TReal>::LookupTableSize][4];
template <class TReal> TReal irtkBSpline<TReal>::LookupTable_I [irtkBSpline<TReal>::LookupTableSize][4];
template <class TReal> TReal irtkBSpline<TReal>::LookupTable_II[irtkBSpline<TReal>::LookupTableSize][4];

// Wether lookup tables of B-spline kernel were initialized
template <class TReal> bool  irtkBSpline<TReal>::_initialized = false;

// Explicit template instantiations
template class irtkBSpline<float>;
template class irtkBSpline<double>;
