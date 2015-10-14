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

#ifndef _IRTKFLOAT_H

#define _IRTKFLOAT_H

#include <limits.h>
#include <float.h>
#include <math.h>
#include <boost/math/special_functions/round.hpp>

// ---------------------------------------------------------------------------
/// Check if floating point value is not a number (NaN)
inline bool IsNaN(double x)
{
#if WINDOWS
  return _isnan(x);
#else
  return ::isnan(x);
#endif
}

// ---------------------------------------------------------------------------
/// Check if floating point value represents infinity
inline bool IsInf(double x)
{
#if WINDOWS
  return !_finite(x);
#else
  return ::isinf(x);
#endif
}

// ---------------------------------------------------------------------------
// "Overwrite" use of round function in math.h by own implementation
using boost::math::iround;

// ---------------------------------------------------------------------------
/// Determine equality of two floating point numbers
inline bool fequal(double a, double b, double tol = 1e-12)
{
  return fabs(a - b) < tol;
}

// ---------------------------------------------------------------------------
/// Sign function - https://en.wikipedia.org/wiki/Sign_function
template <typename T>
int sgn(T val)
{
  return (T(0) < val) - (val < T(0));
}

// ---------------------------------------------------------------------------
/// Round floating-point value
IRTKCU_API inline int round2(double x)
{
  return x > 0 ? int(x + 0.5) : int(x - 0.5);
}

// ---------------------------------------------------------------------------
/// Increment floating-point number by the smallest possible amount such that
/// the resulting number is greater than the original number.
inline double finc(double f)
{
  int e;
  double m = frexp(f, &e);
  return ::ldexp(m + std::numeric_limits<double>::epsilon(), e);
}

// ---------------------------------------------------------------------------
/// Decrement floating-point number by the smallest possible amount such that
/// the resulting number is less than the original number.
inline double fdec(double f)
{
  int e;
  double m = frexp(f, &e);
  return ::ldexp(m - std::numeric_limits<double>::epsilon(), e);
}

// ---------------------------------------------------------------------------
/// Increment floating point number by a given amount, ensuring that the result
/// is not equal f.
///
/// Note that due to roundoff errors, adding a small number to
/// a big number, may result in a number which is yet equal the initial big number.
/// This function adjusts the increment if necessary such that the result is
/// guaranteed to be greater (df > 0) or smaller (df < 0) than f.
/// If df is zero, f remains unchanged.
inline double finc(double f, double df)
{
  if (df == 0) return f;
  double s = f + df;
  if (s == f) {
    if (df < 0) s = fdec(f);
    else        s = finc(f);
  }
  return s;
}

// ---------------------------------------------------------------------------
/// Decrement floating point number by a given amount, ensuring that the result
/// is not equal f.
///
/// Note that due to roundoff errors, subtracting a small number
/// from a big number, may result in a number which is yet equal the initial big
/// number. This function adjusts the decrement if necessary such that the result
/// is guaranteed to be smaller (df > 0) or greater (df < 0) than f.
/// If df is zero, f remains unchanged.
inline double fdec(double f, double df)
{
  if (df == 0) return f;
  double s = f - df;
  if (s == f) {
    if (df < 0) s = finc(f);
    else        s = fdec(f);
  }
  return s;
}


#endif
