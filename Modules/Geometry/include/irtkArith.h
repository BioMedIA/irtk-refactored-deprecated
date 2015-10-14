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

#ifndef IRTKARITH_H

#define IRTKARITH_H


class irtkArith
{
public:
  static int IAbs (int iValue);
  static int ICeil (float fValue);
  static int IFloor (float fValue);
  static int ISign (int iValue);

  static double Abs (double fValue);
  static double ACos (double fValue);
  static double ASin (double fValue);
  static double ATan (double fValue);
  static double ATan2 (double fY, double fX);
  static double Ceil (double fValue);
  static double Cos (double fValue);
  static double Exp (double fValue);
  static double Floor (double fValue);
  static double Log (double fValue);
  static double Pow (double kBase, double kExponent);
  static double Sign (double fValue);
  static double Sin (double fValue);
  static double Sqr (double fValue);
  static double Sqrt (double fValue);
  static double UnitRandom ();  // in [0,1]
  static double SymmetricRandom ();  // in [-1,1]

};

#endif
