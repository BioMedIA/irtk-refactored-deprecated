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

#include <irtkVector4D.h>

template <typename T> irtkVector4D<T> irtkVector4D<T>::operator/(const irtkVector4D<T>& v) const
{
  irtkVector4D<T> val(0, 0, 0, 0);

  if (v._x != 0) {
    val._x = _x/v._x;
  }

  if (v._y != 0) {
    val._y = _y/v._y;
  }

  if (v._z != 0) {
    val._z = _z/v._z;
  }

  if (v._t != 0) {
    val._t = _t/v._t;
  }

  return val;
}

template <typename T> irtkVector4D<T>& irtkVector4D<T>::operator/=(const irtkVector4D<T>& v)
{
  if (v._x != 0) {
    _x /= v._x;
  }

  if (v._y != 0) {
    _y /= v._y;
  }

  if (v._z != 0) {
    _z /= v._z;
  }

  if (v._t != 0) {
    _t /= v._t;
  }

  return *this;
}

template struct irtkVector4D<char>;
template struct irtkVector4D<short>;
template struct irtkVector4D<float>;
template struct irtkVector4D<double>;
