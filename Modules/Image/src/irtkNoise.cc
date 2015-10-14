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

#include <sys/types.h>

#ifndef WIN32
#include <sys/time.h>
#endif

#include <irtkImage.h>

#include <irtkNoise.h>

template <class VoxelType> irtkNoise<VoxelType>::irtkNoise(double Amplitude)
{
#ifndef WIN32
  timeval tv;
  gettimeofday(&tv, NULL);
  _Init = tv.tv_usec;
  _Amplitude = Amplitude;
#else
  cerr << "irtkNoise<VoxelType>::irtkNoise: Not implemented yet for Windows" << endl;
  exit(1);
#endif
}

template class irtkNoise<irtkBytePixel>;
template class irtkNoise<irtkGreyPixel>;
template class irtkNoise<float>;
template class irtkNoise<double>;
