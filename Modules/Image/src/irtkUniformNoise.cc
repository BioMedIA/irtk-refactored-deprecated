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

#include <irtkImage.h>

#include <irtkNoise.h>

#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>

template <class VoxelType> irtkUniformNoise<VoxelType>::irtkUniformNoise(double Amplitude) : irtkNoise<VoxelType>(Amplitude)
{
  boost::mt19937 rng;
  rng.seed(this->_Init);

  boost::uniform_int<> ud(0, 1);
  boost::variate_generator<boost::mt19937&,
			   boost::uniform_int<> > uni_engine(rng, ud); 
  (void) uni_engine();
}

template <class VoxelType> double irtkUniformNoise<VoxelType>::Run(int x, int y, int z, int t)
{
  boost::mt19937 rng;
  rng.seed(this->_Init);

  boost::uniform_int<> ud(0, this->_Amplitude);
  boost::variate_generator<boost::mt19937&,
			   boost::uniform_int<> > uni_engine(rng, ud);
  return double(this->_input->Get(x, y, z, t)) + uni_engine();
}

template class irtkUniformNoise<irtkBytePixel>;
template class irtkUniformNoise<irtkGreyPixel>;
template class irtkUniformNoise<float>;
template class irtkUniformNoise<double>;
