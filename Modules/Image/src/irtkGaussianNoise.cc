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
#include <boost/random/normal_distribution.hpp>


template <class VoxelType> irtkGaussianNoise<VoxelType>::irtkGaussianNoise() : irtkNoise<VoxelType>()
{
  _Mean   = 0;
  _Sigma  = 1;
  _MinVal = VoxelType(MIN_GREY);
  _MaxVal = VoxelType(MAX_GREY);

  long temp = -1 * this->_Init;

  boost::mt19937 rng;
  rng.seed(temp);

  boost::normal_distribution<> nd(0, 1);
  boost::variate_generator<boost::mt19937&,
                           boost::normal_distribution<> > var_nor(rng, nd);
  (void) var_nor();
}

template <class VoxelType> irtkGaussianNoise<VoxelType>::irtkGaussianNoise(double Mean, double Sigma, VoxelType MinVal, VoxelType MaxVal) : irtkNoise<VoxelType>()
{
  this->_Mean   = Mean;
  this->_Sigma  = Sigma;
  this->_MinVal = MinVal;
  this->_MaxVal = MaxVal;

  long temp = -1 * this->_Init;

  boost::mt19937 rng;
  rng.seed(temp);

  boost::normal_distribution<> nd(0, 1);
  boost::variate_generator<boost::mt19937&,
                           boost::normal_distribution<> > var_nor(rng, nd);
  (void) var_nor();
}

template <class VoxelType> double irtkGaussianNoise<VoxelType>::Run(int x, int y, int z, int t)
{
  boost::mt19937 rng;
  rng.seed(this->_Init);

  boost::normal_distribution<> nd(0, 1);
  boost::variate_generator<boost::mt19937&,
                           boost::normal_distribution<> > var_nor(rng, nd);

  double tmp = this->_input->Get(x, y, z, t) + this->_Sigma * var_nor() + this->_Mean;
  if (tmp < this->_MinVal) return this->_MinVal;
  if (tmp > this->_MaxVal) return this->_MaxVal;
  return tmp;
}

template class irtkGaussianNoise<irtkBytePixel>;
template class irtkGaussianNoise<irtkGreyPixel>;
template class irtkGaussianNoise<float>;
template class irtkGaussianNoise<double>;
