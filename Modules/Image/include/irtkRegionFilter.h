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

#ifndef IRTKREGIONFILTER_H

#define IRTKREGIONFILTER_H

#include <irtkImageToImage2.h>

/**
 * Class for extracting regions for images
 *
 */

class irtkRegionFilter : public irtkImageToImage2
{
  irtkImageFilterMacro(irtkRegionFilter);

protected:

	/// Region of interest
	int _i1, _j1, _k1, _l1, _i2, _j2, _k2, _l2; 

public:

  /// Constructor
	irtkRegionFilter();

  /// Destructor
  ~irtkRegionFilter();

  /// Define ROI
  virtual void PutRegion(int, int, int, int, int, int, int, int);
  
  /// Extract ROI
  virtual void Run();

};

inline void irtkRegionFilter::PutRegion(int i1, int j1, int k1, int l1, int i2, int j2, int k2, int l2)
{
	_i1 = i1;
	_j1 = j1;
	_k1 = k1;
	_l1 = l1;
	_i2 = i2;
	_j2 = j2;
	_k2 = k2;
	_l2 = l2;
}

#endif 
