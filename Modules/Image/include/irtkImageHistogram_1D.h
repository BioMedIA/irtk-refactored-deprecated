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

#ifndef _IRTKIMAGEHISTOGRAM_1D_H

#define _IRTKIMAGEHISTOGRAM_1D_H

#include <irtkGenericImage.h>
#include <irtkHistogram_1D.h>

/** Class for 2D histograms.
 *  This class defines and implements 2D histograms.
 */

template <class VoxelType> class irtkImageHistogram_1D : public irtkHistogram_1D<double>
{
protected:
    /// min for equalize
    VoxelType _emin;
    /// max for equalize
    VoxelType _emax;

public:
	/// Evaluate the histogram from a given image with padding value
	virtual void Evaluate(irtkGenericImage<VoxelType> *, double padding = -10000);
	/// Histogram Equalization
	virtual void Equalize(VoxelType min,VoxelType max);
	/// Back project the equalized histogram to image
	virtual void BackProject(irtkGenericImage<VoxelType> *); 
};

#endif
