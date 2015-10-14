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

#ifndef _IRTKIMAGE_H

#define _IRTKIMAGE_H

#include <irtkCommon.h>

// Basic voxel classes
#include <irtkVoxel.h>
#include <irtkVoxelCast.h>

// Basic image class
template <class Type> class irtkGenericImage;

// Includes
#include <irtkGeometry.h>

#include <irtkBaseImage.h>
#include <irtkGenericImage.h>

#if defined(HAS_OPENCV) && !defined(__CUDACC__)
#  include <irtkImageToOpenCv.h>
#endif

// Image iteration helpers
#include <irtkImageRegion.h>
#include <irtkConstImageIterator.h>
#include <irtkConstGenericImageIterator.h>
#include <irtkImageIterator.h>
#include <irtkGenericImageIterator.h>

#endif
