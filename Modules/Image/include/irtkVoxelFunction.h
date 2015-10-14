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

#ifndef _IRTKVOXELFUNCTION_H
#define _IRTKVOXELFUNCTION_H

#include <irtkImageDomain.h>


////////////////////////////////////////////////////////////////////////////////
// Base classes
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Base class for voxel functions
// =============================================================================

/**
 * Base class for voxel functions
 *
 * The basic voxel functions which are part of this library and included by the
 * irtkVoxelFunction.h module demonstrate the implementation and use of the
 * ForEachVoxel function templates. These templates can be well optimized by the
 * compiler and provide the fastest way of iterating over all voxels in a given
 * image region (or the entire image) of one or more images. When processing
 * multiple images at the same time, these must have the same image size.
 *
 * Example usage:
 * \code
 * // Determine intensity range of grayscale image
 * irtkGreyImage image(attr);
 * irtkUnaryVoxelFunction::GetMinMax minmax;
 * ForEachVoxel(image, minmax);
 * double min = minmax.GetMinAsDouble();
 * double max = minmax.GetMaxAsDouble();
 *
 * // Add three images within a given region and store result in another
 * irtkNaryVoxelFunction::Sum sum;
 * irtkGreyImage input1(attr);
 * irtkGreyImage input2(attr);
 * irtkGreyImage input3(attr);
 * irtkGreyImage output(attr);
 * blocked_range3d<int> region(1, 2, 0, attr._y, 0, attr._x);
 * ForEachVoxel(region, input1, input2, input3, output, sum);
 *
 * // Subtract one image from another
 * // Note: Voxel function passed in this case by copy instead of reference!
 * //       This is only possible if it is not a voxel reduction.
 * irtkGreyImage minuend   (attr);
 * irtkGreyImage subtrahend(attr);
 * ParallelForEachVoxel(irtkBinaryVoxelFunction::Sub(), subtrahend, minuend);
 * \endcode
 *
 * If a voxel function needs to consider also the intensities of neighboring
 * voxels, it can access these by manipulating the given pointer to the current
 * voxel to be processed. A useful helper for this is the irtkImageRegion class.
 *
 * \sa irtkImageRegion
 *
 */
struct irtkVoxelFunction
{
  /// Used by ParallelForEachVoxel to determine if voxel function has to
  /// be executed using parallel_reduce
  static bool IsReduction() { return false; }

  /// Split "constructor"
  ///
  /// \note A method is used instead of an actual split constructor such that
  ///       subclasses which are not used for reduction do not need to
  ///       implement such constructor.
  void split(irtkVoxelFunction &) {}

  /// Join results
  void join(irtkVoxelFunction &) {}
};


/**
 * Base class for voxel functions which implement a reduction of voxel values
 *
 * Voxel functions derived from this base class are run by ParallelForEachVoxel
 * using parallel_reduce instead of parallel_for. Typical voxel reduction
 * functions are unary functions computing image statistics such as the min,
 * max, mean, or variance of voxel values within an image region.
 *
 */
struct irtkVoxelReduction : public irtkVoxelFunction
{
  /// Used by ParallelForEachVoxel to determine if voxel function has to
  /// be executed using parallel_reduce
  static bool IsReduction() { return true; }

  /// Split "constructor"
  void split(irtkVoxelFunction &)
  {
    cerr << "irtkVoxelReduction::split must be overriden by each subclass!" << endl;
    cerr << "Otherwise you should use irtkVoxelFunction with parallel_for instead." << endl;
    exit(1);
  }

  /// Join results
  void join(irtkVoxelFunction &)
  {
    cerr << "irtkVoxelReduction::join must be overriden by each subclass!" << endl;
    cerr << "Otherwise you should use irtkVoxelFunction with parallel_for instead." << endl;
    exit(1);
  }
};

// =============================================================================
// Base class for ForEachVoxel body with single voxel function
// =============================================================================

/**
 * Base class for ForEachVoxel template function body with single voxel
 * function for each voxel
 */
template <class VoxelFunc>
struct irtkForEachVoxelBody
{
  // ---------------------------------------------------------------------------
  // Members

  VoxelFunc _VoxelFunc; ///< Functor executed for each voxel
  int       _k, _l;     ///< Indices for fixed dimensions

  // ---------------------------------------------------------------------------
  // Construction

  /// Constructor
  irtkForEachVoxelBody(const VoxelFunc &vf)
  :
    _VoxelFunc(vf), _k(0), _l(0)
  {}

  /// Copy constructor
  irtkForEachVoxelBody(const irtkForEachVoxelBody &o)
  :
    _VoxelFunc(o._VoxelFunc), _k(o._k), _l(o._l)
  {}

  // ---------------------------------------------------------------------------
  // Parallel reduction

  /// Split constructor
  irtkForEachVoxelBody(irtkForEachVoxelBody &o, split s)
  :
    _VoxelFunc(o._VoxelFunc), _k(o._k), _l(o._l)
  {
    _VoxelFunc.split(o._VoxelFunc);
  }

  /// Join results
  void join(irtkForEachVoxelBody &rhs)
  {
    _VoxelFunc.join(rhs._VoxelFunc);
  }
};

// =============================================================================
// Base class for ForEachVoxelIf body with separate inside/outside voxel functions
// =============================================================================

/**
 * Base class for ForEachVoxelIf template function body with separate voxel
 * function for inside and outside voxels
 */
template <class VoxelFunc, class OutsideFunc>
struct irtkForEachVoxelIfBody : public irtkForEachVoxelBody<VoxelFunc>
{
  // ---------------------------------------------------------------------------
  // Members

  OutsideFunc _OutsideFunc; ///< Functor executed for each background voxel

  // ---------------------------------------------------------------------------
  // Construction

  /// Constructor
  irtkForEachVoxelIfBody(const VoxelFunc &vf, const OutsideFunc &of)
  :
    irtkForEachVoxelBody<VoxelFunc>(vf), _OutsideFunc(of)
  {}

  /// Copy constructor
  irtkForEachVoxelIfBody(const irtkForEachVoxelIfBody &o)
  :
    irtkForEachVoxelBody<VoxelFunc>(o), _OutsideFunc(o._OutsideFunc)
  {}

  // ---------------------------------------------------------------------------
  // Parallel reduction

  /// Split constructor
  irtkForEachVoxelIfBody(irtkForEachVoxelIfBody &o, split s)
  :
    irtkForEachVoxelBody<VoxelFunc>(o, s),
    _OutsideFunc(o._OutsideFunc)
  {
    _OutsideFunc.split(o._OutsideFunc);
  }

  /// Join results
  void join(irtkForEachVoxelIfBody &rhs)
  {
    irtkForEachVoxelBody<VoxelFunc>::join(rhs);
    _OutsideFunc.join(rhs._OutsideFunc);
  }
};

// =============================================================================
// Basic voxel functions - also for demonstrating their implementation
// =============================================================================

// Note: In the following, basic and generally useful voxel functions should
//       be defined. More specialized voxel functions which might only be
//       used by a single image filter implementation should be defined in the
//       .cc module where they are actually used.

#include <irtkUnaryVoxelFunction.h>
#include <irtkBinaryVoxelFunction.h>
#include <irtkTernaryVoxelFunction.h>
//#include <irtkQuaternaryVoxelFunction.h>
//#include <irtkQuinternaryVoxelFunction.h>
// ...
#include <irtkNaryVoxelFunction.h> // i.e., NOP voxel function

// =============================================================================
// ForEachVoxel template functions for various number of image arguments
// =============================================================================

#include <irtkForEachUnaryVoxelFunction.h>
#include <irtkForEachBinaryVoxelFunction.h>
#include <irtkForEachTernaryVoxelFunction.h>
#include <irtkForEachQuaternaryVoxelFunction.h>
#include <irtkForEachQuinaryVoxelFunction.h>
#include <irtkForEachSenaryVoxelFunction.h>
#include <irtkForEachSeptenaryVoxelFunction.h>
#include <irtkForEachOctaryVoxelFunction.h>
#include <irtkForEachNonaryVoxelFunction.h>


#endif
