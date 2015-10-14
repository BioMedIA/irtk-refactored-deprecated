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

#ifndef _irtkForEachSeptenaryVoxelFunction_H
#define _irtkForEachSeptenaryVoxelFunction_H

#include <irtkVoxelFunction.h>


inline void _irtkforeachseptenaryvoxelfunction_must_not_be_reduction()
{
  cerr << "(Parallel)ForEachVoxel(If): Voxel reductions must be passed by reference!"
          " Pass voxel functor object(s) as last argument(s) instead of first." << endl;
  exit(1);
}


// =============================================================================
// 7 const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 7 const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
struct irtkSeptenaryForEachVoxelBody_Const : public irtkForEachVoxelBody<VoxelFunc>
{
  const irtkGenericImage<T1> &im1;
  const irtkGenericImage<T2> &im2;
  const irtkGenericImage<T3> &im3;
  const irtkGenericImage<T4> &im4;
  const irtkGenericImage<T5> &im5;
  const irtkGenericImage<T6> &im6;
  const irtkGenericImage<T7> &im7;

  /// Constructor
  irtkSeptenaryForEachVoxelBody_Const(const irtkGenericImage<T1> &im1,
                                      const irtkGenericImage<T2> &im2,
                                      const irtkGenericImage<T3> &im3,
                                      const irtkGenericImage<T4> &im4,
                                      const irtkGenericImage<T5> &im5,
                                      const irtkGenericImage<T6> &im6,
                                      const irtkGenericImage<T7> &im7,
                                      VoxelFunc &vf)
  :
    irtkForEachVoxelBody<VoxelFunc>(vf), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  irtkSeptenaryForEachVoxelBody_Const(const irtkSeptenaryForEachVoxelBody_Const &o)
  :
    irtkForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  irtkSeptenaryForEachVoxelBody_Const(irtkSeptenaryForEachVoxelBody_Const &o, split s)
  :
    irtkForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const irtkImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
    const T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody_Const *>(this)->_VoxelFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 1D image region
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
    const T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody_Const *>(this)->_VoxelFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody_Const *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody_Const *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 7 const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7,
          class VoxelFunc, class OutsideFunc = irtkNaryVoxelFunction::NOP,
          class Domain = irtkImageDomain::Foreground>
struct irtkSeptenaryForEachVoxelIfBody_Const : public irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  const irtkGenericImage<T1> &im1;
  const irtkGenericImage<T2> &im2;
  const irtkGenericImage<T3> &im3;
  const irtkGenericImage<T4> &im4;
  const irtkGenericImage<T5> &im5;
  const irtkGenericImage<T6> &im6;
  const irtkGenericImage<T7> &im7;

  /// Constructor
  irtkSeptenaryForEachVoxelIfBody_Const(const irtkGenericImage<T1> &im1,
                                        const irtkGenericImage<T2> &im2,
                                        const irtkGenericImage<T3> &im3,
                                        const irtkGenericImage<T4> &im4,
                                        const irtkGenericImage<T5> &im5,
                                        const irtkGenericImage<T6> &im6,
                                        const irtkGenericImage<T7> &im7,
                                        VoxelFunc &vf, OutsideFunc &of)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  irtkSeptenaryForEachVoxelIfBody_Const(const irtkSeptenaryForEachVoxelIfBody_Const &o)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  irtkSeptenaryForEachVoxelIfBody_Const(irtkSeptenaryForEachVoxelIfBody_Const &o, split s)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const irtkImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
    const T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      if (Domain::IsInside(im7, i, j, k, l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody_Const *>(this)->_VoxelFunc  (i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody_Const *>(this)->_OutsideFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 1D image region
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
    const T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, idx, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody_Const *>(this)->_VoxelFunc  (im7, idx, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody_Const *>(this)->_OutsideFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, this->_k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody_Const *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody_Const *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody_Const *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody_Const *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
// ForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, const irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, const irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// =============================================================================
// 6 const, 1 non-const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 6 const, 1 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
struct irtkSeptenaryForEachVoxelBody_6Const : public irtkForEachVoxelBody<VoxelFunc>
{
  const irtkGenericImage<T1> &im1;
  const irtkGenericImage<T2> &im2;
  const irtkGenericImage<T3> &im3;
  const irtkGenericImage<T4> &im4;
  const irtkGenericImage<T5> &im5;
  const irtkGenericImage<T6> &im6;
        irtkGenericImage<T7> &im7;

  /// Constructor
  irtkSeptenaryForEachVoxelBody_6Const(const irtkGenericImage<T1> &im1,
                                       const irtkGenericImage<T2> &im2,
                                       const irtkGenericImage<T3> &im3,
                                       const irtkGenericImage<T4> &im4,
                                       const irtkGenericImage<T5> &im5,
                                       const irtkGenericImage<T6> &im6,
                                             irtkGenericImage<T7> &im7,
                                       VoxelFunc &vf)
  :
    irtkForEachVoxelBody<VoxelFunc>(vf), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  irtkSeptenaryForEachVoxelBody_6Const(const irtkSeptenaryForEachVoxelBody_6Const &o)
  :
    irtkForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  irtkSeptenaryForEachVoxelBody_6Const(irtkSeptenaryForEachVoxelBody_6Const &o, split s)
  :
    irtkForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const irtkImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody_6Const *>(this)->_VoxelFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 1D image region
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody_6Const *>(this)->_VoxelFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody_6Const *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody_6Const *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 6 const, 1 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7,
          class VoxelFunc, class OutsideFunc = irtkNaryVoxelFunction::NOP,
          class Domain = irtkImageDomain::Foreground>
struct irtkSeptenaryForEachVoxelIfBody_6Const : public irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  const irtkGenericImage<T1> &im1;
  const irtkGenericImage<T2> &im2;
  const irtkGenericImage<T3> &im3;
  const irtkGenericImage<T4> &im4;
  const irtkGenericImage<T5> &im5;
  const irtkGenericImage<T6> &im6;
        irtkGenericImage<T7> &im7;

  /// Constructor
  irtkSeptenaryForEachVoxelIfBody_6Const(const irtkGenericImage<T1> &im1,
                                         const irtkGenericImage<T2> &im2,
                                         const irtkGenericImage<T3> &im3,
                                         const irtkGenericImage<T4> &im4,
                                         const irtkGenericImage<T5> &im5,
                                         const irtkGenericImage<T6> &im6,
                                               irtkGenericImage<T7> &im7,
                                         VoxelFunc &vf, OutsideFunc &of)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  irtkSeptenaryForEachVoxelIfBody_6Const(const irtkSeptenaryForEachVoxelIfBody_6Const &o)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  irtkSeptenaryForEachVoxelIfBody_6Const(irtkSeptenaryForEachVoxelIfBody_6Const &o, split s)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const irtkImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      if (Domain::IsInside(im7, i, j, k, l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody_6Const *>(this)->_VoxelFunc  (i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody_6Const *>(this)->_OutsideFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 1D image region
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, idx, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody_6Const *>(this)->_VoxelFunc  (im7, idx, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody_6Const *>(this)->_OutsideFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, this->_k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody_6Const *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody_6Const *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody_6Const *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody_6Const *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
// ForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_6Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// =============================================================================
// 5 const, 2 non-const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 5 const, 2 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
struct irtkSeptenaryForEachVoxelBody_5Const : public irtkForEachVoxelBody<VoxelFunc>
{
  const irtkGenericImage<T1> &im1;
  const irtkGenericImage<T2> &im2;
  const irtkGenericImage<T3> &im3;
  const irtkGenericImage<T4> &im4;
  const irtkGenericImage<T5> &im5;
        irtkGenericImage<T6> &im6;
        irtkGenericImage<T7> &im7;

  /// Constructor
  irtkSeptenaryForEachVoxelBody_5Const(const irtkGenericImage<T1> &im1,
                                       const irtkGenericImage<T2> &im2,
                                       const irtkGenericImage<T3> &im3,
                                       const irtkGenericImage<T4> &im4,
                                       const irtkGenericImage<T5> &im5,
                                             irtkGenericImage<T6> &im6,
                                             irtkGenericImage<T7> &im7,
                                       VoxelFunc &vf)
  :
    irtkForEachVoxelBody<VoxelFunc>(vf), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  irtkSeptenaryForEachVoxelBody_5Const(const irtkSeptenaryForEachVoxelBody_5Const &o)
  :
    irtkForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  irtkSeptenaryForEachVoxelBody_5Const(irtkSeptenaryForEachVoxelBody_5Const &o, split s)
  :
    irtkForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const irtkImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody_5Const *>(this)->_VoxelFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 1D image region
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody_5Const *>(this)->_VoxelFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody_5Const *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody_5Const *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 5 const, 2 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7,
          class VoxelFunc, class OutsideFunc = irtkNaryVoxelFunction::NOP,
          class Domain = irtkImageDomain::Foreground>
struct irtkSeptenaryForEachVoxelIfBody_5Const : public irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  const irtkGenericImage<T1> &im1;
  const irtkGenericImage<T2> &im2;
  const irtkGenericImage<T3> &im3;
  const irtkGenericImage<T4> &im4;
  const irtkGenericImage<T5> &im5;
        irtkGenericImage<T6> &im6;
        irtkGenericImage<T7> &im7;

  /// Constructor
  irtkSeptenaryForEachVoxelIfBody_5Const(const irtkGenericImage<T1> &im1,
                                         const irtkGenericImage<T2> &im2,
                                         const irtkGenericImage<T3> &im3,
                                         const irtkGenericImage<T4> &im4,
                                         const irtkGenericImage<T5> &im5,
                                               irtkGenericImage<T6> &im6,
                                               irtkGenericImage<T7> &im7,
                                         VoxelFunc &vf, OutsideFunc &of)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  irtkSeptenaryForEachVoxelIfBody_5Const(const irtkSeptenaryForEachVoxelIfBody_5Const &o)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  irtkSeptenaryForEachVoxelIfBody_5Const(irtkSeptenaryForEachVoxelIfBody_5Const &o, split s)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const irtkImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      if (Domain::IsInside(im7, i, j, k, l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody_5Const *>(this)->_VoxelFunc  (i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody_5Const *>(this)->_OutsideFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 1D image region
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, idx, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody_5Const *>(this)->_VoxelFunc  (im7, idx, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody_5Const *>(this)->_OutsideFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, this->_k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody_5Const *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody_5Const *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody_5Const *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody_5Const *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
// ForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// =============================================================================
// 4 const, 3 non-const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 4 const, 3 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
struct irtkSeptenaryForEachVoxelBody_4Const : public irtkForEachVoxelBody<VoxelFunc>
{
  const irtkGenericImage<T1> &im1;
  const irtkGenericImage<T2> &im2;
  const irtkGenericImage<T3> &im3;
  const irtkGenericImage<T4> &im4;
        irtkGenericImage<T5> &im5;
        irtkGenericImage<T6> &im6;
        irtkGenericImage<T7> &im7;

  /// Constructor
  irtkSeptenaryForEachVoxelBody_4Const(const irtkGenericImage<T1> &im1,
                                       const irtkGenericImage<T2> &im2,
                                       const irtkGenericImage<T3> &im3,
                                       const irtkGenericImage<T4> &im4,
                                             irtkGenericImage<T5> &im5,
                                             irtkGenericImage<T6> &im6,
                                             irtkGenericImage<T7> &im7,
                                       VoxelFunc &vf)
  :
    irtkForEachVoxelBody<VoxelFunc>(vf), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  irtkSeptenaryForEachVoxelBody_4Const(const irtkSeptenaryForEachVoxelBody_4Const &o)
  :
    irtkForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  irtkSeptenaryForEachVoxelBody_4Const(irtkSeptenaryForEachVoxelBody_4Const &o, split s)
  :
    irtkForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const irtkImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody_4Const *>(this)->_VoxelFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 1D image region
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody_4Const *>(this)->_VoxelFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody_4Const *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody_4Const *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 4 const, 3 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7,
          class VoxelFunc, class OutsideFunc = irtkNaryVoxelFunction::NOP,
          class Domain = irtkImageDomain::Foreground>
struct irtkSeptenaryForEachVoxelIfBody_4Const : public irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  const irtkGenericImage<T1> &im1;
  const irtkGenericImage<T2> &im2;
  const irtkGenericImage<T3> &im3;
  const irtkGenericImage<T4> &im4;
        irtkGenericImage<T5> &im5;
        irtkGenericImage<T6> &im6;
        irtkGenericImage<T7> &im7;

  /// Constructor
  irtkSeptenaryForEachVoxelIfBody_4Const(const irtkGenericImage<T1> &im1,
                                         const irtkGenericImage<T2> &im2,
                                         const irtkGenericImage<T3> &im3,
                                         const irtkGenericImage<T4> &im4,
                                               irtkGenericImage<T5> &im5,
                                               irtkGenericImage<T6> &im6,
                                               irtkGenericImage<T7> &im7,
                                         VoxelFunc &vf, OutsideFunc &of)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  irtkSeptenaryForEachVoxelIfBody_4Const(const irtkSeptenaryForEachVoxelIfBody_4Const &o)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  irtkSeptenaryForEachVoxelIfBody_4Const(irtkSeptenaryForEachVoxelIfBody_4Const &o, split s)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const irtkImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      if (Domain::IsInside(im7, i, j, k, l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody_4Const *>(this)->_VoxelFunc  (i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody_4Const *>(this)->_OutsideFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 1D image region
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, idx, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody_4Const *>(this)->_VoxelFunc  (im7, idx, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody_4Const *>(this)->_OutsideFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, this->_k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody_4Const *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody_4Const *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody_4Const *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody_4Const *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
// ForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// =============================================================================
// 3 const, 4 non-const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 3 const, 4 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
struct irtkSeptenaryForEachVoxelBody_3Const : public irtkForEachVoxelBody<VoxelFunc>
{
  const irtkGenericImage<T1> &im1;
  const irtkGenericImage<T2> &im2;
  const irtkGenericImage<T3> &im3;
        irtkGenericImage<T4> &im4;
        irtkGenericImage<T5> &im5;
        irtkGenericImage<T6> &im6;
        irtkGenericImage<T7> &im7;

  /// Constructor
  irtkSeptenaryForEachVoxelBody_3Const(const irtkGenericImage<T1> &im1,
                                       const irtkGenericImage<T2> &im2,
                                       const irtkGenericImage<T3> &im3,
                                             irtkGenericImage<T4> &im4,
                                             irtkGenericImage<T5> &im5,
                                             irtkGenericImage<T6> &im6,
                                             irtkGenericImage<T7> &im7,
                                       VoxelFunc &vf)
  :
    irtkForEachVoxelBody<VoxelFunc>(vf), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  irtkSeptenaryForEachVoxelBody_3Const(const irtkSeptenaryForEachVoxelBody_3Const &o)
  :
    irtkForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  irtkSeptenaryForEachVoxelBody_3Const(irtkSeptenaryForEachVoxelBody_3Const &o, split s)
  :
    irtkForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const irtkImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody_3Const *>(this)->_VoxelFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 1D image region
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody_3Const *>(this)->_VoxelFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody_3Const *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody_3Const *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 3 const, 4 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7,
          class VoxelFunc, class OutsideFunc = irtkNaryVoxelFunction::NOP,
          class Domain = irtkImageDomain::Foreground>
struct irtkSeptenaryForEachVoxelIfBody_3Const : public irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  const irtkGenericImage<T1> &im1;
  const irtkGenericImage<T2> &im2;
  const irtkGenericImage<T3> &im3;
        irtkGenericImage<T4> &im4;
        irtkGenericImage<T5> &im5;
        irtkGenericImage<T6> &im6;
        irtkGenericImage<T7> &im7;

  /// Constructor
  irtkSeptenaryForEachVoxelIfBody_3Const(const irtkGenericImage<T1> &im1,
                                         const irtkGenericImage<T2> &im2,
                                         const irtkGenericImage<T3> &im3,
                                               irtkGenericImage<T4> &im4,
                                               irtkGenericImage<T5> &im5,
                                               irtkGenericImage<T6> &im6,
                                               irtkGenericImage<T7> &im7,
                                         VoxelFunc &vf, OutsideFunc &of)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  irtkSeptenaryForEachVoxelIfBody_3Const(const irtkSeptenaryForEachVoxelIfBody_3Const &o)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  irtkSeptenaryForEachVoxelIfBody_3Const(irtkSeptenaryForEachVoxelIfBody_3Const &o, split s)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const irtkImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      if (Domain::IsInside(im7, i, j, k, l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody_3Const *>(this)->_VoxelFunc  (i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody_3Const *>(this)->_OutsideFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 1D image region
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, idx, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody_3Const *>(this)->_VoxelFunc  (im7, idx, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody_3Const *>(this)->_OutsideFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, this->_k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody_3Const *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody_3Const *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody_3Const *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody_3Const *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
// ForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// =============================================================================
// 2 const, 5 non-const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 2 const, 5 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
struct irtkSeptenaryForEachVoxelBody_2Const : public irtkForEachVoxelBody<VoxelFunc>
{
  const irtkGenericImage<T1> &im1;
  const irtkGenericImage<T2> &im2;
        irtkGenericImage<T3> &im3;
        irtkGenericImage<T4> &im4;
        irtkGenericImage<T5> &im5;
        irtkGenericImage<T6> &im6;
        irtkGenericImage<T7> &im7;

  /// Constructor
  irtkSeptenaryForEachVoxelBody_2Const(const irtkGenericImage<T1> &im1,
                                       const irtkGenericImage<T2> &im2,
                                             irtkGenericImage<T3> &im3,
                                             irtkGenericImage<T4> &im4,
                                             irtkGenericImage<T5> &im5,
                                             irtkGenericImage<T6> &im6,
                                             irtkGenericImage<T7> &im7,
                                       VoxelFunc &vf)
  :
    irtkForEachVoxelBody<VoxelFunc>(vf), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  irtkSeptenaryForEachVoxelBody_2Const(const irtkSeptenaryForEachVoxelBody_2Const &o)
  :
    irtkForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  irtkSeptenaryForEachVoxelBody_2Const(irtkSeptenaryForEachVoxelBody_2Const &o, split s)
  :
    irtkForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const irtkImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody_2Const *>(this)->_VoxelFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 1D image region
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody_2Const *>(this)->_VoxelFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody_2Const *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody_2Const *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 2 const, 5 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7,
          class VoxelFunc, class OutsideFunc = irtkNaryVoxelFunction::NOP,
          class Domain = irtkImageDomain::Foreground>
struct irtkSeptenaryForEachVoxelIfBody_2Const : public irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  const irtkGenericImage<T1> &im1;
  const irtkGenericImage<T2> &im2;
        irtkGenericImage<T3> &im3;
        irtkGenericImage<T4> &im4;
        irtkGenericImage<T5> &im5;
        irtkGenericImage<T6> &im6;
        irtkGenericImage<T7> &im7;

  /// Constructor
  irtkSeptenaryForEachVoxelIfBody_2Const(const irtkGenericImage<T1> &im1,
                                         const irtkGenericImage<T2> &im2,
                                               irtkGenericImage<T3> &im3,
                                               irtkGenericImage<T4> &im4,
                                               irtkGenericImage<T5> &im5,
                                               irtkGenericImage<T6> &im6,
                                               irtkGenericImage<T7> &im7,
                                         VoxelFunc &vf, OutsideFunc &of)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  irtkSeptenaryForEachVoxelIfBody_2Const(const irtkSeptenaryForEachVoxelIfBody_2Const &o)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  irtkSeptenaryForEachVoxelIfBody_2Const(irtkSeptenaryForEachVoxelIfBody_2Const &o, split s)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const irtkImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      if (Domain::IsInside(im7, i, j, k, l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody_2Const *>(this)->_VoxelFunc  (i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody_2Const *>(this)->_OutsideFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 1D image region
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, idx, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody_2Const *>(this)->_VoxelFunc  (im7, idx, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody_2Const *>(this)->_OutsideFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, this->_k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody_2Const *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody_2Const *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody_2Const *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody_2Const *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
// ForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// =============================================================================
// 1 const, 6 non-const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 1 const, 6 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
struct irtkSeptenaryForEachVoxelBody_1Const : public irtkForEachVoxelBody<VoxelFunc>
{
  const irtkGenericImage<T1> &im1;
        irtkGenericImage<T2> &im2;
        irtkGenericImage<T3> &im3;
        irtkGenericImage<T4> &im4;
        irtkGenericImage<T5> &im5;
        irtkGenericImage<T6> &im6;
        irtkGenericImage<T7> &im7;

  /// Constructor
  irtkSeptenaryForEachVoxelBody_1Const(const irtkGenericImage<T1> &im1,
                                             irtkGenericImage<T2> &im2,
                                             irtkGenericImage<T3> &im3,
                                             irtkGenericImage<T4> &im4,
                                             irtkGenericImage<T5> &im5,
                                             irtkGenericImage<T6> &im6,
                                             irtkGenericImage<T7> &im7,
                                       VoxelFunc &vf)
  :
    irtkForEachVoxelBody<VoxelFunc>(vf), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  irtkSeptenaryForEachVoxelBody_1Const(const irtkSeptenaryForEachVoxelBody_1Const &o)
  :
    irtkForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  irtkSeptenaryForEachVoxelBody_1Const(irtkSeptenaryForEachVoxelBody_1Const &o, split s)
  :
    irtkForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const irtkImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody_1Const *>(this)->_VoxelFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 1D image region
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody_1Const *>(this)->_VoxelFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody_1Const *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody_1Const *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 1 const, 6 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7,
          class VoxelFunc, class OutsideFunc = irtkNaryVoxelFunction::NOP,
          class Domain = irtkImageDomain::Foreground>
struct irtkSeptenaryForEachVoxelIfBody_1Const : public irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  const irtkGenericImage<T1> &im1;
        irtkGenericImage<T2> &im2;
        irtkGenericImage<T3> &im3;
        irtkGenericImage<T4> &im4;
        irtkGenericImage<T5> &im5;
        irtkGenericImage<T6> &im6;
        irtkGenericImage<T7> &im7;

  /// Constructor
  irtkSeptenaryForEachVoxelIfBody_1Const(const irtkGenericImage<T1> &im1,
                                               irtkGenericImage<T2> &im2,
                                               irtkGenericImage<T3> &im3,
                                               irtkGenericImage<T4> &im4,
                                               irtkGenericImage<T5> &im5,
                                               irtkGenericImage<T6> &im6,
                                               irtkGenericImage<T7> &im7,
                                         VoxelFunc &vf, OutsideFunc &of)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  irtkSeptenaryForEachVoxelIfBody_1Const(const irtkSeptenaryForEachVoxelIfBody_1Const &o)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  irtkSeptenaryForEachVoxelIfBody_1Const(irtkSeptenaryForEachVoxelIfBody_1Const &o, split s)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const irtkImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      if (Domain::IsInside(im7, i, j, k, l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody_1Const *>(this)->_VoxelFunc  (i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody_1Const *>(this)->_OutsideFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 1D image region
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, idx, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody_1Const *>(this)->_VoxelFunc  (im7, idx, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody_1Const *>(this)->_OutsideFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, this->_k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody_1Const *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody_1Const *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
          T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody_1Const *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody_1Const *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
// ForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// =============================================================================
// 7 non-const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 7 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
struct irtkSeptenaryForEachVoxelBody : public irtkForEachVoxelBody<VoxelFunc>
{
  irtkGenericImage<T1> &im1;
  irtkGenericImage<T2> &im2;
  irtkGenericImage<T3> &im3;
  irtkGenericImage<T4> &im4;
  irtkGenericImage<T5> &im5;
  irtkGenericImage<T6> &im6;
  irtkGenericImage<T7> &im7;

  /// Constructor
  irtkSeptenaryForEachVoxelBody(irtkGenericImage<T1> &im1,
                                irtkGenericImage<T2> &im2,
                                irtkGenericImage<T3> &im3,
                                irtkGenericImage<T4> &im4,
                                irtkGenericImage<T5> &im5,
                                irtkGenericImage<T6> &im6,
                                irtkGenericImage<T7> &im7,
                                VoxelFunc &vf)
  :
    irtkForEachVoxelBody<VoxelFunc>(vf), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  irtkSeptenaryForEachVoxelBody(const irtkSeptenaryForEachVoxelBody &o)
  :
    irtkForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  irtkSeptenaryForEachVoxelBody(irtkSeptenaryForEachVoxelBody &o, split s)
  :
    irtkForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const irtkImageAttributes &attr) const
  {
    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
    T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
    T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
    T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
    T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody *>(this)->_VoxelFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 1D image region
  void operator ()(const blocked_range<int> &re) const
  {
    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
    T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
    T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
    T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
    T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody *>(this)->_VoxelFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
    T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
    T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
    T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSeptenaryForEachVoxelBody *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 7 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7,
          class VoxelFunc, class OutsideFunc = irtkNaryVoxelFunction::NOP,
          class Domain = irtkImageDomain::Foreground>
struct irtkSeptenaryForEachVoxelIfBody : public irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  irtkGenericImage<T1> &im1;
  irtkGenericImage<T2> &im2;
  irtkGenericImage<T3> &im3;
  irtkGenericImage<T4> &im4;
  irtkGenericImage<T5> &im5;
  irtkGenericImage<T6> &im6;
  irtkGenericImage<T7> &im7;

  /// Constructor
  irtkSeptenaryForEachVoxelIfBody(irtkGenericImage<T1> &im1,
                                  irtkGenericImage<T2> &im2,
                                  irtkGenericImage<T3> &im3,
                                  irtkGenericImage<T4> &im4,
                                  irtkGenericImage<T5> &im5,
                                  irtkGenericImage<T6> &im6,
                                  irtkGenericImage<T7> &im7,
                                  VoxelFunc &vf, OutsideFunc &of)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6), im7(im7)
  {}

  /// Copy constructor
  irtkSeptenaryForEachVoxelIfBody(const irtkSeptenaryForEachVoxelIfBody &o)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Split constructor
  irtkSeptenaryForEachVoxelIfBody(irtkSeptenaryForEachVoxelIfBody &o, split s)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6), im7(o.im7)
  {}

  /// Process entire image
  void operator ()(const irtkImageAttributes &attr) const
  {
    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();
    T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels();
    T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels();
    T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels();
    T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels();
    T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6, ++p7) {
      if (Domain::IsInside(im7, i, j, k, l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody *>(this)->_VoxelFunc  (i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody *>(this)->_OutsideFunc(i, j, k, l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 1D image region
  void operator ()(const blocked_range<int> &re) const
  {
    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();
    T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels() + re.begin();
    T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels() + re.begin();
    T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels() + re.begin();
    T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels() + re.begin();
    T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, idx, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody *>(this)->_VoxelFunc  (im7, idx, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody *>(this)->_OutsideFunc(im7, idx, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im7.GetX() - (ei - bi);

    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, this->_k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }

  /// Process 3D image region
  void operator ()(const blocked_range3d<int> &re) const
  {
    const int bi = re.cols ().begin();
    const int bj = re.rows ().begin();
    const int bk = re.pages().begin();
    const int ei = re.cols ().end();
    const int ej = re.rows ().end();
    const int ek = re.pages().end();

    const int s1 =  im7.GetX() - (ei - bi);
    const int s2 = (im7.GetY() - (ej - bj)) * im7.GetX();

    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
    T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
    T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);
    T7 *p7 = im7.IsEmpty() ? NULL : im7.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2, p7 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1, p7 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1, p7 +=  1) {
      if (Domain::IsInside(im7, i, j, k, this->_l, p7)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSeptenaryForEachVoxelIfBody *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
      } else const_cast<irtkSeptenaryForEachVoxelIfBody *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6, p7);
    }
  }
};

// -----------------------------------------------------------------------------
// ForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  if (im7->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  if (im7.GetTSize()) {
    ParallelForEachScalar(im1, im2, im3, im4, im5, im6, im7, vf);
  } else {
    irtkSeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkSeptenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc> body(im1, im2, im3, im4, im5, im6, im7, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range<int> re(0, im7->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
    blocked_range<int> re(0, im7->GetNumberOfVoxels() / im7->GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, irtkGenericImage<T7> *im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, *im7, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range<int> re(0, im7.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  if (im7.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
  } else {
    irtkSeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
    blocked_range<int> re(0, im7.GetNumberOfVoxels() / im7.GetT());
    if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
      parallel_reduce(re, body);
      vf.join(body._VoxelFunc);
      of.join(body._OutsideFunc);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  blocked_range3d<int> re(0, attr._z, 0, attr._y, 0, attr._x);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_reduce(re, body);
    } else {
      parallel_reduce(re, body);
    }
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    if (attr._dt) {
      for (body._l = 0; body._l < attr._t; ++body._l) parallel_for(re, body);
    } else {
      parallel_for(re, body);
    }
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSeptenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, T7, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, im7, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, irtkGenericImage<T7> &im7)
{
  if (VoxelFunc::IsReduction()) _irtkforeachseptenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, im7, vf);
}

#endif
