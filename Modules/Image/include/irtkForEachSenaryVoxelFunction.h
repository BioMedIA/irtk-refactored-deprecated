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

#ifndef _irtkForEachSenaryVoxelFunction_H
#define _irtkForEachSenaryVoxelFunction_H

#include <irtkVoxelFunction.h>


inline void _irtkforeachsenaryvoxelfunction_must_not_be_reduction()
{
  cerr << "(Parallel)ForEachVoxel(If): Voxel reductions must be passed by reference!"
          " Pass voxel functor object(s) as last argument(s) instead of first." << endl;
  exit(1);
}


// =============================================================================
// 6 const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 6 const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
struct irtkSenaryForEachVoxelBody_Const : public irtkForEachVoxelBody<VoxelFunc>
{
  const irtkGenericImage<T1> &im1;
  const irtkGenericImage<T2> &im2;
  const irtkGenericImage<T3> &im3;
  const irtkGenericImage<T4> &im4;
  const irtkGenericImage<T5> &im5;
  const irtkGenericImage<T6> &im6;

  /// Constructor
  irtkSenaryForEachVoxelBody_Const(const irtkGenericImage<T1> &im1,
                                   const irtkGenericImage<T2> &im2,
                                   const irtkGenericImage<T3> &im3,
                                   const irtkGenericImage<T4> &im4,
                                   const irtkGenericImage<T5> &im5,
                                   const irtkGenericImage<T6> &im6,
                                   VoxelFunc &vf)
  :
    irtkForEachVoxelBody<VoxelFunc>(vf), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6)
  {}

  /// Copy constructor
  irtkSenaryForEachVoxelBody_Const(const irtkSenaryForEachVoxelBody_Const &o)
  :
    irtkForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6)
  {}

  /// Split constructor
  irtkSenaryForEachVoxelBody_Const(irtkSenaryForEachVoxelBody_Const &o, split s)
  :
    irtkForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6)
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

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSenaryForEachVoxelBody_Const *>(this)->_VoxelFunc(i, j, k, l, p1, p2, p3, p4, p5, p6);
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

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSenaryForEachVoxelBody_Const *>(this)->_VoxelFunc(im6, idx, p1, p2, p3, p4, p5, p6);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im6.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSenaryForEachVoxelBody_Const *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6);
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

    const int s1 =  im6.GetX() - (ei - bi);
    const int s2 = (im6.GetY() - (ej - bj)) * im6.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSenaryForEachVoxelBody_Const *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 6 const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6,
          class VoxelFunc, class OutsideFunc = irtkNaryVoxelFunction::NOP,
          class Domain = irtkImageDomain::Foreground>
struct irtkSenaryForEachVoxelIfBody_Const : public irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  const irtkGenericImage<T1> &im1;
  const irtkGenericImage<T2> &im2;
  const irtkGenericImage<T3> &im3;
  const irtkGenericImage<T4> &im4;
  const irtkGenericImage<T5> &im5;
  const irtkGenericImage<T6> &im6;

  /// Constructor
  irtkSenaryForEachVoxelIfBody_Const(const irtkGenericImage<T1> &im1,
                                     const irtkGenericImage<T2> &im2,
                                     const irtkGenericImage<T3> &im3,
                                     const irtkGenericImage<T4> &im4,
                                     const irtkGenericImage<T5> &im5,
                                     const irtkGenericImage<T6> &im6,
                                     VoxelFunc &vf, OutsideFunc &of)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6)
  {}

  /// Copy constructor
  irtkSenaryForEachVoxelIfBody_Const(const irtkSenaryForEachVoxelIfBody_Const &o)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6)
  {}

  /// Split constructor
  irtkSenaryForEachVoxelIfBody_Const(irtkSenaryForEachVoxelIfBody_Const &o, split s)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6)
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

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6) {
      if (Domain::IsInside(im6, i, j, k, l, p6)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSenaryForEachVoxelIfBody_Const *>(this)->_VoxelFunc  (i, j, k, l, p1, p2, p3, p4, p5, p6);
      } else const_cast<irtkSenaryForEachVoxelIfBody_Const *>(this)->_OutsideFunc(i, j, k, l, p1, p2, p3, p4, p5, p6);
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

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      if (Domain::IsInside(im6, idx, p6)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSenaryForEachVoxelIfBody_Const *>(this)->_VoxelFunc  (im6, idx, p1, p2, p3, p4, p5, p6);
      } else const_cast<irtkSenaryForEachVoxelIfBody_Const *>(this)->_OutsideFunc(im6, idx, p1, p2, p3, p4, p5, p6);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im6.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      if (Domain::IsInside(im6, i, j, this->_k, this->_l, p6)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSenaryForEachVoxelIfBody_Const *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6);
      } else const_cast<irtkSenaryForEachVoxelIfBody_Const *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6);
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

    const int s1 =  im6.GetX() - (ei - bi);
    const int s2 = (im6.GetY() - (ej - bj)) * im6.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      if (Domain::IsInside(im6, i, j, k, this->_l, p6)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSenaryForEachVoxelIfBody_Const *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2, p3, p4, p5, p6);
      } else const_cast<irtkSenaryForEachVoxelIfBody_Const *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6);
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
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  blocked_range<int> re(0, im6->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  if (im6->GetTSize()) {
    ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  } else {
    irtkSenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
    blocked_range<int> re(0, im6->GetNumberOfVoxels() / im6->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  blocked_range<int> re(0, im6.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  if (im6.GetTSize()) {
    ForEachScalar(im1, im2, im3, im4, im5, im6, vf);
  } else {
    irtkSenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
    blocked_range<int> re(0, im6.GetNumberOfVoxels() / im6.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  blocked_range<int> re(0, im6->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  if (im6->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  } else {
    irtkSenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
    blocked_range<int> re(0, im6->GetNumberOfVoxels() / im6->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  blocked_range<int> re(0, im6.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  if (im6.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
  } else {
    irtkSenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
    blocked_range<int> re(0, im6.GetNumberOfVoxels() / im6.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  blocked_range<int> re(0, im6->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  if (im6->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  } else {
    irtkSenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
    blocked_range<int> re(0, im6->GetNumberOfVoxels() / im6->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
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
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  blocked_range<int> re(0, im6.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  if (im6.GetTSize()) {
    ParallelForEachScalar(im1, im2, im3, im4, im5, im6, vf);
  } else {
    irtkSenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
    blocked_range<int> re(0, im6.GetNumberOfVoxels() / im6.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
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
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  blocked_range<int> re(0, im6->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  if (im6->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  } else {
    irtkSenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
    blocked_range<int> re(0, im6->GetNumberOfVoxels() / im6->GetT());
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
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
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
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, const irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  blocked_range<int> re(0, im6.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  if (im6.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
  } else {
    irtkSenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
    blocked_range<int> re(0, im6.GetNumberOfVoxels() / im6.GetT());
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
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
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
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, const irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// =============================================================================
// 5 const, 1 non-const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 5 const, 1 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
struct irtkSenaryForEachVoxelBody_5Const : public irtkForEachVoxelBody<VoxelFunc>
{
  const irtkGenericImage<T1> &im1;
  const irtkGenericImage<T2> &im2;
  const irtkGenericImage<T3> &im3;
  const irtkGenericImage<T4> &im4;
  const irtkGenericImage<T5> &im5;
        irtkGenericImage<T6> &im6;

  /// Constructor
  irtkSenaryForEachVoxelBody_5Const(const irtkGenericImage<T1> &im1,
                                    const irtkGenericImage<T2> &im2,
                                    const irtkGenericImage<T3> &im3,
                                    const irtkGenericImage<T4> &im4,
                                    const irtkGenericImage<T5> &im5,
                                          irtkGenericImage<T6> &im6,
                                    VoxelFunc &vf)
  :
    irtkForEachVoxelBody<VoxelFunc>(vf), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6)
  {}

  /// Copy constructor
  irtkSenaryForEachVoxelBody_5Const(const irtkSenaryForEachVoxelBody_5Const &o)
  :
    irtkForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6)
  {}

  /// Split constructor
  irtkSenaryForEachVoxelBody_5Const(irtkSenaryForEachVoxelBody_5Const &o, split s)
  :
    irtkForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6)
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

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSenaryForEachVoxelBody_5Const *>(this)->_VoxelFunc(i, j, k, l, p1, p2, p3, p4, p5, p6);
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

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSenaryForEachVoxelBody_5Const *>(this)->_VoxelFunc(im6, idx, p1, p2, p3, p4, p5, p6);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im6.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSenaryForEachVoxelBody_5Const *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6);
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

    const int s1 =  im6.GetX() - (ei - bi);
    const int s2 = (im6.GetY() - (ej - bj)) * im6.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSenaryForEachVoxelBody_5Const *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 5 const, 1 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6,
          class VoxelFunc, class OutsideFunc = irtkNaryVoxelFunction::NOP,
          class Domain = irtkImageDomain::Foreground>
struct irtkSenaryForEachVoxelIfBody_5Const : public irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  const irtkGenericImage<T1> &im1;
  const irtkGenericImage<T2> &im2;
  const irtkGenericImage<T3> &im3;
  const irtkGenericImage<T4> &im4;
  const irtkGenericImage<T5> &im5;
        irtkGenericImage<T6> &im6;

  /// Constructor
  irtkSenaryForEachVoxelIfBody_5Const(const irtkGenericImage<T1> &im1,
                                      const irtkGenericImage<T2> &im2,
                                      const irtkGenericImage<T3> &im3,
                                      const irtkGenericImage<T4> &im4,
                                      const irtkGenericImage<T5> &im5,
                                            irtkGenericImage<T6> &im6,
                                      VoxelFunc &vf, OutsideFunc &of)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6)
  {}

  /// Copy constructor
  irtkSenaryForEachVoxelIfBody_5Const(const irtkSenaryForEachVoxelIfBody_5Const &o)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6)
  {}

  /// Split constructor
  irtkSenaryForEachVoxelIfBody_5Const(irtkSenaryForEachVoxelIfBody_5Const &o, split s)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6)
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

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6) {
      if (Domain::IsInside(im6, i, j, k, l, p6)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSenaryForEachVoxelIfBody_5Const *>(this)->_VoxelFunc  (i, j, k, l, p1, p2, p3, p4, p5, p6);
      } else const_cast<irtkSenaryForEachVoxelIfBody_5Const *>(this)->_OutsideFunc(i, j, k, l, p1, p2, p3, p4, p5, p6);
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

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      if (Domain::IsInside(im6, idx, p6)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSenaryForEachVoxelIfBody_5Const *>(this)->_VoxelFunc  (im6, idx, p1, p2, p3, p4, p5, p6);
      } else const_cast<irtkSenaryForEachVoxelIfBody_5Const *>(this)->_OutsideFunc(im6, idx, p1, p2, p3, p4, p5, p6);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im6.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      if (Domain::IsInside(im6, i, j, this->_k, this->_l, p6)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSenaryForEachVoxelIfBody_5Const *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6);
      } else const_cast<irtkSenaryForEachVoxelIfBody_5Const *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6);
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

    const int s1 =  im6.GetX() - (ei - bi);
    const int s2 = (im6.GetY() - (ej - bj)) * im6.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      if (Domain::IsInside(im6, i, j, k, this->_l, p6)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSenaryForEachVoxelIfBody_5Const *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2, p3, p4, p5, p6);
      } else const_cast<irtkSenaryForEachVoxelIfBody_5Const *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6);
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
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  blocked_range<int> re(0, im6->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  if (im6->GetTSize()) {
    ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  } else {
    irtkSenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
    blocked_range<int> re(0, im6->GetNumberOfVoxels() / im6->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  blocked_range<int> re(0, im6.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  if (im6.GetTSize()) {
    ForEachScalar(im1, im2, im3, im4, im5, im6, vf);
  } else {
    irtkSenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
    blocked_range<int> re(0, im6.GetNumberOfVoxels() / im6.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  blocked_range<int> re(0, im6->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  if (im6->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  } else {
    irtkSenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
    blocked_range<int> re(0, im6->GetNumberOfVoxels() / im6->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  blocked_range<int> re(0, im6.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  if (im6.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
  } else {
    irtkSenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
    blocked_range<int> re(0, im6.GetNumberOfVoxels() / im6.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  blocked_range<int> re(0, im6->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  if (im6->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  } else {
    irtkSenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
    blocked_range<int> re(0, im6->GetNumberOfVoxels() / im6->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
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
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  blocked_range<int> re(0, im6.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  if (im6.GetTSize()) {
    ParallelForEachScalar(im1, im2, im3, im4, im5, im6, vf);
  } else {
    irtkSenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
    blocked_range<int> re(0, im6.GetNumberOfVoxels() / im6.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
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
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  blocked_range<int> re(0, im6->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  if (im6->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  } else {
    irtkSenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
    blocked_range<int> re(0, im6->GetNumberOfVoxels() / im6->GetT());
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
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
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
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, const irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  blocked_range<int> re(0, im6.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  if (im6.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
  } else {
    irtkSenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
    blocked_range<int> re(0, im6.GetNumberOfVoxels() / im6.GetT());
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
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
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
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_5Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, const irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// =============================================================================
// 4 const, 2 non-const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 4 const, 2 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
struct irtkSenaryForEachVoxelBody_4Const : public irtkForEachVoxelBody<VoxelFunc>
{
  const irtkGenericImage<T1> &im1;
  const irtkGenericImage<T2> &im2;
  const irtkGenericImage<T3> &im3;
  const irtkGenericImage<T4> &im4;
        irtkGenericImage<T5> &im5;
        irtkGenericImage<T6> &im6;

  /// Constructor
  irtkSenaryForEachVoxelBody_4Const(const irtkGenericImage<T1> &im1,
                                    const irtkGenericImage<T2> &im2,
                                    const irtkGenericImage<T3> &im3,
                                    const irtkGenericImage<T4> &im4,
                                          irtkGenericImage<T5> &im5,
                                          irtkGenericImage<T6> &im6,
                                    VoxelFunc &vf)
  :
    irtkForEachVoxelBody<VoxelFunc>(vf), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6)
  {}

  /// Copy constructor
  irtkSenaryForEachVoxelBody_4Const(const irtkSenaryForEachVoxelBody_4Const &o)
  :
    irtkForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6)
  {}

  /// Split constructor
  irtkSenaryForEachVoxelBody_4Const(irtkSenaryForEachVoxelBody_4Const &o, split s)
  :
    irtkForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6)
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

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSenaryForEachVoxelBody_4Const *>(this)->_VoxelFunc(i, j, k, l, p1, p2, p3, p4, p5, p6);
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

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSenaryForEachVoxelBody_4Const *>(this)->_VoxelFunc(im6, idx, p1, p2, p3, p4, p5, p6);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im6.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSenaryForEachVoxelBody_4Const *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6);
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

    const int s1 =  im6.GetX() - (ei - bi);
    const int s2 = (im6.GetY() - (ej - bj)) * im6.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSenaryForEachVoxelBody_4Const *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 4 const, 2 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6,
          class VoxelFunc, class OutsideFunc = irtkNaryVoxelFunction::NOP,
          class Domain = irtkImageDomain::Foreground>
struct irtkSenaryForEachVoxelIfBody_4Const : public irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  const irtkGenericImage<T1> &im1;
  const irtkGenericImage<T2> &im2;
  const irtkGenericImage<T3> &im3;
  const irtkGenericImage<T4> &im4;
        irtkGenericImage<T5> &im5;
        irtkGenericImage<T6> &im6;

  /// Constructor
  irtkSenaryForEachVoxelIfBody_4Const(const irtkGenericImage<T1> &im1,
                                      const irtkGenericImage<T2> &im2,
                                      const irtkGenericImage<T3> &im3,
                                      const irtkGenericImage<T4> &im4,
                                            irtkGenericImage<T5> &im5,
                                            irtkGenericImage<T6> &im6,
                                      VoxelFunc &vf, OutsideFunc &of)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6)
  {}

  /// Copy constructor
  irtkSenaryForEachVoxelIfBody_4Const(const irtkSenaryForEachVoxelIfBody_4Const &o)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6)
  {}

  /// Split constructor
  irtkSenaryForEachVoxelIfBody_4Const(irtkSenaryForEachVoxelIfBody_4Const &o, split s)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6)
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

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6) {
      if (Domain::IsInside(im6, i, j, k, l, p6)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSenaryForEachVoxelIfBody_4Const *>(this)->_VoxelFunc  (i, j, k, l, p1, p2, p3, p4, p5, p6);
      } else const_cast<irtkSenaryForEachVoxelIfBody_4Const *>(this)->_OutsideFunc(i, j, k, l, p1, p2, p3, p4, p5, p6);
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

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      if (Domain::IsInside(im6, idx, p6)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSenaryForEachVoxelIfBody_4Const *>(this)->_VoxelFunc  (im6, idx, p1, p2, p3, p4, p5, p6);
      } else const_cast<irtkSenaryForEachVoxelIfBody_4Const *>(this)->_OutsideFunc(im6, idx, p1, p2, p3, p4, p5, p6);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im6.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      if (Domain::IsInside(im6, i, j, this->_k, this->_l, p6)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSenaryForEachVoxelIfBody_4Const *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6);
      } else const_cast<irtkSenaryForEachVoxelIfBody_4Const *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6);
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

    const int s1 =  im6.GetX() - (ei - bi);
    const int s2 = (im6.GetY() - (ej - bj)) * im6.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      if (Domain::IsInside(im6, i, j, k, this->_l, p6)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSenaryForEachVoxelIfBody_4Const *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2, p3, p4, p5, p6);
      } else const_cast<irtkSenaryForEachVoxelIfBody_4Const *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6);
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
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  blocked_range<int> re(0, im6->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  if (im6->GetTSize()) {
    ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  } else {
    irtkSenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
    blocked_range<int> re(0, im6->GetNumberOfVoxels() / im6->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  blocked_range<int> re(0, im6.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  if (im6.GetTSize()) {
    ForEachScalar(im1, im2, im3, im4, im5, im6, vf);
  } else {
    irtkSenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
    blocked_range<int> re(0, im6.GetNumberOfVoxels() / im6.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  blocked_range<int> re(0, im6->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  if (im6->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  } else {
    irtkSenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
    blocked_range<int> re(0, im6->GetNumberOfVoxels() / im6->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  blocked_range<int> re(0, im6.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  if (im6.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
  } else {
    irtkSenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
    blocked_range<int> re(0, im6.GetNumberOfVoxels() / im6.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  blocked_range<int> re(0, im6->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  if (im6->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  } else {
    irtkSenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
    blocked_range<int> re(0, im6->GetNumberOfVoxels() / im6->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
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
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  blocked_range<int> re(0, im6.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  if (im6.GetTSize()) {
    ParallelForEachScalar(im1, im2, im3, im4, im5, im6, vf);
  } else {
    irtkSenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
    blocked_range<int> re(0, im6.GetNumberOfVoxels() / im6.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
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
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  blocked_range<int> re(0, im6->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  if (im6->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  } else {
    irtkSenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
    blocked_range<int> re(0, im6->GetNumberOfVoxels() / im6->GetT());
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
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
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
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, const irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  blocked_range<int> re(0, im6.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  if (im6.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
  } else {
    irtkSenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
    blocked_range<int> re(0, im6.GetNumberOfVoxels() / im6.GetT());
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
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
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
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_4Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, const irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// =============================================================================
// 3 const, 3 non-const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 3 const, 3 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
struct irtkSenaryForEachVoxelBody_3Const : public irtkForEachVoxelBody<VoxelFunc>
{
  const irtkGenericImage<T1> &im1;
  const irtkGenericImage<T2> &im2;
  const irtkGenericImage<T3> &im3;
        irtkGenericImage<T4> &im4;
        irtkGenericImage<T5> &im5;
        irtkGenericImage<T6> &im6;

  /// Constructor
  irtkSenaryForEachVoxelBody_3Const(const irtkGenericImage<T1> &im1,
                                    const irtkGenericImage<T2> &im2,
                                    const irtkGenericImage<T3> &im3,
                                          irtkGenericImage<T4> &im4,
                                          irtkGenericImage<T5> &im5,
                                          irtkGenericImage<T6> &im6,
                                    VoxelFunc &vf)
  :
    irtkForEachVoxelBody<VoxelFunc>(vf), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6)
  {}

  /// Copy constructor
  irtkSenaryForEachVoxelBody_3Const(const irtkSenaryForEachVoxelBody_3Const &o)
  :
    irtkForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6)
  {}

  /// Split constructor
  irtkSenaryForEachVoxelBody_3Const(irtkSenaryForEachVoxelBody_3Const &o, split s)
  :
    irtkForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6)
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

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSenaryForEachVoxelBody_3Const *>(this)->_VoxelFunc(i, j, k, l, p1, p2, p3, p4, p5, p6);
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

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSenaryForEachVoxelBody_3Const *>(this)->_VoxelFunc(im6, idx, p1, p2, p3, p4, p5, p6);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im6.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSenaryForEachVoxelBody_3Const *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6);
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

    const int s1 =  im6.GetX() - (ei - bi);
    const int s2 = (im6.GetY() - (ej - bj)) * im6.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSenaryForEachVoxelBody_3Const *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 3 const, 3 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6,
          class VoxelFunc, class OutsideFunc = irtkNaryVoxelFunction::NOP,
          class Domain = irtkImageDomain::Foreground>
struct irtkSenaryForEachVoxelIfBody_3Const : public irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  const irtkGenericImage<T1> &im1;
  const irtkGenericImage<T2> &im2;
  const irtkGenericImage<T3> &im3;
        irtkGenericImage<T4> &im4;
        irtkGenericImage<T5> &im5;
        irtkGenericImage<T6> &im6;

  /// Constructor
  irtkSenaryForEachVoxelIfBody_3Const(const irtkGenericImage<T1> &im1,
                                      const irtkGenericImage<T2> &im2,
                                      const irtkGenericImage<T3> &im3,
                                            irtkGenericImage<T4> &im4,
                                            irtkGenericImage<T5> &im5,
                                            irtkGenericImage<T6> &im6,
                                      VoxelFunc &vf, OutsideFunc &of)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6)
  {}

  /// Copy constructor
  irtkSenaryForEachVoxelIfBody_3Const(const irtkSenaryForEachVoxelIfBody_3Const &o)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6)
  {}

  /// Split constructor
  irtkSenaryForEachVoxelIfBody_3Const(irtkSenaryForEachVoxelIfBody_3Const &o, split s)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6)
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

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6) {
      if (Domain::IsInside(im6, i, j, k, l, p6)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSenaryForEachVoxelIfBody_3Const *>(this)->_VoxelFunc  (i, j, k, l, p1, p2, p3, p4, p5, p6);
      } else const_cast<irtkSenaryForEachVoxelIfBody_3Const *>(this)->_OutsideFunc(i, j, k, l, p1, p2, p3, p4, p5, p6);
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

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      if (Domain::IsInside(im6, idx, p6)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSenaryForEachVoxelIfBody_3Const *>(this)->_VoxelFunc  (im6, idx, p1, p2, p3, p4, p5, p6);
      } else const_cast<irtkSenaryForEachVoxelIfBody_3Const *>(this)->_OutsideFunc(im6, idx, p1, p2, p3, p4, p5, p6);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im6.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      if (Domain::IsInside(im6, i, j, this->_k, this->_l, p6)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSenaryForEachVoxelIfBody_3Const *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6);
      } else const_cast<irtkSenaryForEachVoxelIfBody_3Const *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6);
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

    const int s1 =  im6.GetX() - (ei - bi);
    const int s2 = (im6.GetY() - (ej - bj)) * im6.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      if (Domain::IsInside(im6, i, j, k, this->_l, p6)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSenaryForEachVoxelIfBody_3Const *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2, p3, p4, p5, p6);
      } else const_cast<irtkSenaryForEachVoxelIfBody_3Const *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6);
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
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  blocked_range<int> re(0, im6->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  if (im6->GetTSize()) {
    ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  } else {
    irtkSenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
    blocked_range<int> re(0, im6->GetNumberOfVoxels() / im6->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  blocked_range<int> re(0, im6.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  if (im6.GetTSize()) {
    ForEachScalar(im1, im2, im3, im4, im5, im6, vf);
  } else {
    irtkSenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
    blocked_range<int> re(0, im6.GetNumberOfVoxels() / im6.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  blocked_range<int> re(0, im6->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  if (im6->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  } else {
    irtkSenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
    blocked_range<int> re(0, im6->GetNumberOfVoxels() / im6->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  blocked_range<int> re(0, im6.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  if (im6.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
  } else {
    irtkSenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
    blocked_range<int> re(0, im6.GetNumberOfVoxels() / im6.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  blocked_range<int> re(0, im6->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  if (im6->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  } else {
    irtkSenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
    blocked_range<int> re(0, im6->GetNumberOfVoxels() / im6->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
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
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  blocked_range<int> re(0, im6.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  if (im6.GetTSize()) {
    ParallelForEachScalar(im1, im2, im3, im4, im5, im6, vf);
  } else {
    irtkSenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
    blocked_range<int> re(0, im6.GetNumberOfVoxels() / im6.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
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
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  blocked_range<int> re(0, im6->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  if (im6->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  } else {
    irtkSenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
    blocked_range<int> re(0, im6->GetNumberOfVoxels() / im6->GetT());
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
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
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
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, const irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  blocked_range<int> re(0, im6.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  if (im6.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
  } else {
    irtkSenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
    blocked_range<int> re(0, im6.GetNumberOfVoxels() / im6.GetT());
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
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
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
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_3Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, const irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// =============================================================================
// 2 const, 4 non-const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 2 const, 4 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
struct irtkSenaryForEachVoxelBody_2Const : public irtkForEachVoxelBody<VoxelFunc>
{
  const irtkGenericImage<T1> &im1;
  const irtkGenericImage<T2> &im2;
        irtkGenericImage<T3> &im3;
        irtkGenericImage<T4> &im4;
        irtkGenericImage<T5> &im5;
        irtkGenericImage<T6> &im6;

  /// Constructor
  irtkSenaryForEachVoxelBody_2Const(const irtkGenericImage<T1> &im1,
                                    const irtkGenericImage<T2> &im2,
                                          irtkGenericImage<T3> &im3,
                                          irtkGenericImage<T4> &im4,
                                          irtkGenericImage<T5> &im5,
                                          irtkGenericImage<T6> &im6,
                                    VoxelFunc &vf)
  :
    irtkForEachVoxelBody<VoxelFunc>(vf), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6)
  {}

  /// Copy constructor
  irtkSenaryForEachVoxelBody_2Const(const irtkSenaryForEachVoxelBody_2Const &o)
  :
    irtkForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6)
  {}

  /// Split constructor
  irtkSenaryForEachVoxelBody_2Const(irtkSenaryForEachVoxelBody_2Const &o, split s)
  :
    irtkForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6)
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

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSenaryForEachVoxelBody_2Const *>(this)->_VoxelFunc(i, j, k, l, p1, p2, p3, p4, p5, p6);
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

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSenaryForEachVoxelBody_2Const *>(this)->_VoxelFunc(im6, idx, p1, p2, p3, p4, p5, p6);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im6.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSenaryForEachVoxelBody_2Const *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6);
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

    const int s1 =  im6.GetX() - (ei - bi);
    const int s2 = (im6.GetY() - (ej - bj)) * im6.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSenaryForEachVoxelBody_2Const *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 2 const, 4 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6,
          class VoxelFunc, class OutsideFunc = irtkNaryVoxelFunction::NOP,
          class Domain = irtkImageDomain::Foreground>
struct irtkSenaryForEachVoxelIfBody_2Const : public irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  const irtkGenericImage<T1> &im1;
  const irtkGenericImage<T2> &im2;
        irtkGenericImage<T3> &im3;
        irtkGenericImage<T4> &im4;
        irtkGenericImage<T5> &im5;
        irtkGenericImage<T6> &im6;

  /// Constructor
  irtkSenaryForEachVoxelIfBody_2Const(const irtkGenericImage<T1> &im1,
                                      const irtkGenericImage<T2> &im2,
                                            irtkGenericImage<T3> &im3,
                                            irtkGenericImage<T4> &im4,
                                            irtkGenericImage<T5> &im5,
                                            irtkGenericImage<T6> &im6,
                                      VoxelFunc &vf, OutsideFunc &of)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6)
  {}

  /// Copy constructor
  irtkSenaryForEachVoxelIfBody_2Const(const irtkSenaryForEachVoxelIfBody_2Const &o)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6)
  {}

  /// Split constructor
  irtkSenaryForEachVoxelIfBody_2Const(irtkSenaryForEachVoxelIfBody_2Const &o, split s)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6)
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

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6) {
      if (Domain::IsInside(im6, i, j, k, l, p6)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSenaryForEachVoxelIfBody_2Const *>(this)->_VoxelFunc  (i, j, k, l, p1, p2, p3, p4, p5, p6);
      } else const_cast<irtkSenaryForEachVoxelIfBody_2Const *>(this)->_OutsideFunc(i, j, k, l, p1, p2, p3, p4, p5, p6);
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

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      if (Domain::IsInside(im6, idx, p6)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSenaryForEachVoxelIfBody_2Const *>(this)->_VoxelFunc  (im6, idx, p1, p2, p3, p4, p5, p6);
      } else const_cast<irtkSenaryForEachVoxelIfBody_2Const *>(this)->_OutsideFunc(im6, idx, p1, p2, p3, p4, p5, p6);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im6.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      if (Domain::IsInside(im6, i, j, this->_k, this->_l, p6)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSenaryForEachVoxelIfBody_2Const *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6);
      } else const_cast<irtkSenaryForEachVoxelIfBody_2Const *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6);
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

    const int s1 =  im6.GetX() - (ei - bi);
    const int s2 = (im6.GetY() - (ej - bj)) * im6.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      if (Domain::IsInside(im6, i, j, k, this->_l, p6)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSenaryForEachVoxelIfBody_2Const *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2, p3, p4, p5, p6);
      } else const_cast<irtkSenaryForEachVoxelIfBody_2Const *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6);
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
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  blocked_range<int> re(0, im6->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  if (im6->GetTSize()) {
    ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  } else {
    irtkSenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
    blocked_range<int> re(0, im6->GetNumberOfVoxels() / im6->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  blocked_range<int> re(0, im6.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  if (im6.GetTSize()) {
    ForEachScalar(im1, im2, im3, im4, im5, im6, vf);
  } else {
    irtkSenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
    blocked_range<int> re(0, im6.GetNumberOfVoxels() / im6.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  blocked_range<int> re(0, im6->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  if (im6->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  } else {
    irtkSenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
    blocked_range<int> re(0, im6->GetNumberOfVoxels() / im6->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  blocked_range<int> re(0, im6.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  if (im6.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
  } else {
    irtkSenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
    blocked_range<int> re(0, im6.GetNumberOfVoxels() / im6.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  blocked_range<int> re(0, im6->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  if (im6->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  } else {
    irtkSenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
    blocked_range<int> re(0, im6->GetNumberOfVoxels() / im6->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
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
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  blocked_range<int> re(0, im6.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  if (im6.GetTSize()) {
    ParallelForEachScalar(im1, im2, im3, im4, im5, im6, vf);
  } else {
    irtkSenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
    blocked_range<int> re(0, im6.GetNumberOfVoxels() / im6.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
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
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  blocked_range<int> re(0, im6->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  if (im6->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  } else {
    irtkSenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
    blocked_range<int> re(0, im6->GetNumberOfVoxels() / im6->GetT());
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
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
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
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  blocked_range<int> re(0, im6.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  if (im6.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
  } else {
    irtkSenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
    blocked_range<int> re(0, im6.GetNumberOfVoxels() / im6.GetT());
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
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
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
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_2Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// =============================================================================
// 1 const, 5 non-const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 1 const, 5 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
struct irtkSenaryForEachVoxelBody_1Const : public irtkForEachVoxelBody<VoxelFunc>
{
  const irtkGenericImage<T1> &im1;
        irtkGenericImage<T2> &im2;
        irtkGenericImage<T3> &im3;
        irtkGenericImage<T4> &im4;
        irtkGenericImage<T5> &im5;
        irtkGenericImage<T6> &im6;

  /// Constructor
  irtkSenaryForEachVoxelBody_1Const(const irtkGenericImage<T1> &im1,
                                          irtkGenericImage<T2> &im2,
                                          irtkGenericImage<T3> &im3,
                                          irtkGenericImage<T4> &im4,
                                          irtkGenericImage<T5> &im5,
                                          irtkGenericImage<T6> &im6,
                                    VoxelFunc &vf)
  :
    irtkForEachVoxelBody<VoxelFunc>(vf), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6)
  {}

  /// Copy constructor
  irtkSenaryForEachVoxelBody_1Const(const irtkSenaryForEachVoxelBody_1Const &o)
  :
    irtkForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6)
  {}

  /// Split constructor
  irtkSenaryForEachVoxelBody_1Const(irtkSenaryForEachVoxelBody_1Const &o, split s)
  :
    irtkForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6)
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

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSenaryForEachVoxelBody_1Const *>(this)->_VoxelFunc(i, j, k, l, p1, p2, p3, p4, p5, p6);
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

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSenaryForEachVoxelBody_1Const *>(this)->_VoxelFunc(im6, idx, p1, p2, p3, p4, p5, p6);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im6.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSenaryForEachVoxelBody_1Const *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6);
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

    const int s1 =  im6.GetX() - (ei - bi);
    const int s2 = (im6.GetY() - (ej - bj)) * im6.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSenaryForEachVoxelBody_1Const *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 1 const, 5 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6,
          class VoxelFunc, class OutsideFunc = irtkNaryVoxelFunction::NOP,
          class Domain = irtkImageDomain::Foreground>
struct irtkSenaryForEachVoxelIfBody_1Const : public irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  const irtkGenericImage<T1> &im1;
        irtkGenericImage<T2> &im2;
        irtkGenericImage<T3> &im3;
        irtkGenericImage<T4> &im4;
        irtkGenericImage<T5> &im5;
        irtkGenericImage<T6> &im6;

  /// Constructor
  irtkSenaryForEachVoxelIfBody_1Const(const irtkGenericImage<T1> &im1,
                                            irtkGenericImage<T2> &im2,
                                            irtkGenericImage<T3> &im3,
                                            irtkGenericImage<T4> &im4,
                                            irtkGenericImage<T5> &im5,
                                            irtkGenericImage<T6> &im6,
                                      VoxelFunc &vf, OutsideFunc &of)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6)
  {}

  /// Copy constructor
  irtkSenaryForEachVoxelIfBody_1Const(const irtkSenaryForEachVoxelIfBody_1Const &o)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6)
  {}

  /// Split constructor
  irtkSenaryForEachVoxelIfBody_1Const(irtkSenaryForEachVoxelIfBody_1Const &o, split s)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6)
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

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6) {
      if (Domain::IsInside(im6, i, j, k, l, p6)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSenaryForEachVoxelIfBody_1Const *>(this)->_VoxelFunc  (i, j, k, l, p1, p2, p3, p4, p5, p6);
      } else const_cast<irtkSenaryForEachVoxelIfBody_1Const *>(this)->_OutsideFunc(i, j, k, l, p1, p2, p3, p4, p5, p6);
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

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      if (Domain::IsInside(im6, idx, p6)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSenaryForEachVoxelIfBody_1Const *>(this)->_VoxelFunc  (im6, idx, p1, p2, p3, p4, p5, p6);
      } else const_cast<irtkSenaryForEachVoxelIfBody_1Const *>(this)->_OutsideFunc(im6, idx, p1, p2, p3, p4, p5, p6);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im6.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      if (Domain::IsInside(im6, i, j, this->_k, this->_l, p6)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSenaryForEachVoxelIfBody_1Const *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6);
      } else const_cast<irtkSenaryForEachVoxelIfBody_1Const *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6);
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

    const int s1 =  im6.GetX() - (ei - bi);
    const int s2 = (im6.GetY() - (ej - bj)) * im6.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
          T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
          T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
          T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
          T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      if (Domain::IsInside(im6, i, j, k, this->_l, p6)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSenaryForEachVoxelIfBody_1Const *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2, p3, p4, p5, p6);
      } else const_cast<irtkSenaryForEachVoxelIfBody_1Const *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6);
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
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  blocked_range<int> re(0, im6->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  if (im6->GetTSize()) {
    ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  } else {
    irtkSenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
    blocked_range<int> re(0, im6->GetNumberOfVoxels() / im6->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  blocked_range<int> re(0, im6.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  if (im6.GetTSize()) {
    ForEachScalar(im1, im2, im3, im4, im5, im6, vf);
  } else {
    irtkSenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
    blocked_range<int> re(0, im6.GetNumberOfVoxels() / im6.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  blocked_range<int> re(0, im6->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  if (im6->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  } else {
    irtkSenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
    blocked_range<int> re(0, im6->GetNumberOfVoxels() / im6->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  blocked_range<int> re(0, im6.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  if (im6.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
  } else {
    irtkSenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
    blocked_range<int> re(0, im6.GetNumberOfVoxels() / im6.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  blocked_range<int> re(0, im6->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  if (im6->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  } else {
    irtkSenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
    blocked_range<int> re(0, im6->GetNumberOfVoxels() / im6->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
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
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  blocked_range<int> re(0, im6.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  if (im6.GetTSize()) {
    ParallelForEachScalar(im1, im2, im3, im4, im5, im6, vf);
  } else {
    irtkSenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
    blocked_range<int> re(0, im6.GetNumberOfVoxels() / im6.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
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
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  blocked_range<int> re(0, im6->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  if (im6->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  } else {
    irtkSenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
    blocked_range<int> re(0, im6->GetNumberOfVoxels() / im6->GetT());
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
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
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
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  blocked_range<int> re(0, im6.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  if (im6.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
  } else {
    irtkSenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
    blocked_range<int> re(0, im6.GetNumberOfVoxels() / im6.GetT());
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
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
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
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody_1Const<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// =============================================================================
// 6 non-const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 6 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
struct irtkSenaryForEachVoxelBody : public irtkForEachVoxelBody<VoxelFunc>
{
  irtkGenericImage<T1> &im1;
  irtkGenericImage<T2> &im2;
  irtkGenericImage<T3> &im3;
  irtkGenericImage<T4> &im4;
  irtkGenericImage<T5> &im5;
  irtkGenericImage<T6> &im6;

  /// Constructor
  irtkSenaryForEachVoxelBody(irtkGenericImage<T1> &im1,
                             irtkGenericImage<T2> &im2,
                             irtkGenericImage<T3> &im3,
                             irtkGenericImage<T4> &im4,
                             irtkGenericImage<T5> &im5,
                             irtkGenericImage<T6> &im6,
                             VoxelFunc &vf)
  :
    irtkForEachVoxelBody<VoxelFunc>(vf), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6)
  {}

  /// Copy constructor
  irtkSenaryForEachVoxelBody(const irtkSenaryForEachVoxelBody &o)
  :
    irtkForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6)
  {}

  /// Split constructor
  irtkSenaryForEachVoxelBody(irtkSenaryForEachVoxelBody &o, split s)
  :
    irtkForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6)
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

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSenaryForEachVoxelBody *>(this)->_VoxelFunc(i, j, k, l, p1, p2, p3, p4, p5, p6);
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

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSenaryForEachVoxelBody *>(this)->_VoxelFunc(im6, idx, p1, p2, p3, p4, p5, p6);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im6.GetX() - (ei - bi);

    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSenaryForEachVoxelBody *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6);
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

    const int s1 =  im6.GetX() - (ei - bi);
    const int s2 = (im6.GetY() - (ej - bj)) * im6.GetX();

    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
    T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
    T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkSenaryForEachVoxelBody *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 6 non-const images
 */
template <class T1, class T2, class T3, class T4, class T5, class T6,
          class VoxelFunc, class OutsideFunc = irtkNaryVoxelFunction::NOP,
          class Domain = irtkImageDomain::Foreground>
struct irtkSenaryForEachVoxelIfBody : public irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  irtkGenericImage<T1> &im1;
  irtkGenericImage<T2> &im2;
  irtkGenericImage<T3> &im3;
  irtkGenericImage<T4> &im4;
  irtkGenericImage<T5> &im5;
  irtkGenericImage<T6> &im6;

  /// Constructor
  irtkSenaryForEachVoxelIfBody(irtkGenericImage<T1> &im1,
                               irtkGenericImage<T2> &im2,
                               irtkGenericImage<T3> &im3,
                               irtkGenericImage<T4> &im4,
                               irtkGenericImage<T5> &im5,
                               irtkGenericImage<T6> &im6,
                               VoxelFunc &vf, OutsideFunc &of)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of), im1(im1), im2(im2), im3(im3), im4(im4), im5(im5), im6(im6)
  {}

  /// Copy constructor
  irtkSenaryForEachVoxelIfBody(const irtkSenaryForEachVoxelIfBody &o)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6)
  {}

  /// Split constructor
  irtkSenaryForEachVoxelIfBody(irtkSenaryForEachVoxelIfBody &o, split s)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2), im3(o.im3), im4(o.im4), im5(o.im5), im6(o.im6)
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

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2, ++p3, ++p4, ++p5, ++p6) {
      if (Domain::IsInside(im6, i, j, k, l, p6)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSenaryForEachVoxelIfBody *>(this)->_VoxelFunc  (i, j, k, l, p1, p2, p3, p4, p5, p6);
      } else const_cast<irtkSenaryForEachVoxelIfBody *>(this)->_OutsideFunc(i, j, k, l, p1, p2, p3, p4, p5, p6);
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

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      if (Domain::IsInside(im6, idx, p6)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSenaryForEachVoxelIfBody *>(this)->_VoxelFunc  (im6, idx, p1, p2, p3, p4, p5, p6);
      } else const_cast<irtkSenaryForEachVoxelIfBody *>(this)->_OutsideFunc(im6, idx, p1, p2, p3, p4, p5, p6);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im6.GetX() - (ei - bi);

    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      if (Domain::IsInside(im6, i, j, this->_k, this->_l, p6)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSenaryForEachVoxelIfBody *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6);
      } else const_cast<irtkSenaryForEachVoxelIfBody *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2, p3, p4, p5, p6);
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

    const int s1 =  im6.GetX() - (ei - bi);
    const int s2 = (im6.GetY() - (ej - bj)) * im6.GetX();

    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);
    T3 *p3 = im3.IsEmpty() ? NULL : im3.GetPointerToVoxels(bi, bj, bk, this->_l);
    T4 *p4 = im4.IsEmpty() ? NULL : im4.GetPointerToVoxels(bi, bj, bk, this->_l);
    T5 *p5 = im5.IsEmpty() ? NULL : im5.GetPointerToVoxels(bi, bj, bk, this->_l);
    T6 *p6 = im6.IsEmpty() ? NULL : im6.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2, p3 += s2, p4 += s2, p5 += s2, p6 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1, p3 += s1, p4 += s1, p5 += s1, p6 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1, p3 +=  1, p4 +=  1, p5 +=  1, p6 +=  1) {
      if (Domain::IsInside(im6, i, j, k, this->_l, p6)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkSenaryForEachVoxelIfBody *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2, p3, p4, p5, p6);
      } else const_cast<irtkSenaryForEachVoxelIfBody *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2, p3, p4, p5, p6);
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
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalar(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  blocked_range<int> re(0, im6->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  if (im6->GetTSize()) {
    ForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  } else {
    irtkSenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
    blocked_range<int> re(0, im6->GetNumberOfVoxels() / im6->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalar(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  blocked_range<int> re(0, im6.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  if (im6.GetTSize()) {
    ForEachScalar(im1, im2, im3, im4, im5, im6, vf);
  } else {
    irtkSenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
    blocked_range<int> re(0, im6.GetNumberOfVoxels() / im6.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  blocked_range<int> re(0, im6->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalarIf(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  if (im6->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  } else {
    irtkSenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
    blocked_range<int> re(0, im6->GetNumberOfVoxels() / im6->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  blocked_range<int> re(0, im6.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalarIf(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  if (im6.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
  } else {
    irtkSenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
    blocked_range<int> re(0, im6.GetNumberOfVoxels() / im6.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalar(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  blocked_range<int> re(0, im6->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  if (im6->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  } else {
    irtkSenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
    blocked_range<int> re(0, im6->GetNumberOfVoxels() / im6->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
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
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, VoxelFunc> body(*im1, *im2, *im3, *im4, *im5, *im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalar(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  blocked_range<int> re(0, im6.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  if (im6.GetTSize()) {
    ParallelForEachScalar(im1, im2, im3, im4, im5, im6, vf);
  } else {
    irtkSenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
    blocked_range<int> re(0, im6.GetNumberOfVoxels() / im6.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
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
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkSenaryForEachVoxelBody<T1, T2, T3, T4, T5, T6, VoxelFunc> body(im1, im2, im3, im4, im5, im6, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  blocked_range<int> re(0, im6->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalarIf(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  if (im6->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  } else {
    irtkSenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
    blocked_range<int> re(0, im6->GetNumberOfVoxels() / im6->GetT());
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
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
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
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, *im3, *im4, *im5, *im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, irtkGenericImage<T3> *im3, irtkGenericImage<T4> *im4, irtkGenericImage<T5> *im5, irtkGenericImage<T6> *im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, *im3, *im4, *im5, *im6, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  blocked_range<int> re(0, im6.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalarIf(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  if (im6.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
  } else {
    irtkSenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
    blocked_range<int> re(0, im6.GetNumberOfVoxels() / im6.GetT());
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
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
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
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf, OutsideFunc &of)
{
  irtkSenaryForEachVoxelIfBody<T1, T2, T3, T4, T5, T6, VoxelFunc, OutsideFunc, Domain> body(im1, im2, im3, im4, im5, im6, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class T3, class T4, class T5, class T6, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, irtkGenericImage<T3> &im3, irtkGenericImage<T4> &im4, irtkGenericImage<T5> &im5, irtkGenericImage<T6> &im6)
{
  if (VoxelFunc::IsReduction()) _irtkforeachsenaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, im3, im4, im5, im6, vf);
}

#endif
