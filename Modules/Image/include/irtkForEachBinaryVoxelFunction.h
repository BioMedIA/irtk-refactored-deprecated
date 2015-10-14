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

#ifndef _irtkForEachBinaryVoxelFunction_H
#define _irtkForEachBinaryVoxelFunction_H

#include <irtkVoxelFunction.h>


inline void _irtkforeachbinaryvoxelfunction_must_not_be_reduction()
{
  cerr << "(Parallel)ForEachVoxel(If): Voxel reductions must be passed by reference!"
          " Pass voxel functor object(s) as last argument(s) instead of first." << endl;
  exit(1);
}


// =============================================================================
// 2 const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 2 const images
 */
template <class T1, class T2, class VoxelFunc>
struct irtkBinaryForEachVoxelBody_Const : public irtkForEachVoxelBody<VoxelFunc>
{
  const irtkGenericImage<T1> &im1;
  const irtkGenericImage<T2> &im2;

  /// Constructor
  irtkBinaryForEachVoxelBody_Const(const irtkGenericImage<T1> &im1,
                                   const irtkGenericImage<T2> &im2,
                                   VoxelFunc &vf)
  :
    irtkForEachVoxelBody<VoxelFunc>(vf), im1(im1), im2(im2)
  {}

  /// Copy constructor
  irtkBinaryForEachVoxelBody_Const(const irtkBinaryForEachVoxelBody_Const &o)
  :
    irtkForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2)
  {}

  /// Split constructor
  irtkBinaryForEachVoxelBody_Const(irtkBinaryForEachVoxelBody_Const &o, split s)
  :
    irtkForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2)
  {}

  /// Process entire image
  void operator ()(const irtkImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkBinaryForEachVoxelBody_Const *>(this)->_VoxelFunc(i, j, k, l, p1, p2);
    }
  }

  /// Process 1D image region
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkBinaryForEachVoxelBody_Const *>(this)->_VoxelFunc(im2, idx, p1, p2);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im2.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkBinaryForEachVoxelBody_Const *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2);
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

    const int s1 =  im2.GetX() - (ei - bi);
    const int s2 = (im2.GetY() - (ej - bj)) * im2.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkBinaryForEachVoxelBody_Const *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 2 const images
 */
template <class T1, class T2,
          class VoxelFunc, class OutsideFunc = irtkNaryVoxelFunction::NOP,
          class Domain = irtkImageDomain::Foreground>
struct irtkBinaryForEachVoxelIfBody_Const : public irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  const irtkGenericImage<T1> &im1;
  const irtkGenericImage<T2> &im2;

  /// Constructor
  irtkBinaryForEachVoxelIfBody_Const(const irtkGenericImage<T1> &im1,
                                     const irtkGenericImage<T2> &im2,
                                     VoxelFunc &vf, OutsideFunc &of)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of), im1(im1), im2(im2)
  {}

  /// Copy constructor
  irtkBinaryForEachVoxelIfBody_Const(const irtkBinaryForEachVoxelIfBody_Const &o)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2)
  {}

  /// Split constructor
  irtkBinaryForEachVoxelIfBody_Const(irtkBinaryForEachVoxelIfBody_Const &o, split s)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2)
  {}

  /// Process entire image
  void operator ()(const irtkImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2) {
      if (Domain::IsInside(im2, i, j, k, l, p2)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkBinaryForEachVoxelIfBody_Const *>(this)->_VoxelFunc  (i, j, k, l, p1, p2);
      } else const_cast<irtkBinaryForEachVoxelIfBody_Const *>(this)->_OutsideFunc(i, j, k, l, p1, p2);
    }
  }

  /// Process 1D image region
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1) {
      if (Domain::IsInside(im2, idx, p2)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkBinaryForEachVoxelIfBody_Const *>(this)->_VoxelFunc  (im2, idx, p1, p2);
      } else const_cast<irtkBinaryForEachVoxelIfBody_Const *>(this)->_OutsideFunc(im2, idx, p1, p2);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im2.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1) {
      if (Domain::IsInside(im2, i, j, this->_k, this->_l, p2)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkBinaryForEachVoxelIfBody_Const *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2);
      } else const_cast<irtkBinaryForEachVoxelIfBody_Const *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2);
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

    const int s1 =  im2.GetX() - (ei - bi);
    const int s2 = (im2.GetY() - (ej - bj)) * im2.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    const T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1) {
      if (Domain::IsInside(im2, i, j, k, this->_l, p2)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkBinaryForEachVoxelIfBody_Const *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2);
      } else const_cast<irtkBinaryForEachVoxelIfBody_Const *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2);
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
template <class T1, class T2, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_Const<T1, T2, VoxelFunc> body(*im1, *im2, vf);
  blocked_range<int> re(0, im2->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  if (im2->GetTSize()) {
    ForEachScalar(*im1, *im2, vf);
  } else {
    irtkBinaryForEachVoxelBody_Const<T1, T2, VoxelFunc> body(*im1, *im2, vf);
    blocked_range<int> re(0, im2->GetNumberOfVoxels() / im2->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_Const<T1, T2, VoxelFunc> body(*im1, *im2, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_Const<T1, T2, VoxelFunc> body(*im1, *im2, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_Const<T1, T2, VoxelFunc> body(*im1, *im2, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_Const<T1, T2, VoxelFunc> body(*im1, *im2, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_Const<T1, T2, VoxelFunc> body(im1, im2, vf);
  blocked_range<int> re(0, im2.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  if (im2.GetTSize()) {
    ForEachScalar(im1, im2, vf);
  } else {
    irtkBinaryForEachVoxelBody_Const<T1, T2, VoxelFunc> body(im1, im2, vf);
    blocked_range<int> re(0, im2.GetNumberOfVoxels() / im2.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_Const<T1, T2, VoxelFunc> body(im1, im2, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_Const<T1, T2, VoxelFunc> body(im1, im2, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_Const<T1, T2, VoxelFunc> body(im1, im2, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_Const<T1, T2, VoxelFunc> body(im1, im2, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
  blocked_range<int> re(0, im2->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  if (im2->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, vf, of);
  } else {
    irtkBinaryForEachVoxelIfBody_Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
    blocked_range<int> re(0, im2->GetNumberOfVoxels() / im2->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
  blocked_range<int> re(0, im2.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  if (im2.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, vf, of);
  } else {
    irtkBinaryForEachVoxelIfBody_Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
    blocked_range<int> re(0, im2.GetNumberOfVoxels() / im2.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_Const<T1, T2, VoxelFunc> body(*im1, *im2, vf);
  blocked_range<int> re(0, im2->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  if (im2->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, vf);
  } else {
    irtkBinaryForEachVoxelBody_Const<T1, T2, VoxelFunc> body(*im1, *im2, vf);
    blocked_range<int> re(0, im2->GetNumberOfVoxels() / im2->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_Const<T1, T2, VoxelFunc> body(*im1, *im2, vf);
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
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_Const<T1, T2, VoxelFunc> body(*im1, *im2, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_Const<T1, T2, VoxelFunc> body(*im1, *im2, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_Const<T1, T2, VoxelFunc> body(*im1, *im2, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_Const<T1, T2, VoxelFunc> body(im1, im2, vf);
  blocked_range<int> re(0, im2.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  if (im2.GetTSize()) {
    ParallelForEachScalar(im1, im2, vf);
  } else {
    irtkBinaryForEachVoxelBody_Const<T1, T2, VoxelFunc> body(im1, im2, vf);
    blocked_range<int> re(0, im2.GetNumberOfVoxels() / im2.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_Const<T1, T2, VoxelFunc> body(im1, im2, vf);
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
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_Const<T1, T2, VoxelFunc> body(im1, im2, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_Const<T1, T2, VoxelFunc> body(im1, im2, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_Const<T1, T2, VoxelFunc> body(im1, im2, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
  blocked_range<int> re(0, im2->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  if (im2->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, vf, of);
  } else {
    irtkBinaryForEachVoxelIfBody_Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
    blocked_range<int> re(0, im2->GetNumberOfVoxels() / im2->GetT());
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
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
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
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, const irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
  blocked_range<int> re(0, im2.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  if (im2.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, vf, of);
  } else {
    irtkBinaryForEachVoxelIfBody_Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
    blocked_range<int> re(0, im2.GetNumberOfVoxels() / im2.GetT());
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
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
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
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, const irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, vf);
}

// =============================================================================
// 1 const, 1 non-const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 1 const, 1 non-const images
 */
template <class T1, class T2, class VoxelFunc>
struct irtkBinaryForEachVoxelBody_1Const : public irtkForEachVoxelBody<VoxelFunc>
{
  const irtkGenericImage<T1> &im1;
        irtkGenericImage<T2> &im2;

  /// Constructor
  irtkBinaryForEachVoxelBody_1Const(const irtkGenericImage<T1> &im1,
                                          irtkGenericImage<T2> &im2,
                                    VoxelFunc &vf)
  :
    irtkForEachVoxelBody<VoxelFunc>(vf), im1(im1), im2(im2)
  {}

  /// Copy constructor
  irtkBinaryForEachVoxelBody_1Const(const irtkBinaryForEachVoxelBody_1Const &o)
  :
    irtkForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2)
  {}

  /// Split constructor
  irtkBinaryForEachVoxelBody_1Const(irtkBinaryForEachVoxelBody_1Const &o, split s)
  :
    irtkForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2)
  {}

  /// Process entire image
  void operator ()(const irtkImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkBinaryForEachVoxelBody_1Const *>(this)->_VoxelFunc(i, j, k, l, p1, p2);
    }
  }

  /// Process 1D image region
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkBinaryForEachVoxelBody_1Const *>(this)->_VoxelFunc(im2, idx, p1, p2);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im2.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkBinaryForEachVoxelBody_1Const *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2);
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

    const int s1 =  im2.GetX() - (ei - bi);
    const int s2 = (im2.GetY() - (ej - bj)) * im2.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkBinaryForEachVoxelBody_1Const *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 1 const, 1 non-const images
 */
template <class T1, class T2,
          class VoxelFunc, class OutsideFunc = irtkNaryVoxelFunction::NOP,
          class Domain = irtkImageDomain::Foreground>
struct irtkBinaryForEachVoxelIfBody_1Const : public irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  const irtkGenericImage<T1> &im1;
        irtkGenericImage<T2> &im2;

  /// Constructor
  irtkBinaryForEachVoxelIfBody_1Const(const irtkGenericImage<T1> &im1,
                                            irtkGenericImage<T2> &im2,
                                      VoxelFunc &vf, OutsideFunc &of)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of), im1(im1), im2(im2)
  {}

  /// Copy constructor
  irtkBinaryForEachVoxelIfBody_1Const(const irtkBinaryForEachVoxelIfBody_1Const &o)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2)
  {}

  /// Split constructor
  irtkBinaryForEachVoxelIfBody_1Const(irtkBinaryForEachVoxelIfBody_1Const &o, split s)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2)
  {}

  /// Process entire image
  void operator ()(const irtkImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2) {
      if (Domain::IsInside(im2, i, j, k, l, p2)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkBinaryForEachVoxelIfBody_1Const *>(this)->_VoxelFunc  (i, j, k, l, p1, p2);
      } else const_cast<irtkBinaryForEachVoxelIfBody_1Const *>(this)->_OutsideFunc(i, j, k, l, p1, p2);
    }
  }

  /// Process 1D image region
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1) {
      if (Domain::IsInside(im2, idx, p2)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkBinaryForEachVoxelIfBody_1Const *>(this)->_VoxelFunc  (im2, idx, p1, p2);
      } else const_cast<irtkBinaryForEachVoxelIfBody_1Const *>(this)->_OutsideFunc(im2, idx, p1, p2);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im2.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1) {
      if (Domain::IsInside(im2, i, j, this->_k, this->_l, p2)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkBinaryForEachVoxelIfBody_1Const *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2);
      } else const_cast<irtkBinaryForEachVoxelIfBody_1Const *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2);
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

    const int s1 =  im2.GetX() - (ei - bi);
    const int s2 = (im2.GetY() - (ej - bj)) * im2.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
          T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1) {
      if (Domain::IsInside(im2, i, j, k, this->_l, p2)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkBinaryForEachVoxelIfBody_1Const *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2);
      } else const_cast<irtkBinaryForEachVoxelIfBody_1Const *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2);
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
template <class T1, class T2, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_1Const<T1, T2, VoxelFunc> body(*im1, *im2, vf);
  blocked_range<int> re(0, im2->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  if (im2->GetTSize()) {
    ForEachScalar(*im1, *im2, vf);
  } else {
    irtkBinaryForEachVoxelBody_1Const<T1, T2, VoxelFunc> body(*im1, *im2, vf);
    blocked_range<int> re(0, im2->GetNumberOfVoxels() / im2->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_1Const<T1, T2, VoxelFunc> body(*im1, *im2, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_1Const<T1, T2, VoxelFunc> body(*im1, *im2, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_1Const<T1, T2, VoxelFunc> body(*im1, *im2, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_1Const<T1, T2, VoxelFunc> body(*im1, *im2, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_1Const<T1, T2, VoxelFunc> body(im1, im2, vf);
  blocked_range<int> re(0, im2.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  if (im2.GetTSize()) {
    ForEachScalar(im1, im2, vf);
  } else {
    irtkBinaryForEachVoxelBody_1Const<T1, T2, VoxelFunc> body(im1, im2, vf);
    blocked_range<int> re(0, im2.GetNumberOfVoxels() / im2.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_1Const<T1, T2, VoxelFunc> body(im1, im2, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_1Const<T1, T2, VoxelFunc> body(im1, im2, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_1Const<T1, T2, VoxelFunc> body(im1, im2, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_1Const<T1, T2, VoxelFunc> body(im1, im2, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_1Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
  blocked_range<int> re(0, im2->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  if (im2->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, vf, of);
  } else {
    irtkBinaryForEachVoxelIfBody_1Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
    blocked_range<int> re(0, im2->GetNumberOfVoxels() / im2->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_1Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_1Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_1Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_1Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_1Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
  blocked_range<int> re(0, im2.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  if (im2.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, vf, of);
  } else {
    irtkBinaryForEachVoxelIfBody_1Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
    blocked_range<int> re(0, im2.GetNumberOfVoxels() / im2.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_1Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_1Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_1Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_1Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_1Const<T1, T2, VoxelFunc> body(*im1, *im2, vf);
  blocked_range<int> re(0, im2->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  if (im2->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, vf);
  } else {
    irtkBinaryForEachVoxelBody_1Const<T1, T2, VoxelFunc> body(*im1, *im2, vf);
    blocked_range<int> re(0, im2->GetNumberOfVoxels() / im2->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_1Const<T1, T2, VoxelFunc> body(*im1, *im2, vf);
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
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_1Const<T1, T2, VoxelFunc> body(*im1, *im2, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_1Const<T1, T2, VoxelFunc> body(*im1, *im2, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_1Const<T1, T2, VoxelFunc> body(*im1, *im2, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_1Const<T1, T2, VoxelFunc> body(im1, im2, vf);
  blocked_range<int> re(0, im2.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  if (im2.GetTSize()) {
    ParallelForEachScalar(im1, im2, vf);
  } else {
    irtkBinaryForEachVoxelBody_1Const<T1, T2, VoxelFunc> body(im1, im2, vf);
    blocked_range<int> re(0, im2.GetNumberOfVoxels() / im2.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_1Const<T1, T2, VoxelFunc> body(im1, im2, vf);
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
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_1Const<T1, T2, VoxelFunc> body(im1, im2, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_1Const<T1, T2, VoxelFunc> body(im1, im2, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody_1Const<T1, T2, VoxelFunc> body(im1, im2, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_1Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
  blocked_range<int> re(0, im2->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  if (im2->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, vf, of);
  } else {
    irtkBinaryForEachVoxelIfBody_1Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
    blocked_range<int> re(0, im2->GetNumberOfVoxels() / im2->GetT());
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
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_1Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
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
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_1Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_1Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_1Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_1Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
  blocked_range<int> re(0, im2.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  if (im2.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, vf, of);
  } else {
    irtkBinaryForEachVoxelIfBody_1Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
    blocked_range<int> re(0, im2.GetNumberOfVoxels() / im2.GetT());
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
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_1Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
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
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_1Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_1Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody_1Const<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, vf);
}

// =============================================================================
// 2 non-const images
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 2 non-const images
 */
template <class T1, class T2, class VoxelFunc>
struct irtkBinaryForEachVoxelBody : public irtkForEachVoxelBody<VoxelFunc>
{
  irtkGenericImage<T1> &im1;
  irtkGenericImage<T2> &im2;

  /// Constructor
  irtkBinaryForEachVoxelBody(irtkGenericImage<T1> &im1,
                             irtkGenericImage<T2> &im2,
                             VoxelFunc &vf)
  :
    irtkForEachVoxelBody<VoxelFunc>(vf), im1(im1), im2(im2)
  {}

  /// Copy constructor
  irtkBinaryForEachVoxelBody(const irtkBinaryForEachVoxelBody &o)
  :
    irtkForEachVoxelBody<VoxelFunc>(o), im1(o.im1), im2(o.im2)
  {}

  /// Split constructor
  irtkBinaryForEachVoxelBody(irtkBinaryForEachVoxelBody &o, split s)
  :
    irtkForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1), im2(o.im2)
  {}

  /// Process entire image
  void operator ()(const irtkImageAttributes &attr) const
  {
    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkBinaryForEachVoxelBody *>(this)->_VoxelFunc(i, j, k, l, p1, p2);
    }
  }

  /// Process 1D image region
  void operator ()(const blocked_range<int> &re) const
  {
    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkBinaryForEachVoxelBody *>(this)->_VoxelFunc(im2, idx, p1, p2);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im2.GetX() - (ei - bi);

    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkBinaryForEachVoxelBody *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1, p2);
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

    const int s1 =  im2.GetX() - (ei - bi);
    const int s2 = (im2.GetY() - (ej - bj)) * im2.GetX();

    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkBinaryForEachVoxelBody *>(this)->_VoxelFunc(i, j, k, this->_l, p1, p2);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 2 non-const images
 */
template <class T1, class T2,
          class VoxelFunc, class OutsideFunc = irtkNaryVoxelFunction::NOP,
          class Domain = irtkImageDomain::Foreground>
struct irtkBinaryForEachVoxelIfBody : public irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  irtkGenericImage<T1> &im1;
  irtkGenericImage<T2> &im2;

  /// Constructor
  irtkBinaryForEachVoxelIfBody(irtkGenericImage<T1> &im1,
                               irtkGenericImage<T2> &im2,
                               VoxelFunc &vf, OutsideFunc &of)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of), im1(im1), im2(im2)
  {}

  /// Copy constructor
  irtkBinaryForEachVoxelIfBody(const irtkBinaryForEachVoxelIfBody &o)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1), im2(o.im2)
  {}

  /// Split constructor
  irtkBinaryForEachVoxelIfBody(irtkBinaryForEachVoxelIfBody &o, split s)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1), im2(o.im2)
  {}

  /// Process entire image
  void operator ()(const irtkImageAttributes &attr) const
  {
    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1, ++p2) {
      if (Domain::IsInside(im2, i, j, k, l, p2)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkBinaryForEachVoxelIfBody *>(this)->_VoxelFunc  (i, j, k, l, p1, p2);
      } else const_cast<irtkBinaryForEachVoxelIfBody *>(this)->_OutsideFunc(i, j, k, l, p1, p2);
    }
  }

  /// Process 1D image region
  void operator ()(const blocked_range<int> &re) const
  {
    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1, p2 +=  1) {
      if (Domain::IsInside(im2, idx, p2)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkBinaryForEachVoxelIfBody *>(this)->_VoxelFunc  (im2, idx, p1, p2);
      } else const_cast<irtkBinaryForEachVoxelIfBody *>(this)->_OutsideFunc(im2, idx, p1, p2);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im2.GetX() - (ei - bi);

    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1) {
      if (Domain::IsInside(im2, i, j, this->_k, this->_l, p2)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkBinaryForEachVoxelIfBody *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1, p2);
      } else const_cast<irtkBinaryForEachVoxelIfBody *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1, p2);
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

    const int s1 =  im2.GetX() - (ei - bi);
    const int s2 = (im2.GetY() - (ej - bj)) * im2.GetX();

    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);
    T2 *p2 = im2.IsEmpty() ? NULL : im2.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2, p2 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1, p2 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1, p2 +=  1) {
      if (Domain::IsInside(im2, i, j, k, this->_l, p2)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkBinaryForEachVoxelIfBody *>(this)->_VoxelFunc  (i, j, k, this->_l, p1, p2);
      } else const_cast<irtkBinaryForEachVoxelIfBody *>(this)->_OutsideFunc(i, j, k, this->_l, p1, p2);
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
template <class T1, class T2, class VoxelFunc>
void ForEachScalar(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody<T1, T2, VoxelFunc> body(*im1, *im2, vf);
  blocked_range<int> re(0, im2->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  if (im2->GetTSize()) {
    ForEachScalar(*im1, *im2, vf);
  } else {
    irtkBinaryForEachVoxelBody<T1, T2, VoxelFunc> body(*im1, *im2, vf);
    blocked_range<int> re(0, im2->GetNumberOfVoxels() / im2->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody<T1, T2, VoxelFunc> body(*im1, *im2, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody<T1, T2, VoxelFunc> body(*im1, *im2, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody<T1, T2, VoxelFunc> body(*im1, *im2, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody<T1, T2, VoxelFunc> body(*im1, *im2, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, *im2, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachScalar(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody<T1, T2, VoxelFunc> body(im1, im2, vf);
  blocked_range<int> re(0, im2.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  if (im2.GetTSize()) {
    ForEachScalar(im1, im2, vf);
  } else {
    irtkBinaryForEachVoxelBody<T1, T2, VoxelFunc> body(im1, im2, vf);
    blocked_range<int> re(0, im2.GetNumberOfVoxels() / im2.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody<T1, T2, VoxelFunc> body(im1, im2, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody<T1, T2, VoxelFunc> body(im1, im2, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody<T1, T2, VoxelFunc> body(im1, im2, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody<T1, T2, VoxelFunc> body(im1, im2, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
  blocked_range<int> re(0, im2->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachScalarIf(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  if (im2->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, *im2, vf, of);
  } else {
    irtkBinaryForEachVoxelIfBody<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
    blocked_range<int> re(0, im2->GetNumberOfVoxels() / im2->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, *im2, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
  blocked_range<int> re(0, im2.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachScalarIf(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  if (im2.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, im2, vf, of);
  } else {
    irtkBinaryForEachVoxelIfBody<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
    blocked_range<int> re(0, im2.GetNumberOfVoxels() / im2.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachScalar(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody<T1, T2, VoxelFunc> body(*im1, *im2, vf);
  blocked_range<int> re(0, im2->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  if (im2->GetTSize()) {
    ParallelForEachScalar(*im1, *im2, vf);
  } else {
    irtkBinaryForEachVoxelBody<T1, T2, VoxelFunc> body(*im1, *im2, vf);
    blocked_range<int> re(0, im2->GetNumberOfVoxels() / im2->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody<T1, T2, VoxelFunc> body(*im1, *im2, vf);
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
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody<T1, T2, VoxelFunc> body(*im1, *im2, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody<T1, T2, VoxelFunc> body(*im1, *im2, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody<T1, T2, VoxelFunc> body(*im1, *im2, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, *im2, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachScalar(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody<T1, T2, VoxelFunc> body(im1, im2, vf);
  blocked_range<int> re(0, im2.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  if (im2.GetTSize()) {
    ParallelForEachScalar(im1, im2, vf);
  } else {
    irtkBinaryForEachVoxelBody<T1, T2, VoxelFunc> body(im1, im2, vf);
    blocked_range<int> re(0, im2.GetNumberOfVoxels() / im2.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody<T1, T2, VoxelFunc> body(im1, im2, vf);
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
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody<T1, T2, VoxelFunc> body(im1, im2, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody<T1, T2, VoxelFunc> body(im1, im2, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkBinaryForEachVoxelBody<T1, T2, VoxelFunc> body(im1, im2, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class T2, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
  blocked_range<int> re(0, im2->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachScalarIf(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  if (im2->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, *im2, vf, of);
  } else {
    irtkBinaryForEachVoxelIfBody<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
    blocked_range<int> re(0, im2->GetNumberOfVoxels() / im2->GetT());
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
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
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
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody<T1, T2, VoxelFunc, OutsideFunc, Domain> body(*im1, *im2, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, irtkGenericImage<T2> *im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, *im2, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
  blocked_range<int> re(0, im2.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachScalarIf(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  if (im2.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, im2, vf, of);
  } else {
    irtkBinaryForEachVoxelIfBody<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
    blocked_range<int> re(0, im2.GetNumberOfVoxels() / im2.GetT());
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
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
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
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf, OutsideFunc &of)
{
  irtkBinaryForEachVoxelIfBody<T1, T2, VoxelFunc, OutsideFunc, Domain> body(im1, im2, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, im2, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class T2, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, irtkGenericImage<T2> &im2)
{
  if (VoxelFunc::IsReduction()) _irtkforeachbinaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, im2, vf);
}

#endif
