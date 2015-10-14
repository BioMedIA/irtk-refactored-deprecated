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

#ifndef _irtkForEachUnaryVoxelFunction_H
#define _irtkForEachUnaryVoxelFunction_H

#include <irtkVoxelFunction.h>


inline void _irtkforeachunaryvoxelfunction_must_not_be_reduction()
{
  cerr << "(Parallel)ForEachVoxel(If): Voxel reductions must be passed by reference!"
          " Pass voxel functor object(s) as last argument(s) instead of first." << endl;
  exit(1);
}


// =============================================================================
// 1 const image
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 1 const image
 */
template <class T1, class VoxelFunc>
struct irtkUnaryForEachVoxelBody_Const : public irtkForEachVoxelBody<VoxelFunc>
{
  const irtkGenericImage<T1> &im1;

  /// Constructor
  irtkUnaryForEachVoxelBody_Const(const irtkGenericImage<T1> &im1,
                                  VoxelFunc &vf)
  :
    irtkForEachVoxelBody<VoxelFunc>(vf), im1(im1)
  {}

  /// Copy constructor
  irtkUnaryForEachVoxelBody_Const(const irtkUnaryForEachVoxelBody_Const &o)
  :
    irtkForEachVoxelBody<VoxelFunc>(o), im1(o.im1)
  {}

  /// Split constructor
  irtkUnaryForEachVoxelBody_Const(irtkUnaryForEachVoxelBody_Const &o, split s)
  :
    irtkForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1)
  {}

  /// Process entire image
  void operator ()(const irtkImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkUnaryForEachVoxelBody_Const *>(this)->_VoxelFunc(i, j, k, l, p1);
    }
  }

  /// Process 1D image region
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkUnaryForEachVoxelBody_Const *>(this)->_VoxelFunc(im1, idx, p1);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im1.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkUnaryForEachVoxelBody_Const *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1);
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

    const int s1 =  im1.GetX() - (ei - bi);
    const int s2 = (im1.GetY() - (ej - bj)) * im1.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkUnaryForEachVoxelBody_Const *>(this)->_VoxelFunc(i, j, k, this->_l, p1);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 1 const image
 */
template <class T1,
          class VoxelFunc, class OutsideFunc = irtkNaryVoxelFunction::NOP,
          class Domain = irtkImageDomain::Foreground>
struct irtkUnaryForEachVoxelIfBody_Const : public irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  const irtkGenericImage<T1> &im1;

  /// Constructor
  irtkUnaryForEachVoxelIfBody_Const(const irtkGenericImage<T1> &im1,
                                    VoxelFunc &vf, OutsideFunc &of)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of), im1(im1)
  {}

  /// Copy constructor
  irtkUnaryForEachVoxelIfBody_Const(const irtkUnaryForEachVoxelIfBody_Const &o)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1)
  {}

  /// Split constructor
  irtkUnaryForEachVoxelIfBody_Const(irtkUnaryForEachVoxelIfBody_Const &o, split s)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1)
  {}

  /// Process entire image
  void operator ()(const irtkImageAttributes &attr) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1) {
      if (Domain::IsInside(im1, i, j, k, l, p1)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkUnaryForEachVoxelIfBody_Const *>(this)->_VoxelFunc  (i, j, k, l, p1);
      } else const_cast<irtkUnaryForEachVoxelIfBody_Const *>(this)->_OutsideFunc(i, j, k, l, p1);
    }
  }

  /// Process 1D image region
  void operator ()(const blocked_range<int> &re) const
  {
    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1) {
      if (Domain::IsInside(im1, idx, p1)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkUnaryForEachVoxelIfBody_Const *>(this)->_VoxelFunc  (im1, idx, p1);
      } else const_cast<irtkUnaryForEachVoxelIfBody_Const *>(this)->_OutsideFunc(im1, idx, p1);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im1.GetX() - (ei - bi);

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1) {
      if (Domain::IsInside(im1, i, j, this->_k, this->_l, p1)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkUnaryForEachVoxelIfBody_Const *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1);
      } else const_cast<irtkUnaryForEachVoxelIfBody_Const *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1);
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

    const int s1 =  im1.GetX() - (ei - bi);
    const int s2 = (im1.GetY() - (ej - bj)) * im1.GetX();

    const T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1) {
      if (Domain::IsInside(im1, i, j, k, this->_l, p1)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkUnaryForEachVoxelIfBody_Const *>(this)->_VoxelFunc  (i, j, k, this->_l, p1);
      } else const_cast<irtkUnaryForEachVoxelIfBody_Const *>(this)->_OutsideFunc(i, j, k, this->_l, p1);
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
template <class T1, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody_Const<T1, VoxelFunc> body(*im1, vf);
  blocked_range<int> re(0, im1->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  if (im1->GetTSize()) {
    ForEachScalar(*im1, vf);
  } else {
    irtkUnaryForEachVoxelBody_Const<T1, VoxelFunc> body(*im1, vf);
    blocked_range<int> re(0, im1->GetNumberOfVoxels() / im1->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody_Const<T1, VoxelFunc> body(*im1, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody_Const<T1, VoxelFunc> body(*im1, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody_Const<T1, VoxelFunc> body(*im1, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody_Const<T1, VoxelFunc> body(*im1, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachScalar(const irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody_Const<T1, VoxelFunc> body(im1, vf);
  blocked_range<int> re(0, im1.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(const irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  if (im1.GetTSize()) {
    ForEachScalar(im1, vf);
  } else {
    irtkUnaryForEachVoxelBody_Const<T1, VoxelFunc> body(im1, vf);
    blocked_range<int> re(0, im1.GetNumberOfVoxels() / im1.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody_Const<T1, VoxelFunc> body(im1, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody_Const<T1, VoxelFunc> body(im1, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody_Const<T1, VoxelFunc> body(im1, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody_Const<T1, VoxelFunc> body(im1, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody_Const<T1, VoxelFunc, OutsideFunc, Domain> body(*im1, vf, of);
  blocked_range<int> re(0, im1->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, VoxelFunc &vf, OutsideFunc &of)
{
  if (im1->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, vf, of);
  } else {
    irtkUnaryForEachVoxelIfBody_Const<T1, VoxelFunc, OutsideFunc, Domain> body(*im1, vf, of);
    blocked_range<int> re(0, im1->GetNumberOfVoxels() / im1->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody_Const<T1, VoxelFunc, OutsideFunc, Domain> body(*im1, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody_Const<T1, VoxelFunc, OutsideFunc, Domain> body(*im1, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody_Const<T1, VoxelFunc, OutsideFunc, Domain> body(*im1, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody_Const<T1, VoxelFunc, OutsideFunc, Domain> body(*im1, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody_Const<T1, VoxelFunc, OutsideFunc, Domain> body(im1, vf, of);
  blocked_range<int> re(0, im1.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachScalarIf(const irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, VoxelFunc &vf, OutsideFunc &of)
{
  if (im1.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, vf, of);
  } else {
    irtkUnaryForEachVoxelIfBody_Const<T1, VoxelFunc, OutsideFunc, Domain> body(im1, vf, of);
    blocked_range<int> re(0, im1.GetNumberOfVoxels() / im1.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(const irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody_Const<T1, VoxelFunc, OutsideFunc, Domain> body(im1, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody_Const<T1, VoxelFunc, OutsideFunc, Domain> body(im1, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody_Const<T1, VoxelFunc, OutsideFunc, Domain> body(im1, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody_Const<T1, VoxelFunc, OutsideFunc, Domain> body(im1, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody_Const<T1, VoxelFunc> body(*im1, vf);
  blocked_range<int> re(0, im1->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  if (im1->GetTSize()) {
    ParallelForEachScalar(*im1, vf);
  } else {
    irtkUnaryForEachVoxelBody_Const<T1, VoxelFunc> body(*im1, vf);
    blocked_range<int> re(0, im1->GetNumberOfVoxels() / im1->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody_Const<T1, VoxelFunc> body(*im1, vf);
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
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody_Const<T1, VoxelFunc> body(*im1, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody_Const<T1, VoxelFunc> body(*im1, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody_Const<T1, VoxelFunc> body(*im1, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachScalar(const irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody_Const<T1, VoxelFunc> body(im1, vf);
  blocked_range<int> re(0, im1.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(const irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  if (im1.GetTSize()) {
    ParallelForEachScalar(im1, vf);
  } else {
    irtkUnaryForEachVoxelBody_Const<T1, VoxelFunc> body(im1, vf);
    blocked_range<int> re(0, im1.GetNumberOfVoxels() / im1.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody_Const<T1, VoxelFunc> body(im1, vf);
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
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody_Const<T1, VoxelFunc> body(im1, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody_Const<T1, VoxelFunc> body(im1, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody_Const<T1, VoxelFunc> body(im1, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody_Const<T1, VoxelFunc, OutsideFunc, Domain> body(*im1, vf, of);
  blocked_range<int> re(0, im1->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, VoxelFunc &vf, OutsideFunc &of)
{
  if (im1->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, vf, of);
  } else {
    irtkUnaryForEachVoxelIfBody_Const<T1, VoxelFunc, OutsideFunc, Domain> body(*im1, vf, of);
    blocked_range<int> re(0, im1->GetNumberOfVoxels() / im1->GetT());
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
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody_Const<T1, VoxelFunc, OutsideFunc, Domain> body(*im1, vf, of);
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
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody_Const<T1, VoxelFunc, OutsideFunc, Domain> body(*im1, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody_Const<T1, VoxelFunc, OutsideFunc, Domain> body(*im1, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody_Const<T1, VoxelFunc, OutsideFunc, Domain> body(*im1, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody_Const<T1, VoxelFunc, OutsideFunc, Domain> body(im1, vf, of);
  blocked_range<int> re(0, im1.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachScalarIf(const irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, VoxelFunc &vf, OutsideFunc &of)
{
  if (im1.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, vf, of);
  } else {
    irtkUnaryForEachVoxelIfBody_Const<T1, VoxelFunc, OutsideFunc, Domain> body(im1, vf, of);
    blocked_range<int> re(0, im1.GetNumberOfVoxels() / im1.GetT());
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
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody_Const<T1, VoxelFunc, OutsideFunc, Domain> body(im1, vf, of);
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
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody_Const<T1, VoxelFunc, OutsideFunc, Domain> body(im1, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, const irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody_Const<T1, VoxelFunc, OutsideFunc, Domain> body(im1, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody_Const<T1, VoxelFunc, OutsideFunc, Domain> body(im1, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, const irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, vf);
}

// =============================================================================
// 1 non-const image
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for voxel function of 1 non-const image
 */
template <class T1, class VoxelFunc>
struct irtkUnaryForEachVoxelBody : public irtkForEachVoxelBody<VoxelFunc>
{
  irtkGenericImage<T1> &im1;

  /// Constructor
  irtkUnaryForEachVoxelBody(irtkGenericImage<T1> &im1,
                            VoxelFunc &vf)
  :
    irtkForEachVoxelBody<VoxelFunc>(vf), im1(im1)
  {}

  /// Copy constructor
  irtkUnaryForEachVoxelBody(const irtkUnaryForEachVoxelBody &o)
  :
    irtkForEachVoxelBody<VoxelFunc>(o), im1(o.im1)
  {}

  /// Split constructor
  irtkUnaryForEachVoxelBody(irtkUnaryForEachVoxelBody &o, split s)
  :
    irtkForEachVoxelBody<VoxelFunc>(o, s), im1(o.im1)
  {}

  /// Process entire image
  void operator ()(const irtkImageAttributes &attr) const
  {
    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkUnaryForEachVoxelBody *>(this)->_VoxelFunc(i, j, k, l, p1);
    }
  }

  /// Process 1D image region
  void operator ()(const blocked_range<int> &re) const
  {
    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkUnaryForEachVoxelBody *>(this)->_VoxelFunc(im1, idx, p1);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im1.GetX() - (ei - bi);

    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkUnaryForEachVoxelBody *>(this)->_VoxelFunc(i, j, this->_k, this->_l, p1);
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

    const int s1 =  im1.GetX() - (ei - bi);
    const int s2 = (im1.GetY() - (ej - bj)) * im1.GetX();

    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1) {
      // const_cast such that voxel functions need only implement
      // non-const operator() which is required for parallel_reduce
      const_cast<irtkUnaryForEachVoxelBody *>(this)->_VoxelFunc(i, j, k, this->_l, p1);
    }
  }
};

// -----------------------------------------------------------------------------
/**
 * ForEachVoxel body for inside and outside unary voxel function of 1 non-const image
 */
template <class T1,
          class VoxelFunc, class OutsideFunc = irtkNaryVoxelFunction::NOP,
          class Domain = irtkImageDomain::Foreground>
struct irtkUnaryForEachVoxelIfBody : public irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>
{
  irtkGenericImage<T1> &im1;

  /// Constructor
  irtkUnaryForEachVoxelIfBody(irtkGenericImage<T1> &im1,
                              VoxelFunc &vf, OutsideFunc &of)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(vf, of), im1(im1)
  {}

  /// Copy constructor
  irtkUnaryForEachVoxelIfBody(const irtkUnaryForEachVoxelIfBody &o)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o), im1(o.im1)
  {}

  /// Split constructor
  irtkUnaryForEachVoxelIfBody(irtkUnaryForEachVoxelIfBody &o, split s)
  :
    irtkForEachVoxelIfBody<VoxelFunc, OutsideFunc>(o, s), im1(o.im1)
  {}

  /// Process entire image
  void operator ()(const irtkImageAttributes &attr) const
  {
    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels();

    const int T = (attr._dt ? attr._t : 1);

    for (int l = 0; l < T;       ++l)
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i, ++p1) {
      if (Domain::IsInside(im1, i, j, k, l, p1)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkUnaryForEachVoxelIfBody *>(this)->_VoxelFunc  (i, j, k, l, p1);
      } else const_cast<irtkUnaryForEachVoxelIfBody *>(this)->_OutsideFunc(i, j, k, l, p1);
    }
  }

  /// Process 1D image region
  void operator ()(const blocked_range<int> &re) const
  {
    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels() + re.begin();

    for (int idx = re.begin(); idx < re.end(); ++idx, p1 +=  1) {
      if (Domain::IsInside(im1, idx, p1)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkUnaryForEachVoxelIfBody *>(this)->_VoxelFunc  (im1, idx, p1);
      } else const_cast<irtkUnaryForEachVoxelIfBody *>(this)->_OutsideFunc(im1, idx, p1);
    }
  }

  /// Process 2D image region
  void operator ()(const blocked_range2d<int> &re) const
  {
    const int bi = re.cols().begin();
    const int bj = re.rows().begin();
    const int ei = re.cols().end();
    const int ej = re.rows().end();

    const int s1 = im1.GetX() - (ei - bi);

    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, this->_k, this->_l);

    for (int j = bj; j < ej; ++j, p1 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1) {
      if (Domain::IsInside(im1, i, j, this->_k, this->_l, p1)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkUnaryForEachVoxelIfBody *>(this)->_VoxelFunc  (i, j, this->_k, this->_l, p1);
      } else const_cast<irtkUnaryForEachVoxelIfBody *>(this)->_OutsideFunc(i, j, this->_k, this->_l, p1);
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

    const int s1 =  im1.GetX() - (ei - bi);
    const int s2 = (im1.GetY() - (ej - bj)) * im1.GetX();

    T1 *p1 = im1.IsEmpty() ? NULL : im1.GetPointerToVoxels(bi, bj, bk, this->_l);

    for (int k = bk; k < ek; ++k, p1 += s2)
    for (int j = bj; j < ej; ++j, p1 += s1)
    for (int i = bi; i < ei; ++i, p1 +=  1) {
      if (Domain::IsInside(im1, i, j, k, this->_l, p1)) {
             // const_cast such that voxel functions need only implement
             // non-const operator() which is required for parallel_reduce
             const_cast<irtkUnaryForEachVoxelIfBody *>(this)->_VoxelFunc  (i, j, k, this->_l, p1);
      } else const_cast<irtkUnaryForEachVoxelIfBody *>(this)->_OutsideFunc(i, j, k, this->_l, p1);
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
template <class T1, class VoxelFunc>
void ForEachScalar(irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody<T1, VoxelFunc> body(*im1, vf);
  blocked_range<int> re(0, im1->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(*im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  if (im1->GetTSize()) {
    ForEachScalar(*im1, vf);
  } else {
    irtkUnaryForEachVoxelBody<T1, VoxelFunc> body(*im1, vf);
    blocked_range<int> re(0, im1->GetNumberOfVoxels() / im1->GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(*im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody<T1, VoxelFunc> body(*im1, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, *im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody<T1, VoxelFunc> body(*im1, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody<T1, VoxelFunc> body(*im1, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody<T1, VoxelFunc> body(*im1, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, *im1, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachScalar(irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody<T1, VoxelFunc> body(im1, vf);
  blocked_range<int> re(0, im1.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachScalar(VoxelFunc vf, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachScalar(im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  if (im1.GetTSize()) {
    ForEachScalar(im1, vf);
  } else {
    irtkUnaryForEachVoxelBody<T1, VoxelFunc> body(im1, vf);
    blocked_range<int> re(0, im1.GetNumberOfVoxels() / im1.GetT());
    body(re);
    vf.join(body._VoxelFunc);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody<T1, VoxelFunc> body(im1, vf);
  body(attr);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(attr, im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(const blocked_range<int> &re, irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody<T1, VoxelFunc> body(im1, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody<T1, VoxelFunc> body(im1, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody<T1, VoxelFunc> body(im1, vf);
  body(re);
  vf.join(body._VoxelFunc);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxel(re, im1, vf);
}

// -----------------------------------------------------------------------------
// ForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(irtkGenericImage<T1> *im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody<T1, VoxelFunc, OutsideFunc, Domain> body(*im1, vf, of);
  blocked_range<int> re(0, im1->GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachScalarIf(irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(*im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(*im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(irtkGenericImage<T1> *im1, VoxelFunc &vf, OutsideFunc &of)
{
  if (im1->GetTSize()) {
    ForEachScalarIf<Domain>(*im1, vf, of);
  } else {
    irtkUnaryForEachVoxelIfBody<T1, VoxelFunc, OutsideFunc, Domain> body(*im1, vf, of);
    blocked_range<int> re(0, im1->GetNumberOfVoxels() / im1->GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(*im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody<T1, VoxelFunc, OutsideFunc, Domain> body(*im1, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, *im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, *im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, irtkGenericImage<T1> *im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody<T1, VoxelFunc, OutsideFunc, Domain> body(*im1, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody<T1, VoxelFunc, OutsideFunc, Domain> body(*im1, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody<T1, VoxelFunc, OutsideFunc, Domain> body(*im1, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(*im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, *im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, *im1, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(irtkGenericImage<T1> &im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody<T1, VoxelFunc, OutsideFunc, Domain> body(im1, vf, of);
  blocked_range<int> re(0, im1.GetNumberOfVoxels());
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachScalarIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachScalarIf(irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachScalarIf<Domain>(im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachScalarIf(VoxelFunc vf, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachScalarIf<Domain>(im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(irtkGenericImage<T1> &im1, VoxelFunc &vf, OutsideFunc &of)
{
  if (im1.GetTSize()) {
    ForEachVoxelIf<Domain>(im1, vf, of);
  } else {
    irtkUnaryForEachVoxelIfBody<T1, VoxelFunc, OutsideFunc, Domain> body(im1, vf, of);
    blocked_range<int> re(0, im1.GetNumberOfVoxels() / im1.GetT());
    body(re);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody<T1, VoxelFunc, OutsideFunc, Domain> body(im1, vf, of);
  body(attr);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(attr, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(attr, im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range<int> &re, irtkGenericImage<T1> &im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody<T1, VoxelFunc, OutsideFunc, Domain> body(im1, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody<T1, VoxelFunc, OutsideFunc, Domain> body(im1, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody<T1, VoxelFunc, OutsideFunc, Domain> body(im1, vf, of);
  body(re);
  vf.join(body._VoxelFunc);
  of.join(body._OutsideFunc);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ForEachVoxelIf<Domain>(re, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ForEachVoxelIf<Domain>(re, im1, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxel
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachScalar(irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody<T1, VoxelFunc> body(*im1, vf);
  blocked_range<int> re(0, im1->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(*im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  if (im1->GetTSize()) {
    ParallelForEachScalar(*im1, vf);
  } else {
    irtkUnaryForEachVoxelBody<T1, VoxelFunc> body(*im1, vf);
    blocked_range<int> re(0, im1->GetNumberOfVoxels() / im1->GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(*im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody<T1, VoxelFunc> body(*im1, vf);
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
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, *im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody<T1, VoxelFunc> body(*im1, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody<T1, VoxelFunc> body(*im1, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody<T1, VoxelFunc> body(*im1, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, *im1, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachScalar(irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody<T1, VoxelFunc> body(im1, vf);
  blocked_range<int> re(0, im1.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachScalar(VoxelFunc vf, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalar(im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  if (im1.GetTSize()) {
    ParallelForEachScalar(im1, vf);
  } else {
    irtkUnaryForEachVoxelBody<T1, VoxelFunc> body(im1, vf);
    blocked_range<int> re(0, im1.GetNumberOfVoxels() / im1.GetT());
    if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
    else                            parallel_for   (re, body);
  }
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody<T1, VoxelFunc> body(im1, vf);
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
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(attr, im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range<int> &re, irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody<T1, VoxelFunc> body(im1, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range<int> &re, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody<T1, VoxelFunc> body(im1, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, vf);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkUnaryForEachVoxelBody<T1, VoxelFunc> body(im1, vf);
  if (VoxelFunc::IsReduction()) { parallel_reduce(re, body); vf.join(body._VoxelFunc); }
  else                            parallel_for   (re, body);
}

// -----------------------------------------------------------------------------
template <class T1, class VoxelFunc>
void ParallelForEachVoxel(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxel(re, im1, vf);
}

// -----------------------------------------------------------------------------
// ParallelForEachVoxelIf
// -----------------------------------------------------------------------------

//
// Image arguments by pointer
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(irtkGenericImage<T1> *im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody<T1, VoxelFunc, OutsideFunc, Domain> body(*im1, vf, of);
  blocked_range<int> re(0, im1->GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachScalarIf(irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(*im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(*im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(irtkGenericImage<T1> *im1, VoxelFunc &vf, OutsideFunc &of)
{
  if (im1->GetTSize()) {
    ParallelForEachVoxelIf<Domain>(*im1, vf, of);
  } else {
    irtkUnaryForEachVoxelIfBody<T1, VoxelFunc, OutsideFunc, Domain> body(*im1, vf, of);
    blocked_range<int> re(0, im1->GetNumberOfVoxels() / im1->GetT());
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
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(*im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(*im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody<T1, VoxelFunc, OutsideFunc, Domain> body(*im1, vf, of);
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
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, *im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, *im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, irtkGenericImage<T1> *im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody<T1, VoxelFunc, OutsideFunc, Domain> body(*im1, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody<T1, VoxelFunc, OutsideFunc, Domain> body(*im1, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody<T1, VoxelFunc, OutsideFunc, Domain> body(*im1, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> *im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, *im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> *im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, *im1, vf);
}

//
// Image arguments by reference
//

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(irtkGenericImage<T1> &im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody<T1, VoxelFunc, OutsideFunc, Domain> body(im1, vf, of);
  blocked_range<int> re(0, im1.GetNumberOfVoxels());
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachScalarIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachScalarIf(irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachScalarIf<Domain>(im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachScalarIf(VoxelFunc vf, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachScalarIf<Domain>(im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(irtkGenericImage<T1> &im1, VoxelFunc &vf, OutsideFunc &of)
{
  if (im1.GetTSize()) {
    ParallelForEachVoxelIf<Domain>(im1, vf, of);
  } else {
    irtkUnaryForEachVoxelIfBody<T1, VoxelFunc, OutsideFunc, Domain> body(im1, vf, of);
    blocked_range<int> re(0, im1.GetNumberOfVoxels() / im1.GetT());
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
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody<T1, VoxelFunc, OutsideFunc, Domain> body(im1, vf, of);
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
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const irtkImageAttributes &attr, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(const irtkImageAttributes &attr, irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(attr, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const irtkImageAttributes &attr, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(attr, im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, irtkGenericImage<T1> &im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody<T1, VoxelFunc, OutsideFunc, Domain> body(im1, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range<int> &re, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range<int> &re, irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range<int> &re, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody<T1, VoxelFunc, OutsideFunc, Domain> body(im1, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range2d<int> &re, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range2d<int> &re, irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range2d<int> &re, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, vf);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, VoxelFunc &vf, OutsideFunc &of)
{
  irtkUnaryForEachVoxelIfBody<T1, VoxelFunc, OutsideFunc, Domain> body(im1, vf, of);
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) {
    parallel_reduce(re, body);
    vf.join(body._VoxelFunc);
    of.join(body._OutsideFunc);
  } else {
    parallel_for(re, body);
  }
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc, class OutsideFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, OutsideFunc of, const blocked_range3d<int> &re, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction() || OutsideFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(const blocked_range3d<int> &re, irtkGenericImage<T1> &im1, VoxelFunc &vf)
{
  irtkNaryVoxelFunction::NOP of;
  ParallelForEachVoxelIf<Domain>(re, im1, vf, of);
}

// -----------------------------------------------------------------------------
template <class Domain, class T1, class VoxelFunc>
void ParallelForEachVoxelIf(VoxelFunc vf, const blocked_range3d<int> &re, irtkGenericImage<T1> &im1)
{
  if (VoxelFunc::IsReduction()) _irtkforeachunaryvoxelfunction_must_not_be_reduction();
  ParallelForEachVoxelIf<Domain>(re, im1, vf);
}

#endif
