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

#ifndef IRTK_NULLPTR_H
#define IRTK_NULLPTR_H

// LLVM defines a nullptr macro in cstddef to __get_nullptr_t
#if !defined(nullptr) && (!defined(HAVE_nullptr) || HAVE_nullptr == 0)


// keywords required if nullptr used in CUDA code (compiled by nvcc)
#ifdef __CUDACC__
#  include <cuda.h>
#  define _NULLPTR_DECL __host__ __device__
#else
#  define _NULLPTR_DECL
#endif


namespace irtk {

const                                    // this is a const object...
class nullptr_t {
public:
  template<class T>                      // convertible to any type
  _NULLPTR_DECL operator T    *() const  // of null non-member pointer...
  { return 0; }
  template<class C, class T>             // or any type of null
  _NULLPTR_DECL operator T C::*() const  // member pointer...
  { return 0; }
  
private:
  _NULLPTR_DECL void operator&() const;  // whose address can't be taken
} nullptr = {};

} // namespace irtk


namespace std { using ::irtk::nullptr_t; }
using ::irtk::nullptr;


#undef _NULLPTR_DECL
#endif // emulation


#endif
