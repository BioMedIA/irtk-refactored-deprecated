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

#ifndef _IRTKZLIB_H
#define _IRTKZLIB_H


#ifdef HAS_ZLIB
#  include <zlib.h>
#  define irtkGzFile gzFile
#else
#  define irtkGzFile FILE*
#endif


/// Read bytes from compressed file
int ReadCompressed(irtkGzFile, char *, long, long);


#endif
