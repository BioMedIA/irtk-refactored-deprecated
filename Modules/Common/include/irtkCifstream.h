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

#ifndef _IRTKCIFSTREAM_H
#define _IRTKCIFSTREAM_H

#include <irtkObject.h>
#ifdef HAS_ZLIB
#include <zlib.h>
#endif

/**
 * Class for reading (compressed) file streams.
 *
 * This class defines and implements functions for reading compressed file
 * streams. The file streams can be either uncompressed or compressed.
 */

class irtkCifstream : public irtkObject
{
  irtkObjectMacro(irtkCifstream);

  /// File pointer to potentially compressed file
#ifdef HAS_ZLIB
  gzFile _File;
#else
  FILE  *_File;
#endif

  /// Flag indicating whether file bytes are swapped
  irtkPublicAttributeMacro(bool, Swapped);

public:

  /// Constructor
  irtkCifstream(const char * = NULL);

  /// Destructor
  ~irtkCifstream();

  /// Read n data as array (possibly compressed) from offset
  bool Read(char *, long, long);

  /// Read n data as array of char (possibly compressed) from offset
  bool ReadAsChar(char *, long, long = -1);

  /// Read n data as array of unsigned char (possibly compressed) from offset
  bool ReadAsUChar(unsigned char *, long, long = -1);

  /// Read n data as array of short (possibly compressed) from offset
  bool ReadAsShort(short *, long, long = -1);

  /// Read n data as array of unsigned short (possibly compressed) from offset
  bool ReadAsUShort(unsigned short *, long, long = -1);

  /// Read n data as array of short (possibly compressed) from offset
  bool ReadAsInt(int *, long, long = -1);

  /// Read n data as array of unsigned short (possibly compressed) from offset
  bool ReadAsUInt(unsigned int *, long, long = -1);

  /// Read n data as array of short  compressed) from offset
  bool ReadAsFloat(float *, long, long offset = -1);

  /// Read n data as array of unsigned short (possibly compressed) from offset
  bool ReadAsDouble(double *, long, long = -1);

  /// Read n data as string (possibly compressed) from current position
  bool ReadAsString(char *, long, long = -1);

  /// Open file
  void Open(const char *);

  /// Close file
  void Close();

  /// Go to position in file
  void Seek(long);

  /// Current position in file
  long Tell() const;

  /// Returns whether file is swapped
  /// \deprecated Use Swapped() instead.
  int IsSwapped() const;

  /// Sets whether file is swapped
  /// \deprecated Use Swapped(bool) instead.
  void IsSwapped(int);

};


#endif
