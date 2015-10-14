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

#include <irtkCifstream.h>

// -----------------------------------------------------------------------------
irtkCifstream::irtkCifstream(const char *fname)
:
  _File(NULL), _Swapped(true)
{
#ifdef WORDS_BIGENDIAN
  _Swapped = false;
#endif
  if (fname) Open(fname);
}

// -----------------------------------------------------------------------------
irtkCifstream::~irtkCifstream()
{
  Close();
}

// -----------------------------------------------------------------------------
void irtkCifstream::Open(const char *fname)
{
#ifdef HAS_ZLIB
  _File = gzopen(fname, "rb");
#else
  _File = fopen(fname, "rb");
#endif
  if (_File == NULL) {
    cerr << "cifstream::Open: Cannot open file " << fname << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
void irtkCifstream::Close()
{
  if (_File != NULL) {
#ifdef HAS_ZLIB
    gzclose(_File);
#else
    fclose(_File);
#endif
    _File = NULL;
  }
}

// -----------------------------------------------------------------------------
int irtkCifstream::IsSwapped() const
{
  return static_cast<int>(_Swapped);
}

// -----------------------------------------------------------------------------
void irtkCifstream::IsSwapped(int swapped)
{
  _Swapped = static_cast<bool>(swapped);
}

// -----------------------------------------------------------------------------
long irtkCifstream::Tell() const
{
#ifdef HAS_ZLIB
  return gztell(_File);
#else
  return ftell(_File);
#endif
}

// -----------------------------------------------------------------------------
void irtkCifstream::Seek(long offset)
{
#ifdef HAS_ZLIB
  gzseek(_File, offset, SEEK_SET);
#else
  fseek(_File, offset, SEEK_SET);
#endif
}

// -----------------------------------------------------------------------------
bool irtkCifstream::Read(char *mem, long start, long num)
{
#ifdef HAS_ZLIB
  if (start != -1) gzseek(_File, start, SEEK_SET);
  return (gzread(_File, mem, num) == num);
#else
  if (start != -1) fseek(_File, start, SEEK_SET);
  return (fread(mem, num, 1, _File) == 1);
#endif
}

// -----------------------------------------------------------------------------
bool irtkCifstream::ReadAsChar(char *data, long length, long offset)
{
  return Read((char *)data, offset, length * sizeof(char));
}

// -----------------------------------------------------------------------------
bool irtkCifstream::ReadAsUChar(unsigned char *data, long length, long offset)
{
  return Read((char *)data, offset, length * sizeof(unsigned char));
}

// -----------------------------------------------------------------------------
bool irtkCifstream::ReadAsShort(short *data, long length, long offset)
{
  if (!Read((char *)data, offset, length * sizeof(short))) return false;
  if (_Swapped) swap16((char *)data, (char *)data, length);
  return true;
}

// -----------------------------------------------------------------------------
bool irtkCifstream::ReadAsUShort(unsigned short *data, long length, long offset)
{
  if (!Read((char *)data, offset, length * sizeof(unsigned short))) return false;
  if (_Swapped) swap16((char *)data, (char *)data, length);
  return true;
}

// -----------------------------------------------------------------------------
bool irtkCifstream::ReadAsInt(int *data, long length, long offset)
{
  if (!Read((char *)data, offset, length * sizeof(int))) return false;
  if (_Swapped) swap32((char *)data, (char *)data, length);
  return true;
}

// -----------------------------------------------------------------------------
bool irtkCifstream::ReadAsUInt(unsigned int *data, long length, long offset)
{
  if (!Read((char *)data, offset, length * sizeof(unsigned int))) return false;
  if (_Swapped) swap32((char *)data, (char *)data, length);
  return true;
}

// -----------------------------------------------------------------------------
bool irtkCifstream::ReadAsFloat(float *data, long length, long offset)
{
  if (!Read((char *)data, offset, length * sizeof(float))) return false;
  if (_Swapped) swap32((char *)data, (char *)data, length);
  return true;
}

// -----------------------------------------------------------------------------
bool irtkCifstream::ReadAsDouble(double *data, long length, long offset)
{
  if (!Read((char *)data, offset, length * sizeof(double))) return false;
  if (_Swapped) swap64((char *)data, (char *)data, length);
  return true;
}

// -----------------------------------------------------------------------------
bool irtkCifstream::ReadAsString(char *data, long length, long offset)
{
  // Read string
#ifdef HAS_ZLIB
  if (offset != -1) gzseek(_File, offset, SEEK_SET);
  if (gzgets(_File, data, length) == Z_NULL) return false;
#else
  if (offset != -1) fseek(_File, offset, SEEK_SET);
  if (fgets(data, length, _File) == NULL) return false;
#endif

  // Discard end-of-line character(s)
  const size_t len = strlen(data);
  if (len > 0 && (data[len-1] == '\n' || data[len-1] == '\r')) {
    data[len-1] = '\0';
    if (len > 1 && data[len-2] == '\r') {
      data[len-2] = '\0';
    }
  }
  return true;
}

