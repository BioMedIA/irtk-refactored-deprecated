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

#include <irtkCofstream.h>

// -----------------------------------------------------------------------------
irtkCofstream::irtkCofstream(const char *fname)
{
  _File = NULL;
#ifdef HAS_ZLIB
  _ZFile = NULL;
#endif
#ifndef WORDS_BIGENDIAN
  _Swapped = true;
#else
  _Swapped = false;
#endif
  _Compressed = false;
  if (fname) Open(fname);
}

// -----------------------------------------------------------------------------
irtkCofstream::~irtkCofstream()
{
  Close();
}

// -----------------------------------------------------------------------------
void irtkCofstream::Open(const char *fname)
{
  size_t len = strlen(fname);
  if (len == 0) {
    cerr << "irtkCofstream::Open: Filename is empty!" << endl;
    exit(1);
  }
  if (len > 3 && (strncmp(fname + len-3, ".gz", 3) == 0 || strncmp(fname + len-3, ".GZ", 3) == 0)) {
#ifdef HAS_ZLIB
    _ZFile = gzopen(fname, "wb");
    if (_ZFile == NULL) {
      cerr << "cofstream::Open: Cannot open file " << fname << endl;
      exit(1);
    }
    _Compressed = true;
#else
    cerr << "cofstream::Open: Cannot write compressed file without zlib" << endl;
    exit(1);
#endif
  } else {
    _File = fopen(fname, "wb");
    if (_File == NULL) {
      cerr << "cofstream::Open: Cannot open file " << fname << endl;
      exit(1);
    }
    _Compressed = false;
  }
}

// -----------------------------------------------------------------------------
void irtkCofstream::Close()
{
#ifdef HAS_ZLIB
  if (_ZFile) {
    gzclose(_ZFile);
    _ZFile = NULL;
  }
#endif
  if (_File) {
    fclose(_File);
    _File = NULL;
  }
}

// -----------------------------------------------------------------------------
int irtkCofstream::IsCompressed() const
{
  return static_cast<int>(_Compressed);
}

// -----------------------------------------------------------------------------
void irtkCofstream::IsCompressed(int compressed)
{
  _Compressed = static_cast<bool>(compressed);
}

// -----------------------------------------------------------------------------
int irtkCofstream::IsSwapped() const
{
  return static_cast<int>(_Swapped);
}

// -----------------------------------------------------------------------------
void irtkCofstream::IsSwapped(int swapped)
{
  _Swapped = static_cast<bool>(swapped);
}

// -----------------------------------------------------------------------------
bool irtkCofstream::Write(const char *data, long offset, long length)
{
#ifdef HAS_ZLIB
  if (_Compressed) {
    if (offset != -1) {
      if (gztell(_ZFile) > offset) {
        cerr << "Error: Writing compressed files only supports forward seek (pos="
             << gztell(_ZFile) << ", offset=" << offset << ")" << endl;
        exit(1);
      }
      gzseek(_ZFile, offset, SEEK_SET);
    }
    return (gzwrite(_ZFile, data, length) == length);
  }
#endif
  if (offset != -1) fseek(_File, offset, SEEK_SET);
  return (fwrite(data, length, 1, _File) == 1);
}

// -----------------------------------------------------------------------------
bool irtkCofstream::WriteAsChar(char data, long offset)
{
  return Write((char *)&data, offset, sizeof(char));
}

// -----------------------------------------------------------------------------
bool irtkCofstream::WriteAsChar(const char *data, long length, long offset)
{
  return Write((char *)data, offset, length*sizeof(char));
}

// -----------------------------------------------------------------------------
bool irtkCofstream::WriteAsUChar(unsigned char data, long offset)
{
  return Write((char *)&data, offset, sizeof(unsigned char));
}

// -----------------------------------------------------------------------------
bool irtkCofstream::WriteAsUChar(const unsigned char *data, long length, long offset)
{
  return Write((char *)data, offset, length*sizeof(unsigned char));
}

// -----------------------------------------------------------------------------
bool irtkCofstream::WriteAsShort(short data, long offset)
{
  if (_Swapped) swap16((char *)&data, (char *)&data, 1);
  return Write((char *)&data, offset, sizeof(short));
}

// -----------------------------------------------------------------------------
bool irtkCofstream::WriteAsShort(const short *data, long length, long offset)
{
  if (_Swapped) swap16((char *)data, (char *)data, length);
  if (!Write((char *)data, offset, length*sizeof(short))) return false;
  if (_Swapped) swap16((char *)data, (char *)data, length);
  return true;
}

// -----------------------------------------------------------------------------
bool irtkCofstream::WriteAsUShort(unsigned short data, long offset)
{
  if (_Swapped) swap16((char *)&data, (char *)&data, 1);
  return Write((char *)&data, offset, sizeof(unsigned short));
}

// -----------------------------------------------------------------------------
bool irtkCofstream::WriteAsUShort(const unsigned short *data, long length, long offset)
{
  if (_Swapped) swap16((char *)data, (char *)data, length);
  if (!Write((char *)data, offset, length*sizeof(unsigned short))) return false;
  if (_Swapped) swap16((char *)data, (char *)data, length);
  return true;
}

// -----------------------------------------------------------------------------
bool irtkCofstream::WriteAsInt(int data, long offset)
{
  if (_Swapped) swap32((char *)&data, (char *)&data, 1);
  return Write((char *)&data, offset, sizeof(int));
}

// -----------------------------------------------------------------------------
bool irtkCofstream::WriteAsInt(const int *data, long length, long offset)
{
  if (_Swapped) swap32((char *)data, (char *)data, length);
  if (!Write((char *)data, offset, length*sizeof(int))) return false;
  if (_Swapped) swap32((char *)data, (char *)data, length);
  return true;
}

// -----------------------------------------------------------------------------
bool irtkCofstream::WriteAsUInt(unsigned int data, long offset)
{
  if (_Swapped) swap32((char *)&data, (char *)&data, 1);
  return Write((char *)&data, offset, sizeof(unsigned int));
}

// -----------------------------------------------------------------------------
bool irtkCofstream::WriteAsUInt(const unsigned int *data, long length, long offset)
{
  if (_Swapped) swap32((char *)data, (char *)data, length);
  if (!Write((char *)data, offset, length*sizeof(unsigned int))) return false;
  if (_Swapped) swap32((char *)data, (char *)data, length);
  return true;
}

// -----------------------------------------------------------------------------
bool irtkCofstream::WriteAsFloat(float data, long offset)
{
  if (_Swapped) swap32((char *)&data, (char *)&data, 1);
  return Write((char *)&data, offset, sizeof(float));
}

// -----------------------------------------------------------------------------
bool irtkCofstream::WriteAsFloat(const float *data, long length, long offset)
{
  if (_Swapped) swap32((char *)data, (char *)data, length);
  if (!Write((char *)data, offset, length*sizeof(float))) return false;
  if (_Swapped) swap32((char *)data, (char *)data, length);
  return true;
}

// -----------------------------------------------------------------------------
bool irtkCofstream::WriteAsDouble(const double *data, long length, long offset)
{
  if (_Swapped) swap64((char *)data, (char *)data, length);
  if (!Write((char *)data, offset, length*sizeof(double))) return false;
  if (_Swapped) swap64((char *)data, (char *)data, length);
  return true;
}

// -----------------------------------------------------------------------------
bool irtkCofstream::WriteAsDouble(double data, long offset)
{
  if (_Swapped) swap64((char *)&data, (char *)&data, 1);
  return Write((char *)&data, offset, sizeof(double));
}

// -----------------------------------------------------------------------------
bool irtkCofstream::WriteAsString(const char *data, long offset)
{
#ifdef HAS_ZLIB
  if (_Compressed) {
    if (offset != -1) gzseek(_ZFile, offset, SEEK_SET);
    return (gzputs(_ZFile, data) > 0);
  }
#endif
  if (offset != -1) fseek(_File, offset, SEEK_SET);
  return (fputs(data, _File) != EOF);
}
