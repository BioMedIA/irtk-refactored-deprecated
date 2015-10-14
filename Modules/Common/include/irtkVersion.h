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

#ifndef _IRTKVERSION_H
#define _IRTKVERSION_H

#include <irtkObject.h>


// -----------------------------------------------------------------------------
/// Software version object
class irtkVersion : public irtkObject
{
  irtkObjectMacro(irtkVersion);

  irtkPublicAttributeMacro(unsigned int, Major); ///< Major version number
  irtkPublicAttributeMacro(unsigned int, Minor); ///< Minor version number
  irtkPublicAttributeMacro(unsigned int, Patch); ///< Patch number
 
public:

  /// Constructor
  irtkVersion(unsigned int major = 0u,
              unsigned int minor = 0u,
              unsigned int patch = 0u);

  /// Constructor
  irtkVersion(int major, int minor = 0, int patch = 0);

  /// Constructor
  irtkVersion(const char *);
 
  /// Copy constructor
  irtkVersion(const irtkVersion &);
 
  /// Assignment operator
  irtkVersion &operator =(const irtkVersion &);
 
  /// Whether this version is valid
  operator bool() const;

  /// Equality operator
  bool operator ==(const irtkVersion &) const;
 
  /// Inequality operator
  bool operator !=(const irtkVersion &) const;

  /// Less operator
  bool operator <(const irtkVersion &) const;

  /// Greater operator
  bool operator >(const irtkVersion &) const;

  /// Less or equal operator
  bool operator <=(const irtkVersion &) const;

  /// Greater or equal operator
  bool operator >=(const irtkVersion &) const;

  /// Get version information as string
  string ToString() const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

#include <irtkCommon.h>

// -----------------------------------------------------------------------------
/// Write version to output stream
inline ostream &operator <<(ostream &os, const irtkVersion &version)
{
  os << version.Major() << "." << version.Minor();
  if (version.Patch() != 0u) os << "." << version.Patch();
  return os;
}

// -----------------------------------------------------------------------------
/// Read version from input stream
inline istream &operator >>(istream &is, irtkVersion &version)
{
  string token;
  char   c;
  while (is.get(c) && (isdigit(c) || c == '.')) token.push_back(c);
  if (is.bad()) return is;
  if (is.eof()) is.clear(ios::eofbit); // reset failbit
  istringstream ss(token);
  unsigned int number[3] = {0u, 0u, 0u};
  int          i         = -1;
  while (getline(ss, token, '.')) {
    if (++i == 3 || !FromString(token.c_str(), number[i])) {
      i = -1;
      break;
    }
  }
  if (i == -1) is.setstate(ios::failbit);
  else {
    version.Major(number[0]);
    version.Minor(number[1]);
    version.Patch(number[2]);
  }
  return is;
}

// -----------------------------------------------------------------------------
inline string irtkVersion::ToString() const
{
  ostringstream ss;
  ss << (*this);
  return ss.str();
}

// -----------------------------------------------------------------------------
inline irtkVersion::irtkVersion(unsigned int major,
                                unsigned int minor,
                                unsigned int patch)
:
  _Major(major), _Minor(minor), _Patch(patch)
{
}

// -----------------------------------------------------------------------------
inline irtkVersion::irtkVersion(int major, int minor, int patch)
:
  _Major(major >= 0 ? static_cast<unsigned int>(major) : 0u),
  _Minor(minor >= 0 ? static_cast<unsigned int>(minor) : 0u),
  _Patch(patch >= 0 ? static_cast<unsigned int>(patch) : 0u)
{
}

// -----------------------------------------------------------------------------
inline irtkVersion::irtkVersion(const char *str)
:
  _Major(0u), _Minor(0u), _Patch(0u)
{
  FromString(str, *this);
}

// -----------------------------------------------------------------------------
inline irtkVersion::irtkVersion(const irtkVersion &other)
:
  _Major(other._Major), _Minor(other._Minor), _Patch(other._Patch)
{
}

// -----------------------------------------------------------------------------
inline irtkVersion &irtkVersion::operator =(const irtkVersion &other)
{
  _Major = other._Major;
  _Minor = other._Minor;
  _Patch = other._Patch;
  return *this;
}

// -----------------------------------------------------------------------------
inline irtkVersion::operator bool() const
{
  return (_Major != 0u || _Minor != 0u || _Patch != 0u);
}

// -----------------------------------------------------------------------------
inline bool irtkVersion::operator ==(const irtkVersion &rhs) const
{
  return _Major == rhs._Major && _Minor == rhs._Minor && _Patch == rhs._Patch;
}

// -----------------------------------------------------------------------------
inline bool irtkVersion::operator !=(const irtkVersion &rhs) const
{
  return !(*this == rhs);
}

// -----------------------------------------------------------------------------
// Note: Invalid/development branch version 0.0.0 is greater than any other
inline bool irtkVersion::operator <(const irtkVersion &rhs) const
{
  return bool(*this) && (!rhs || (_Major < rhs._Major || (_Major == rhs._Major && (_Minor < rhs._Minor || (_Minor == rhs._Minor && _Patch < rhs._Patch)))));
}

// -----------------------------------------------------------------------------
inline bool irtkVersion::operator <=(const irtkVersion &rhs) const
{
  return (*this == rhs) || (*this < rhs);
}

// -----------------------------------------------------------------------------
// Note: Invalid/development branch version 0.0.0 is greater than any other
inline bool irtkVersion::operator >(const irtkVersion &rhs) const
{
  return !(*this <= rhs);
}

// -----------------------------------------------------------------------------
inline bool irtkVersion::operator >=(const irtkVersion &rhs) const
{
  return (*this == rhs) || (*this > rhs);
}


#endif
