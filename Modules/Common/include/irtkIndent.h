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

#ifndef _IRTKINDENT_H

#define _IRTKINDENT_H

/**
 * Auxiliary class for output indentation.
 */
class irtkIndent
{
protected:

  /// Number of indentation levels
  int _N;

  /// Number of whitespace characters per indentation level
  int _NumberOfSpaces;

public:

  /// Constructor
  irtkIndent(int n = 0, int s = 2) : _N(n), _NumberOfSpaces(s) {}

  /// Copy constructor
  irtkIndent(const irtkIndent &o) : _N(o._N), _NumberOfSpaces(o._NumberOfSpaces) {}

  /// Assignment operator
  irtkIndent &operator= (int n)
  {
    _N = n;
    return *this;
  }

  /// Assignment operator
  irtkIndent &operator= (const irtkIndent &indent)
  {
    _N              = indent._N;
    _NumberOfSpaces = indent._NumberOfSpaces;
    return *this;
  }

  /// Pre-decrement operator
  irtkIndent &operator-- ()
  {
    _N--;
    return *this;
  }

  /// Pre-increment operator
  irtkIndent &operator++ ()
  {
    _N++;
    return *this;
  }

  /// Post-decrement operator
  irtkIndent operator-- (int)
  {
    irtkIndent pre(*this);
    --(*this);
    return pre;
  }

  /// Post-increment operator
  irtkIndent operator++ (int)
  {
    irtkIndent pre(*this);
    ++(*this);
    return pre;
  }

  /// Add indentation to this indentation
  irtkIndent &operator+= (const irtkIndent &indent)
  {
    _N += indent._N;
    return *this;
  }

  /// Add two indentations
  irtkIndent operator+ (const irtkIndent &indent) const
  {
    return irtkIndent(_N + indent._N);
  }

  /// Subtract indentation from this indentation
  irtkIndent &operator-= (const irtkIndent &indent)
  {
    _N -= indent._N;
    return *this;
  }

  /// Subtract two indentations from another
  irtkIndent operator- (const irtkIndent &indent) const
  {
    return irtkIndent(_N + indent._N);
  }

  /// Get indentation level
  int Level() const
  {
    return _N;
  }

  /// Set number of indentation whitespace characters per level
  void SpacesPerLevel(int n)
  {
    _NumberOfSpaces = n;
  }

  /// Get number of indentation whitespace characters per level
  int SpacesPerLevel() const
  {
    return _NumberOfSpaces;
  }

  /// Number of space characters
  int Spaces() const
  {
    return _N * _NumberOfSpaces;
  }
};

// ---------------------------------------------------------------------------
// Streaming operator
inline ostream &operator<< (ostream &os, const irtkIndent &indent)
{
  for (int i = 0; i < indent.Spaces(); i++) os << " ";
  return os;
}


#endif
