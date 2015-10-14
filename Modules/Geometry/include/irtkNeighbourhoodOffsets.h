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

#ifndef IRTKNEIGHBOURHOODOFFSETS_H_
#define IRTKNEIGHBOURHOODOFFSETS_H_

#include <irtkBaseImage.h>

typedef enum {CONNECTIVITY_04,
	      CONNECTIVITY_06,
	      CONNECTIVITY_18,
	      CONNECTIVITY_26
} irtkConnectivityType;

class irtkNeighbourhoodOffsets : public irtkObject
{
  irtkObjectMacro(irtkNeighbourhoodOffsets);

  // Set to maximum possible size.
  int _Offsets[26];

  int _Size;

  irtkConnectivityType _Connectivity;

public:
  //
  // Constructors and destructor
  //

  /// Constructor
  irtkNeighbourhoodOffsets();

  /// Constructor with image and connectivity specified
  irtkNeighbourhoodOffsets(irtkBaseImage*, irtkConnectivityType);

  /// Initializer with image and connectivity specified
  void Initialize(irtkBaseImage*, irtkConnectivityType);

  /// Initializer with slice dimensions and connectivity specified
  void Initialize(int, int, irtkConnectivityType);

  /// Default destructor
  virtual ~irtkNeighbourhoodOffsets(void);

  SetMacro(Connectivity, irtkConnectivityType);

  GetMacro(Connectivity, irtkConnectivityType);

  GetMacro(Size, int);

  //
  // Operators
  //

  int  operator()(int) const;

};

inline int irtkNeighbourhoodOffsets::operator()(int pos) const
{
#ifdef NO_BOUNDS
  return this->_Offsets[pos];
#else
  if (pos >= 0 && pos < this->_Size){
    return this->_Offsets[pos];
  } else {
    cout << "irtkNeighbourhoodOffsets::operator(): parameter out of range\n";
    return 0;
  }

#endif
}

#endif /* IRTKNEIGHBOURHOODOFFSETS_H_ */
