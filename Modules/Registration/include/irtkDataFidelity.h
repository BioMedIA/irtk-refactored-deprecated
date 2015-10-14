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

#ifndef _IRTKDATAFIDELITY_H
#define _IRTKDATAFIDELITY_H

#include <irtkEnergyTerm.h>


/**
 * Base class for energy terms measuring the amount of data fidelity
 *
 * Lower data fidelity corresponds to a better alignment of the registered
 * data. The optimizer should minimize the data fidelity.
 */
class irtkDataFidelity : public irtkEnergyTerm
{
  irtkAbstractMacro(irtkDataFidelity);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
protected:

  /// Constructor
  irtkDataFidelity(const char * = "", double = 1.0);

  /// Copy constructor
  irtkDataFidelity(const irtkDataFidelity &);

  /// Assignment operator
  irtkDataFidelity &operator =(const irtkDataFidelity &);

public:

  /// Destructor
  virtual ~irtkDataFidelity();

  // ---------------------------------------------------------------------------
  // Parameters

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

};


#endif
