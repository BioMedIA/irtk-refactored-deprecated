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

#ifndef _IRTKMAXSTEPLINESEARCH_H
#define _IRTKMAXSTEPLINESEARCH_H

#include <irtkInexactLineSearch.h>


/**
 * Dummy line search which always takes the maximum step
 *
 * This line search implements the LS_None line search strategy.
 */
class irtkMaxStepLineSearch : public irtkInexactLineSearch
{
  irtkLineSearchMacro(irtkMaxStepLineSearch, LS_None);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  irtkMaxStepLineSearch(irtkObjectiveFunction * = NULL);

  /// Copy constructor
  irtkMaxStepLineSearch(const irtkMaxStepLineSearch &);

  /// Assignment operator
  irtkMaxStepLineSearch &operator =(const irtkMaxStepLineSearch &);

  /// Destructor
  virtual ~irtkMaxStepLineSearch();

  // ---------------------------------------------------------------------------
  // Optimization

  /// Make optimal step along search direction
  virtual double Run();

};


#endif
