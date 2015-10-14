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

#ifndef _IRTKOBSERVER_H
#define _IRTKOBSERVER_H

#include <set>
using namespace std;

#include <irtkObject.h>
#include <irtkEvent.h>

class irtkObservable;


/**
 * Observer of an observable object
 */
class irtkObserver : public irtkObject
{
  irtkAbstractMacro(irtkObserver);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
protected:

  /// Default constructor
  irtkObserver();

  /// Copy constructor
  irtkObserver(const irtkObserver &);

  /// Assignment operator
  irtkObserver &operator =(const irtkObserver &);

  /// Destructor
  virtual ~irtkObserver();

  // ---------------------------------------------------------------------------
  // Observation
public:

  /// Stop observing any of the currently monitored observables
  void ClearObservables();

  /// Receives event messages from observed objects
  virtual void HandleEvent(irtkObservable *, irtkEvent, const void * = NULL) = 0;

  // ---------------------------------------------------------------------------
  // Attributes
private:

  set<irtkObservable *> _Observable; ///< Observed objects

  void StartObserving(irtkObservable *);
  void StopObserving (irtkObservable *);

  friend class irtkObservable;
};


#endif
