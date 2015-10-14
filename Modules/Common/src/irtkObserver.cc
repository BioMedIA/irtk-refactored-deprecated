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

#include <irtkCommon.h>
#include <irtkObserver.h>
#include <irtkObservable.h>


// -----------------------------------------------------------------------------
irtkObserver::irtkObserver()
:
  irtkObject()
{
}

// -----------------------------------------------------------------------------
irtkObserver::irtkObserver(const irtkObserver &other)
:
  irtkObject(other)
{
  set<irtkObservable *>::iterator it = other._Observable.begin();
  while (it != other._Observable.end()) {
    (*it)->AddObserver(*this);
    ++it;
  }
}

// -----------------------------------------------------------------------------
irtkObserver &irtkObserver::operator =(const irtkObserver &other)
{
  set<irtkObservable *>::iterator it = other._Observable.begin();
  while (it != other._Observable.end()) {
    (*it)->AddObserver(*this);
    ++it;
  }
  return *this;
}

// -----------------------------------------------------------------------------
void irtkObserver::ClearObservables()
{
  while (!_Observable.empty()) {
    irtkObservable *observable = *_Observable.begin();
    // Calls this->StopObserving(observable) which removes the observable
    // and thus would invalidate any iterators to _Observable
    observable->DeleteObserver(*this);
  }
}

// -----------------------------------------------------------------------------
irtkObserver::~irtkObserver()
{
  ClearObservables();
}

// -----------------------------------------------------------------------------
void irtkObserver::StartObserving(irtkObservable *obj)
{
  if (_Observable.find(obj) == _Observable.end()) {
    _Observable.insert(obj);
    this->HandleEvent(obj, RegisteredEvent);
  }
}

// -----------------------------------------------------------------------------
void irtkObserver::StopObserving(irtkObservable *obj)
{
  if (_Observable.find(obj) != _Observable.end()) {
    _Observable.erase(obj);
    this->HandleEvent(obj, UnregisteredEvent);
  }
}
