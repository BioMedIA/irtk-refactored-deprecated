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

#ifndef _IRTKOBSERVABLE_H
#define _IRTKOBSERVABLE_H

#include <irtkObject.h>
#include <irtkObserver.h>
#include <irtkEvent.h>


/**
 * Base class for an observable object
 *
 * Any object which should be observable must derive from irtkObservable
 * instead of irtkObject. A client which wants to observe events emitted
 * by an observable must use an instance of irtkObservable and register
 * callback functions at this instance for the events of interest.
 *
 * \attention Adding and removing observers is not thread-safe!
 */
class irtkObservable : public irtkObject
{
  irtkAbstractMacro(irtkObservable);

protected:

  /// Default constructor
  irtkObservable();

  /// Copy constructor
  irtkObservable(const irtkObservable &);

  /// Assignment operator
  irtkObservable &operator =(const irtkObservable &);

public:

  /// Destructor
  virtual ~irtkObservable();

  /// Number of current observers
  int NumberOfObservers() const;

  /// Add observer
  void AddObserver(irtkObserver &);

  /// Delete observer
  void DeleteObserver(irtkObserver &);

  /// Delete all observers
  void ClearObservers();

  /// Broadcast event to observers
  void Broadcast(irtkEvent, const void * = NULL);

  /// Notify all observers about given event if this object has changed
  void NotifyObservers(irtkEvent, const void * = NULL);

  /// Set whether object has changed and should notify observers upon request
  virtual void Changed(bool);

  /// Whether this object has changed and will notify observers upon request
  bool Changed() const;

private:

  set<irtkObserver *> _Observer; ///< Registered observers
  bool                _Changed;  ///< Whether this object has changed

};

////////////////////////////////////////////////////////////////////////////////
// Auxiliary macro for broadcasting log messages
////////////////////////////////////////////////////////////////////////////////

/// Broadcast a LogEvent message
#define IRTK_LOG_EVENT(msg) \
  do { \
    ostringstream ss; \
    ss << msg; \
    Broadcast(LogEvent, ss.str().c_str()); \
  } while(false)

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
inline irtkObservable::irtkObservable()
:
  irtkObject(),
  _Changed(false)
{
}

// -----------------------------------------------------------------------------
inline irtkObservable::irtkObservable(const irtkObservable &other)
:
  irtkObject(other),
  _Changed(false)
  // Do not copy observers!
{
}

// -----------------------------------------------------------------------------
inline irtkObservable &irtkObservable::operator =(const irtkObservable &)
{
  _Changed = true;
  // Do not copy observers!
  return *this;
}

// -----------------------------------------------------------------------------
inline void irtkObservable::ClearObservers()
{
  for (set<irtkObserver *>::iterator it = _Observer.begin(); it != _Observer.end(); ++it) {
    (*it)->StopObserving(this);
  }
  _Observer.clear();
}

// -----------------------------------------------------------------------------
inline irtkObservable::~irtkObservable()
{
  ClearObservers();
}

// =============================================================================
// Add/Remove observer
// =============================================================================

// -----------------------------------------------------------------------------
inline int irtkObservable::NumberOfObservers() const
{
  return static_cast<int>(_Observer.size());
}

// -----------------------------------------------------------------------------
inline void irtkObservable::AddObserver(irtkObserver &observer)
{
  _Observer.insert(&observer);
  observer.StartObserving(this);
}

// -----------------------------------------------------------------------------
inline void irtkObservable::DeleteObserver(irtkObserver &observer)
{
  observer.StopObserving(this);
  _Observer.erase(&observer);
}

// =============================================================================
// Events
// =============================================================================

// -----------------------------------------------------------------------------
inline void irtkObservable::Broadcast(irtkEvent event, const void *data)
{
  for (set<irtkObserver *>::iterator it = _Observer.begin(); it != _Observer.end(); ++it) {
    (*it)->HandleEvent(this, event, data);
  }
}

// -----------------------------------------------------------------------------
inline void irtkObservable::NotifyObservers(irtkEvent event, const void *data)
{
  if (_Changed) {
    Broadcast(event, data);
    _Changed = false;
  }
}

// -----------------------------------------------------------------------------
inline void irtkObservable::Changed(bool changed)
{
  _Changed = changed;
}

// -----------------------------------------------------------------------------
inline bool irtkObservable::Changed() const
{
  return _Changed;
}


#endif
