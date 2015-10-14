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

#ifndef _IRTKEVENTDELEGATE_H
#define _IRTKEVENTDELEGATE_H

#include <irtkObserver.h>
#include "FastDelegate.h"
#include <map>


/**
 * Invokes a callback function for observed events
 *
 * This class can be registered as observer of an observable object derived
 * from irtkObservable. It receives the event notifications from the observed
 * object(s) and dispatches the corresponding registered callback (member)
 * function (if any). It therefore uses instances of fastdelegate::FastDelegate.
 *
 * \note Do not derive from this class. Instead add instances of it as attributes
 *       to classes which observe other objects and handle the events.
 *
 * \code
 * irtkEventDelegate delegate;
 * // Bind non-member function to any event
 * delegate.Bind(&function);
 * // Bind static member function to EndEvent
 * delegate.Bind(EndEvent, &function);
 * // Bind non-static member function to IterationEvent
 * delegate.Bind(IterationEvent, MakeDelegate(&obj, &method));
 * \endcode
 */
class irtkEventDelegate : public irtkObserver
{
  irtkObjectMacro(irtkEventDelegate);

  // ---------------------------------------------------------------------------
  // Types

  typedef fastdelegate::FastDelegate0<>                                          Delegate0;
  typedef fastdelegate::FastDelegate1<irtkEvent>                                 Delegate1;
  typedef fastdelegate::FastDelegate2<irtkEvent, const void *>                   Delegate2;
  typedef fastdelegate::FastDelegate3<irtkObservable *, irtkEvent, const void *> Delegate3;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Copy constructor
  /// \note Intentionally not implemented
  irtkEventDelegate(const irtkEventDelegate &);

  /// Assignment operator
  /// \note Intentionally not implemented
  irtkEventDelegate &operator =(const irtkEventDelegate &);

public:

  /// Default constructor
  irtkEventDelegate();

  /// Destructor
  ~irtkEventDelegate();

  // ---------------------------------------------------------------------------
  // Bind callback functions

  /// Bind callback (member) function for any event with no arguments
  void Bind(Delegate0);

  /// Bind callback (member) function for any event with one argument
  void Bind(Delegate1);

  /// Bind callback (member) function for any event with two arguments
  void Bind(Delegate2);

  /// Bind callback (member) function for any event with four arguments
  void Bind(Delegate3);

  /// Bind callback (member) function to specified event
  void Bind(irtkEvent, Delegate0);

  /// Bind callback (member) function to specified event
  void Bind(irtkEvent, Delegate1);

  /// Bind callback (member) function to specified event
  void Bind(irtkEvent, Delegate2);

  /// Bind callback (member) function to specified event
  void Bind(irtkEvent, Delegate3);

  // ---------------------------------------------------------------------------
  // Event forwarding

  /// Forward event notification to registered callback (member) function
  void HandleEvent(irtkObservable *, irtkEvent, const void * = NULL);

  // ---------------------------------------------------------------------------
  // Attributes
private:

  /// Auxiliary structure to store disperse callback delegates
  struct Callback
  {
    int                           _N;
    fastdelegate::DelegateMemento _Memento;
    Callback() : _N(0), _Memento() {}
    Callback(fastdelegate::DelegateMemento memento, int n) : _N(n), _Memento(memento) {}
    Callback(const Callback &o) : _N(o._N), _Memento(o._Memento) {}
  };

  void Forward(Callback &func, irtkObservable *obj, irtkEvent event, const void *data);

  std::map<irtkEvent, Callback> _SingleEventDelegate; ///< Event delegates
  Callback                      _AnyEventDelegate;    ///< AnyEvent delegate
};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline irtkEventDelegate::irtkEventDelegate()
:
  irtkObserver()
{
}

// -----------------------------------------------------------------------------
inline irtkEventDelegate::~irtkEventDelegate()
{
}

// -----------------------------------------------------------------------------
inline void irtkEventDelegate::Bind(Delegate0 delegate)
{
  _AnyEventDelegate = Callback(delegate.GetMemento(), 0);
}

// -----------------------------------------------------------------------------
inline void irtkEventDelegate::Bind(Delegate1 delegate)
{
  _AnyEventDelegate = Callback(delegate.GetMemento(), 1);
}

// -----------------------------------------------------------------------------
inline void irtkEventDelegate::Bind(Delegate2 delegate)
{
  _AnyEventDelegate = Callback(delegate.GetMemento(), 2);
}

// -----------------------------------------------------------------------------
inline void irtkEventDelegate::Bind(Delegate3 delegate)
{
  _AnyEventDelegate = Callback(delegate.GetMemento(), 3);
}

// -----------------------------------------------------------------------------
inline void irtkEventDelegate::Bind(irtkEvent event, Delegate0 delegate)
{
  _SingleEventDelegate[event] = Callback(delegate.GetMemento(), 0);
}

// -----------------------------------------------------------------------------
inline void irtkEventDelegate::Bind(irtkEvent event, Delegate1 delegate)
{
  _SingleEventDelegate[event] = Callback(delegate.GetMemento(), 1);
}

// -----------------------------------------------------------------------------
inline void irtkEventDelegate::Bind(irtkEvent event, Delegate2 delegate)
{
  _SingleEventDelegate[event] = Callback(delegate.GetMemento(), 2);
}

// -----------------------------------------------------------------------------
inline void irtkEventDelegate::Bind(irtkEvent event, Delegate3 delegate)
{
  _SingleEventDelegate[event] = Callback(delegate.GetMemento(), 3);
}

// -----------------------------------------------------------------------------
inline void irtkEventDelegate::Forward(Callback &func, irtkObservable *obj, irtkEvent event, const void *data)
{
  switch (func._N) {
    case 0: {
      Delegate0 delegate;
      delegate.SetMemento(func._Memento);
      if (delegate) delegate();
      break;
    }
    case 1: {
      Delegate1 delegate;
      delegate.SetMemento(func._Memento);
      if (delegate) delegate(event);
      break;
    }
    case 2: {
      Delegate2 delegate;
      delegate.SetMemento(func._Memento);
      if (delegate) delegate(event, data);
      break;
    }
    case 3: {
      Delegate3 delegate;
      delegate.SetMemento(func._Memento);
      if (delegate) delegate(obj, event, data);
      break;
    }
  }
}

// -----------------------------------------------------------------------------
inline void irtkEventDelegate::HandleEvent(irtkObservable *obj, irtkEvent event, const void *data)
{
  Forward(_AnyEventDelegate, obj, event, data);
  std::map<irtkEvent, Callback>::iterator it = _SingleEventDelegate.find(event);
  if (it != _SingleEventDelegate.end()) Forward(it->second, obj, event, data);
}


#endif
