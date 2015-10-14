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

#ifndef _IRTKOBJECT_H
#define _IRTKOBJECT_H

#include <irtkCxxLib.h>


// =============================================================================
// Basic object interface
// =============================================================================

/// Ordered list of parameter name/value pairs
typedef vector<pair<string, string> >     irtkParameterList;
typedef irtkParameterList::iterator       irtkParameterIterator;
typedef irtkParameterList::const_iterator irtkParameterConstIterator;

/**
 * Base class for all IRTK object classes
 *
 * \note This base class must be a virtual interface without any own data
 *       members to not change the type size of subclasses! Derive another
 *       intermediate abstract base class from it to add data members shared
 *       by a more specific class of objects.
 */
class irtkObject
{
public:

  /// Get name of this class type
  inline static const char *NameOfType() { return "irtkObject"; }

  /// Get name of class, which this object is an instance of
  virtual const char *NameOfClass() const = 0;

  /// Destructor
  virtual ~irtkObject() {}

  /// Set parameter value from string
  virtual bool Set(const char *, const char *) { return false; }

  /// Get parameter name/value pairs
  virtual irtkParameterList Parameter() const { return irtkParameterList(); }

  /// Set parameters from name/value pairs
  inline void Parameter(const irtkParameterList &param) {
    for (irtkParameterConstIterator it = param.begin(); it != param.end(); ++it) {
      this->Set(it->first.c_str(), it->second.c_str());
    }
  }
};

// =============================================================================
// Auxiliary functions for subclass implementation
// =============================================================================

// -----------------------------------------------------------------------------
/// Find parameter in parameters list
inline irtkParameterConstIterator Find(const irtkParameterList &params, string name)
{
  irtkParameterConstIterator it = params.begin();
  while (it != params.end() && it->first != name) ++it;
  return it;
}

// -----------------------------------------------------------------------------
/// Find parameter in parameters list
inline irtkParameterIterator Find(irtkParameterList &params, string name)
{
  irtkParameterIterator it = params.begin();
  while (it != params.end() && it->first != name) ++it;
  return it;
}

// -----------------------------------------------------------------------------
/// Whether parameter is in parameters list
inline bool Contains(const irtkParameterList &params, string name)
{
  return Find(params, name) != params.end();
}

// -----------------------------------------------------------------------------
/// Get parameter value from parameters list
inline string Get(const irtkParameterList &params, string name)
{
  irtkParameterConstIterator pos = Find(params, name);
  if (pos == params.end()) return string("");
  return pos->second;
}

// -----------------------------------------------------------------------------
/// Insert/replace value into/in parameters list
template <class T>
inline irtkParameterList &Insert(irtkParameterList &params, string name, T value)
{
  irtkParameterIterator pos = Find(params, name);
  if (pos == params.end()) {
    params.push_back(make_pair(name, ToString(value)));
  } else {
    pos->second = ToString(value);
  }
  return params;
}

// -----------------------------------------------------------------------------
/// Insert/replace string value into/in parameters list
template <>
inline irtkParameterList &Insert(irtkParameterList &params, string name, const char *value)
{
  irtkParameterIterator pos = Find(params, name);
  if (pos == params.end()) {
    params.push_back(make_pair(name, string(value)));
  } else {
    pos->second = value;
  }
  return params;
}

// -----------------------------------------------------------------------------
/// Insert/replace string value into/in parameters list
template <>
inline irtkParameterList &Insert(irtkParameterList &params, string name, string value)
{
  irtkParameterIterator pos = Find(params, name);
  if (pos == params.end()) {
    params.push_back(make_pair(name, value));
  } else {
    pos->second = value;
  }
  return params;
}

// -----------------------------------------------------------------------------
/// Insert/replace values into/in parameters list
inline irtkParameterList &Insert(irtkParameterList       &params,
                                 const irtkParameterList &other,
                                 const char              *prefix = NULL)
{
  if (prefix) {
    string name;
    for (irtkParameterConstIterator it = other.begin(); it != other.end(); ++it) {
      name    = it->first;
      name[0] = ::tolower(name[0]);
      Insert(params, string(prefix) + " " + name, it->second);
    }
  } else {
    for (irtkParameterConstIterator it = other.begin(); it != other.end(); ++it) {
      Insert(params, it->first, it->second);
    }
  }
  return params;
}

// -----------------------------------------------------------------------------
/// Remove parameter from parameters list
inline irtkParameterList &Remove(irtkParameterList &params, string name)
{
  irtkParameterIterator pos = Find(params, name);
  if (pos != params.end()) params.erase(pos);
  return params;
}

// =============================================================================
// Auxiliary macros for subclass implementation
// =============================================================================

// -----------------------------------------------------------------------------
/// Declare abstract base class derived from irtkObject
#define irtkAbstractMacro(name)                                                \
  public:                                                                      \
    /** Get name of this class type */                                         \
    inline static const char *NameOfType() { return #name; }                   \
    /** Get name of class, which this object is an instance of */              \
    virtual const char *NameOfClass() const = 0;                               \
  private:                                                                     \
    static void _irtkAbstractMacro_needs_trailing_semicolon()

// -----------------------------------------------------------------------------
/// Declare class derived from irtkObject
#define irtkObjectMacro(name)                                                  \
  public:                                                                      \
    /** Get name of this class type */                                         \
    inline static  const char *NameOfType() { return #name; }                  \
    /** Get name of class, which this object is an instance of */              \
    inline virtual const char *NameOfClass() const { return #name; }           \
  private:                                                                     \
    static void _irtkObjectMacro_needs_trailing_semicolon()

// -----------------------------------------------------------------------------
/// Declare class of mutable objects, i.e., ones which define their own
/// NameOfClass implementation that returns a different type identifier
/// depending on the state of the object.
#define irtkMutableObjectMacro(name)                                           \
  public:                                                                      \
    /** Get name of this class type */                                         \
    inline static const char *NameOfType() { return #name; }                   \
    /** Get name of class, which this object is an instance of */              \
    virtual const char *NameOfClass() const;                                   \
  private:                                                                     \
    static void _irtkMutableObjectMacro_needs_trailing_semicolon()

// -----------------------------------------------------------------------------
/// Define setter for class member variable
/// \sa irtkPublicAttributeMacro
#define irtkSetMacro(name, type)                                               \
    virtual void Set##name(type arg) { this->_##name = arg; }                  \
    static void _irtkSetMacro_##name##_needs_trailing_semicolon()

// -----------------------------------------------------------------------------
/// Define getter for class member variable
/// \sa irtkPublicAttributeMacro, irtkReadOnlyAttributeMacro
#define irtkGetMacro(name, type)                                               \
    type Get##name() const { return this->_##name; }                           \
    static void _irtkGetMacro_##name##_needs_trailing_semicolon()

// -----------------------------------------------------------------------------
/// Define VTK-like On/Off setter for boolean class member variable
/// \sa irtkPublicAttributeMacro
#define irtkOnOffMacro(name)                                                   \
    virtual void name##On()  { this->_##name = true;  }                        \
    virtual void name##Off() { this->_##name = false; }                        \
    static void _irtkOnOffMacro_##name##_needs_trailing_semicolon()

// -----------------------------------------------------------------------------
/// Define read-only class attribute and corresponding accessors
#define irtkDefineReadOnlyAttributeMacro(typedecl, type, name)                 \
  protected:                                                                   \
    typedecl _##name;                                                          \
  public:                                                                      \
    /** Get value of _##name attribute */                                      \
    inline type &name() { return _##name; }                                    \
    /** Get value of _##name attribute */                                      \
    inline const type &name() const { return _##name; }                        \
  private:                                                                     \
    /* require developer to end macro with a semicolon */                      \
    static void _irtkAttributeMacro_##name##_needs_trailing_semicolon()

// -----------------------------------------------------------------------------
/// Define class attribute and corresponding accessors
#define irtkDefineAttributeMacro(access, typedecl, type, name)                 \
  protected:                                                                   \
    typedecl _##name;                                                          \
  access:                                                                      \
    /** Set value of _##name attribute */                                      \
    inline virtual void name(type arg) { _##name = arg; }                      \
    /** Get value of _##name attribute */                                      \
    inline type &name() { return _##name; }                                    \
    /** Get value of _##name attribute */                                      \
    inline const type &name() const { return _##name; }                        \
  private:                                                                     \
    /* require developer to end macro with a semicolon */                      \
    static void _irtkAttributeMacro_##name##_needs_trailing_semicolon()

// -----------------------------------------------------------------------------
/// Define public attribute
#define irtkPublicAttributeMacro(type, name)                                   \
  irtkDefineAttributeMacro(public, type, type, name)
/// Define public read-only attribute
#define irtkReadOnlyAttributeMacro(type, name)                                 \
  irtkDefineReadOnlyAttributeMacro(type, type, name)
/// Define public mutable attribute
#define irtkPublicMutableAttributeMacro(type, name)                            \
  irtkDefineAttributeMacro(public, mutable type, type, name)

// -----------------------------------------------------------------------------
/// Define protected attribute
#define irtkAttributeMacro(type, name)                                         \
  irtkDefineAttributeMacro(protected, type, type, name)
/// Define protected mutable attribute
#define irtkMutableAttributeMacro(type, name)                                  \
  irtkDefineAttributeMacro(protected, mutable type, type, name)

// -----------------------------------------------------------------------------
/// Define pointer to aggregate (cf. UML aggregation) and corresponding accessors
#define irtkDefineAggregateMacro(access, type, name)                           \
  protected:                                                                   \
    type *_##name;                                                             \
  access:                                                                      \
    /** Set value of _##name attribute */                                      \
    inline virtual void name(type *arg) { _##name = arg; }                     \
    /** Get value of _##name attribute */                                      \
    inline type *name() const { return _##name; }                              \
  private:                                                                     \
    /* require developer to end macro with a semicolon */                      \
    static void _irtkAggregateMacro_##name##_needs_trailing_semicolon()

// -----------------------------------------------------------------------------
/// Define pointer to aggregate (cf. UML aggregation) and corresponding accessors
#define irtkDefineReadOnlyAggregateMacro(access, type, name)                   \
  protected:                                                                   \
    type *_##name;                                                             \
  access:                                                                      \
    /** Get value of _##name attribute */                                      \
    inline type *name() const { return _##name; }                              \
  private:                                                                     \
    /* require developer to end macro with a semicolon */                      \
    static void _irtkReadOnlyAggregateMacro_##name##_needs_trailing_semicolon()

// -----------------------------------------------------------------------------
/// Define public pointer to aggregate (cf. UML aggregation)
#define irtkPublicAggregateMacro(type, name)                                   \
  irtkDefineAggregateMacro(public, type, name)
/// Define public read-only pointer to aggregate (cf. UML aggregation)
#define irtkReadOnlyAggregateMacro(type, name)                                 \
  irtkDefineReadOnlyAggregateMacro(public, type, name)

// -----------------------------------------------------------------------------
/// Define protected pointer to aggregate (cf. UML aggregation)
#define irtkAggregateMacro(type, name)                                         \
  irtkDefineAggregateMacro(protected, type, name)

// -----------------------------------------------------------------------------
/// Define pointer to component (cf. UML composition) and corresponding accessors
#define irtkDefineComponentMacro(access, type, name)                           \
  protected:                                                                   \
    type *_##name;                                                             \
  access:                                                                      \
    /** Set pointer to _##name attribute */                                    \
    inline virtual void name(type *arg) { delete _##name; _##name = arg; }     \
    /** Get pointer to _##name attribute */                                    \
    inline type *name() { return _##name; }                                    \
    /** Get const pointer to _##name attribute */                              \
    inline const type *name() const { return _##name; }                        \
  private:                                                                     \
    /* require developer to end macro with a semicolon */                      \
    static void _irtkComponentMacro_##name##_needs_trailing_semicolon()

// -----------------------------------------------------------------------------
/// Define public pointer to component (cf. UML composition)
#define irtkPublicComponentMacro(type, name)                                   \
  irtkDefineComponentMacro(public, type, name)
/// Define public read-only pointer to component (cf. UML composition)
#define irtkReadOnlyComponentMacro(type, name)                                 \
  irtkDefineReadOnlyAggregateMacro(public, type, name)

// -----------------------------------------------------------------------------
/// Define protected pointer to component (cf. UML composition)
#define irtkComponentMacro(type, name)                                         \
  irtkDefineComponentMacro(protected, type, name)

// -----------------------------------------------------------------------------
/// \deprecated Use irtkSetMacro or irtkPublicAttributeMacro et al. instead.
#define SetMacro(name, type) void Set##name(type arg) { this->_##name = arg; }

// -----------------------------------------------------------------------------
/// \deprecated Use irtkGetMacro or irtkReadOnlyAttributeMacro et al. instead.
#define GetMacro(name, type) irtkGetMacro(name, type);


#endif

