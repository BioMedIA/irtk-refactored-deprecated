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

#ifndef _IRTKPOINTSET_H
#define _IRTKPOINTSET_H

#define POINTSET_SIZE 4096

#include <irtkPoint.h>
#ifdef HAS_VTK
class vtkAbstractArray;
#endif


/**

  Point set class.

*/

class irtkPointSet : public irtkObject
{
  irtkObjectMacro(irtkPointSet);

protected:

  /// Allocated size of irtkPointSet
  int _m;

  /// Actual size of irtkPointSet
  int _n;

  /// Pointer to Points
  irtkPoint *_data;

public:

  //
  // Constructors and destructor
  //

  /// Default constructor
  irtkPointSet(int = 0);

  /// Copy constructor
  irtkPointSet(const irtkPointSet &);

  /// Destructor
  virtual ~irtkPointSet();

  //
  // Size get acccessor and clearing of irtkPointSet
  //

  /// Resizes container so it contains n elements
  void Resize(int, const irtkPoint & = irtkPoint());

  /// Request a change in capacity
  void Reserve(int);

  /// Return size of allocated storage capacity
  int Capacity() const;

  /// Set container size (and capacity)
  ///
  /// The result of this function is equivalent to
  /// \code
  /// Resize(n);
  /// ShrinkToFit();
  /// \endcode
  /// but with only one reallocation operation. It sets both the size and
  /// the capacity of the point set container.
  void Size(int);

  /// Access function for size
  int Size() const;

  /// Requests the container to reduce its capacity to fit its size
  void ShrinkToFit();

  /// Clearing of irtkPointSet
  void Clear();

  //
  // Operators for access
  //

  /// Get n-th point
  irtkPoint &GetPoint(int);

  /// Get n-th point
  const irtkPoint &GetPoint(int) const;

  /// Get n-th point
  void GetPoint(int, irtkPoint &) const;

  /// Get n-th point
  void GetPoint(int, double [3]) const;

  /// Set n-th point
  void SetPoint(int, const irtkPoint &);

  /// Set n-th point
  void SetPoint(int, const double [3]);

  /// Operator for Point put access
  irtkPoint &operator()(int);

  /// Operator for Point get access
  const irtkPoint &operator()(int) const;

  /// Operator
  irtkPointSet operator()(int, int) const;

  //
  // Operators for irtkPointSet arithmetic
  //

  // irtkPointSet operation for =
  irtkPointSet &operator=(const irtkPointSet &);

  // irtkPointSet operator for Point adding
  irtkPointSet& operator+=(const irtkPoint&);

  // irtkPointSet operator for Point substraction
  irtkPointSet& operator-=(const irtkPoint&);

  // irtkPointSet operator for irtkPointSet adding
  irtkPointSet& operator+=(const irtkPointSet&);

  // irtkPointSet operator for irtkPointSet substraction
  irtkPointSet& operator-=(const irtkPointSet&);

  /// Centre of gravity
  virtual irtkPoint CenterOfGravity() const;

  /// Closest point to given point
  virtual irtkPoint ClosestPoint(irtkPoint&);

  /// Point set distance to given point
  virtual double PointDistance(irtkPoint&);

  /// Bounding box
  virtual void BoundingBox(irtkPoint &, irtkPoint &) const;

  //
  // Operators for irtkPointSet comparison
  //

  /// Test for equality
  bool operator==(const irtkPointSet &) const;

  /// Test for inequality
  bool operator!=(const irtkPointSet &) const;

  //
  // Explicit methods to add or delete points
  //

  /// Adding of a Point to Pointset
  void Add(const irtkPoint&);

  /// Deleting of a Point from Pointset
  void Del(const irtkPoint&);

  /// Adding of a irtkPointSet to Pointset
  void Add(const irtkPointSet &);

  /// Deleting of a irtkPointSet from Pointset
  void Del(const irtkPointSet &);

  /// Adding of a Point to Pointset
  void Add(double *);

  /// Deleting of a Point from Pointset
  void Del(double *);

  //
  // irtkPointSet in- and output
  //

  /// Interface to output stream
  friend ostream& operator<< (ostream&, const irtkPointSet&);

  /// Interface to input stream
  friend istream& operator>> (istream&, irtkPointSet&);

  /// Read pointset from file
  void Read(const char *);

  /// Write pointset to file
  void Write(const char *) const;

  /// Read pointset from file in VTK format
  void ReadVTK(const char *);

  /// Add pointset from file in VTK format
  void AddVTK(const char *);

  /// Write pointset to file in VTK format
#ifdef HAS_VTK
  void WriteVTK(const char *, vtkAbstractArray * = NULL) const;
#else
  void WriteVTK(const char *, void * = NULL) const;
#endif

  //
  // Misc. functions
  //

  /// Computes the standard deviation ellipsoid about the centroid of a point set.
  irtkPoint StandardDeviationEllipsoid() const;

  /// Tests if a point is inside the polygon defined by the point set
  int IsInside(double, double) const;

};

inline irtkPointSet::irtkPointSet(int n)
:
  _m(n), _n(n), _data(Allocate<irtkPoint>(n))
{
}

inline irtkPointSet::irtkPointSet(const irtkPointSet &pset)
:
  irtkObject(pset),
  _m(0), _n(0), _data(NULL)
{
  (*this) = pset;
}

inline irtkPointSet::~irtkPointSet()
{
  Clear();
}

inline irtkPoint& irtkPointSet::operator()(int j)
{
#ifdef NO_BOUNDS
  return _data[j];
#else
  if ((j >= 0) && (j < _n)) {
    return _data[j];
  } else {
    cerr << "irtkPointSet::operator(int) parameter out of range\n";
    exit(1);
  }
#endif
}

inline const irtkPoint &irtkPointSet::operator()(int j) const
{
#ifdef NO_BOUNDS
  return _data[j];
#else
  if ((j >= 0) && (j < _n)) {
    return _data[j];
  } else {
    cerr << "irtkPointSet::operator(int) parameter out of range\n";
    exit(1);
  }
#endif
}

inline irtkPoint &irtkPointSet::GetPoint(int i)
{
  return this->operator()(i);
}

inline const irtkPoint &irtkPointSet::GetPoint(int i) const
{
  return this->operator()(i);
}

inline void irtkPointSet::GetPoint(int i, irtkPoint &p) const
{
  p = this->operator()(i);
}

inline void irtkPointSet::GetPoint(int i, double p[3]) const
{
  const irtkPoint &pt = this->operator()(i);
  p[0] = pt._x, p[1] = pt._y, p[2] = pt._z;
}

inline void irtkPointSet::SetPoint(int i, const irtkPoint &p)
{
  this->operator()(i) = p;
}

inline void irtkPointSet::SetPoint(int i, const double p[3])
{
  irtkPoint &pt = this->operator()(i);
  pt._x = p[0], pt._y = p[1], pt._z = p[2];
}

inline irtkPointSet& irtkPointSet::operator+=(const irtkPoint &p)
{
  this->Add(p);
  return *this;
}

inline irtkPointSet& irtkPointSet::operator-=(const irtkPoint &p)
{
  this->Del(p);
  return *this;
}

inline irtkPointSet& irtkPointSet::operator+=(const irtkPointSet &pset)
{
  this->Add(pset);
  return *this;
}

inline irtkPointSet& irtkPointSet::operator-=(const irtkPointSet &pset)
{
  this->Del(pset);
  return *this;
}

inline irtkPointSet irtkPointSet::operator()(int j, int k) const
{
  irtkPointSet pset;
#ifdef NO_BOUNDS
  for (int i = j; i < k; j++) {
    pset += _data[i];
  }
#else
  if ((j >= 0) && (k < _n)) {
    for (int i = j; i < k; j++) {
      pset += _data[i];
    }
  } else {
    cerr << "irtkPointSet::operator(int, int) parameter out of range\n";
    exit(1);
  }
#endif
  return pset;
}

inline void irtkPointSet::Resize(int n, const irtkPoint &p)
{
  int m = _m;
  while (n > m) m += POINTSET_SIZE;
  this->Reserve(m);
  _n = n;
  for (int i = _n; i < _m; ++i) _data[i] = p;
}

inline void irtkPointSet::Reserve(int m)
{
  if (_m < m) {
    _m = m;
    irtkPoint *new_data = Allocate<irtkPoint>(m);
    for (int i = 0; i < _n; ++i) new_data[i] = _data[i];
    Deallocate(_data);
    _data = new_data;
  }
}

inline int irtkPointSet::Capacity() const
{
  return _m;
}

inline void irtkPointSet::Size(int n)
{
  if (_m != n) {
    irtkPoint *new_data = Allocate<irtkPoint>(n);
    for (int i = 0; i < _n; ++i) new_data[i] = _data[i];
    Deallocate(_data);
    _data = new_data;
    _m = _n = n;
  }
}

inline int irtkPointSet::Size() const
{
  return _n;
}

inline void irtkPointSet::ShrinkToFit()
{
  if (_n < _m) {
    irtkPoint *new_data = Allocate<irtkPoint>(_n);
    for (int i = 0; i < _n; ++i) new_data[i] = _data[i];
    Deallocate(_data);
    _data = new_data;
    _m    = _n;
  }
}

inline void irtkPointSet::Add(double *p)
{
  this->Add(irtkPoint(p[0], p[1], p[2]));
}

inline void irtkPointSet::Del(double *p)
{
  this->Del(irtkPoint(p[0], p[1], p[2]));
}

#endif

