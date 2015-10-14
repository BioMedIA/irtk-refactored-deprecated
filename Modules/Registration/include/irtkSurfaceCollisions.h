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

#ifndef _IRTKSURFACECOLLISIONS_H
#define _IRTKSURFACECOLLISIONS_H

#include <irtkObject.h>
#include <vtkSmartPointer.h>

class vtkPolyData;
class vtkDataArray;


/**
 * Auxiliary class used to find self-collisions of a triangulated surface mesh
 *
 * Instances of this class detect different types of self-collisions such as
 * in particular self-intersections between non-adjacent as well as adjacent
 * triangular faces and a list of non-adjacent faces which are about to collide,
 * i.e., very close to each other. They are used to impose either hard or soft
 * non-self-intersection constraints on a deformable surface.
 */
class irtkSurfaceCollisions : public irtkObject
{
  irtkObjectMacro(irtkSurfaceCollisions);

  // ---------------------------------------------------------------------------
  // Types
public:

  /// Type of self-collision
  enum CollisionType
  {
    NoCollision,          ///< No collision with another triangle within search range
    Collision,            ///< Unknown or mixed collision types
    FrontfaceCollision,   ///< Second triangle is within minimum distance towards the front of first triangle
    BackfaceCollision,    ///< Second triangle is within minimum distance towards the back of first triangle
    Intersection,         ///< Unknown or mixed intersection types
    SelfIntersection,     ///< Non-adjacent triangles intersect each other
    AdjacentIntersection, ///< Adjacent triangles intersect each other
    Ambiguous             ///< Both collisions and self-intersections found
  };

  /// Whether a given collision type indicates a near miss collision
  static bool IsCollision(CollisionType);

  /// Whether a given collision type indicates a self-intersection
  static bool IsIntersection(CollisionType);

  /// Structure storing information about detected self-intersection
  struct IntersectionInfo
  {
    int  _CellId;   ///< ID of other triangle
    bool _Adjacent; ///< Whether triangles are adjacent

    IntersectionInfo(int cellId = -1, bool adj = false)
    :
      _CellId(cellId), _Adjacent(adj)
    {}

    IntersectionInfo(const IntersectionInfo &other)
    :
      _CellId(other._CellId), _Adjacent(other._Adjacent)
    {}

    IntersectionInfo &operator =(const IntersectionInfo &other)
    {
      _CellId   = other._CellId;
      _Adjacent = other._Adjacent;
      return *this;
    }

    bool operator <(const IntersectionInfo &rhs) const
    {
      return _CellId < rhs._CellId;
    }
  };

  /// Structure storing information about detected collision
  struct CollisionInfo
  {
    int           _CellId;    ///< ID of other triangle
    double        _Point1[3]; ///< Point in this  triangle which is closest
    double        _Point2[3]; ///< Point in other triangle which is closest
    double        _Distance;  ///< Distance between closest points
    CollisionType _Type;      ///< Type of collision

    CollisionInfo(int cellId = -1)
    :
      _CellId(cellId),
      _Distance(numeric_limits<double>::quiet_NaN()),
      _Type(Collision)
    {
      memset(_Point1, 0, 3 * sizeof(double));
      memset(_Point2, 0, 3 * sizeof(double));
    }

    CollisionInfo(int cellId, const double p1[3], const double p2[3], CollisionType type = Collision)
    :
      _CellId(cellId),
      _Type(type)
    {
      memcpy(_Point1, p1, 3 * sizeof(double));
      memcpy(_Point2, p2, 3 * sizeof(double));
      const double a = _Point2[0] - _Point1[0];
      const double b = _Point2[1] - _Point1[1];
      const double c = _Point2[2] - _Point1[2];
      _Distance = sqrt(a * a + b * b + c * c);
    }

    CollisionInfo(int cellId, const double p1[3], const double p2[3], double dist, CollisionType type = Collision)
    :
      _CellId(cellId),
      _Distance(dist),
      _Type(type)
    {
      memcpy(_Point1, p1, 3 * sizeof(double));
      memcpy(_Point2, p2, 3 * sizeof(double));
    }

    CollisionInfo(const CollisionInfo &other)
    :
      _CellId(other._CellId),
      _Distance(other._Distance),
      _Type(other._Type)
    {
      memcpy(_Point1, other._Point1, 3 * sizeof(double));
      memcpy(_Point2, other._Point2, 3 * sizeof(double));
    }

    CollisionInfo &operator =(const CollisionInfo &other)
    {
      _CellId   = other._CellId;
      _Distance = other._Distance;
      _Type     = other._Type;
      memcpy(_Point1, other._Point1, 3 * sizeof(double));
      memcpy(_Point2, other._Point2, 3 * sizeof(double));
      return *this;
    }

    bool operator <(const CollisionInfo &rhs) const
    {
      return _CellId < rhs._CellId;
    }
  };

  // ---------------------------------------------------------------------------
  // Attributes
private:

  /// Triangulated input surface mesh
  irtkPublicAttributeMacro(vtkSmartPointer<vtkPolyData>, Input);

  /// Annotated output surface mesh
  irtkReadOnlyAttributeMacro(vtkSmartPointer<vtkPolyData>, Output);

  /// Minimum search radius around triangle center used to determine the set
  /// of nearby triangles to be tested for collisions with this triangle
  irtkPublicAttributeMacro(double, MinSearchRadius);

  /// Maximum search radius around triangle center used to determine the set
  /// of nearby triangles to be tested for collisions with this triangle
  irtkPublicAttributeMacro(double, MaxSearchRadius);

  /// Minimum required distance between non-adjacent triangular faces
  irtkPublicAttributeMacro(double, MinDistance);

  /// Maximum angle between face normal and center to closest point vector
  /// required for collision to be detected
  irtkPublicAttributeMacro(double, MaxAngle);

  /// Whether to perform intersection tests between adjacent triangles
  irtkPublicAttributeMacro(bool, AdjacentIntersectionTest);

  /// Whether to perform intersection tests between non-adjacent triangles
  irtkPublicAttributeMacro(bool, NonAdjacentIntersectionTest);

  /// Whether to detect collisions with the front of a face
  irtkPublicAttributeMacro(bool, FrontfaceCollisionTest);

  /// Whether to detect collisions with the back of a face
  irtkPublicAttributeMacro(bool, BackfaceCollisionTest);

  /// IDs of cells for which self-intersections were found
  irtkReadOnlyAttributeMacro(set<int>, IntersectionCells);

  /// IDs of cells for which collisions were found
  irtkReadOnlyAttributeMacro(set<int>, CollisionCells);

  /// Total number of self-intersections
  irtkReadOnlyAttributeMacro(int, NumberOfIntersections);

  /// Total number of collisions
  irtkReadOnlyAttributeMacro(int, NumberOfCollisions);

  /// Found self-intersections per face
  vector<set<IntersectionInfo> > _Intersections;

  /// Found collisions per face
  vector<set<CollisionInfo> > _Collisions;

  // ---------------------------------------------------------------------------
  // Construction/destruction
private:

  /// Copy constructor -- not implemented
  irtkSurfaceCollisions(const irtkSurfaceCollisions &);

  /// Assignment operator -- not implemented
  irtkSurfaceCollisions &operator =(const irtkSurfaceCollisions &);

public:

  /// Constructor
  irtkSurfaceCollisions();

  /// Destructor
  virtual ~irtkSurfaceCollisions();

  // ---------------------------------------------------------------------------
  // Execution
protected:

  /// Initialize filter execution
  void Initialize();

  /// Finalize filter execution
  void Finalize();

public:

  /// Detect self-collisions of input surface
  void Run();

  // ---------------------------------------------------------------------------
  // Output
public:

  /// Get cell data array storing radii of bounding spheres
  vtkDataArray *GetRadiusArray() const;

  /// Get cell data array storing center points of bounding spheres
  vtkDataArray *GetCenterArray() const;

  /// Get cell data array storing type of cell collision
  vtkDataArray *GetCollisionTypeArray() const;

  /// Set of intersections of other faces with the specified cell
  const set<IntersectionInfo> &Intersections(int) const;

  /// Set of collisions of other faces with the specified cell
  const set<CollisionInfo> &Collisions(int) const;

  /// Number of self-intersections with specified cell
  int NumberOfIntersections(int) const;

  /// Number of collisions with specified cell
  int NumberOfCollisions(int) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline bool irtkSurfaceCollisions::IsCollision(CollisionType type)
{
  return Collision <= type && type < Intersection;
}

// -----------------------------------------------------------------------------
inline bool irtkSurfaceCollisions::IsIntersection(CollisionType type)
{
  return Intersection <= type && type < Ambiguous;
}

// -----------------------------------------------------------------------------
inline const set<irtkSurfaceCollisions::IntersectionInfo> &
irtkSurfaceCollisions::Intersections(int cellId) const
{
  return _Intersections[cellId];
}

// -----------------------------------------------------------------------------
inline const set<irtkSurfaceCollisions::CollisionInfo> &
irtkSurfaceCollisions::Collisions(int cellId) const
{
  return _Collisions[cellId];
}

// -----------------------------------------------------------------------------
inline int irtkSurfaceCollisions::NumberOfIntersections(int cellId) const
{
  return static_cast<int>(Intersections(cellId).size());
}

// -----------------------------------------------------------------------------
inline int irtkSurfaceCollisions::NumberOfCollisions(int cellId) const
{
  return static_cast<int>(Collisions(cellId).size());
}


#endif
