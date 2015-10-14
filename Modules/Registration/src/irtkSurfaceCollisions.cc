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

#include <irtkSurfaceCollisions.h>


#include <irtkTriangle.h>
#include <irtkPolyDataUtils.h>
#include <irtkPointLocator.h>

#include <vtkMath.h>
#include <vtkPlane.h>
#include <vtkTriangle.h>

#include <vtkPolyData.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkUnsignedCharArray.h>

#include <vtkAbstractPointLocator.h>
#include <vtkOctreePointLocator.h>
#include <vtkIntersectionPolyDataFilter.h>

using namespace irtk::polydata;


// =============================================================================
// Auxiliary functors
// =============================================================================

namespace irtkSurfaceCollisionsUtils {


// -----------------------------------------------------------------------------
/// Compute center and radius of each bounding spheres
class ComputeBoundingSpheres
{
  vtkPolyData  *_Surface;
  vtkDataArray *_Center;
  vtkDataArray *_Radius;

  ComputeBoundingSpheres(vtkPolyData *surface, vtkDataArray *center, vtkDataArray *radius)
  :
    _Surface(surface), _Center(center), _Radius(radius)
  {}

public:

  /// Compute bounding spheres of specified range of cells
  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    vtkIdType npts, *pts;
    double a[3], b[3], c[3], origin[3], radius;

    for (vtkIdType cellId = re.begin(); cellId != re.end(); ++cellId) {
      // Get triangle vertices
      _Surface->GetCellPoints(cellId, npts, pts);
      irtkAssert(npts == 3, "surface is triangular mesh");

      // Get triangle vertex positions
      _Surface->GetPoint(pts[0], a);
      _Surface->GetPoint(pts[1], b);
      _Surface->GetPoint(pts[2], c);

      // Get center of bounding sphere
      vtkTriangle::TriangleCenter(a, b, c, origin);
      _Center->SetTuple(cellId, origin);

      // Compute radius of bounding sphere
      radius = sqrt(max(max(vtkMath::Distance2BetweenPoints(a, origin),
                            vtkMath::Distance2BetweenPoints(b, origin)),
                            vtkMath::Distance2BetweenPoints(c, origin)));
      _Radius->SetTuple1(cellId, radius);
    }
  }

  /// Compute bounding spheres of surface faces
  static void Run(vtkPolyData *surface, vtkDataArray *center, vtkDataArray *radius)
  {
    center->SetNumberOfComponents(3);
    center->SetNumberOfTuples(surface->GetNumberOfCells());
    radius->SetNumberOfComponents(1);
    radius->SetNumberOfTuples(surface->GetNumberOfCells());
    if (surface->GetNumberOfCells() == 0) return;
    ComputeBoundingSpheres eval(surface, center, radius);
    parallel_for(blocked_range<vtkIdType>(0, surface->GetNumberOfCells()), eval);
  }
};

// -----------------------------------------------------------------------------
/// Check for collisions such as self-intersections and faces too close
class FindCollisions
{
  typedef irtkSurfaceCollisions::CollisionType    CollisionType;
  typedef irtkSurfaceCollisions::IntersectionInfo IntersectionInfo;
  typedef irtkSurfaceCollisions::CollisionInfo    CollisionInfo;

  irtkSurfaceCollisions          *_Filter;
  vtkAbstractPointLocator        *_Locator;
  vector<set<IntersectionInfo> > *_Intersections;
  vector<set<CollisionInfo> >    *_Collisions;
  double                          _MaxRadius;
  double                          _MinDistance;
  double                          _MinAngleCos;

  /// Check whether point c is on the "left" of the line defined by a and b
  inline int IsLeft(const double a[2], const double b[2], const double c[2]) const
  {
     return sgn((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]));
  }

public:

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    const double TOL = 1e-6; // Tolerance for point coordinate equality check

    vtkPolyData  *surface   = _Filter->Output();
    vtkDataArray *center    = _Filter->GetCenterArray();
    vtkDataArray *radius    = _Filter->GetRadiusArray();
    vtkDataArray *coll_type = _Filter->GetCollisionTypeArray();

    const double R = _MaxRadius + 1.1 * _Filter->MinDistance();

    double         tri1[3][3], tri2[3][3], tri1_2D[3][2], tri2_2D[3][2];
    double         n1[3], n2[3], p1[3], p2[3], r1, c1[3], v[3], search_radius, dot;
    int            tri12[3], i1, i2, shared_vertex1, shared_vertex2, coplanar, s1, s2;
    vtkIdType      npts, *pts1, *pts2, *cells, cellId1, cellId2;
    unsigned short ncells;
    CollisionInfo  collision;
    CollisionType  type;

    vtkSmartPointer<vtkIdList> ptIds   = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();

    for (cellId1 = re.begin(); cellId1 != re.end(); ++cellId1) {

      surface->GetCellPoints(cellId1, npts, pts1);
      irtkAssert(npts == 3, "surface is triangular mesh");

      // No collision if none detected
      set<IntersectionInfo> &intersections = (*_Intersections)[cellId1];
      set<CollisionInfo>    &collisions    = (*_Collisions   )[cellId1];
      type = CollisionType::NoCollision;

      // Get vertices and normal of this triangle
      surface->GetPoint(pts1[0], tri1[0]);
      surface->GetPoint(pts1[1], tri1[1]);
      surface->GetPoint(pts1[2], tri1[2]);
      vtkTriangle::ComputeNormal(tri1[0], tri1[1], tri1[2], n1);

      // Get bounding sphere
      center->GetTuple(cellId1, c1);
      r1 = radius->GetComponent(cellId1, 0);

      // Find other triangles within search radius
      search_radius = min(max(_Filter->MinSearchRadius(), r1 + R), _Filter->MaxSearchRadius());
      _Locator->FindPointsWithinRadius(search_radius, c1, ptIds);
      cellIds->Reset();
      for (vtkIdType i = 0; i < ptIds->GetNumberOfIds(); ++i) {
        surface->GetPointCells(ptIds->GetId(i), ncells, cells);
        for (unsigned short j = 0; j < ncells; ++j) {
          if (cells[j] != cellId1) cellIds->InsertUniqueId(cells[j]);
        }
      }

      // Check for collisions between this triangle and the found nearby triangles
      for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {
        cellId2 = cellIds->GetId(i);

        // Get vertex positions of nearby candidate triangle
        surface->GetCellPoints(cellId2, npts, pts2);
        surface->GetPoint(pts2[0], tri2[0]);
        surface->GetPoint(pts2[1], tri2[1]);
        surface->GetPoint(pts2[2], tri2[2]);

        // Get corresponding indices of shared vertices
        for (i1 = 0; i1 < 3; ++i1) {
          tri12[i1] = -1;
          for (i2 = 0; i2 < 3; ++i2) {
            if (pts1[i1] == pts2[i2]) {
              tri12[i1] = i2;
              break;
            }
          }
        }

        // Determine whether triangles share one vertex or even an edge
        shared_vertex1 = shared_vertex2 = -1;
        for (i1 = 0; i1 < 3; ++i1) {
          if (tri12[i1] != -1) {
            shared_vertex1 = i1;
            for (i2 = i1 + 1; i2 < 3; ++i2) {
              if (tri12[i2] != -1) {
                shared_vertex2 = i2;
                break;
              }
            }
            break;
          }
        }

        // Triangles with a shared edge can only overlap, but not intersect
        if (shared_vertex2 != -1) {
          if (_Filter->AdjacentIntersectionTest()) {
            coplanar = 0;
            for (i2 = 0; i2 < 3; ++i2) {
              if (i2 != tri12[shared_vertex1] && i2 != tri12[shared_vertex2]) {
                coplanar = int(vtkPlane::DistanceToPlane(tri2[i2], tri1[0], n1) < TOL);
                break;
              }
            }
            if (coplanar) {
              // Project triangles into their common plane
              vtkTriangle::ProjectTo2D(tri1   [0], tri1   [1], tri1   [2],
                                       tri1_2D[0], tri1_2D[1], tri1_2D[2]);
              vtkTriangle::ProjectTo2D(tri2   [0], tri2   [1], tri2   [2],
                                       tri2_2D[0], tri2_2D[1], tri2_2D[2]);
              // Check on which side of the shared edge the non-shared vertices are
              s1 = 0;
              for (i1 = 0; i1 < 3; ++i1) {
                if (i1 != shared_vertex1 && i1 != shared_vertex2) {
                  s1 = IsLeft(tri1_2D[shared_vertex1], tri1_2D[shared_vertex2], tri1_2D[i1]);
                  break;
                }
              }
              // Note: re-uses i2 index found outside of this if-block
              s2 = IsLeft(tri1_2D[shared_vertex1], tri1_2D[shared_vertex2], tri2_2D[i2]);
              // If they are on the same side of the edge, the triangles must overlap
              if (s1 == s2) {
                intersections.insert(IntersectionInfo(cellId2, true));
                if (_Filter->IsCollision(type)) {
                  type = CollisionType::Ambiguous;
                } else if (type == CollisionType::SelfIntersection) {
                  type = CollisionType::Intersection;
                } else {
                  type = CollisionType::AdjacentIntersection;
                }
              }
            }
          }
        }
        // If triangles share a single vertex, use different triangle/triangle
        // intersection check which gives us the points forming the line/point
        // of intersection, but does not handle the coplanar case itself
        else if (shared_vertex1 != -1) {
          if (_Filter->AdjacentIntersectionTest()) {
            if (vtkIntersectionPolyDataFilter::TriangleTriangleIntersection(
                  tri1[0], tri1[1], tri1[2],
                  tri2[0], tri2[1], tri2[2], coplanar, p1, p2)) {
              // Ignore valid intersection of single shared vertex
              if (!fequal(p1[0], p2[0], TOL) ||
                  !fequal(p1[1], p2[1], TOL) ||
                  !fequal(p1[2], p2[2], TOL)) {
                intersections.insert(IntersectionInfo(cellId2, true));
                if (_Filter->IsCollision(type)) {
                  type = CollisionType::Ambiguous;
                } else if (type == CollisionType::SelfIntersection) {
                  type = CollisionType::Intersection;
                } else {
                  type = CollisionType::AdjacentIntersection;
                }
              }
            } else if (coplanar) {
              // TODO: Either one triangle fully contained within the other
              //       or one edge of the first triangle intersects an edge of
              //       the second triangle.
            }
          }
        }
        // In case of non-adjacent triangles, use fast intersection test which
        // also checks for overlap of coplanar triangles
        else if (_Filter->NonAdjacentIntersectionTest() &&
                 irtkTriangle::TriangleTriangleIntersection(tri1[0], tri1[1], tri1[2],
                                                            tri2[0], tri2[1], tri2[2])) {
          intersections.insert(IntersectionInfo(cellId2, false));
          if (_Filter->IsCollision(type)) {
            type = CollisionType::Ambiguous;
          } else if (type == CollisionType::AdjacentIntersection) {
            type = CollisionType::Intersection;
          } else {
            type = CollisionType::SelfIntersection;
          }
        }
        // If self-intersection check of non-adjacent triangles disabled or negative,
        // check for near miss collision if minimum distance set
        else if (_MinDistance > .0) {
          vtkTriangle::ComputeNormal(tri2[0], tri2[1], tri2[2], n2);
          collision._Distance = irtkTriangle::DistanceBetweenTriangles(
                                    tri1[0], tri1[1], tri1[2], n1,
                                    tri2[0], tri2[1], tri2[2], n2,
                                    collision._Point1, collision._Point2);
          if (collision._Distance < _MinDistance) {
            vtkMath::Subtract(collision._Point2, c1, v);
            vtkMath::Normalize(v);
            dot = vtkMath::Dot(v, n1);
            if (fabs(dot) >= _MinAngleCos) {
              if (dot < .0) {
                if (_Filter->BackfaceCollisionTest()) {
                  collision._Type = CollisionType::BackfaceCollision;
                } else {
                  collision._Type = CollisionType::NoCollision;
                }
              } else {
                if (_Filter->FrontfaceCollisionTest()) {
                  collision._Type = CollisionType::FrontfaceCollision;
                } else {
                  collision._Type = CollisionType::NoCollision;
                }
              }
              if (collision._Type != CollisionType::NoCollision) {
                collision._CellId = cellId2;
                collisions.insert(collision);
                if (_Filter->IsIntersection(type)) {
                  type = CollisionType::Ambiguous;
                } else if (type != CollisionType::NoCollision) {
                  if (type != collision._Type) type = CollisionType::Collision;
                } else {
                  type = collision._Type;
                }
              }
            }
          }
        }
      }

      // Set collision type of this cell
      coll_type->SetComponent(cellId1, 0, type);
    }
  }

  /// Find collision and self-intersections
  static void Run(irtkSurfaceCollisions          *filter,
                  vector<set<IntersectionInfo> > &intersections,
                  vector<set<CollisionInfo   > > &collisions)
  {
    vtkDataArray *radius = filter->GetRadiusArray();
    vtkSmartPointer<vtkAbstractPointLocator> locator;
    locator = vtkSmartPointer<vtkOctreePointLocator>::New();
    locator->SetDataSet(filter->Output());
    locator->BuildLocator();
    FindCollisions body;
    body._Filter        = filter;
    body._Locator       = locator;
    body._Intersections = &intersections;
    body._Collisions    = &collisions;
    body._MaxRadius     = radius->GetRange(0)[1];
    body._MinDistance   = filter->MinDistance();
    body._MinAngleCos   = cos(filter->MaxAngle() * M_PI / 180.0);
    if (!filter->BackfaceCollisionTest() && !filter->FrontfaceCollisionTest()) {
      body._MinDistance = .0;
    }
    blocked_range<vtkIdType> cellIds(0, filter->Output()->GetNumberOfCells());
    parallel_for(cellIds, body);
  }
};


} // namespace irtkSurfaceCollisionsUtils
using namespace irtkSurfaceCollisionsUtils;

// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkSurfaceCollisions::irtkSurfaceCollisions()
:
  _MinSearchRadius(.0),
  _MaxSearchRadius(numeric_limits<double>::infinity()),
  _MinDistance(.0),
  _MaxAngle(90.0),
  _AdjacentIntersectionTest(false),
  _NonAdjacentIntersectionTest(false),
  _FrontfaceCollisionTest(false),
  _BackfaceCollisionTest(false)
{
}

// -----------------------------------------------------------------------------
irtkSurfaceCollisions::~irtkSurfaceCollisions()
{
}

// =============================================================================
// Output attributes
// =============================================================================

// -----------------------------------------------------------------------------
vtkDataArray *irtkSurfaceCollisions::GetCenterArray() const
{
  return _Output->GetCellData()->GetArray("BoundingSphereCenter");
}

// -----------------------------------------------------------------------------
vtkDataArray *irtkSurfaceCollisions::GetRadiusArray() const
{
  return _Output->GetCellData()->GetArray("BoundingSphereRadius");
}

// -----------------------------------------------------------------------------
vtkDataArray *irtkSurfaceCollisions::GetCollisionTypeArray() const
{
  return _Output->GetCellData()->GetArray("CollisionType");
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void irtkSurfaceCollisions::Initialize()
{
  _NumberOfIntersections = 0;
  _Intersections.clear();
  _IntersectionCells.clear();

  _NumberOfCollisions = 0;
  _Collisions.clear();
  _CollisionCells.clear();

  if (!_Input) {
    cerr << "irtkSurfaceCollisions::Initialize: Missing input surface" << endl;
    exit(1);
  }

  if (_MinDistance     < .0) _MinDistance     = .0;
  if (_MinSearchRadius < .0) _MinSearchRadius = .0;
  if (_MinDistance <= .0) {
    if (!_AdjacentIntersectionTest && !_NonAdjacentIntersectionTest) {
      _AdjacentIntersectionTest = _NonAdjacentIntersectionTest = true;
    }
  } else {
    if (!_BackfaceCollisionTest && !_FrontfaceCollisionTest) {
      _FrontfaceCollisionTest = _BackfaceCollisionTest = true;
    }
  }

  if (IsTriangularMesh(_Input)) {
    _Output = _Input->NewInstance();
    _Output->ShallowCopy(_Input);
  } else {
    _Output = Triangulate(_Input);
  }

  vtkSmartPointer<vtkDataArray> center = vtkSmartPointer<vtkFloatArray>::New();
  vtkSmartPointer<vtkDataArray> radius = vtkSmartPointer<vtkFloatArray>::New();
  vtkSmartPointer<vtkDataArray> types  = vtkSmartPointer<vtkUnsignedCharArray>::New();

  center->SetName("BoundingSphereCenter");
  radius->SetName("BoundingSphereRadius");
  types ->SetName("CollisionType");

  ComputeBoundingSpheres::Run(_Output, center, radius);
  types->SetNumberOfComponents(1);
  types->SetNumberOfTuples(_Output->GetNumberOfCells());

  _Output->GetCellData()->AddArray(center);
  _Output->GetCellData()->AddArray(radius);
  _Output->GetCellData()->AddArray(types);

  _Intersections.resize(_Output->GetNumberOfCells());
  _Collisions   .resize(_Output->GetNumberOfCells());
}

// -----------------------------------------------------------------------------
void irtkSurfaceCollisions::Finalize()
{
  // Count total number of intersections
  for (size_t i = 0; i < _Intersections.size(); ++i) {
    if (!_Intersections[i].empty()) {
      _NumberOfIntersections += static_cast<int>(_Intersections[i].size());
      _IntersectionCells.insert(i);
    }
  }
  _NumberOfIntersections /= 2;

  // Count total number of collisions
  for (size_t i = 0; i < _Collisions.size(); ++i) {
    if (!_Collisions[i].empty()) {
      _NumberOfCollisions += static_cast<int>(_Collisions[i].size());
      _CollisionCells.insert(i);
    }
  }
  _NumberOfCollisions /= 2;
}

// -----------------------------------------------------------------------------
void irtkSurfaceCollisions::Run()
{
  Initialize();
  if (_Output->GetNumberOfCells() > 0) {
    FindCollisions::Run(this, _Intersections, _Collisions);
  }
  Finalize();
}
