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

#include <irtkVolumeParameterizer.h>

#include <irtkHarmonicVolumeParameterizer.h>
#include <irtkAsConformalAsPossibleVolumeParameterizer.h>

#include <irtkPolyDataUtils.h>

#include <vtkPointData.h>
#include <vtkModifiedBSPTree.h>



namespace irtk { namespace polydata {

// =============================================================================
// Factory method
// =============================================================================

// -----------------------------------------------------------------------------
irtkVolumeParameterizer *irtkVolumeParameterizer::New(MapType map)
{
  switch (map) {
    case DefaultMap:  return new irtkHarmonicVolumeParameterizer;
    case HarmonicMap: return new irtkHarmonicVolumeParameterizer;
    case ACAPMap:     return new irtkAsConformalAsPossibleVolumeParameterizer;
    default:
      cerr << "irtkVolumeParameterizer::New: Unknown map type/method: " << map << endl;
      exit(1);
  }
}

// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void irtkVolumeParameterizer::Copy(const irtkVolumeParameterizer &other)
{
  _Input                  = other._Input;
  _CoordsName             = other._CoordsName;
  _NumberOfPoints         = other._NumberOfPoints;
  _NumberOfBoundaryPoints = other._NumberOfBoundaryPoints;
  _NumberOfInteriorPoints = other._NumberOfInteriorPoints;
  // TODO
}

// -----------------------------------------------------------------------------
irtkVolumeParameterizer::irtkVolumeParameterizer()
:
  _CoordsName("coords"),
  _NumberOfPoints(0),
  _NumberOfBoundaryPoints(0),
  _NumberOfInteriorPoints(0)
{
}

// -----------------------------------------------------------------------------
irtkVolumeParameterizer::irtkVolumeParameterizer(const irtkVolumeParameterizer &other)
{
  Copy(other);
}

// -----------------------------------------------------------------------------
irtkVolumeParameterizer &irtkVolumeParameterizer::operator =(const irtkVolumeParameterizer &other)
{
  irtkObject::operator =(other);
  Copy(other);
  return *this;
}

// -----------------------------------------------------------------------------
irtkVolumeParameterizer::~irtkVolumeParameterizer()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void irtkVolumeParameterizer::Run()
{
  this->Initialize();
  this->Parameterize();
  this->Finalize();
}

// -----------------------------------------------------------------------------
void irtkVolumeParameterizer::Initialize()
{
  // Check input
  if (!_Input) {
    cerr << "irtkVolumeParameterizer::Initialize: Missing input point set" << endl;
    exit(1);
  }
  if (_Input->GetNumberOfCells() == 0) {
    cerr << "irtkVolumeParameterizer::Initialize: Input has no cells" << endl;
    exit(1);
  }

  _Coords = GetArrayByCaseInsensitiveName(_Input->GetPointData(), _CoordsName.c_str());
  if (!_Coords) {
    cerr << "irtkVolumeParameterizer::Initialize: Missing input point coordinates array named: "
         << _CoordsName << " (case insensitive)" << endl;
    exit(1);
  }
  if (_Coords->GetNumberOfTuples() != _Input->GetNumberOfPoints()) {
    cerr << "irtkVolumeParameterizer::Initialize: Input point coordinates array "
         << _CoordsName << " must have one tuple for each input point" << endl;
    exit(1);
  }
  if (_Coords->GetNumberOfComponents() < 3) {
    cerr << "irtkVolumeParameterizer::Initialize: Input point coordinates array "
         << _CoordsName << " must have at least 3 components" << endl;
    exit(1);
  }

  // Tetrahedralize interior of input point set
  _Volume = Tetrahedralize(_Input);

  // Keep only first three output domain coordinates
  int loc;
  _Coords = GetArrayByCaseInsensitiveName(_Volume->GetPointData(), _CoordsName.c_str(), &loc);
  if (_Coords->GetNumberOfComponents() > 3) {
    vtkSmartPointer<vtkDataArray> coords = _Coords->NewInstance();
    coords->SetName(_Coords->GetName());
    coords->SetNumberOfComponents(3);
    coords->SetNumberOfTuples(_Coords->GetNumberOfTuples());
    for (vtkIdType ptId = 0; ptId < _Coords->GetNumberOfTuples(); ++ptId) {
      for (int j = 0; j < 3; ++j) {
        coords->SetComponent(ptId, j, _Coords->GetComponent(ptId, j));
      }
    }
    _Volume->GetPointData()->RemoveArray(loc);
    loc = _Volume->GetPointData()->AddArray(coords);
    _Coords = coords;
  }

  // Extract surface of tetrahedral mesh
  _Surface = DataSetSurface(_Volume, true);

  // Make shallow copy of tetrahedral mesh
  _Output = _Volume->NewInstance();
  _Output->ShallowCopy(_Volume);

  // Set read-only attributes of number of interior/boundary points
  _NumberOfPoints         = static_cast<int>(_Volume ->GetNumberOfPoints());
  _NumberOfBoundaryPoints = static_cast<int>(_Surface->GetNumberOfPoints());
  _NumberOfInteriorPoints = _NumberOfPoints - _NumberOfBoundaryPoints;

  // Initialize boundary point mask
  _IsBoundaryPoint.clear();
  _IsBoundaryPoint.resize(_NumberOfPoints, false);
  vtkDataArray *origPtIds = _Surface->GetPointData()->GetArray("vtkOriginalPointIds");
  for (vtkIdType ptId = 0; ptId < _Surface->GetNumberOfPoints(); ++ptId) {
    _IsBoundaryPoint[static_cast<vtkIdType>(origPtIds->GetComponent(ptId, 0))] = true;
  }
}

// -----------------------------------------------------------------------------
void irtkVolumeParameterizer::Finalize()
{
  // Replace point coordinates of output tetrahedral mesh
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(_Volume->GetNumberOfPoints());
  for (vtkIdType ptId = 0; ptId < _Volume->GetNumberOfPoints(); ++ptId) {
    points->SetPoint(ptId, _Coords->GetTuple(ptId));
  }
  _Output->SetPoints(points);

  // Build cell locator used to map any point in the input domain (cf. MapPoint)
  _Locator = vtkSmartPointer<vtkModifiedBSPTree>::New();
  _Locator->SetDataSet(_Volume);
  _Locator->BuildLocator();
}


} } // namespace irtk::polydata
