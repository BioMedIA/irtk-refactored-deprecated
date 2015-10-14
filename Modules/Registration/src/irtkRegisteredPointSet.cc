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

#include <irtkRegisteredPointSet.h>

#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkIdTypeArray.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyDataNormals.h>

#include <irtkLinearInterpolateImageFunction.hxx>
#include <irtkPolyDataUtils.h>

using namespace irtk::polydata;


// =============================================================================
// Copy points to irtkPointSet structure
// =============================================================================

namespace irtkRegisteredPointSetUtils {

/// Type of interpolator used to interpolate dense displacement field
typedef irtkGenericLinearInterpolateImageFunction<irtkGenericImage<double> > DisplacementInterpolator;

// -----------------------------------------------------------------------------
/// Copy VTK points to IRTK point set
class CopyVtkPointsToIrtkPointSet
{
  vtkPoints    *_Input;
  irtkPointSet *_Output;

public:

  CopyVtkPointsToIrtkPointSet(vtkPoints *in, irtkPointSet &out)
  :
    _Input(in), _Output(&out)
  {}

  void operator ()(const blocked_range<int> &re) const
  {
    double p[3];
    for (int i = re.begin(); i != re.end(); ++i) {
      _Input ->GetPoint(i, p);
      _Output->SetPoint(i, p);
    }
  }
};

// -----------------------------------------------------------------------------
/// Update points of extracted point set surface
class UpdateSurfacePoints
{
  vtkPoints      *_Points;
  vtkPolyData    *_Surface;
  vtkIdTypeArray *_PtIds;

public:

  UpdateSurfacePoints(vtkPoints *points, vtkPolyData *surface, vtkIdTypeArray *ptIds)
  :
    _Points(points), _Surface(surface), _PtIds(ptIds)
  {}

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    double p[3];
    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId) {
      _Points->GetPoint(_PtIds->GetComponent(ptId, 0), p);
      _Surface->GetPoints()->SetPoint(ptId, p);
    }
  }
};

// -----------------------------------------------------------------------------
/// Update point data of extracted point set surface
class UpdateSurfacePointData
{
  vtkPointData *_PointData;
  vtkPolyData  *_Surface;
  vtkDataArray *_PtIds;
  int           _MaxNumberOfComponents;

public:

  UpdateSurfacePointData(vtkPointData *pd, vtkPolyData *surface, vtkIdTypeArray *ptIds)
  :
    _PointData(pd), _Surface(surface), _PtIds(ptIds)
  {
    irtkAssert(_PointData->GetNumberOfArrays() == _Surface->GetPointData()->GetNumberOfArrays(),
               "surface has expected number of point data arrays");
    _MaxNumberOfComponents = 0;
    for (vtkIdType i = 0; i < _PointData->GetNumberOfArrays(); ++i) {
      int n = _PointData->GetArray(i)->GetNumberOfComponents();
      irtkAssert(_Surface->GetPointData()->GetArray(i)->GetNumberOfComponents() == n,
                 "surface arrays have identical number of components");
      if (n > _MaxNumberOfComponents) _MaxNumberOfComponents = n;
    }
  }

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    vtkIdType origPtId;
    double *tuple = new double[_MaxNumberOfComponents];
    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId) {
      origPtId = static_cast<vtkIdType>(_PtIds->GetComponent(ptId, 0));
      for (vtkIdType i = 0; i < _PointData->GetNumberOfArrays(); ++i) {
        _PointData->GetArray(i)->GetTuple(origPtId, tuple);
        _Surface->GetPointData()->GetArray(i)->SetTuple(ptId, tuple);
      }
    }
    delete[] tuple;
  }
};

// -----------------------------------------------------------------------------
/// Transform each point by applying the transformation evaluated at this point
struct TransformPoint
{
  vtkPoints                *_InputPoints;
  vtkPoints                *_OutputPoints;
  const irtkTransformation *_Transformation;
  double                    _t, _t0;

  void operator ()(const blocked_range<vtkIdType> &ids) const
  {
    double p[3];
    for (vtkIdType id = ids.begin(); id != ids.end(); ++id) {
      _InputPoints->GetPoint(id, p);
      _Transformation->Transform(p[0], p[1], p[2], _t, _t0);
      _OutputPoints->SetPoint(id, p);
    }
  }
};

// -----------------------------------------------------------------------------
/// Transform each point by applying the displacement interpolated at this point
struct TransformPointUsingInterpolatedDisplacements
{
  vtkPoints                      *_InputPoints;
  vtkPoints                      *_OutputPoints;
  const DisplacementInterpolator *_Displacement;

  void operator ()(const blocked_range<vtkIdType> &ids) const
  {
    double p[3], d[3];
    for (vtkIdType id = ids.begin(); id != ids.end(); ++id) {
      _InputPoints->GetPoint(id, p);
      memcpy(d, p, 3 * sizeof(double));
      _Displacement->WorldToImage(d[0], d[1], d[2]);
      _Displacement->Evaluate(d, d[0], d[1], d[2]);
      p[0] += d[0], p[1] += d[1], p[2] += d[2];
      _OutputPoints->SetPoint(id, p);
    }
  }
};

// -----------------------------------------------------------------------------
/// Rescale points
struct RescalePoints
{
  vtkPoints *_Points;
  double     _Slope;
  double     _Intercept;

  RescalePoints(vtkPoints *points, double m, double t)
  :
    _Points(points), _Slope(m), _Intercept(t)
  {}

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    double p[3];
    for (vtkIdType i = re.begin(); i != re.end(); ++i) {
      _Points->GetPoint(i, p);
      p[0] = p[0] * _Slope + _Intercept;
      p[1] = p[1] * _Slope + _Intercept;
      p[2] = p[2] * _Slope + _Intercept;
      _Points->SetPoint(i, p);
    }
  }
};

// -----------------------------------------------------------------------------
/// Rescale tuples of given data array
struct RescaleData
{
  vtkDataArray *_Data;
  double        _Slope;
  double        _Intercept;

  RescaleData(vtkDataArray *data, double m, double t)
  :
    _Data(data), _Slope(m), _Intercept(t)
  {}

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    double *v = new double[_Data->GetNumberOfComponents()];
    for (vtkIdType i = re.begin(); i != re.end(); ++i) {
      _Data->GetTuple(i, v);
      for (int c = 0; c < _Data->GetNumberOfComponents(); ++c) {
        v[c] = v[c] * _Slope + _Intercept;
      }
      _Data->SetTuple(i, v);
    }
    delete[] v;
  }
};


} // namespace irtkRegisteredPointSetUtils
using namespace irtkRegisteredPointSetUtils;

// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
irtkRegisteredPointSet::irtkRegisteredPointSet(vtkPointSet              *data,
                                               const irtkTransformation *t)
:
  _InputPointSet       (data),
  _InputSurfacePoints  (NULL),
  _InputTime           (1.0),
  _Time                (.0),
  _Transformation      (t),
  _CopyAll             (true),
  _SelfUpdate          (true),
  _UpdateSurfaceNormals(false),
  _ExternalDisplacement(NULL),
  _Displacement        (NULL)
{
}

// -----------------------------------------------------------------------------
void irtkRegisteredPointSet::Copy(const irtkRegisteredPointSet &other)
{
  if (_InputSurfacePoints != &_InputPoints) Delete(_InputSurfacePoints);
  Delete(_Displacement);

  _InputPointSet        = other._InputPointSet;
  _InputSurface         = other._InputSurface;
  _InputPoints          = other._InputPoints;
  _InputTime            = other._InputTime;
  _Time                 = other._Time;
  _Transformation       = other._Transformation;
  _CopyAll              = other._CopyAll;
  _EdgeTable            = other._EdgeTable;
  _SurfaceEdgeTable     = other._SurfaceEdgeTable;
  _SelfUpdate           = other._SelfUpdate;
  _UpdateSurfaceNormals = other._UpdateSurfaceNormals;
  _Domain               = other._Domain;
  _ExternalDisplacement = other._ExternalDisplacement;

  if (other._InputSurfacePoints) {
    if (other._InputSurfacePoints == &other._InputPoints) {
      _InputSurfacePoints = &_InputPoints;
    } else {
      _InputSurfacePoints = new irtkPointSet(*other._InputSurfacePoints);
    }
  }
  if (other._OutputPointSet) {
    _OutputPointSet = other._OutputPointSet->NewInstance();
    _OutputPointSet->DeepCopy(other._OutputPointSet);
  } else {
    _OutputPointSet = NULL;
  }
  if (other._OutputSurface) {
    _OutputSurface = other._OutputSurface->NewInstance();
    _OutputSurface->DeepCopy(other._OutputSurface);
  } else {
    _OutputSurface = NULL;
  }
}

// -----------------------------------------------------------------------------
irtkRegisteredPointSet::irtkRegisteredPointSet(const irtkRegisteredPointSet &other)
:
  irtkObject(other),
  _InputSurfacePoints(NULL),
  _Displacement(NULL)
{
  Copy(other);
}

// -----------------------------------------------------------------------------
irtkRegisteredPointSet &irtkRegisteredPointSet::operator =(const irtkRegisteredPointSet &other)
{
  if (this != &other) {
    irtkObject::operator =(other);
    Copy(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
irtkRegisteredPointSet::~irtkRegisteredPointSet()
{
  if (_InputSurfacePoints != &_InputPoints) Delete(_InputSurfacePoints);
  Delete(_Displacement);
}

// -----------------------------------------------------------------------------
void irtkRegisteredPointSet::Initialize(bool deep_copy_points, bool init_edge_tables)
{
  // Check input
  if (!_InputPointSet) {
    cerr << "irtkRegisteredPointSet::Initialize: Missing input point set!" << endl;
    exit(1);
  }

  // Extract input point set surface
  _InputSurface = vtkPolyData::SafeDownCast(_InputPointSet);
  if (!_InputSurface) _InputSurface = DataSetSurface(_InputPointSet, true, true);
  _InputSurface->BuildLinks();

  // Copy input points to IRTK structure needed by the energy terms to compute
  // the gradient w.r.t. the transformation parameters
  // irtkTransformation::ParametricGradient and thus used by energy terms
  this->GetInputPoints(_InputPoints);
  if (_InputSurface == _InputPointSet) {
    if (_InputSurfacePoints != &_InputPoints) {
      delete _InputSurfacePoints;
      _InputSurfacePoints = &_InputPoints;
    }
  } else {
    if (!_InputSurfacePoints || _InputSurfacePoints == &_InputPoints) {
      _InputSurfacePoints = new irtkPointSet;
    }
    this->GetInputSurfacePoints(*_InputSurfacePoints);
  }

  // Make shallow copy of input point set and its surface
  //
  // Note that the vtkPoints of the point set are replaced by a new instance
  // upon the first Update of a transformed point set. If no transformation is
  // set or if the point set is never updated, no copy of input points is made.
  // This is overridden by the deep_copy_points argument.
  _OutputPointSet = _InputPointSet->NewInstance();
  _OutputPointSet->ShallowCopy(_InputPointSet);
  if (deep_copy_points) {
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->DeepCopy(_InputPointSet->GetPoints());
    _OutputPointSet->SetPoints(points);
  }

  _OutputSurface = vtkPolyData::SafeDownCast(_OutputPointSet);
  if (!_OutputSurface) {
    _OutputSurface = _InputSurface->NewInstance();
    _OutputSurface->ShallowCopy(_InputSurface);
    // ShallowCopy also copies the previously built cells and links
    if (deep_copy_points) {
      vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
      points->DeepCopy(_InputSurface->GetPoints());
      _OutputSurface->SetPoints(points);
    }
  }

  // Reset point and cell data
  if (!_CopyAll || !_PointDataToCopy.empty()) {
    _OutputPointSet->GetPointData()->Initialize();
    _OutputSurface ->GetPointData()->Initialize();
  }
  if (!_CopyAll) {
    _OutputPointSet->GetCellData()->Initialize();
    _OutputSurface ->GetCellData()->Initialize();
  }

  // Mark input surface normals as invalid
  _UpdateSurfaceNormals = true;

  // Initialize edge tables
  if (init_edge_tables) {
    _EdgeTable.Initialize(_InputPointSet);
    if (_InputPointSet != _InputSurface) {
      _SurfaceEdgeTable.Initialize(_InputSurface);
    }
  } else {
    _EdgeTable.Clear();
    _SurfaceEdgeTable.Clear();
  }
}

// -----------------------------------------------------------------------------
void irtkRegisteredPointSet::PointsChanged()
{
  _UpdateSurfaceNormals = true;
}

// =============================================================================
// Copy points to irtkPointSet structure
// =============================================================================

// -----------------------------------------------------------------------------
vtkIdTypeArray *irtkRegisteredPointSet::OriginalSurfacePointIds() const
{
  return vtkIdTypeArray::SafeDownCast(_InputSurface->GetPointData()->GetArray("vtkOriginalPointIds"));
}

// -----------------------------------------------------------------------------
vtkIdTypeArray *irtkRegisteredPointSet::OriginalSurfaceCellIds() const
{
  return vtkIdTypeArray::SafeDownCast(_InputSurface->GetCellData()->GetArray("vtkOriginalCellIds"));
}

// -----------------------------------------------------------------------------
void irtkRegisteredPointSet::GetInputPoints(irtkPointSet &pset) const
{
  const int n = this->NumberOfPoints();
  pset.Reserve(n);
  pset.Resize(n);
  pset.ShrinkToFit();
  CopyVtkPointsToIrtkPointSet copy(_InputPointSet->GetPoints(), pset);
  parallel_for(blocked_range<int>(0, pset.Size()), copy);
}

// -----------------------------------------------------------------------------
void irtkRegisteredPointSet::GetInputSurfacePoints(irtkPointSet &pset) const
{
  const int n = this->NumberOfSurfacePoints();
  pset.Reserve(n);
  pset.Resize (n);
  pset.ShrinkToFit();
  CopyVtkPointsToIrtkPointSet copy(_InputSurface->GetPoints(), pset);
  parallel_for(blocked_range<int>(0, pset.Size()), copy);
}

// -----------------------------------------------------------------------------
void irtkRegisteredPointSet::GetPoints(irtkPointSet &pset) const
{
  const int n = this->NumberOfPoints();
  pset.Reserve(n);
  pset.Resize(n);
  pset.ShrinkToFit();
  CopyVtkPointsToIrtkPointSet copy(_OutputPointSet->GetPoints(), pset);
  parallel_for(blocked_range<int>(0, pset.Size()), copy);
}

// -----------------------------------------------------------------------------
void irtkRegisteredPointSet::GetSurfacePoints(irtkPointSet &pset) const
{
  const int n = this->NumberOfSurfacePoints();
  pset.Reserve(n);
  pset.Resize(n);
  pset.ShrinkToFit();
  CopyVtkPointsToIrtkPointSet copy(_OutputSurface->GetPoints(), pset);
  parallel_for(blocked_range<int>(0, pset.Size()), copy);
}

// -----------------------------------------------------------------------------
vtkDataArray *irtkRegisteredPointSet::SurfaceNormals() const
{
  if (_OutputSurface->GetPointData()->GetNormals() == NULL || _UpdateSurfaceNormals) {
    vtkSmartPointer<vtkPolyDataNormals> filter = vtkSmartPointer<vtkPolyDataNormals>::New();
    SetVTKInput(filter, _OutputSurface);
    filter->SplittingOff();
    filter->ComputePointNormalsOn();
    filter->ComputeCellNormalsOff();
    filter->ConsistencyOn();
    filter->AutoOrientNormalsOff();
    filter->FlipNormalsOff();
    filter->Update();
    vtkDataArray *normals = filter->GetOutput()->GetPointData()->GetNormals();
    _OutputSurface->GetPointData()->SetNormals(normals);
  }
  return _OutputSurface->GetPointData()->GetNormals();
}

// -----------------------------------------------------------------------------
const irtkRegisteredPointSet::EdgeTable *irtkRegisteredPointSet::Edges() const
{
  if (_EdgeTable.Rows() == 0) {
    _EdgeTable.Initialize(_InputPointSet);
  }
  return &_EdgeTable;
}

// -----------------------------------------------------------------------------
const irtkRegisteredPointSet::EdgeTable *irtkRegisteredPointSet::SurfaceEdges() const
{
  if (_InputSurface == _InputPointSet) return Edges();
  if (_SurfaceEdgeTable.Rows() == 0) {
    _SurfaceEdgeTable.Initialize(_InputSurface);
  }
  return &_SurfaceEdgeTable;
}

// =============================================================================
// Update
// =============================================================================

// -----------------------------------------------------------------------------
void irtkRegisteredPointSet::Update(bool force)
{
  // Do nothing if self-update is disabled (i.e., external process is
  // responsible for update of registered data)
  if (!force && !_SelfUpdate) return;
  IRTK_START_TIMING();

  // ---------------------------------------------------------------------------
  // Transform points
  if (_Transformation) {
    // Update points of output point set
    vtkPoints *inputPoints  = _InputPointSet ->GetPoints();
    vtkPoints *outputPoints = _OutputPointSet->GetPoints();
    const vtkIdType npoints = inputPoints->GetNumberOfPoints();
    if (outputPoints == inputPoints) {
      vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
      _OutputPointSet->SetPoints(points);
      outputPoints = points;
    }
    outputPoints->SetNumberOfPoints(npoints);
    if (_ExternalDisplacement) {
      IRTK_START_TIMING();
      DisplacementInterpolator disp;
      disp.Input(_ExternalDisplacement);
      disp.Initialize();
      TransformPointUsingInterpolatedDisplacements transform;
      transform._InputPoints  = inputPoints;
      transform._OutputPoints = outputPoints;
      transform._Displacement = &disp;
      parallel_for(blocked_range<vtkIdType>(0, npoints), transform);
      IRTK_DEBUG_TIMING(7, "transforming points");
    } else if (_Transformation->RequiresCachingOfDisplacements() && _Domain) {
      IRTK_START_TIMING();
      if (!_Displacement) _Displacement = new irtkGenericImage<double>();
      _Displacement->Initialize(_Domain, 3);
      _Transformation->Displacement(*_Displacement, _Time, _InputTime);
      DisplacementInterpolator disp;
      disp.Input(_Displacement);
      disp.Initialize();
      IRTK_DEBUG_TIMING(7, "caching displacements");
      IRTK_RESET_TIMING();
      TransformPointUsingInterpolatedDisplacements transform;
      transform._InputPoints  = inputPoints;
      transform._OutputPoints = outputPoints;
      transform._Displacement = &disp;
      parallel_for(blocked_range<vtkIdType>(0, npoints), transform);
      IRTK_DEBUG_TIMING(7, "transforming points");
    } else {
      IRTK_START_TIMING();
      TransformPoint transform;
      transform._InputPoints    = inputPoints;
      transform._OutputPoints   = outputPoints;
      transform._Transformation = _Transformation;
      transform._t              = _Time;
      transform._t0             = _InputTime;
      parallel_for(blocked_range<vtkIdType>(0, npoints), transform);
      IRTK_DEBUG_TIMING(7, "transforming points");
    }
    // Update points of output point set surface
    if (_OutputSurface->GetPoints() != outputPoints) {
      if (_OutputSurface->GetPoints() == _InputSurface->GetPoints()) {
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        _OutputSurface->SetPoints(points);
      }
      const vtkIdType nspoints = _InputSurface->GetNumberOfPoints();
      _OutputSurface->GetPoints()->SetNumberOfPoints(nspoints);
      UpdateSurfacePoints update(outputPoints, _OutputSurface, OriginalSurfacePointIds());
      parallel_for(blocked_range<vtkIdType>(0, nspoints), update);
    }
    // Invalidate cached attributes which are recomputed on demand only
    PointsChanged();
  }

  // ---------------------------------------------------------------------------
  // Copy and optionally normalize point data
  if (!_PointDataToCopy.empty()) {
    IRTK_START_TIMING();
    // Update point data of output point set
    vtkPointData * const inputPD  = _InputPointSet ->GetPointData();
    vtkPointData * const outputPD = _OutputPointSet->GetPointData();
    vtkDataArray *input_array;
    vtkSmartPointer<vtkDataArray> output_array;
    vector<ScalingFunction>::iterator it;
    for (it = _PointDataToCopy.begin(); it != _PointDataToCopy.end(); ++it) {
      if (it->_Slope == .0) continue;
      if (it->_InputIndex < 0 || it->_InputIndex >= inputPD->GetNumberOfArrays()) {
        cerr << "irtkRegisteredPointSet::Update: Invalid point data index: " << it->_InputIndex << endl;
        exit(1);
      }
      int attr = inputPD->IsArrayAnAttribute(it->_InputIndex);
      input_array = inputPD->GetArray(it->_InputIndex);
      if (0 <= it->_OutputIndex && it->_OutputIndex < outputPD->GetNumberOfArrays()) {
        output_array = outputPD->GetArray(it->_OutputIndex);
      } else {
        if (it->_Slope == 1.0 && it->_Intercept == .0) {
          output_array = input_array->NewInstance();
        } else {
          output_array = vtkSmartPointer<vtkFloatArray>::New();
        }
        output_array->SetName(input_array->GetName());
        it->_OutputIndex = -1;
      }
      output_array->DeepCopy(input_array);
      if (attr < 0) {
        parallel_for(blocked_range<vtkIdType>(0, outputPD->GetNumberOfTuples()),
                     RescaleData(output_array, it->_Slope, it->_Intercept));
      }
      if (it->_OutputIndex == -1) {
        it->_OutputIndex = outputPD->AddArray(output_array);
        if (attr >= 0) outputPD->SetActiveAttribute(it->_OutputIndex, attr);
      }
    }
    // Update point data of output point set surface
    if (_OutputSurface->GetPointData() != outputPD) {
      UpdateSurfacePointData update(outputPD, _OutputSurface, OriginalSurfacePointIds());
      parallel_for(blocked_range<vtkIdType>(0, _OutputSurface->GetNumberOfPoints()), update);
    }
    IRTK_DEBUG_TIMING(7, "copying and rescaling point data");
  }

  IRTK_DEBUG_TIMING(2, "update of" << (_Transformation ? " moving " : " fixed ")
                                   << (IsSurface() ? "surface" : "point set"));
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
const char *irtkRegisteredPointSet::DefaultExtension() const
{
  return irtk::polydata::DefaultExtension(_InputPointSet);
}

// -----------------------------------------------------------------------------
void irtkRegisteredPointSet::Write(const char       *fname,
                                   vtkAbstractArray *pointdata,
                                   vtkAbstractArray *celldata) const
{
  Write(fname, &pointdata, 1, &celldata, 1);
}

// -----------------------------------------------------------------------------
void irtkRegisteredPointSet::Write(const char        *fname,
                                   vtkAbstractArray **pointdata, int npointdata,
                                   vtkAbstractArray **celldata,  int ncelldata) const
{
  // Add point data
  int *pointdataidx = NULL;
  if (npointdata > 0) {
    pointdataidx = new int[npointdata];
    for (int i = 0; i < npointdata; ++i) {
      pointdataidx[i] = _OutputPointSet->GetPointData()->AddArray(pointdata[i]);
    }
  }
  // Add cell data
  int *celldataidx = NULL;
  if (ncelldata > 0) {
    celldataidx = new int[ncelldata];
    for (int i = 0; i < ncelldata; ++i) {
      celldataidx[i] = _OutputPointSet->GetCellData()->AddArray(celldata[i]);
    }
  }
  // Write data set with additional point and cell data
  WritePointSet(fname, _OutputPointSet);
  // Remove point data
  for (int i = 0; i < npointdata; ++i) {
    _OutputPointSet->GetPointData()->RemoveArray(pointdataidx[i]);
  }
  delete[] pointdataidx;
  // Remove cell data
  for (int i = 0; i < ncelldata; ++i) {
    _OutputPointSet->GetCellData()->RemoveArray(celldataidx[i]);
  }
  delete[] celldataidx;
}
