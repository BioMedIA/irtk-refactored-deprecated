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

#include <irtkGeometry.h>

#ifdef HAS_VTK

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkPolyDataReader.h>
#include <vtkXMLGenericDataObjectReader.h>
#include <vtkOBJReader.h>
#include <vtkPLYReader.h>
#include <vtkSTLReader.h>
#include <vtkPolyDataWriter.h>

#endif

int intersection(double x1, double y1, double x2, double y2,
                 double x3, double y3, double x4, double y4)
{
  double a = (x4 - x3)*(y1 - y3) - (y4 - y3)*(x1 - x3);
  double b = (x2 - x1)*(y1 - y3) - (y2 - y1)*(x1 - x3);
  if ((y4 - y3)*(x2 - x1) - (x4 - x3)*(y2 - y1) != 0) {
    a /= (y4 - y3)*(x2 - x1) - (x4 - x3)*(y2 - y1);
    b /= (y4 - y3)*(x2 - x1) - (x4 - x3)*(y2 - y1);
    if ((a >= 0) && (a < 1) && (b >= 0) && (b < 1)) {
      return 1;
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}

void irtkPointSet::Clear()
{
  Deallocate(_data);
  _m = 0;
  _n = 0;
}

void irtkPointSet::Add(const irtkPoint &p)
{
  if (_n + 1 <= _m) {
    // There is still enough memory left, so just add the point
    _data[_n] = p;
    _n++;
    return;
  }
  // There is not enough memory left, so allocate new point list and copy
  _m += POINTSET_SIZE;
  irtkPoint *new_data = Allocate<irtkPoint>(_m);
  for (int i = 0; i < _n; i++) {
    new_data[i] = _data[i];
  }
  new_data[_n] = p;
  Deallocate(_data);
  _data = new_data;
  _n++;
}

void irtkPointSet::Del(const irtkPoint &p)
{
  irtkPoint *new_data = Allocate<irtkPoint>(_m);

  int new_n = 0;
  for (int i = 0; i < _n; ++i) {
    if (_data[i] == p) continue;
    new_data[new_n] = _data[i];
    ++new_n;
  }

  Deallocate(_data);
  _data = new_data;
  _n    = new_n;

  if ((new_n % POINTSET_SIZE) == 0) {
    _m = (new_n / POINTSET_SIZE) * POINTSET_SIZE;
    new_data = Allocate<irtkPoint>(_m);
    for (int i = 0; i < _n; ++i) {
      new_data[i] = _data[i];
    }
    Deallocate(_data);
    _data = new_data;
  }
}

void irtkPointSet::Add(const irtkPointSet &pset)
{
  for (int i = 0; i < pset.Size(); i++) this->Add(pset(i));
}

void irtkPointSet::Del(const irtkPointSet &pset)
{
  for (int i = 0; i < pset.Size(); i++) this->Del(pset(i));
}

irtkPointSet& irtkPointSet::operator=(const irtkPointSet &pset)
{
  Reserve(pset._n);
  _n = pset._n;
  for (int i = 0; i < _n; ++i) {
    _data[i] = pset._data[i];
  }
  return *this;
}

irtkPoint irtkPointSet::CenterOfGravity() const
{
  int i;
  irtkPoint p;

  if (this->Size() == 0) {
    cerr << "irtkPointSet::CenterOfGravity(): No points in point set" << endl;
    return p;
  }
  for (i = 0; i < this->Size(); i++) {
    p += (*this)(i);
  }
  return p / (double)this->Size();
}

void irtkPointSet::BoundingBox(irtkPoint &p1, irtkPoint &p2) const
{
  int i;

  if (this->Size() == 0) {
    p1 = irtkPoint();
    p2 = irtkPoint();
    return;
  } else {
    p1 = (*this)(0);
    p2 = (*this)(0);
  }
  for (i = 1; i < this->Size(); i++) {
    p1._x = p1._x < (*this)(i)._x ? p1._x : (*this)(i)._x;
    p1._y = p1._y < (*this)(i)._y ? p1._y : (*this)(i)._y;
    p1._z = p1._z < (*this)(i)._z ? p1._z : (*this)(i)._z;
    p2._x = p2._x > (*this)(i)._x ? p2._x : (*this)(i)._x;
    p2._y = p2._y > (*this)(i)._y ? p2._y : (*this)(i)._y;
    p2._z = p2._z > (*this)(i)._z ? p2._z : (*this)(i)._z;
  }
}

bool irtkPointSet::operator ==(const irtkPointSet &rhs) const
{
  if (_n != rhs._n) return false;
  for (int i = 0; i < _n; ++i) {
    if (_data[i] != rhs._data[i]) return false;
  }
  return true;
}

bool irtkPointSet::operator !=(const irtkPointSet &rhs) const
{
  return !(*this == rhs);
}

ostream& operator<< (ostream& os, const irtkPointSet &pset)
{
  int i;

  os << "irtkPointSet " << pset._n << endl;
  os.setf(ios::right);
  os.setf(ios::fixed);
  os.precision(10);
  for (i = 0; i < pset._n; i++) {
    os << setw(15) << pset._data[i] << endl;
  }
  os.precision(6);
  os.unsetf(ios::right);
  os.unsetf(ios::fixed);
  return os;
}

istream& operator>> (istream& is, irtkPointSet &pset)
{
  int i, size;
  char buffer[255];
  irtkPoint p;

  // Read header
  is >> buffer;
  if (strcmp(buffer, "irtkPointSet") != 0) {
    cerr << "Can't read irtkPointSet file: " << buffer << endl;
    exit(1);
  }

  // Read size
  is >> size;

  // Read irtkPointSet
  for (i = 0; i < size; i++) {
    is >> p;
    pset.Add(p);
  }
  return is;
}

void irtkPointSet::Read(const char *filename)
{
  // Read VTK file if extension matches known format supported by VTK
  const string ext = Extension(filename);
  if (ext == ".vtk" ||
     (ext.length() == 4 && ext.substr(0, 3) == ".vt") ||
      ext == ".obj" || ext == ".ply" || ext == ".stl") {
    ReadVTK(filename);
    return;
  }

  // Clear pointset first
  this->Clear();

  // Open file stream
  ifstream from(filename);

  // Check whether file opened ok
  if (!from) {
    cerr << "irtkPointSet::Read: Can't open file " << filename << endl;
    exit(1);
  }

  // Read irtkPointSet
  from >> *this;
}

void irtkPointSet::Write(const char *filename) const
{
  // Open file stream
  ofstream to(filename);

  // Check whether file opened ok
  if (!to) {
    cerr << "irtkPointSet::Write: Can't open file " << filename << endl;
    exit(1);
  }

  // Write irtkPointSet
  to << *this;
}

void irtkPointSet::AddVTK(const char *fname)
{
#ifdef HAS_VTK

  vtkSmartPointer<vtkPointSet> pointset;
  const string ext = Extension(fname);
  if (ext == ".obj") {
    vtkSmartPointer<vtkOBJReader> reader;
    reader = vtkSmartPointer<vtkOBJReader>::New();
    reader->SetFileName(fname);
    reader->Update();
    pointset = reader->GetOutput();
  } else if (ext == ".ply") {
    vtkSmartPointer<vtkPLYReader> reader;
    reader = vtkSmartPointer<vtkPLYReader>::New();
    reader->SetFileName(fname);
    reader->Update();
    pointset = reader->GetOutput();
  } else if (ext == ".stl") {
    vtkSmartPointer<vtkSTLReader> reader;
    reader = vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(fname);
    reader->Update();
    pointset = reader->GetOutput();
  } else if (ext.length() == 4 && ext.substr(0, 3) == ".vt" && ext[3] != 'k') {
    vtkSmartPointer<vtkXMLGenericDataObjectReader> reader;
    reader = vtkSmartPointer<vtkXMLGenericDataObjectReader>::New();
    reader->SetFileName(fname);
    reader->Update();
    pointset = vtkPointSet::SafeDownCast(reader->GetOutput());
  } else {
    vtkSmartPointer<vtkPolyDataReader> reader;
    reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(fname);
    reader->Update();
    pointset = reader->GetOutput();
  }

  double p[3];
  vtkPoints * const points = pointset->GetPoints();
  this->Reserve(this->Size() + points->GetNumberOfPoints());
  for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++) {
    points->GetPoint(i, p);
    this->Add(irtkPoint(p));
  }

#else
  cerr << "irtkPointSet::ReadVTK: Must be compiled with VTK enabled" << endl;
  exit(1);
#endif
}

void irtkPointSet::ReadVTK(const char *filename)
{
  this->Clear();
  AddVTK(filename);
}

#ifdef HAS_VTK
void irtkPointSet::WriteVTK(const char *filename, vtkAbstractArray *data) const
{

  int i;

  // Create data
  vtkPolyData *dataset = vtkPolyData::New();

  // Create point set
  vtkPoints *points = vtkPoints::New();

  // Convert point set
  for (i = 0; i < this->Size(); i++) {
    double p[3];
    irtkPoint point = this->operator()(i);
    p[0] = point._x;
    p[1] = point._y;
    p[2] = point._z;
    points->InsertPoint(i, p);
  }
  dataset->SetPoints(points);
  if (data) dataset->GetPointData()->AddArray(data);

  // Create writer
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetFileName(filename);
#if VTK_MAJOR_VERSION >= 6
  writer->SetInputData(dataset);
#else
  writer->SetInput(dataset);
#endif
  writer ->Update();
  writer ->Delete();
  dataset->Delete();
  points ->Delete();
}
#else
void irtkPointSet::WriteVTK(const char *filename, void *data) const
{
  cerr << "irtkPointSet::WriteVTK: Must be compiled with VTK enabled" << endl;
  exit(1);
}
#endif

irtkPoint irtkPointSet::ClosestPoint(irtkPoint& p_input){
	int i,j;
	double mindistance = 100000;
	double tmpdistance;
	irtkPoint e;

	if (this->Size() == 0) {
		cerr << "irtkPointSet::ClosestPoint(): No points in pointset" << endl;
		return e;
	}
	j = 0;
	for (i = 0; i < this->Size(); i++) {
		tmpdistance = _data[i].Distance(p_input);
		if(tmpdistance < mindistance){
			mindistance = tmpdistance;
			j = i;
		}
	}
	e = _data[j];
	return e ;
}

double irtkPointSet::PointDistance(irtkPoint& p_input){
	int i;
	double mindistance = 100000;
	double tmpdistance;

	if (this->Size() == 0) {
		cerr << "irtkPointSet::PointDistant(): No points in pointset" << endl;
		return 0;
	}
	for (i = 0; i < this->Size(); i++) {
		tmpdistance = _data[i].Distance(p_input);
		if(tmpdistance < mindistance){
			mindistance = tmpdistance;
		}
	}
	return mindistance ;
}

irtkPoint irtkPointSet::StandardDeviationEllipsoid() const
{
  int i;
  irtkPoint p;
  irtkPoint e;

  if (this->Size() == 0) {
    cerr << "irtkPointSet::StandardDeviationEllipsoid(): No points in pointset" << endl;
    return e;
  }
  p = this->CenterOfGravity();

  for (i = 0; i < this->Size(); i++) {
    e._x += ((*this)(i)._x-p._x)*((*this)(i)._x-p._x);
    e._y += ((*this)(i)._y-p._y)*((*this)(i)._y-p._y);
    e._z += ((*this)(i)._z-p._z)*((*this)(i)._z-p._z);
  }
  e._x = sqrt(e._x / ((double)this->Size()-1));
  e._y = sqrt(e._y / ((double)this->Size()-1));
  e._z = sqrt(e._z / ((double)this->Size()-1));
  return e ;
}

int irtkPointSet::IsInside(double x, double y) const
{
  int i;

  // compute no of points
  if (_n == 0) return false;

  // compute centre
  double cx = 0;
  double cy = 0;
  for (i = 0; i < _n; ++i) {
    cx += _data[i]._x;
    cy += _data[i]._y;
  }
  cx /= _n;
  cy /= _n;

  // compute a point outside the polygon.
  double ox = 0, oy = 0;
  for (i = 0; i < _n; ++i) {
    double tmp;

    tmp = _data[i]._x-cx;
    if (tmp<0) tmp = -tmp;
    if (tmp>ox) ox = tmp;

    tmp = _data[i]._y-cy;
    if (tmp<0) tmp = -tmp;
    if (tmp>oy) oy = tmp;
  }
  ox = cx + ox + oy + 1;
  oy = cy + ox + oy + 1;

  // count crossings.
  int crossings = 0;
  for (i = 0; i < _n; ++i) {
    crossings += intersection(_data[i]._x, _data[i]._y, _data[(i+1)%_n]._x,
                              _data[(i+1)%_n]._y, ox, oy, x, y);
  }

  // inside iff there was an odd number of crossings.
  return crossings % 2 != 0;
}

