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

#include <irtkDataOp.h>

#include <irtkImage.h>
#ifdef HAS_VTK
#include <irtkVTK.h>
#include <vtkSmartPointer.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetWriter.h>
#include <vtkXMLGenericDataObjectReader.h>
#include <vtkXMLDataSetWriter.h>
#include <vtkImageData.h>
#include <vtkPointSet.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#endif

#include <memory>

namespace irtk { namespace data {


// -----------------------------------------------------------------------------
DataFileType FileType(const char *name)
{
  string ext = Extension(name);
  if (ext == ".vtk") return LEGACY_VTK;
  if (ext.length() == 4 && ext[0] == '.' && ext[1] == 'v') return XML_VTK;
  return IMAGE;
}

// -----------------------------------------------------------------------------
#ifdef HAS_VTK
int Read(const char *name, double *&data, int *dtype, irtkImageAttributes *attr,
         vtkSmartPointer<vtkDataSet> *dataset, const char *scalar_name)
#else
int Read(const char *name, double *&data, int *dtype, irtkImageAttributes *attr)
#endif
{
  int n = 0;
  data = NULL;
  if (attr) *attr = irtkImageAttributes();
  int ftype = FileType(name);
  switch (ftype) {
#ifdef HAS_VTK
  if (dataset) *dataset = NULL;
    case LEGACY_VTK:
    case XML_VTK: {
      vtkSmartPointer<vtkDataArray> scalars;
      if (ftype == LEGACY_VTK) {
        vtkSmartPointer<vtkDataSetReader> reader;
        reader = vtkSmartPointer<vtkDataSetReader>::New();
        reader->SetFileName(name);
        reader->Update();
        vtkDataSet *output = reader->GetOutput();
        if (output) {
          if (scalar_name) {
            scalars = output->GetPointData()->GetArray(scalar_name);
          } else {
            scalars = output->GetPointData()->GetScalars();
          }
        }
        if (dataset) *dataset = output;
      } else {
        vtkSmartPointer<vtkXMLGenericDataObjectReader> reader;
        reader = vtkSmartPointer<vtkXMLGenericDataObjectReader>::New();
        reader->SetFileName(name);
        reader->Update();
        vtkDataSet *output = vtkDataSet::SafeDownCast(reader->GetOutput());
        if (output) {
          if (scalar_name) {
            scalars = output->GetPointData()->GetArray(scalar_name);
          } else {
            scalars = output->GetPointData()->GetScalars();
          }
        }
        if (dataset) *dataset = output;
      }
      if (!scalars) {
        cerr << "Failed to read VTK dataset! Type is either not supported or dataset has no scalar point data." << endl;
        cerr << "Use -scalars option to specify the name of a point data array to use instead." << endl;
        exit(1);
      }
      if (dtype) *dtype = FromVTKDataType(scalars->GetDataType());
      n = static_cast<int>(scalars->GetNumberOfTuples()) * scalars->GetNumberOfComponents();
      if (n == 0) {
        cerr << "VTK dataset has empty scalar point data!" << endl;
        exit(1);
      }
      data = new double[n];
      double *p = data;
      for (vtkIdType i = 0; i < scalars->GetNumberOfTuples(); ++i) {
        for (int j = 0; j < scalars->GetNumberOfComponents(); ++j, ++p) {
          *p = scalars->GetComponent(i, j);
        }
      }
    } break;
#else // HAS_VTK
    case LEGACY_VTK:
    case XML_VTK:
      cerr << "Cannot read VTK files when compiled without VTK library!" << endl;
      exit(1);
#endif // HAS_VTK
    case IMAGE: {
      std::unique_ptr<irtkBaseImage> image(irtkBaseImage::New(name));
      if (attr) *attr = image->Attributes();
      if (dtype) *dtype = image->GetDataType();
      n = image->NumberOfVoxels();
      data = new double[n];
      for (int i = 0; i < n; ++i) data[i] = image->GetAsDouble(i);
    } break;
    default:
      cerr << "Unsupported input data file: " << name << endl;
      exit(1);
  }
  return n;
}

// -----------------------------------------------------------------------------
void Write::Process(int n, double *data, bool *)
{
  if (_FileName.empty()) {
    cerr << "Output file name not set!" << endl;
    exit(1);
  }
  int type = FileType(_FileName.c_str());
  switch (type) {
#ifdef HAS_VTK
    case LEGACY_VTK:
    case XML_VTK: {
      if (!_DataSet) {
        cerr << "Cannot write data sequence to VTK file when input was not a VTK file itself!" << endl;
        exit(1);
      }
      vtkDataArray *input_scalars;
      if (_ArrayName.empty()) input_scalars = _DataSet->GetPointData()->GetScalars();
      else                    input_scalars = _DataSet->GetPointData()->GetArray(_ArrayName.c_str());
      if (input_scalars == NULL) {
        cerr << "Invalid output array name " << _ArrayName << endl;
        exit(1);
      }
      vtkIdType ntuples = input_scalars->GetNumberOfTuples();
      int       m       = input_scalars->GetNumberOfComponents();
      if (n != static_cast<int>(ntuples * m)) {
        cerr << "Cannot write data sequence to file! Length of data sequence changed." << endl;
        exit(1);
      }
      vtkSmartPointer<vtkDataArray> output_scalars = NewVTKDataArray(ToVTKDataType(_DataType));
      if (!_ArrayName.empty()) output_scalars->SetName(_ArrayName.c_str());
      output_scalars->SetNumberOfComponents(m);
      output_scalars->SetNumberOfTuples(ntuples);
      for (vtkIdType i = 0; i < ntuples; ++i) {
        for (int j = 0; j < m; ++j, ++data) {
          output_scalars->SetComponent(i, j, *data);
        }
      }
      if (input_scalars == _DataSet->GetPointData()->GetScalars()) {
        _DataSet->GetPointData()->SetScalars(output_scalars);
      } else {
        _DataSet->GetPointData()->RemoveArray(_ArrayName.c_str());
        _DataSet->GetPointData()->AddArray(output_scalars);
      }
      if (type == LEGACY_VTK) {
        vtkSmartPointer<vtkDataSetWriter> writer;
        writer = vtkSmartPointer<vtkDataSetWriter>::New();
        writer->SetFileName(_FileName.c_str());
        SetVTKInput(writer, _DataSet);
        writer->Write();
      } else {
        vtkSmartPointer<vtkXMLDataSetWriter> writer;
        writer = vtkSmartPointer<vtkXMLDataSetWriter>::New();
        writer->SetFileName(_FileName.c_str());
        SetVTKInput(writer, _DataSet);
        writer->Write();
      }
    } break;
#else
    case LEGACY_VTK:
    case XML_VTK:
      cerr << "Cannot write data sequence to VTK file when compiled without VTK" << endl;
      exit(1);
#endif
    case IMAGE: {
      if (_Attributes.NumberOfLatticePoints() == 0) {
        cerr << "Cannot write data sequence to image file when input was not an image file itself!" << endl;
        exit(1);
      }
      if (_Attributes.NumberOfLatticePoints() != n) {
        cerr << "Cannot write data sequence to file! Length of data sequence changed." << endl;
        exit(1);
      }
      std::unique_ptr<irtkImage> image(irtkImage::New(_DataType));
      image->Initialize(_Attributes);
      for (int i = 0; i < n; ++i) image->PutAsDouble(i, data[i]);
      image->Write(_FileName.c_str());
    } break;
    default:
      cerr << "Unknown output file type: " << _FileName << endl;
      exit(1);
  }
}


} } // namespace irtk::data

