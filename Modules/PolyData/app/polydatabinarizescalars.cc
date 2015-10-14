#if (defined HAS_VTK)

#include <irtkImage.h>
//#include <nr.h>

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>

char *input_name = NULL, *output_name = NULL;
char *scalar_name = NULL;

void usage()
{
  cerr << "Usage: polydatabinarizescalars [surfaceIn] [surfaceOut] [thresh] [val1] [val2] <options>" << endl;
  cerr << "" << endl;
  cerr << "Options: " << endl;
  cerr << "-name : Name of scalars to use." << endl;
  cerr << "" << endl;

  exit(1);
}

int main(int argc, char **argv)
{
  bool ok;
  double val;
  int i, count;
  int noOfPoints;

  double val1, val2, thresh;

  if (argc < 6){
    usage();
  }

  // Parse image
  input_name  = argv[1];
  argc--;
  argv++;
  output_name  = argv[1];
  argc--;
  argv++;
  thresh = atof(argv[1]);
  argc--;
  argv++;
  val1 = atof(argv[1]);
  argc--;
  argv++;
  val2 = atof(argv[1]);
  argc--;
  argv++;

  // Parse remaining arguments
  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-name") == 0)) {
      argc--;
      argv++;
      scalar_name = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if (!ok){
      cerr << "Cannot parse argument " << argv[1] << endl;
      usage();
    }
  }

  vtkPolyData *polys = vtkPolyData::New();

  // Read surface
  vtkPolyDataReader *surface_reader = vtkPolyDataReader::New();
  surface_reader->SetFileName(input_name);
  surface_reader->Modified();
  surface_reader->Update();
  polys = surface_reader->GetOutput();

  noOfPoints= polys->GetNumberOfPoints();

//  vtkFloatArray *scalars = vtkFloatArray::New();
//  scalars = (vtkFloatArray*) polys->GetPointData()->GetScalars();
  vtkFloatArray *scalars;

  if (scalar_name != NULL){
    scalars = (vtkFloatArray*) polys->GetPointData()->GetArray(scalar_name);
    if (scalars == NULL){
      cerr << "Cannot retrieve scalars : " << scalar_name << endl;
      exit(1);
    }
  } else {
    scalars = (vtkFloatArray*) polys->GetPointData()->GetScalars();
    cerr << "Using scalars :  " << scalars->GetName() << endl;
  }

  count = 0;

  for (i = 0; i < noOfPoints; ++i){
    val = scalars->GetTuple1(i);
    val = val > thresh ? val2 : val1 ;
    scalars->SetTuple1(i, val);
  }
  scalars->Modified();

  polys->GetPointData()->AddArray(scalars);
  polys->Modified();

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
#if VTK_MAJOR_VERSION >= 6
  writer->SetInputData(polys);
#else
  writer->SetInput(polys);
#endif
  writer->SetFileName(output_name);
  writer->SetFileTypeToBinary();
  writer->Write();

  return 0;
}


#else
int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif

