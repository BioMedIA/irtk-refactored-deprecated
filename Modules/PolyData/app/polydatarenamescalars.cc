/*
 * polydatarenamescalars.cc
 *
 *  Created on: Jan 11, 2012
 *      Author: paul
 */

#ifdef HAS_VTK

#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>

#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>

void usage(){
  cerr << "polydatarenamescalars [input] [output] [oldName] [newName]" << endl;
  cerr << "" << endl;
  exit(1);
}

char *input_name = NULL;
char *output_name = NULL;
char *scalar_name_old = NULL;
char *scalar_name_new = NULL;

int main(int argc, char **argv)
{

  bool ok;

  if (argc < 5){
    usage();
  }

  input_name  = argv[1];
  argc--;
  argv++;
  output_name  = argv[1];
  argc--;
  argv++;
  scalar_name_old = argv[1];
  argc--;
  argv++;
  scalar_name_new = argv[1];
  argc--;
  argv++;


  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-someOption") == 0)){
      argc--;
      argv++;
      // Do stuff and maybe
//      argc--;
//      argv++;
      ok = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Read surface
  vtkPolyDataReader *surface_reader = vtkPolyDataReader::New();
  surface_reader->SetFileName(input_name);
  surface_reader->Modified();
  surface_reader->Update();

  vtkPolyData *surface = surface_reader->GetOutput();


  if (surface->GetPointData()->GetNumberOfArrays() < 1){
  	cout << "Surface has no arrays." << endl;
  	exit(1);
  }


  vtkFloatArray *scalars;
  int ind;

  scalars = static_cast<vtkFloatArray*> (surface->GetPointData()->GetArray(scalar_name_old, ind));

  if (ind == -1 || scalars == NULL){
  	cerr << "Scalars unavailable with name " << scalar_name_old << endl;
  	exit(1);
  }

  scalars->SetName(scalar_name_new);
  surface->GetPointData()->SetActiveScalars(scalar_name_new);


  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
#if VTK_MAJOR_VERSION >= 6
  writer->SetInputData(surface);
#else
  writer->SetInput(surface);
#endif
  writer->SetFileName(output_name);
  writer->SetFileTypeToBinary();
  writer->Write();


}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif

