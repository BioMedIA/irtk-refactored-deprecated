#ifdef HAS_VTK

#include <irtkImage.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>

#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>

void usage(){
  cerr << "polydatadeletescalars [input] [output] <-name arrayName | -all>" << endl;
  cerr << "" << endl;
  exit(1);
}

char *input_name = NULL;
char *output_name = NULL;
char *scalar_name = NULL;

int main(int argc, char **argv)
{

  bool ok, deleteAll = false;

  if (argc < 4){
    usage();
  }

  input_name  = argv[1];
  argc--;
  argv++;
  output_name  = argv[1];
  argc--;
  argv++;

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-name") == 0)){
      argc--;
      argv++;
      scalar_name  = argv[1];
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-all") == 0)){
      argc--;
      argv++;
      deleteAll = true;
      ok = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  if (deleteAll == false && scalar_name == NULL){
  	cerr << "No array(s) specified." << endl;
  	usage();
  }

  if (deleteAll == true){
  	cout << "Deleting all arrays." << endl;
  	scalar_name = NULL;
  }


  // Read surface
  vtkPolyDataReader *surface_reader = vtkPolyDataReader::New();
  surface_reader->SetFileName(input_name);
  surface_reader->Modified();
  surface_reader->Update();

  vtkPolyData *surface = surface_reader->GetOutput();

  int i, noOfArrays;
  int success = false;

  noOfArrays = surface->GetPointData()->GetNumberOfArrays();

  if (noOfArrays < 1){
  	cout << "Surface has no arrays." << endl;
  	exit(1);
  }

  if (deleteAll == true){
  	while(surface->GetPointData()->GetNumberOfArrays() > 0){
  		vtkFloatArray *currArray;
  		currArray = (vtkFloatArray*) surface->GetPointData()->GetArray(0);
			cout << "Deleting array : " << currArray->GetName() << endl;
  		surface->GetPointData()->RemoveArray(currArray->GetName());
  		success = true;
  	}
  } else {
  	for (i = 0; i < noOfArrays; i++){
  		vtkFloatArray *currArray;
  		currArray = (vtkFloatArray*) surface->GetPointData()->GetArray(i);

  		if (strcmp(currArray->GetName(), scalar_name) == 0){
  			cout << "Deleting array : " << currArray->GetName() << endl;
  			surface->GetPointData()->RemoveArray(scalar_name);
  			success = true;
  			break;
  		}
  	}
  }

  if (success) {
    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    cout << "Writing surface to file: " << output_name << endl;
#if VTK_MAJOR_VERSION >= 6
    writer->SetInputData(surface);
#else
    writer->SetInput(surface);
#endif
    writer->SetFileName(output_name);
    writer->SetFileTypeToBinary();
    writer->Write();
  } else {
  	if (deleteAll == false){
  		cerr << "No such scalars : " << scalar_name << endl;
  	}
  }
}

#else
int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
