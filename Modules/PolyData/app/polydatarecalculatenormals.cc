#if (defined HAS_VTK)


#include <irtkImage.h>

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataNormals.h>

char *in_name  = NULL, *out_name = NULL;

void usage()
{
  cerr << "Usage: polydatarecalculatenormals [input] [output] <options>\n" << endl;
  cerr << " " << endl;
  cerr << "Options: " << endl;
  cerr << " " << endl;
  cerr << "-split : Allow sharp edges to be split. Not set by default." << endl;
  cerr << "-auto  : Try and guess the correct direction for all normals so " << endl;
  cerr << "         that they all point `outward'.  Only works if the surface" << endl;
  cerr << "         is completely closed (no boundary) and is Hausdorf." << endl;
  cerr << " " << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  bool ok;
  int splitting = false;
  int autoOrient = false;

  if (argc < 3){
    usage();
  }

  // Parse source and target point lists
  in_name  = argv[1];
  argc--;
  argv++;
  out_name = argv[1];
  argc--;
  argv++;

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-split") == 0)){
      argc--;
      argv++;
      splitting = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-auto") == 0)){
      argc--;
      argv++;
      autoOrient = true;
      ok = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  // Read the input surface.
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(in_name);
  reader->Update();

  vtkPolyData *input = vtkPolyData::New();
  input = reader->GetOutput();

  vtkPolyDataNormals *normalsFilter = vtkPolyDataNormals::New();

#if VTK_MAJOR_VERSION >= 6
  normalsFilter->SetInputData(input);
#else
  normalsFilter->SetInput(input);
#endif

  if (splitting == true){
    cerr << "Setting splitting to on" << endl;
    normalsFilter->SplittingOn();
  } else {
    normalsFilter->SplittingOff();
  }

  if (autoOrient == true){
    cout << "Setting auto orientation estimation to on." << endl;
    normalsFilter->AutoOrientNormalsOn();
  } else {
    normalsFilter->AutoOrientNormalsOff();
  }

  normalsFilter->SetConsistency(1);

  normalsFilter->Modified();
  normalsFilter->Update();

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetFileName(out_name);
#if VTK_MAJOR_VERSION >= 6
  writer->SetInputData(normalsFilter->GetOutput());
#else
  writer->SetInput(normalsFilter->GetOutput());
#endif
  writer->Modified();
  writer->SetFileTypeToBinary();
  writer->Write();

}

#else

#include <irtkImage.h>

int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
