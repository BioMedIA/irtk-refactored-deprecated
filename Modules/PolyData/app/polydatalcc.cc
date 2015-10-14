#ifdef HAS_VTK

#include <irtkImage.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkCleanPolyData.h>


char *input_name = NULL;
char *output_name = NULL;

void usage(){
  cerr << "polydatalcc [input] [output] <-options>"<<endl;
  cerr << "Keep only largest connected component of the input polydata." << endl;
  cerr << "Options:" << endl;
  cerr << "-number [n] : number of pieces to return, starting with the largest." << endl;
  cerr << "-info       : Only give information without writing output.  A name for " << endl;
  cerr << "              The output is still required." << endl;
  cerr << "  " << endl;
  exit(1);
}

int main(int argc, char **argv ){

  if (argc <3 ){
    usage();
  }

  int i, noOfPoints, regionPointCount;
  bool ok;
  int number = -1;
  int noOfExtractedRegions;
  bool writeOutput = true;
  char buffer[256];
  int dotPos;


  input_name = argv[1];
  argv++;
  argc--;
  output_name = argv[1];
  argv++;
  argc--;

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-number") == 0)){
      argc--;
      argv++;
      number = atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-info") == 0)){
      argc--;
      argv++;
      writeOutput = false;
      ok = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  vtkPolyData* input;
  vtkPolyDataReader* reader;
  vtkPolyDataConnectivityFilter* connFilter;
  vtkCleanPolyData* cleaner;

  cout << "Reading file " << input_name << endl;

  input = vtkPolyData::New();
  reader    = vtkPolyDataReader::New();
  connFilter = vtkPolyDataConnectivityFilter::New();
  cleaner   = vtkCleanPolyData::New();


  reader->SetFileName(input_name);
  reader->Update();
  input = reader->GetOutput();

  noOfPoints = input->GetNumberOfPoints();
  cout << "Total points : " << noOfPoints << endl;

#if VTK_MAJOR_VERSION >= 6
  connFilter->SetInputData(input);
#else
  connFilter->SetInput(input);
#endif
  connFilter->SetExtractionModeToAllRegions();
  connFilter->Update();
  connFilter->SetExtractionModeToSpecifiedRegions();

  noOfExtractedRegions = connFilter->GetNumberOfExtractedRegions();
  cout << "Extracted " << noOfExtractedRegions << " regions." << endl;

  if (number > 0){
    if (number > noOfExtractedRegions){
      number = noOfExtractedRegions;
      cout << "Reset number of regions requested to " << number << endl;
    }
  } else {
    // Default.
    number = 1;
  }

  vtkPolyData *surface;

  for (i = 0; i < noOfExtractedRegions; ++i){
    connFilter->DeleteSpecifiedRegion(i);
  }
  connFilter->Update();

  // Assuming '.vtk' suffix
  dotPos = strlen(output_name) - 4;
  output_name[dotPos] = 0;

  for (i = 0; i < number; ++i){
    connFilter->AddSpecifiedRegion(i);
    connFilter->Update();

    vtkCleanPolyData *cleaner = vtkCleanPolyData::New();
#if VTK_MAJOR_VERSION >= 6
    cleaner->SetInputData(connFilter->GetOutput());
#else
    cleaner->SetInput(connFilter->GetOutput());
#endif
    cleaner->Update();

    regionPointCount = cleaner->GetOutput()->GetNumberOfPoints();
    cout << "Region : " << i + 1 << " : " << regionPointCount << " = ";
    cout << 100.0 * regionPointCount  / ((double) noOfPoints) << "%" << endl;


    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    if (writeOutput == true){

      if (number > 1){
        sprintf(buffer, "%s_%d.vtk", output_name, i+1);
      } else {
        sprintf(buffer, "%s.vtk", output_name);
      }

      writer->SetFileName(buffer);
#if VTK_MAJOR_VERSION >= 6
      writer->SetInputData(cleaner->GetOutput());
#else
      writer->SetInput(cleaner->GetOutput());
#endif
      writer->Write();
    }

    connFilter->DeleteSpecifiedRegion(i);
    connFilter->Update();
  }

  return 0;
}

#else
int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif
