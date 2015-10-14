/*
 * volumeToMesh.cc
 *
 *  Created on: 15 Oct 2013
 *      Author: rob
 */

#include <irtkVolumeToMesh.h>
#include <irtkImage.h>
#include <irtkResampling.h>
#include <vtkPolyData.h>


char *out_name = NULL, *in_name = NULL,* meshName = NULL;
int ok=false,iterations,maxIterations=100,minIterations=100,levels=1,refineMaxIterations=0,
lowPassIterations=25;
double epsilon=0.0001,smoothing=0,edgeLength=2,refinement=0.5,passBand=0.5,resampleResolution=1.0,boundaryThreshold=0.5;
bool screenshotsOn=false,selfIntersectionOn=false,lowPassOn=false,refineOn=false,resampleOn=false,
closeBinaryVolume=false, gaussianInterpolationOn=false;

enum InitialMesh {input, hull, sphere};

void usage()
{
  cerr << "Usage: [volume_in] [mesh_out] <options> \n" << endl;
  cerr << "<-smoothing double>      Value from 0 to 1 (default 0.5)" << endl;
  cerr << "<-edge_length double>    Set the desired edge length in mm"   << endl;
  cerr << "<-max_iterations int>    Set the maximum number of iterations"   << endl;
  cerr << "<-min_iterations int>    Set the minimum number of iterations"   << endl;
  cerr << "<-epsilon double>        Set the minimum average error improvement (default 0.0001)"   << endl;
  cerr << "<-intersect>             Disallow mesh from self-intersecting"   << endl;
  cerr << "<-low_pass>              Enable low pass filtering"   << endl;
  cerr << "<-low_pass [iterations] [band]>  Enable low pass filtering"   << endl;
  cerr << "  and, set number of iterations and pass band value (0.0-2.0)"<< endl;
  cerr << "<-refine  [iterations]>  Set number of refinement iterations"   << endl;
  cerr << "<-initial_mesh file.vtk|hull|sphere>  Set initial mesh to deform" << endl;
  cerr << "<-screenshots dir>       Save a screenshot after each iteration" << endl;
  cerr << "<-resample resolution>   New input voxel size in mm" << endl;
  cerr << "<-threshold double>      Set isovalue for surface" << endl;
  cerr << "<-close>                 Close binary volume" << endl;
  cerr << "<-gaussian>              Interpolate distance field using a gaussian (default linear)" << endl;
  exit(1);

}
int main(int argc, char **argv){

  InitialMesh initialMesh=sphere;

  if (argc < 3) {
        usage();
        exit(1);
    }


  in_name = argv[1];
  argv++;
  argc--;

  out_name = argv[1];
  argv++;
  argc--;

  char * screenshotDir = "";

  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-smoothing") == 0)){
      argc--;
      argv++;
      smoothing = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-edge_length") == 0)){
      argc--;
      argv++;
      edgeLength = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-max_iterations") == 0)){
      argc--;
      argv++;
      maxIterations= atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-min_iterations") == 0)){
      argc--;
      argv++;
      minIterations= atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-screenshots") == 0)){
      argc--;
      argv++;
      screenshotDir= argv[1];
      screenshotsOn=true;
      argc--;
      argv++;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-initial_mesh") == 0)){
      argc--;
      argv++;

      if (strcmp(argv[1], "hull") == 0){
        initialMesh=hull;
      }

      else if (strcmp(argv[1], "sphere") == 0){
        initialMesh=sphere;
      }

      else{
        initialMesh=input;
        meshName= argv[1];
      }
      argc--;
      argv++;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-epsilon") == 0)){
      argc--;
      argv++;
      epsilon=atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-intersect") == 0)){
      argc--;
      argv++;
      selfIntersectionOn=true;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-low_pass") == 0)){
      argc--;
      argv++;
      lowPassOn=true;



      if (argc > 2 && atoi(argv[1])!=0){
        lowPassIterations=atoi(argv[1]);
        argc--;
        argv++;
        passBand=atof(argv[1]);
        argc--;
        argv++;
      }

      ok = true;
    }


    if ((ok == false) && (strcmp(argv[1], "-refine") == 0)){
      refineOn=true;
      argc--;
      argv++;
      refineMaxIterations=atoi(argv[1]);
      argc--;
      argv++;
      ok = true;
    }


    if ((ok == false) && (strcmp(argv[1], "-resample") == 0)){
      argc--;
      argv++;
      resampleResolution=atof(argv[1]);
      resampleOn=true;
      argc--;
      argv++;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-threshold") == 0)){
      argc--;
      argv++;
      boundaryThreshold=atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-close") == 0)){
      argc--;
      argv++;
      closeBinaryVolume=true;
      argc--;
      argv++;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-gaussian") == 0)){
      argc--;
      argv++;
      gaussianInterpolationOn=true;
      argc--;
      argv++;
      ok = true;
    }


    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  irtkRealImage image;
  image.Read(in_name);

  irtkVolumeToMesh v2m;

  vtkSmartPointer<vtkPolyData> outMesh
    = vtkSmartPointer<vtkPolyData>::New();


  if(resampleOn){//resample input volume
    cout << "resampling...";
    cout.flush();
    irtkImageFunction *interpolator = new irtkLinearInterpolateImageFunction;
    irtkResampling<irtkRealPixel> resampling(resampleResolution, resampleResolution, resampleResolution);
    resampling.SetInput(&image);
    resampling.SetOutput(&image);
    resampling.SetInterpolator(interpolator);
    resampling.Run();
    cout << "done" << endl;
  }
  v2m.SetBoundaryThreshold(boundaryThreshold);
  v2m.SetMaxIterations(maxIterations);
  v2m.SetMinIterations(minIterations);
  v2m.SetEdgeLength(edgeLength);
  v2m.SetNumberOfLevels(levels);
  v2m.SetImprovementThreshold(epsilon);
  if (closeBinaryVolume){
   v2m.CloseBinaryVolumeOn();
  }

  v2m.SetInput(&image);
  if(smoothing>0 && smoothing<=1){
     v2m.SmoothingOn();
     v2m.SetSmoothingValue(smoothing);
  }
  else { v2m.SmoothingOff(); }

  if (lowPassOn){
    v2m.LowPassFilterOn();
    v2m.SetLowPassFilterBand(passBand);
    v2m.SetLowPassFilterIterations(lowPassIterations);
  }

  if (refineOn){
    v2m.FinalRefinementOn();
    v2m.SetRefinementIterations(refineMaxIterations);
  }

  if (selfIntersectionOn) { v2m.SelfIntersectionOn(); }
  else { v2m.SelfIntersectionOff(); }

  if (initialMesh==sphere) { v2m.SetInitialMeshToSphere(); }
  else if (initialMesh==hull) { v2m.SetInitialMeshToHull(); }
  else if (initialMesh==input) {
    vtkSmartPointer<vtkPolyDataReader> polyReader = vtkPolyDataReader::New();
    polyReader->SetFileName(meshName);
    polyReader->Update();
    v2m.SetInitialMesh(polyReader->GetOutput());
  }

  if (screenshotsOn) {
    v2m.ScreenShotsOn();
    v2m.SetScreenShotDirectory(screenshotDir);
  }

  if (gaussianInterpolationOn){
    v2m.GaussianInterpolationOn();
  }


  outMesh=v2m.GetOuput();

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  SetVTKInput(writer, outMesh);
  double vertex [3];
  outMesh->GetPoint(0,vertex);
  cout << vertex[0] << "," << vertex[1] <<  "," << vertex[2] << endl;


  writer->SetFileName(out_name);
  writer->SetFileTypeToBinary();
  int success = writer->Write();
  if (success == 0) exit(1);
  cout << "done" << endl;
  return 0;
}
