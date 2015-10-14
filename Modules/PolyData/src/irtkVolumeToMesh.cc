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


#include<irtkVolumeToMesh.h>
#include<math.h>
#include<irtkErosion.h>
#include<irtkDilation.h>

#include <vtkMath.h>
#include <vtkTriangle.h>
#include <vtkFloatArray.h>
#include <vtkArray.h>
#include <vtkDataArray.h>
#include <vtkCellArray.h>
#include <vtkPlatonicSolidSource.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCleanPolyData.h>
#include <vtkTriangleFilter.h>
#include <irtkEuclideanDistanceTransform.h>
#include <irtkParticleSampleSphere.h>
#include <vtkPointSet.h>

#include <vtkSmoothPolyDataFilter.h>
#include <irtkImage.h>
#include <irtkImageFunction.h>
#include <irtkBaseImage.h>
#include <vtkButterflySubdivisionFilter.h>
#include <math.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkDoubleArray.h>
#include <vtkDecimatePro.h>
#include <vtkQuadricDecimation.h>
#include <vtkOBBTree.h>
#include <vtkModifiedBSPTree.h>
#include <vtkGraphicsFactory.h>
#include <vtkOpenGLRenderer.h>
#include <vtkOpenGLRenderWindow.h>
#include <vtkDelaunay3D.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkHull.h>

#include <vtkPolyDataToImageStencil.h>
#include <vtkImageStencil.h>
#include <vtkStructuredPoints.h>
#include <vtkImageData.h>
#include <irtkResampling.h>
#include <vtkMarchingCubes.h>
#include <vtkImageResample.h>
#include <vtkReverseSense.h>

//Constructor and Destructor

irtkVolumeToMesh::irtkVolumeToMesh(){


  _selfIntersectBSPTree = vtkSmartPointer<vtkModifiedBSPTree>::New();
  _selfIntersectBSPTree->SetMaxLevel(100);
  _selfIntersectBSPTree->SetNumberOfCellsPerNode(1);
  _selfIntersectBSPTree->AutomaticOn();
  _intersectionOn=true;

  _improvementThreshold=0.0001;
  _boundaryThreshold=0.5;
  _maxIterations=200;
  _subdivisions=1;
  _averageEdgeLengthMM=0.5;
  _maxStepSizeMM=1;
  _selfIntersectionMM=3;

  _levels=1;
  _smoothing=0.5;
  _minNumPoints=1000;
  _multiResolutionOn=false;
  _normalFilter = vtkSmartPointer<vtkPolyDataNormals>::New();
  _remeshFilter = new irtkRemesher();
  _smoothingFilter = new irtkPolyDataSmoothing();
  _screenShotsOn=false;
  _screenShotDir="./";
  _selfIntersectionOn=true;
  _lowPassFilterOn=false;
  _finalRefinementOn=false;
  _lowPassFilterBand=0.5;
  _lowPassFilterIterations=100;
  _refinementSmoothing=0.5;
  _smoothingOn=true;
  _lowPassFilter = vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
  _flipNormals=false;
  _refinementIterations=0;
  _closeBinaryVolume=false;
  _gaussianInterpolationOn=false;
}





irtkVolumeToMesh::~irtkVolumeToMesh(){
  delete _remeshFilter;
  delete _smoothingFilter;
  delete _interpolator;
}



//Get and Set Methods
void irtkVolumeToMesh::SetMaxIterations(int maxIterations){
  _maxIterations=maxIterations;
}

void irtkVolumeToMesh::SetMinIterations(int minIterations){
  _minIterations=minIterations;
}


void irtkVolumeToMesh::SetEdgeLength(double averageEdgeLength){
  _averageEdgeLengthMM=averageEdgeLength;
}

void irtkVolumeToMesh::SetImprovementThreshold(double improvementThreshold){
  _improvementThreshold=improvementThreshold;
}

void irtkVolumeToMesh::SetNumberOfLevels(int levels){
  _levels=levels;
}

void irtkVolumeToMesh::ScreenShotsOn(){
  _screenShotsOn=true;
}


void irtkVolumeToMesh::GaussianInterpolationOn(){
  _gaussianInterpolationOn = true;
}

void irtkVolumeToMesh::GaussianInterpolationOff(){
  _gaussianInterpolationOn = false;
}

void irtkVolumeToMesh::ScreenShotsOff(){
  _screenShotsOn=false;
}

void irtkVolumeToMesh::SetScreenShotDirectory(char * screenShotDir){
  _screenShotDir=screenShotDir;
}

void irtkVolumeToMesh::SelfIntersectionOn(){
  _selfIntersectionOn=true;
}
void irtkVolumeToMesh::SelfIntersectionOff(){
  _selfIntersectionOn=false;
}


void irtkVolumeToMesh::LowPassFilterOn(){
  _lowPassFilterOn=true;
}
void irtkVolumeToMesh::LowPassFilterOff(){
  _lowPassFilterOn=false;
}

void irtkVolumeToMesh::SetLowPassFilterBand(double lowPassFilterBand){
  _lowPassFilterBand=lowPassFilterBand;
}

void irtkVolumeToMesh::SetLowPassFilterIterations(int lowPassFilterIterations){
  _lowPassFilterIterations=lowPassFilterIterations;
}

void irtkVolumeToMesh::FinalRefinementOn(){
  _finalRefinementOn=true;
}

void irtkVolumeToMesh::FinalRefinementOff(){
  _finalRefinementOn=false;
}

void irtkVolumeToMesh::SetSmoothingValue(double smoothing){
  _smoothing=smoothing;
}


void irtkVolumeToMesh::SmoothingOn(){
  _smoothingOn=true;
}

void irtkVolumeToMesh::SmoothingOff(){
  _smoothingOn=false;
}

void irtkVolumeToMesh::SetRefinementIterations(int iterations){
  _refinementIterations=iterations;
}

void irtkVolumeToMesh::SetBoundaryThreshold(double threshold){
  _boundaryThreshold=threshold;
}

void irtkVolumeToMesh::CloseBinaryVolumeOn(){
  _closeBinaryVolume=true;
}

void irtkVolumeToMesh::CloseBinaryVolumeOff(){
  _closeBinaryVolume=false;
}

//Initialization Methods

void irtkVolumeToMesh::Initialize(){
  InitializeDistanceField();
  InitializeFilters();
  InitializeMesh();
  _selfIntersectionVoxels=_selfIntersectionMM/_resolution;
  _maxStepSizeVoxels=_maxStepSizeMM/_resolution;

}


void irtkVolumeToMesh::InitializeFilters(){
  _normalFilter->SplittingOff();
  _normalFilter->AutoOrientNormalsOff() ;
  _normalFilter->ConsistencyOn();
  _remeshFilter->SetMaxEdgeLength(2*_averageEdgeLengthMM/_resolution);
  _remeshFilter->SetMinEdgeLength(0.5*_averageEdgeLengthMM/_resolution);
  _smoothingFilter->SetLambda(_smoothing);
  _smoothingFilter->SetSigma(0);
  if (_gaussianInterpolationOn){
    _interpolator = new irtkGaussianInterpolateImageFunction(_resolution);
  }
  else{
    _interpolator = new irtkLinearInterpolateImageFunction();
  }
  _interpolator->SetInput(&_distanceField);
  _interpolator->Initialize();
  _lowPassFilter->SetPassBand(_lowPassFilterBand);
  _lowPassFilter->SetNumberOfIterations(_lowPassFilterIterations);
  _lowPassFilter->FeatureEdgeSmoothingOff();


}


void  irtkVolumeToMesh::CalculateCentreOfMass(){
  int x,y,z,i,j,k,numVoxelsInside;

  x=_binaryVolume.GetX();
  y=_binaryVolume.GetY();
  z=_binaryVolume.GetZ();

  //determine com
  _centreOfMass[0]=0;
  _centreOfMass[1]=0;
  _centreOfMass[2]=0;
  numVoxelsInside=0;
  for (i=0;i<x;i++){
    for (j=0;j<y;j++){
      for (k=0;k<z;k++){
        if (_binaryVolume(i,j,k)==1){
          numVoxelsInside++;
          _centreOfMass[0]+=i;
          _centreOfMass[1]+=j;
          _centreOfMass[2]+=k;
        }
      }
    }
  }
  _centreOfMass[0]/=numVoxelsInside;
  _centreOfMass[1]/=numVoxelsInside;
  _centreOfMass[2]/=numVoxelsInside;

}

void InitializeParameters(){

}

void irtkVolumeToMesh::InitializeMesh(){
  if (_initialMesh==hull) { InitializeMeshAsConvexHull(); }
  if (_initialMesh==sphere) { InitializeMeshAsBoundingSphere(); }
  if (_initialMesh==input) {   TransformPointsToImageCoordinates(); }
}


void irtkVolumeToMesh::InitializeMeshAsBoundingSphere(){
  _mesh=vtkSmartPointer<vtkPolyData>::New();
  irtkParticleSampleSphere pss;
  int x,y,z,i,j,k,numVoxelsInside;

  int minSamples,subdivisions,minSubdivisions;
  double boundingRadius,distance,furthestDistance,surfaceArea, phiAngle, thetaAngle;
  int numPoints,samples;
  vtkSmartPointer<vtkPoints> points;
  double point [3], furthestPoint [3];

  x=_binaryVolume.GetX();
  y=_binaryVolume.GetY();
  z=_binaryVolume.GetZ();


  furthestDistance=0;
  for (i=0;i<x;i++){
    for (j=0;j<y;j++){
      for (k=0;k<z;k++){
        if (_binaryVolume(i,j,k)==1){
          distance = sqrt(
          pow((double)i-_centreOfMass[0],2.0) +
          pow((double)j-_centreOfMass[1],2.0) +
          pow((double)k-_centreOfMass[2],2.0));
          if (distance > furthestDistance){
            furthestPoint[0]=i;
            furthestPoint[1]=j;
            furthestPoint[2]=k;
            furthestDistance=distance;
          }
        }
      }
    }
  }



  _binaryVolume.ImageToWorld(_centreOfMass[0],_centreOfMass[1],_centreOfMass[2]);
  _binaryVolume.ImageToWorld(furthestPoint[0],furthestPoint[1],furthestPoint[2]);
  boundingRadius = sqrt(
      pow(furthestPoint[0]-_centreOfMass[0],2) +
      pow(furthestPoint[1]-_centreOfMass[1],2) +
      pow(furthestPoint[2]-_centreOfMass[2],2));


  cout << "centre of mass: " << _centreOfMass[0] << "," << _centreOfMass[1] << "," << _centreOfMass[2] << endl;
  cout << "furthestPoint: " << furthestPoint[0] << "," << furthestPoint[1] << "," << furthestPoint[2] << endl;

  cout << "radius: " << boundingRadius << endl;

  surfaceArea=4*M_PI*boundingRadius*boundingRadius;
  cout << "surfaceArea: " << surfaceArea << endl;
  //F=2*SA/(EL^2) (number of faces based on equilateral triangles)
  //V−E+F=2(1-g) (Euler's formula)
  //g=1 (genus: number of handles)
  //V=2+E-F
  //2E=3F (Each face has 3 half edges)
  //V=2+F/2;
  //V=2+SA/(EL^2);

  numPoints=(int)(0.5+2+surfaceArea/(_averageEdgeLengthMM*_averageEdgeLengthMM));
  cout << numPoints << " points" << endl;
  //work out number of subdivisions and number of initial points


  minSamples=25; //max samples with be 4*minSamples

  //calc the max number of subdivisions s based on the minimum sample number v
  // F=f*4^s=f*2^2s
  // F=2*(V-2)
  // 2(V-2)=2(v-2)*2^2s
  // 2^2s=(V-2)/(v-2)
  // s=log2((V-2)/(v-2))/2
  // take floor
  subdivisions=(int)(log2((numPoints-2)/(minSamples-2))/2);
  //calc initial samples based on
  //2^2s=(V-2)/(v-2)
  //v=((V-2)/2^2s)+2


  samples=(int)(0.5+(((double)numPoints-2.0)/pow(4.0,subdivisions))+2);
  //recalc number of points
  //V=(2^2s)*(v-2)
  numPoints=pow(4.0,subdivisions)*(minSamples-2);
  pss.SetNumberOfPoints(samples);


  if (_levels>1){
    subdivisions-=_levels-1;
    _averageEdgeLengthMM*=pow(2.0,_levels-1);
  }

  pss.SetNumberOfSubdivisions(subdivisions);
  pss.SphericalCoordsOn();
  _mesh=pss.GetSphere();




  numPoints=_mesh->GetNumberOfPoints();

  points=_mesh->GetPoints();
  //Scale, Translate vertices and convert back to image coords
  for (i=0; i<numPoints; i++){
    points->GetPoint(i,point);
    //cout << point[0] << "," << point[1] << "," << point[2] << endl;
    for (j=0; j<3; j++)  point[j]=point[j]*boundingRadius+_centreOfMass[j];
    //_binaryVolume.WorldToImage(point[0],point[1],point[2]);
    points->SetPoint(i,point);
    //cout << point[0] << "," << point[1] << "," << point[2] << endl;
  }


  _binaryVolume.WorldToImage(_centreOfMass[0],_centreOfMass[1],_centreOfMass[2]);

  points->GetPoint(0,point);
  _mesh->SetPoints(points);
  cout << "scaled" << boundingRadius << endl;
  TransformPointsToImageCoordinates();
}


void irtkVolumeToMesh::InitializeMeshAsConvexHull(){
  irtkPolyDataSmoothing smoothingFilter;
  vtkSmartPointer<vtkMarchingCubes> mcubes
    = vtkSmartPointer<vtkMarchingCubes>::New();
  vtkSmartPointer<vtkImageData> vtkImage;

  vtkSmartPointer<vtkPolyData> hull = GetHull();

  vtkImage=Mesh2VtkMask(hull);
  mcubes->SetValue(0, 0.5);
  mcubes->ComputeNormalsOn();
  mcubes->ComputeGradientsOff();
  SetVTKInput(mcubes, vtkImage);
  mcubes->Update();
  _mesh = mcubes->GetOutput();
  //smooth hull

  smoothingFilter.SetLambda(1);
  smoothingFilter.SetSigma(0);
  for (int i=0;i<20;i++){
    smoothingFilter.SetInput(_mesh);
    smoothingFilter.Update();
    _mesh=smoothingFilter.GetOutput();
  }

  ComputeMeshNormals();
  _mesh->BuildCells();

/*
  // Write result

  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(_mesh);
  writer->SetFileName("hull.vtk");
  writer->SetFileTypeToBinary();
  writer->Write();
  writer->Delete();

  exit(1);

*/


}




void irtkVolumeToMesh::InitializeDistanceField(){
  irtkRealImage input;
  _edt = new irtkEuclideanDistanceTransform<irtkRealPixel>
    (irtkEuclideanDistanceTransform<irtkRealPixel>::irtkDistanceTransform3D);


  // Threshold image
  irtkRealImage insideDistanceField;
  int x,y,z;

  input=(irtkRealImage)_binaryVolume;
  _distanceField=input;

  cout << "Calculating distance map...";
  cout.flush();
  // Calculate EDT
  _edt->SetInput(&input);
  _edt->SetOutput(&_distanceField);
  _edt->Run();

  //get inverse
  for (z = 0; z < input.GetZ(); z++) {
    for (y = 0; y < input.GetY(); y++) {
      for (x = 0; x < input.GetX(); x++) {
         input(x, y, z)=input(x, y, z)==0;
        }
      }
    }

  _edt->SetInput(&input);
  _edt->SetOutput(&insideDistanceField);
  _edt->Run();

  for (z = 0; z < input.GetZ(); z++) {
    for (y = 0; y < input.GetY(); y++) {
      for (x = 0; x < input.GetX(); x++) {
        _distanceField(x, y, z)  = sqrt(_distanceField(x, y, z)) - sqrt(insideDistanceField(x, y, z));
      }
    }
  }

  //_distanceField.Write("dist_field.nii.gz");
  //_binaryVolume.Write("binary_volume.nii.gz");
  cout << "done" << endl;
}


//convert binary image to vtk point field.
//step 1 erode
//step 2 difference

vtkSmartPointer<vtkPolyData> irtkVolumeToMesh::BinaryVolume2VTKSurfacePoints(){
  int i,j,k,x,y,z;
  irtkGreyImage surfaceVoxels;
  irtkErosion<irtkGreyPixel> erosion;
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
  double point [3];

  //Establish surface voxels
  erosion.SetConnectivity(CONNECTIVITY_26);

  surfaceVoxels = _binaryVolume;

  erosion.SetInput(&surfaceVoxels);
  erosion.SetOutput(&surfaceVoxels);
  erosion.Run();

  x=_binaryVolume.GetX();
  y=_binaryVolume.GetY();
  z=_binaryVolume.GetZ();
  for (i=0;i<x;i++){
    for (j=0;j<y;j++){
      for (k=0;k<z;k++){
        if (_binaryVolume(i,j,k) && ~surfaceVoxels(i,j,k)){
          points->InsertNextPoint(i,j,k);
        }
      }
    }
  }

  polydata->SetPoints(points);

  return polydata;
}

vtkSmartPointer<vtkImageData> irtkVolumeToMesh::Mesh2VtkMask(vtkSmartPointer<vtkPolyData> polydata){
  vtkSmartPointer<vtkTriangleFilter> triFilter = vtkTriangleFilter::New();
  irtkGreyImage binaryVolume;
  irtkGreyPixel *irtkImagePtr;
  double resampleFactor,resampleFactorX,resampleFactorY,resampleFactorZ;

  int x,y,z,numVoxels;
  int vtkX,vtkY,vtkZ;
  int i, j, k;
  SetVTKInput(triFilter, polydata);
  triFilter->Update();
  polydata = triFilter->GetOutput();

  x=_binaryVolume.GetX();
  y=_binaryVolume.GetY();
  z=_binaryVolume.GetZ();


  vtkSmartPointer<vtkImageData> vtkImage
    = vtkSmartPointer<vtkImageData>::New();

  resampleFactor=_averageEdgeLengthMM/_resolution;
  vtkX=ceil((int)(x/resampleFactor));
  vtkY=ceil((int)(y/resampleFactor));
  vtkZ=ceil((int)(z/resampleFactor));
  numVoxels=vtkX*vtkY*vtkZ;
  vtkImage->SetDimensions(vtkX,vtkY,vtkZ);

  resampleFactorX=(float)x/(float)vtkX;
  resampleFactorY=(float)y/(float)vtkY;
  resampleFactorZ=(float)z/(float)vtkZ;

  vtkImage->SetSpacing(resampleFactorX,resampleFactorY,resampleFactorZ);
  vtkImage->SetOrigin(0,0,0);
#if VTK_MAJOR_VERSION >= 6
  vtkImage->AllocateScalars(VTK_SHORT, 1);
#else
  vtkImage->SetScalarTypeToShort();
  vtkImage->SetNumberOfScalarComponents(1);
  vtkImage->AllocateScalars();
#endif

  short *vtkImagePtr= (short *)(vtkImage->GetScalarPointer());

  for(i = 0; i < numVoxels; i++) {
    *vtkImagePtr=0;
    vtkImagePtr++;
  }

  vtkSmartPointer<vtkPolyDataToImageStencil> dataToStencil
    = vtkPolyDataToImageStencil::New();
  dataToStencil->SetTolerance(0);
  SetVTKInput(dataToStencil, polydata);
  dataToStencil->SetInformationInput(vtkImage);

  vtkSmartPointer<vtkImageStencil> stencil = vtkImageStencil::New();
  SetVTKInput(stencil, vtkImage);
#if VTK_MAJOR_VERSION >= 6
  stencil->SetStencilData(dataToStencil->GetOutput());
#else
  stencil->SetStencil(dataToStencil->GetOutput());
#endif
  stencil->SetBackgroundValue(1);
  stencil->Update();
  vtkImage = stencil->GetOutput();

  return vtkImage;
}



vtkSmartPointer<vtkPolyData> irtkVolumeToMesh::GetHull(){



  int levels = 3;
  vtkSmartPointer<vtkHull> hullFilter;
  vtkSmartPointer<vtkDelaunay3D> delaunay;
  vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter;
  vtkSmartPointer<vtkPolyData> surfacePoints;
  vtkSmartPointer<vtkPolyData> hull;
  hullFilter= vtkSmartPointer<vtkHull>::New();
  delaunay = vtkSmartPointer<vtkDelaunay3D>::New();
  surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();

  surfacePoints=BinaryVolume2VTKSurfacePoints();
  SetVTKInput(hullFilter, surfacePoints);
  hullFilter->AddRecursiveSpherePlanes(levels);
  delaunay->SetInputConnection(hullFilter->GetOutputPort());
  surfaceFilter->SetInputConnection(delaunay->GetOutputPort());
  surfaceFilter->Update();
  hull=surfaceFilter->GetOutput();

  return hull;
}


void irtkVolumeToMesh::SetInitialMesh (vtkSmartPointer<vtkPolyData> polydata){
  _mesh=polydata;
  _initialMesh=input;
}

void irtkVolumeToMesh::SetInitialMeshToSphere(){
  _initialMesh=sphere;
}
void irtkVolumeToMesh::SetInitialMeshToHull(){
  _initialMesh=hull;
}


void irtkVolumeToMesh::SetInput (irtkRealImage *in){
  //TODO resize non-isotropic volumes
  irtkMatrix i2w;
  cout << _boundaryThreshold << endl;
  cout << "Setting Input" << endl;
  _binaryVolume=*in;

  int x,y,z;
  for (z = 0; z < in->GetZ(); z++) {
    for (y = 0; y < in->GetY(); y++) {
      for (x = 0; x < in->GetX(); x++) {
        _binaryVolume(x, y, z) = in->Get(x,y,z)>_boundaryThreshold;
      }
    }
  }
  if(_closeBinaryVolume){
    cout << "Closing mask...";
    cout.flush();
    ///Close mask
    irtkErosion<irtkGreyPixel> erosion;
    erosion.SetConnectivity(CONNECTIVITY_26);
    irtkDilation<irtkGreyPixel> dilation;
    dilation.SetConnectivity(CONNECTIVITY_26);

    dilation.SetInput(&_binaryVolume);
    dilation.SetOutput(&_binaryVolume);
    dilation.Run();

    erosion.SetInput(&_binaryVolume);
    erosion.SetOutput(&_binaryVolume);
    erosion.Run();
    cout << "done" << endl;
  }
  _resolution=_binaryVolume.GetXSize();
  i2w = _binaryVolume.GetImageToWorldMatrix();
  _flipNormals=i2w(0,0)*i2w(1,1)*i2w(2,2)>0;
  CalculateCentreOfMass();
}

bool irtkVolumeToMesh::Intersect(double vector [], double point [] ){
  //check not outside of volume!

  double normVector [3];
  int i,x,y,z;
  double nf;
  bool searchBackwards,noIntersection, outOfBounds,outOfBounds_;

  int pointRounded [3];

  x=_binaryVolume.GetX();
  y=_binaryVolume.GetY();
  z=_binaryVolume.GetZ();

  for (i=0;i<3;i++)  pointRounded[i]=floor(point[i]+0.5);

  //if starting point is inside volume then return true
  //if ( pointRounded[0]>=0 && pointRounded[0]<x && pointRounded[1]>=0 && pointRounded[1]<y && pointRounded[2]>=0 && pointRounded[2]<z ){
  //  inBounds=false;
    if (_binaryVolume(pointRounded[0],pointRounded[1],pointRounded[2])) //{
      return true;





  int pointRounded_ [3];
  double vectorMag [3],searchPoint [3],searchPoint_ [3];




  for (i=0;i<3;i++)  vectorMag[i]=abs(vector[i]);

  if (vectorMag[0] > vectorMag[1] ){
    if (vectorMag[0] > vectorMag[2] )  nf=vectorMag[0];
    else nf=vectorMag[2];
  }
  else if(vectorMag[1] > vectorMag[2])
    nf=vectorMag[1];
  else
    nf=vectorMag[2];

  for (i=0;i<3;i++)  normVector[i]=vector[i]/nf;


  //step along normal direction
  for (i=0;i<3;i++){
    searchPoint[i]=point[i]+normVector[i];
    searchPoint_[i]=point[i]-normVector[i];
  }






  noIntersection=true;
  outOfBounds=false;
  outOfBounds_=false;

  //while no intersection and outside

  while (noIntersection){


    for (i=0;i<3;i++)  pointRounded[i]=floor(searchPoint[i]+0.5);
    if (pointRounded[0]>=x || pointRounded[0]<0 ||
        pointRounded[1]>=y || pointRounded[1]<0 ||
        pointRounded[2]>=z || pointRounded[2]<0){
      break; //out of bounds
    }
    else{
      if(_binaryVolume(pointRounded[0],pointRounded[1],pointRounded[2])){ //intersection detected!
        noIntersection=false;
        //cout << "binVolume" << _binaryVolume(0,0,0) << endl;
        //cout << "intersection" << endl;
      }
      else{
         //step forwards along normal  direction
        //cout << "step" << endl;
        for (i=0;i<3;i++) searchPoint[i]=searchPoint[i]+normVector[i];
      }
    }


    if (!outOfBounds_){// check backwards direction
      for (i=0;i<3;i++)  pointRounded_[i]=floor(searchPoint_[i]+0.5);
      if (pointRounded_[0]>=x || pointRounded_[0]<0 ||
          pointRounded_[1]>=y || pointRounded_[1]<0 ||
          pointRounded_[2]>=z || pointRounded_[2]<0){
        outOfBounds_=true;
      }
      else{
        if(_binaryVolume(pointRounded_[0],pointRounded_[1],pointRounded_[2])){ //intersection detected!
          break;
          //cout << "intersection" << endl;
        }
        else{
           //step backwards along normal  direction
          //cout << "step" << endl;
          for (i=0;i<3;i++)  searchPoint_[i]=searchPoint_[i]-normVector[i];
        }
      }
    }
  }


  return !noIntersection;

}




double irtkVolumeToMesh::GetVector(double p1 [],double p2 [], double vector []){
  double nf=0;
  int j;
  for(j=0;j<3;j++){
    vector[j]=p2[j]-p1[j]; //head towards com
    nf+=pow(vector[j],2);
  }
  nf=sqrt(nf);
  for(j=0;j<3;j++) vector[j]=vector[j]/nf;
  return nf;
}



void printVector(double v []){
  for (int i = 0; i < 3; ++i){
    cout << v[i];
    if (i != 2) { cout << ","; }
  }
  cout << endl;
}


void copy3DVector(double*v, double*v2 ){
  for (int i = 0; i < 3; i++){
    v2[i]=v[i];
  }
}

void invert3DVector(double*v, double*v2 ){
  for (int i = 0; i < 3; ++i){
    v2[i]=-v[i];
  }
}

void normalize3DVector(double*v, double*v2 ){
  double nf = .0;
  for (int i = 0; i < 3; ++i){
    nf+=pow(v[i],2);
  }
  nf=sqrt(nf);
  for (int i = 0; i < 3; ++i){
    v2[i]=v2[i]/nf;
  }

}



void irtkVolumeToMesh::ComputeMeshNormals(){
  vtkSmartPointer<vtkPolyDataNormals> normalFilter = vtkSmartPointer<vtkPolyDataNormals>::New();
  SetVTKInput(normalFilter, _mesh);
  normalFilter->SplittingOff();
  normalFilter->AutoOrientNormalsOff() ;
  normalFilter->ConsistencyOn();
  normalFilter->Update();
  _mesh = normalFilter->GetOutput();
}


void irtkVolumeToMesh::SmoothMesh(){
  irtkPolyDataSmoothing smoothingFilter;
  smoothingFilter.SetLambda(_smoothing);
  smoothingFilter.SetSigma(0);
  smoothingFilter.SetInput(_mesh);
  smoothingFilter.Update();
  _mesh = NULL;
  _mesh = vtkSmartPointer<vtkPolyData>::New();
  _mesh->DeepCopy(smoothingFilter.GetOutput());
   //cout << "smoothing mesh: " << _smoothing << endl;
}


void irtkVolumeToMesh::TransformPointsToWorldCoordinates(){
  int ptIdx,numPoints;
  vtkSmartPointer<vtkPoints> points;
  vtkSmartPointer<vtkReverseSense> reverseCells;
  double vertex [3];

  points=_mesh->GetPoints();
  numPoints=_mesh->GetNumberOfPoints();

  for(ptIdx=0;ptIdx<numPoints;ptIdx++){
    points->GetPoint(ptIdx,vertex);
    _binaryVolume.ImageToWorld(vertex[0],vertex[1],vertex[2]);
    points->SetPoint(ptIdx,vertex);
  }
  _mesh->SetPoints(points);


  if (_flipNormals) {
    reverseCells=  vtkSmartPointer<vtkReverseSense>::New();
    reverseCells->ReverseCellsOn();
    reverseCells->ReverseNormalsOff();
    SetVTKInput(reverseCells, _mesh);
    reverseCells->Update();
    _mesh=reverseCells->GetOutput();
  }
  ComputeMeshNormals();

}

void irtkVolumeToMesh::TransformPointsToImageCoordinates(){
  int ptIdx,numPoints;
  vtkSmartPointer<vtkPoints> points;
  vtkSmartPointer<vtkReverseSense> reverseCells;
  double vertex [3];

  points=_mesh->GetPoints();
  numPoints=_mesh->GetNumberOfPoints();

  for(ptIdx=0;ptIdx<numPoints;ptIdx++){
    points->GetPoint(ptIdx,vertex);
    _binaryVolume.WorldToImage(vertex[0],vertex[1],vertex[2]);
    points->SetPoint(ptIdx,vertex);
  }
  _mesh->SetPoints(points);

  if (_flipNormals) {
    reverseCells=  vtkSmartPointer<vtkReverseSense>::New();
    reverseCells->ReverseCellsOn();
    reverseCells->ReverseNormalsOff();
    SetVTKInput(reverseCells, _mesh);
    reverseCells->Update();
    _mesh=reverseCells->GetOutput();
  }
  ComputeMeshNormals();
}

//freeze

bool irtkVolumeToMesh::SelfIntersect(double vertex [3], double vector [3], double * stepsize){

  int i, subId=-1;
  double point0 [3],point1 [3], x [3],pcoords[3];
  double t,tol;
  vtkIdType cellId;
  for (i=0; i<3;i++){
     point1[i]=vertex[i]-_selfIntersectionVoxels*vector[i];
     point0[i]=vertex[i]-0.0001*vector[i];
     pcoords[i]=0;
     x[i]=0;
   }
   //outward intersection, possibly twistes do nothing
   t=-1;
   tol=0.00001;
   _selfIntersectBSPTree->IntersectWithLine(point0,point1,tol,t, x, pcoords, subId, cellId);

   if (t>0){
     *stepsize=0;
     return true;
   }

   //inward intersection
   for (i=0; i<3;i++){
     point1[i]=vertex[i]+_selfIntersectionVoxels*vector[i];
     point0[i]=vertex[i]+0.0001*vector[i];
   }
   t=-1;
   _selfIntersectBSPTree->IntersectWithLine(point0,point1,tol,t, x, pcoords, subId, cellId);
   if( t>0){
     *stepsize=0;
     return true;
   }
   return false;
}


//constant repulsion
/*
bool irtkVolumeToMesh::SelfIntersect(double vertex [3], double vector [3], double * stepsize){
  int i, subId=-1;;
  double point [3],point0 [3],point1 [3], x [3],pcoords[3];
  double squareDistance,t,tol;
  vtkIdType cellId;
  for (i=0; i<3;i++){
     point1[i]=vertex[i]-_selfIntersectionVoxels*vector[i];
     point0[i]=vertex[i]-0.0001*vector[i];
   }
   //outward intersection, possibly twistes do nothing
   t=-1;
   tol=0.00001;
   _selfIntersectBSPTree->IntersectWithLine(point0,point1,tol,t, x, pcoords, subId, cellId);

   if (t>0){
     *stepsize=0;
     return true;
   }

   //inward intersection
   for (i=0; i<3;i++){
     point1[i]=vertex[i]+_selfIntersectionVoxels*vector[i];
     point0[i]=vertex[i]+0.0001*vector[i];
   }
   t=-1;
   _selfIntersectBSPTree->IntersectWithLine(point0,point1,tol,t, x, pcoords, subId, cellId);
   if( t>0){
     if (*stepsize>-0.25*_maxStepSizeVoxels) {
       *stepsize=-0.25*_maxStepSize;
     }
     return true;
   }
   else{
     return false;
   }
}
*/
//Progressoin repulsion
/*
bool irtkVolumeToMesh::SelfIntersect(double vertex [3], double vector [3], double * stepsize){
  int i, subId=-1;;
  double point [3],point0 [3],point1 [3], x [3],pcoords[3];
  double distance,projectionLength,t,tol,maxstep;
  vtkIdType cellId;

  projectionLength=_selfIntersectionVoxels;
  for (i=0; i<3;i++){
     point1[i]=vertex[i]-projectionLength*vector[i];
     point0[i]=vertex[i]-0.0001*vector[i];
   }
   //outward intersection, possibly gone to far already self-intersecting!
   t=-1;
   tol=0.00001;
   _selfIntersectBSPTree->IntersectWithLine(point0,point1,tol,t, x, pcoords, subId, cellId);

   if (t>0){
     *stepsize=0;
     return true;
   }

   //inward self intersection
   projectionLength=_selfIntersectionVoxels+_maxStepSizeVoxels;
   for (i=0; i<3;i++){
     point1[i]=vertex[i]+projectionLength*vector[i];
     point0[i]=vertex[i]+0.0001*vector[i];
   }
   t=-1;
   _selfIntersectBSPTree->IntersectWithLine(point0,point1,tol,t, x, pcoords, subId, cellId);
   if( t>0){//reset gap to _selfIntersectionVoxels by moving half of
     distance=projectionLength*t;

     if(distance<=_selfIntersectionVoxels){//too close :/
       *stepsize=-(_selfIntersectionVoxels-distance)/2; //step halfway to safety
     }
     else{ //in the danger zone, only allow half the step to the boundary
       maxstep=(distance-_selfIntersectionVoxels)/2;
       if (*stepsize>maxstep) {
         *stepsize=maxstep;
       }
     }
     return true;
   }

   return false;
}
*/

//Progressoin repulsion 2
/*
bool irtkVolumeToMesh::SelfIntersect(double vertex [3], double vector [3], double * stepsize){
  int i, subId=-1;;
  double point [3],point0 [3],point1 [3], x [3],pcoords[3];
  double distance,projectionLength,t,tol,maxstep;
  vtkIdType cellId;

  projectionLength=_selfIntersectionVoxels;
  for (i=0; i<3;i++){
     point1[i]=vertex[i]-projectionLength*vector[i];
     point0[i]=vertex[i]-0.0001*vector[i];
   }
   //outward intersection, possibly gone to far already self-intersecting!, do nothing
   t=-1;
   tol=0.00001;
   _selfIntersectBSPTree->IntersectWithLine(point0,point1,tol,t, x, pcoords, subId, cellId);

   if (t>0){
     *stepsize=0;
     return true;
   }

   //inward self intersection
   projectionLength=_selfIntersectionVoxels+_maxStepSizeVoxels;
   for (i=0; i<3;i++){
     point1[i]=vertex[i]+projectionLength*vector[i];
     point0[i]=vertex[i]+0.0001*vector[i];
   }
   t=-1;
   _selfIntersectBSPTree->IntersectWithLine(point0,point1,tol,t, x, pcoords, subId, cellId);
   if( t>0){//reset gap to _selfIntersectionVoxels by moving half of
     distance=projectionLength*t;

     if(distance<=_selfIntersectionVoxels){//too close :/
       *stepsize=-(_selfIntersectionVoxels-distance)/2; //step to safety
     }
//     else{ //in the danger zone, only allow half the step to the boundary
//       maxstep=(distance-_selfIntersectionVoxels)/2;
//       if ( *stepsize>maxstep) { //*stepsize > _maxStepSizeVoxels &&
//         *stepsize=maxstep;
//       }
//     }

     return true;
   }


   return false;
}
*/


void irtkVolumeToMesh::BuildSelfIntersectionTree(){
  _selfIntersectBSPTree->SetDataSet(_mesh);
  _selfIntersectBSPTree->BuildLocator();
}

bool irtkVolumeToMesh::InsideVolume(double vertex [3]){
  //check if outside volume
  bool inside;
  int x,y,z;
  x=_distanceField.GetX()-1;
  y=_distanceField.GetY()-1;
  z=_distanceField.GetZ()-1;
  inside=(vertex[0]>=0 && vertex[0]<=x &&
          vertex[1]>=0 && vertex[1]<=y &&
          vertex[2]>=0 && vertex[2]<=z);
  return inside;
}

double irtkVolumeToMesh::UpdatePointPositions(){
  bool update;
  int i,ptId;
  int numPoints;
  char ch;
  double normal [3],direction [3],vertex [3];
  double stepSize,cumulativeDistance;

  vtkSmartPointer<vtkDataArray> normals;
  vtkSmartPointer<vtkPoints> points, newPoints;


  normals=_mesh->GetPointData()->GetNormals();
  points=_mesh->GetPoints();
  numPoints=_mesh->GetNumberOfPoints();

  //copy old points to new points
  points=_mesh->GetPoints();
  newPoints=vtkSmartPointer<vtkPoints>::New();
  newPoints->DeepCopy(points);

  cumulativeDistance=0;


  if (_selfIntersectionOn) { BuildSelfIntersectionTree(); }
  for(ptId=0;ptId<numPoints;ptId++){

    points->GetPoint(ptId,vertex); //Get vertex position
    normals->GetTuple(ptId,normal); //Get normal vector


    if (!InsideVolume(vertex)){
      stepSize=_maxStepSizeVoxels;
      cumulativeDistance+=GetVector(vertex,_centreOfMass,direction);
    }
    else{
      stepSize=_interpolator->Evaluate(vertex[0], vertex[1], vertex[2]);
      cumulativeDistance+=abs(stepSize);
      if (stepSize > _maxStepSizeVoxels) stepSize = _maxStepSizeVoxels;
      if (stepSize < -_maxStepSizeVoxels) stepSize = -_maxStepSizeVoxels;
      copy3DVector(normal,direction);
      if (_intersectionOn && ! Intersect(normal,vertex)==1){
        stepSize=0;
      }

    }

    if(_selfIntersectionOn && stepSize>=0 ){ //modify stepsize if in the proxmity itself
      SelfIntersect(vertex,normal,&stepSize);
    }

    for (i=0; i<3;i++){ vertex[i]=vertex[i]+direction[i]*stepSize; }
    newPoints->SetPoint(ptId,vertex);

  }

  _mesh->SetPoints(newPoints);
  return cumulativeDistance/(double)numPoints;
}



void irtkVolumeToMesh::LowPassFilter(){
  SetVTKInput(_lowPassFilter, _mesh);
  _lowPassFilter->Update();
  _mesh=_lowPassFilter->GetOutput();
}



void irtkVolumeToMesh::Subdivide(){
  vtkSmartPointer<vtkButterflySubdivisionFilter> subdivisionFilter;
  subdivisionFilter = vtkSmartPointer<vtkButterflySubdivisionFilter>::New();
  SetVTKInput(subdivisionFilter, _mesh);
  subdivisionFilter->Update();
  _mesh=subdivisionFilter->GetOutput();
  _averageEdgeLengthMM=_averageEdgeLengthMM/2;
  _remeshFilter->SetMaxEdgeLength((_averageEdgeLengthMM/_resolution)*2);
  _remeshFilter->SetMinEdgeLength(_averageEdgeLengthMM/_resolution/2);
  _smoothing=_smoothing*2;
}

void irtkVolumeToMesh::Screenshot(string fileName){
  // Visualize
  const char * fileNameChars=fileName.c_str();


    vtkSmartPointer<vtkPolyDataMapper> mapper =
      vtkSmartPointer<vtkPolyDataMapper>::New();
  SetVTKInput(mapper, _mesh);

    vtkSmartPointer<vtkActor> actor =
      vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    vtkSmartPointer<vtkRenderer> renderer =
      vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow =
      vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    renderWindow->SetAlphaBitPlanes(1); //enable usage of alpha channel
    renderWindow->SetOffScreenRendering(1);
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
      vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    renderer->AddActor(actor);
    //renderer->SetBackground(1,1,1); // Background color white


    renderWindow->Render();

    // Screenshot
    vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter =
      vtkSmartPointer<vtkWindowToImageFilter>::New();
    windowToImageFilter->SetInput(renderWindow);
    windowToImageFilter->SetMagnification(3); //set the resolution of the output image (3 times the current resolution of vtk render window)
    windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
    windowToImageFilter->Update();

    vtkSmartPointer<vtkPNGWriter> writer =
      vtkSmartPointer<vtkPNGWriter>::New();
    writer->SetFileName(fileNameChars);
    writer->SetInputConnection(windowToImageFilter->GetOutputPort());
    writer->Write();


}

/*
void irtkVolumeToMesh::Screenshot(string fileName){
  const char * fileNameChars=fileName.c_str();



  vtkSmartPointer<vtkGraphicsFactory> graphics_factory =
    vtkSmartPointer<vtkGraphicsFactory>::New();
  graphics_factory->SetOffScreenOnlyMode( 1);
  //graphics_factory->SetUseMesaClasses( 1 );


    vtkSmartPointer<vtkPolyDataMapper> mapper =
      vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInput(_mesh);

    vtkSmartPointer<vtkActor> actor =
      vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);


    // A renderer and render window
    vtkSmartPointer<vtkRenderer> renderer =
      vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow =
      vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->SetOffScreenRendering(1);
    renderWindow->AddRenderer(renderer);


    // Add the actors to the scene
    renderer->AddActor(actor);
    renderer->SetBackground(1,1,1); // Background color white
    renderWindow->Render();

    vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter =
      vtkSmartPointer<vtkWindowToImageFilter>::New();
    windowToImageFilter->SetInput(renderWindow);
    windowToImageFilter->SetMagnification(3); //set the resolution of the output image (3 times the current resolution of vtk render window)
    windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
    windowToImageFilter->Update();

    vtkSmartPointer<vtkPNGWriter> writer =
      vtkSmartPointer<vtkPNGWriter>::New();
    writer->SetFileName(fileNameChars);
    writer->SetInputConnection(windowToImageFilter->GetOutputPort());
    writer->Write();
}
*/

void irtkVolumeToMesh::SetInputAsTarget(){
  //TODO set binary image to binary image!
}

void irtkVolumeToMesh::DeformMeshToTarget(){
  // no self intersections
  //coarse tree for intersections with mesh
  int i,j,iteration, improvementIterations=10;
  double improvement,oldDistanceError,distanceError;
  double distanceErrors [improvementIterations];

  ostringstream os;
  string screenshotFileName;

  distanceError=0;
  oldDistanceError=DBL_MAX;
  improvement=DBL_MAX;

  iteration=0;
  while( iteration < _minIterations ||  (iteration < _maxIterations && improvement > _improvementThreshold) ){
    iteration++;
    cout << "iteration " << iteration << endl;
    ComputeMeshNormals();

    distanceError=UpdatePointPositions();


    //cout << "point positions updated" << endl;
    if (_smoothingOn) { SmoothMesh(); } //     cout << "smoothing done" << endl; }

    if (_lowPassFilterOn) { LowPassFilter(); } //    cout << "low pass filter done" << endl; }

    Remesh();

    if (_screenShotsOn){
      os << _screenShotDir << iteration << ".png";
      screenshotFileName = os.str();
      os.clear();
      os.str(std::string());
      Screenshot(screenshotFileName);
    }

    //oldDistanceError=distanceErrors[(iteration+1)%improvementIterations];
    distanceErrors[iteration%improvementIterations]=distanceError;
    cout << "distanceError: " << distanceError << endl;


    if (iteration>=improvementIterations){
      improvement=-distanceError;
      for (i=1;i<improvementIterations/2;i++){
        improvement-=distanceErrors[(iteration-i)%improvementIterations];
      }
      for (;i<improvementIterations;i++){
        improvement+=distanceErrors[(iteration-i)%improvementIterations];
      }
      improvement/=improvementIterations;
      //cout << "improvement: " << improvement << endl;
    }

  }

}
void irtkVolumeToMesh::Remesh(){
  irtkRemesher remeshFilter;
  remeshFilter.SetMaxEdgeLength(2*_averageEdgeLengthMM/_resolution);
  remeshFilter.SetMinEdgeLength(0.8*_averageEdgeLengthMM/_resolution);
  vtkSmartPointer<vtkPolyData> remeshed = remeshFilter.Remesh(_mesh);

  _mesh = NULL;
  _mesh = vtkSmartPointer<vtkPolyData>::New();
  _mesh->DeepCopy(remeshed);
  _mesh = remeshed;
}

void irtkVolumeToMesh::FinalRefinement(){
  bool lowPassFilterOn=_lowPassFilterOn;
  bool lowPassFilterIterations=_lowPassFilterIterations;
  bool lowPassFilterBand=_lowPassFilterBand;
  int maxIterations=_maxIterations;
  _maxIterations=_refinementIterations;

  _lowPassFilterOn=true;
  _lowPassFilterIterations=25;
  _lowPassFilterBand=0.5;
  _smoothingOn=false;

  // set laplacian smoothing off

  cout << "-------------" << endl;
  cout << "Refining Mesh" << endl;
  cout << "-------------" << endl;
  DeformMeshToTarget();
  //restore filter settings
  //_smoothingFilter->SetLambda(_smoothing);


  _lowPassFilterOn=lowPassFilterOn;
  _lowPassFilterIterations=lowPassFilterIterations;
  _lowPassFilterBand=lowPassFilterBand;
  _smoothingOn=true;
  _maxIterations=maxIterations;

}


void irtkVolumeToMesh::PrintSettings(){
  cout << "SETTINGS: " << endl;
  cout << "Edge length: " << _averageEdgeLengthMM << endl;
  cout << "Smoothing: " << _smoothing << endl;

  cout << "Maximum iterations: " << _maxIterations << endl;
  cout << "Epsilon: " << _improvementThreshold << endl;

  if (_selfIntersectionOn){
    cout << "Self intersection test: On" << endl;
  }
  else{
    cout << "Self intersection test: Off" << endl;
  }
  if (_lowPassFilterOn){
    cout << "Low pass filter On: " << _lowPassFilterIterations << " iterations, filtering band: " <<
        _lowPassFilterBand << endl;
  }

  if (_finalRefinementOn){
    cout << "Final refinement: On, smoothing: " << _refinementSmoothing << endl;
  }

}




vtkSmartPointer<vtkPolyData>  irtkVolumeToMesh::GetOuput(){
  int iteration,level;
  double stepSize,maxStepSize,distanceError,nf;
  ostringstream os;
  string screenshotFileName;
  Initialize();
  PrintSettings();
  cout << "------------------------" << endl;
  cout << "Deforming mesh to target" << endl;
  cout << "------------------------" << endl;

  DeformMeshToTarget();
  if (_finalRefinementOn){ FinalRefinement(); }

  TransformPointsToWorldCoordinates();
  return _mesh;
}
