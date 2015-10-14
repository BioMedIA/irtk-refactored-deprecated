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

#include <irtkParticleSampleSphere.h>
#include <vtkMath.h>
#include <vtkTriangle.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkPlatonicSolidSource.h>
#include <vtkPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyDataWriter.h>
#include <vtkSphereSource.h>
#include <irtkImage.h>
#include <irtkMatrix.h>
#include <vtkDelaunay3D.h>
#include <vtkSmartPointer.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkGeometryFilter.h>
#include <vtkPointSet.h>
#include <vtkLinearSubdivisionFilter.h>
#include <irtkSphereQuadrisect.h>
#include <math.h>
#include <vtkMaskPoints.h>


irtkParticleSampleSphere::irtkParticleSampleSphere(){
  _maxIterations=1000;
  _energyTolerance=1E-5;
  _minStepSize=1E-14;
  _maxStepSize=0.1;
  _numPoints=100;
  _subdivisions=1;
  _sphericalCoordsOn=false;
}

irtkParticleSampleSphere::~irtkParticleSampleSphere(){
}


void irtkParticleSampleSphere::Initialize(){
}


void irtkParticleSampleSphere::SphericalCoordsOn(){
  _sphericalCoordsOn=true;
}

void irtkParticleSampleSphere::SphericalCoordsOff(){
  _sphericalCoordsOn=false;
}


void irtkParticleSampleSphere::SetNumberOfPoints(int numPoints){
  _numPoints=numPoints;
}

void irtkParticleSampleSphere::SetNumberOfSubdivisions(int subdivisions){
  _subdivisions=subdivisions;
}


struct MyComparator
{
    const vector<double> & value_vector;

    MyComparator(const vector<double> & val_vec):
        value_vector(val_vec) {}

    bool operator()(int i1, int i2)
    {
        return value_vector[i1] > value_vector[i2];
    }
};

void irtkParticleSampleSphere::AddSphericalCoordScalars(vtkSmartPointer<vtkPolyData> mesh){
  //TODO fix this mess!
  int i;
  int numPoints;
  double phiAngle, thetaAngle;
  double xAngle, yAngle, zAngle;
  double point [3];
  double M_PI2=2*M_PI;
  double M_PI_by_4=M_PI/4;

  vtkSmartPointer<vtkFloatArray> phiScalars,thetaScalars;
  vtkSmartPointer<vtkFloatArray> xScalars,yScalars,zScalars;
  vtkSmartPointer<vtkUnsignedCharArray> colorScalars;
  vtkSmartPointer<vtkPoints> points;

  mesh->Modified();
  points = mesh->GetPoints();
  numPoints= mesh->GetNumberOfPoints();

  //Spherical coordinates
  //_________________________
  //Use rotation from each positive axis pole as coordinate system

  xScalars = vtkSmartPointer<vtkFloatArray>::New();
  xScalars->SetNumberOfComponents(1);
  xScalars->SetNumberOfTuples(numPoints);

  yScalars = vtkSmartPointer<vtkFloatArray>::New();
  yScalars->SetNumberOfComponents(1);
  yScalars->SetNumberOfTuples(numPoints);

  zScalars = vtkSmartPointer<vtkFloatArray>::New();
  zScalars->SetNumberOfComponents(1);
  zScalars->SetNumberOfTuples(numPoints);

  colorScalars = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colorScalars->SetNumberOfComponents(3);
  colorScalars->SetNumberOfTuples(numPoints);

  for (i=0; i<numPoints; i++){
    mesh->GetPoint(i,point);
    xAngle=acos(point[0]);
    yAngle=acos(point[1]);
    zAngle=acos(point[2]);

    xScalars->SetTuple1(i,xAngle);
    yScalars->SetTuple1(i,yAngle);
    zScalars->SetTuple1(i,zAngle);

    xAngle=255*xAngle/M_PI;
    yAngle=255*yAngle/M_PI;
    zAngle=255*zAngle/M_PI;

    colorScalars->SetTuple3(i,xAngle,yAngle,zAngle);
  }

  xScalars->Modified();
  yScalars->Modified();
  zScalars->Modified();
  colorScalars->Modified();


  xScalars->SetName("xAngle");
  yScalars->SetName("yAngle");
  zScalars->SetName("zAngle");
  colorScalars->SetName("colors");

  mesh->GetPointData()->AddArray(xScalars);
  mesh->GetPointData()->AddArray(yScalars);
  mesh->GetPointData()->AddArray(zScalars);
  mesh->GetPointData()->AddArray(colorScalars);
}


void irtkParticleSampleSphere::OptimizeParticlePositions(vtkSmartPointer<vtkPoints> points){

  // Compute geodesic distances between the points
  int i,j,k,numPoints,s,iteration;

  double geoDist,totalEnergy,energy,energyNew,energyOld,totalEnergyOld,dEnergy,dot;
  double djx,djy,djz,dj,dx,dy,dz,nx,ny,nz,tx,ty,tz,xNew,yNew,zNew;
  double normalize;
  double * stepSize;
  bool updating,adaptiveStepSize;
  s=1;
  adaptiveStepSize=true;
  double p1 [3];
  double p2 [3];
  double previousEnergies [4];

  numPoints=points->GetNumberOfPoints();
  irtkMatrix geodesicMatrix(numPoints,numPoints);
  irtkMatrix dotMatrix(numPoints,numPoints);
  irtkMatrix energyMatrix(numPoints,numPoints);
  totalEnergy=0;

  stepSize = new double [numPoints];

  vector<double> energyI;
  vector<double> sortIdx;
  energyI.resize(numPoints);
  sortIdx.resize(numPoints);

  for (i=0;i<numPoints;i++){
	    energyI[i]=0;
    points->GetPoint(i,p1);

    for (j=0;j<numPoints;j++){

      if (i != j){
        points->GetPoint(j,p2);
        //cout << p1[0] <<" " << p1[1] <<" " << p1[2] << endl;
        //cout << p2[0] <<" " << p2[1] <<" " << p2[2] << endl;
        dot=p1[0]*p2[0]+p1[1]*p2[1]+p1[2]*p2[2];


        if(dot<-1) dot=-1;
        if(dot>1)  dot=1;
        dotMatrix.Put(i,j,dot);
        geoDist=acos(dot);


        geodesicMatrix.Put(i,j,geoDist);
        energy=1/(pow(geoDist,s)+DBL_EPSILON);
        energyMatrix.Put(i,j,energy);
        energyI[i]+=energy;
      }
      else
        geodesicMatrix.Put(i,j,DBL_MAX);
    }
   // cout << energyI[i] << endl;
    totalEnergy+=energyI[i];


  }

  //for (i=0;i<numPoints;i++){
   // for (j=0;j<numPoints;j++){
     //cout << geodesicMatrix.Get(i,j) << ",";
   // }
   // cout << endl;
  //}





  // Iteratively optimize the position of the particles along the negative
  // gradient of the energy functional using an adaptive Gauss Seidel
  // update scheme

  //initialize step_size
  for(i=0;i<numPoints;i++)  stepSize[i]=_maxStepSize; //maxStepSize;



  dEnergy=DBL_MAX;
  iteration=0;

  while ( iteration<_maxIterations && dEnergy>_energyTolerance){
	  //iteration=maxIterations;
    // Sort the particles according to their energyI contribution
    for(i=0;i<numPoints;i++)  sortIdx[i]=i;

    sort(sortIdx.begin(), sortIdx.end(), MyComparator(energyI));

    // Update the position of individual particles

    for(k=0;k<numPoints;k++) {
      i=sortIdx[k];
      //cout << energy[i] << endl;

	  //compute geodesic distances to other particles
	  // idx_j=true(N,1); % particle indices, except the current one
	  // idx_j(j)=false;
      // DOTj=DOT(idx_j,j); %just dotproduct to other particles
      // GDj=GD(idx_j,j); %just distances to other particles

	  // geodesicMatrix.Get(i,j);


      //calculate gradient for the ith particle

	  //dVj=bsxfun(@times,s./(sqrt(1-DOTj.^2)+eps),V(idx_j,:));
	  //dVj=bsxfun(@rdivide,dVj,(GDj.^(s+1))+eps);

	  points->GetPoint(i,p1);
	  for(j=0;j<numPoints;j++){
	    if (i != j){
          points->GetPoint(j,p2);
          dj=s/(sqrt(1-pow(dotMatrix(i,j),2))+DBL_EPSILON);
          dj/=pow(geodesicMatrix(i,j),s+1)+DBL_EPSILON;
          djx=dj*p2[0];
          djy=dj*p2[1];
          djz=dj*p2[2];


		  if (djx*djx+djy*djy+djz*djz<100000000){ //particle isn't too close
			dx+=djx;
			dy+=djy;
			dz+=djz;
		  }
		}
	  }

	   //cout << "gradient:" << dx << "," << dy << "," << dz << endl;


	   //Only retain the tangential component of the gradient

	   dot=dx*p1[0]+dy*p1[1]+dz*p1[2];
	   nx=p1[0]*dot;
	   ny=p1[1]*dot;
	   nz=p1[2]*dot;
     tx=dx-nx;
     ty=dy-ny;
     tz=dz-nz;

     //cout << "tangent:" << tx << "," << ty << "," << tz << endl;

	   // Adaptively update position of the i-th particle


       energyOld=energyI[i];
       updating=true;



       while (updating){


         // Update the position
         xNew=p1[0]-stepSize[i]*tx;  // V(j,:)-a(j)*dVj_t;
         yNew=p1[1]-stepSize[i]*ty;
         zNew=p1[2]-stepSize[i]*tz;
         // Constrain the point to surface of the sphere
         normalize=sqrt(xNew*xNew+yNew*yNew+zNew*zNew);




         xNew/=normalize; //		   Vj_new=Vj_new/norm(Vj_new);
         yNew/=normalize;
         zNew/=normalize;

         //cout << "old pos:" << p1[0] << "," << p1[1] << "," << p1[2] << endl;
         //cout << "new pos:" << xNew << "," << yNew << "," << zNew << endl;

         // Compute new energyI
         //DOTj=sum(bsxfun(@times,V,Vj_new),2);
         //DOTj(DOTj<-1)=-1; DOTj(DOTj>1)=1;
         //GDj=acos(DOTj);
         //GDj(j)=Inf;
         //Ue_ij_j=1./((GDj.^s)+eps);
         energyNew=0;
         for(j=0;j<numPoints;j++){
           if (i != j){
             points->GetPoint(j,p2);
             dot=xNew*p2[0]+yNew*p2[1]+zNew*p2[2];
             if(dot<-1) dot=-1;
             if(dot>1)  dot=1;
             geoDist=acos(dot);
             energyNew+=1/(pow(geoDist,s)+DBL_EPSILON);
           }
         }

         //cout << "energy:" << energyNew << "," << energyOld << "," << endl;

         if (energyNew<=energyOld){ //if energyI in system decreases
           //keep point
           points->SetPoint(i,xNew,yNew,zNew);
           //cout << "update!" << endl;
           //update dot products, geodesics and energyI
           energyI[i]=0;
           for (j=0;j<numPoints;j++){
             if (i != j){

               points->GetPoint(j,p2);
               dot=xNew*p2[0]+yNew*p2[1]+zNew*p2[2];
               if(dot<-1) dot=-1;
               if(dot>1)  dot=1;
               dotMatrix.Put(i,j,dot);
               dotMatrix.Put(j,i,dot);
               geoDist=acos(dot);
               geodesicMatrix.Put(i,j,geoDist);
               geodesicMatrix.Put(j,i,geoDist);
               energy=1/(pow(geoDist,s)+DBL_EPSILON);
               energyI[j]-=energyMatrix.Get(i,j);
               energyMatrix.Put(i,j,energy);
               energyMatrix.Put(j,i,energy);
               energyI[j]+=energy;
               energyI[i]+=energy;
               //cout << "new,old:" << energy << energyMatrix.Get(i,j) <<   endl;
               //cout << "new:" << xNew << "," << yNew << "," << zNew << endl;
               //cout << "old:" << p1[0] << "," << p1[1] << "," << p1[2] << endl;
             }
           }
           stepSize[i]*=2.5;
           if (stepSize[i]>_maxStepSize) stepSize[i]=_maxStepSize;


           updating=false;
           //return;
         }
      else{
        if (stepSize[i]>_minStepSize){
          stepSize[i]/=2.5;
          //cout << "newStep:" << stepSize[i] << endl;
            if (stepSize[i]<_minStepSize) stepSize[i]=_minStepSize; //try smaller step size
          }
        else {
          updating=false;
          //cout << "noUpdate!" << endl;
          //cout << "Stepsize:" << stepSize[i] << endl;
        }
      }
    }

       //cout << "point" << i << endl;
  }


       // Evaluate the total energyI of the system
       totalEnergyOld=totalEnergy;
       totalEnergy=0;
       for(i=0;i<numPoints;i++){
         totalEnergy+=energyI[i];
       }

       previousEnergies[iteration % 4]=totalEnergy;


       // Progress update
       if (iteration==1){
           cout << "Iteration \tEnergy Score\n" << endl;
           cout << iteration << "\t\t" << totalEnergy  << endl;
           }
       else if ( (iteration % 1)==0) cout << iteration << "\t\t" <<totalEnergy  <<  endl;

       // Change in energyI
       if (iteration>=10) dEnergy=previousEnergies[(iteration-3) % 4]-totalEnergy;
       //cout << previousEnergies[(iteration-3) % 4] << endl;
       // Reset the step sizes every 20 iterations ?
       //if (iteration % 20 == 0) {
         for(i=0;i<numPoints;i++) stepSize[i]=_maxStepSize;
       //}

     iteration=iteration+1;

   }

  delete [] stepSize;
}

/*
vtkSmartPointer<vtkPoints> irtkParticleSampleSphere::StratifiedSample(){
  vtkSmartPointer<vtkPoints> surfacePoints
    = vtkSmartPointer<vtkPoints>::New();
  int i,j,numDomains,numDomainsSq,R,idxJ;
  double ds,randomNumber,longitude,latitude;
  double *x,*y,*z;

  // Sample the unfolded right cylinder

      // Partition the [-1,1]x[0,2*pi] domain into ceil(sqrt(N))^2 subdomains
      // and then draw a random sample for each
      cout << numDomains << endl;
      numDomains=ceil(sqrt(_numPoints));
      cout << numDomains << endl;
      numDomainsSq=numDomains*numDomains;
      ds=2/(double)numDomains;

      //get random x and y
      //x = new double [numDomains*numDomains];
      double *u = new double[numDomainsSq]();
      double *v = new double[numDomainsSq]();
      double *change = new double[numDomainsSq]();

      change[0]=-1+ds/2;
      for (i=1;i<numDomains;i++){
        change[i]=change[i-1]+ds;
      }

      for (j=0;j<numDomains;j++){
        idxJ=j*numDomains;
        for (i=0;i<numDomains;i++){
          randomNumber = ((float)rand()/(float)RAND_MAX)-0.5;
          u[i+idxJ]=ds*randomNumber+change[i];
          randomNumber = ((float)rand()/(float)RAND_MAX)-0.5;
          v[i+idxJ]=ds*randomNumber+change[j];
        }
      }

  // Remove excess samples
  //R=numDomainsSq-numPoints;
  //if (R>0){
  //    idx=rand() % numDomainsSq;
  //    u(idx)=[];
  //   v(idx)=[];
  //}


  x = new double[numDomainsSq]();
  y = new double[numDomainsSq]();
  z = new double[numDomainsSq]();

  for (i=0;i<numDomainsSq;i++){
    // Calc longitude and latitude
    longitude=(u[i]+1)*M_PI;
    latitude=acos(z[i]);
    // Convert spherical to rectangular co-ords
    x[i]=cos(longitude)*sin(latitude);
    y[i]=sin(longitude)*sin(latitude);
  }

  surfacePoints->SetNumberOfPoints(numDomainsSq);
  for (i=0;i<numDomainsSq;i++){
    surfacePoints->InsertPoint(i,x[i],y[i],z[i]);
  }


  delete[] x;
  delete[] y;
  delete[] z;


  delete[] u;
  delete[] v;
  delete[] change;

  return surfacePoints;
}
*/


vtkSmartPointer<vtkPoints> irtkParticleSampleSphere::StratifiedSample(){

  vtkSmartPointer<vtkSphereSource> sphereSource =
      vtkSmartPointer<vtkSphereSource>::New();
  sphereSource->SetPhiResolution(100);
  sphereSource->SetThetaResolution(100);
  sphereSource->SetRadius(1);
  sphereSource->Update();
  vtkSmartPointer<vtkPolyData> sphere = sphereSource->GetOutput();

  vtkSmartPointer<vtkMaskPoints> masker =
      vtkSmartPointer<vtkMaskPoints>::New();
  masker->SetOnRatio(1);
  masker->SetMaximumNumberOfPoints(_numPoints);
  masker->RandomModeOn();
  masker->SetRandomMode(2); //spatially stratified
  SetVTKInput(masker, sphere);
  masker->Update();

  return masker->GetOutput()->GetPoints();
}


vtkSmartPointer<vtkPolyData> irtkParticleSampleSphere::GetSphere(){
  vtkSmartPointer<vtkPoints> points;
  vtkSmartPointer<vtkPolyData> polydata;
  vtkSmartPointer<vtkDelaunay3D> delaunay3D;
  vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter;
  vtkSmartPointer<vtkLinearSubdivisionFilter> subdivisionFilter;
  irtkSphereQuadrisect sq;

  points=StratifiedSample();
  OptimizeParticlePositions(points);


  polydata = vtkSmartPointer<vtkPolyData>::New();
  delaunay3D = vtkSmartPointer<vtkDelaunay3D>::New();
  surfaceFilter=vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
  subdivisionFilter = vtkSmartPointer<vtkLinearSubdivisionFilter>::New();

  polydata->SetPoints(points);
  SetVTKInput(delaunay3D, polydata);
  surfaceFilter->SetInputConnection(delaunay3D->GetOutputPort());
  surfaceFilter->Update();
  polydata=surfaceFilter->GetOutput();

  if (_subdivisions>0)
    polydata=sq.Quadrisect(polydata,_subdivisions);



  if (_sphericalCoordsOn){
    AddSphericalCoordScalars(polydata);
  }
  return polydata;
}

