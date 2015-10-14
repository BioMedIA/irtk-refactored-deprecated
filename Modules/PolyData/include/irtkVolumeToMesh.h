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

#ifndef IRTK_VOLUMETOMESH_H_
#define IRTK_VOLUMETOMESH_H_

#include <irtkImage.h>
#include <vtkPolyData.h>
#include <irtkEuclideanDistanceTransform.h>
#include <vtkPolyDataNormals.h>
#include <irtkRemesher.h>
#include <irtkPolyDataSmoothing.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <irtkImageFunction.h>
#include <vtkModifiedBSPTree.h>
#include <vtkOBBTree.h>

using namespace irtk::polydata;


class irtkVolumeToMesh{

  enum InitialMesh {input, hull, sphere};

private:
  bool _multiResolutionOn,_screenShotsOn,_selfIntersectionOn,_intersectionOn,
  _meshInitialized, _lowPassFilterOn,_finalRefinementOn,_smoothingOn,
  _flipNormals,_closeBinaryVolume,_gaussianInterpolationOn;
  char * _screenShotDir;
  double _smoothing,_averageEdgeLengthMM,_boundaryThreshold,_resolution,_maxStepSize,
  _maxStepSizeMM,_maxStepSizeVoxels,_selfIntersectionVoxels,_selfIntersectionMM,
  _boundingRadius,_improvementThreshold,_lowPassFilterBand,_refinementSmoothing;
  double _centreOfMass [3];
  int _maxIterations,_minIterations,_levels,_selfIntersectionTreeLeavesPerNode,_minNumPoints,_subdivisions,
  _lowPassFilterIterations,_refinementIterations;
  int _boundingBoxP0 [3], _boundingBoxP1 [3];
  // vtkPolyData mesh;
  irtkGreyImage _binaryVolume;
  irtkRealImage _distanceField;
  vtkSmartPointer<vtkPolyData> _mesh;
  irtkEuclideanDistanceTransform<irtkRealPixel> *_edt;
  vtkSmartPointer<vtkPolyDataNormals> _normalFilter;
  irtkRemesher* _remeshFilter;
  irtkPolyDataSmoothing* _smoothingFilter;
  vtkSmartPointer<vtkWindowedSincPolyDataFilter> _lowPassFilter;
  irtkImageFunction *_interpolator;

  vtkSmartPointer<vtkModifiedBSPTree> _selfIntersectBSPTree;

  InitialMesh _initialMesh;

  void Initialize();
  void InitializeParameters();
  void InitializeMesh();
  void InitializeMeshAsBoundingSphere();
  void InitializeMeshAsConvexHull();
  void InitializeDistanceField();
  void InitializeFilters();

  bool Intersect(double [], double []);
  double GetVector(double [],double [],double []);
  void GetAdjacency(vtkSmartPointer<vtkCellArray>, vtkSmartPointer<vtkIdList> *);
  void ComputeMeshNormals();
  void SmoothMesh();
  void TransformPointsToWorldCoordinates();
  void TransformPointsToImageCoordinates();
  double UpdatePointPositions();
  void Subdivide();
  void Screenshot(string);
  void GetAdjs(vtkSmartPointer<vtkIdList> adj []);
  void InterpolateVelocities(double**, bool  [],  vtkSmartPointer<vtkIdList>  []);
  void AdaptivelySmoothMesh();
  void Remesh();
  void LowPassFilter();
  vtkSmartPointer<vtkPolyData> GetHull();
  void DeformMeshToTarget();
  bool SelfIntersect(double [3], double [3], double *);
  double Repulsion(double [3], double [3]);
  void SetInputAsTarget();
  void BuildSelfIntersectionTree();
  void CalculateCentreOfMass();
  bool InsideVolume(double [3]);
  void FinalRefinement();
  void PrintSettings();

  vtkSmartPointer<vtkPolyData> BinaryVolume2VTKSurfacePoints();
  vtkSmartPointer<vtkImageData> Mesh2VtkMask(vtkSmartPointer<vtkPolyData>);
public:

// Constructor
  irtkVolumeToMesh();

/// Deconstuctor
  ~irtkVolumeToMesh();




/// Set input image for filter
  void SetInput (irtkRealImage *);
/// Set output image for filter
 void SetMaxIterations(int);
 void SetMinIterations(int);
 void SetEdgeLength(double);
 void SetNumberOfLevels(int);


 void SetScreenShotDirectory(char *);
 void SetInitialMeshToSphere();
 void SetInitialMeshToHull();
 void SetInitialMesh(vtkSmartPointer<vtkPolyData>);

 void SetEdgeSmoothing(double);
 void SetImprovementThreshold(double);

 void ScreenShotsOn();
 void ScreenShotsOff();
 void SelfIntersectionOn();
 void SelfIntersectionOff();
 void LowPassFilterOn();
 void LowPassFilterOff();
 void FinalRefinementOn();
 void FinalRefinementOff();
 void SmoothingOn();
 void SmoothingOff();
 void GaussianInterpolationOn();
 void GaussianInterpolationOff();


 void SetLowPassFilterBand(double);
 void SetLowPassFilterIterations(int);


 void SetSmoothingValue(double);
 void SetRefinementIterations(int);
 void SetBoundaryThreshold(double);
 void CloseBinaryVolumeOn();
 void CloseBinaryVolumeOff();

 vtkSmartPointer<vtkPolyData>  GetOuput();


};



#endif /* VOLUMETOMESH_H_ */
