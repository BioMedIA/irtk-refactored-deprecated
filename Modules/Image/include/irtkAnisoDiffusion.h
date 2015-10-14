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

#ifndef _ANISODIFFUSION_H
#define _ANISODIFFUSION_H

#include <irtkImageToImage.h>

/**
 * Class for anisotopic diffusion filtering
 */

template <class VoxelType>
class anisoDiffusion : public irtkImageToImage<VoxelType>
{
  irtkImageFilterMacro(anisoDiffusion);

public:

  /// Constructor
  anisoDiffusion();
  
  /// Destructor
  ~anisoDiffusion();

  /// Run anisotropic diffusion filtering
  virtual void Run();
  
  virtual void Run_4D_semiImplicit();
  virtual void Run_3D_semiImplicit();
  virtual void Run_4D_Explicit();
  virtual void Run_3D_Explicit();
  
  /// Parameters
  float ax; //caracteristic gray level variations considered significant between to neighbors voxels along direction x
  float ay; //caracteristic gray level variations considered significant between to neighbors voxels along direction y
  float az; //caracteristic gray level variations considered significant between to neighbors voxels along direction z
  float at; //caracteristic gray level variations considered significant between to neighbors voxels along time
  float dx; //distance between to neighbors voxels along direction x
  float dy; //distance between to neighbors voxels along direction y
  float dz; //distance between to neighbors voxels along direction z
  float dt; //distance between to neighbors voxels along time
  float dTau;         //virtual time step (for the diffusion pde)
  int ITERATIONS_NB; //number of virtual time iterations 
  bool TimeDependent;  //1 -> anisotrop filtering along the time / 0 -> otherwise
  bool SemiImplicit; //1 -> The scheme will be semi implicit (ADI) / 0 -> explicit scheme
};


void TridiagonalSolveFloat(const float *, const float *, float *, float *, float *, int);

#endif
