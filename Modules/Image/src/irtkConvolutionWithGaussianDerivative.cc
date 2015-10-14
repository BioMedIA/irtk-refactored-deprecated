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

#include <irtkImage.h>

#include <irtkConvolution.h>
#include <irtkScalarFunctionToImage.h>
#include <irtkScalarGaussianDx.h>
#include <irtkScalarGaussianDy.h>
#include <irtkScalarGaussianDz.h>
#include <irtkConvolutionWithGaussianDerivative.h>

template <class VoxelType> irtkConvolutionWithGaussianDerivative<VoxelType>::irtkConvolutionWithGaussianDerivative(double Sigma)
{
  _Sigma = Sigma;
}

template <class VoxelType> irtkConvolutionWithGaussianDerivative<VoxelType>::~irtkConvolutionWithGaussianDerivative()
{}

template <class VoxelType> void irtkConvolutionWithGaussianDerivative<VoxelType>::Ix()
{
  double xsize, ysize, zsize;
  
  // Do the initial set up
  this->Initialize();

  // Get voxel dimensions
  this->_input->GetPixelSize(&xsize, &ysize, &zsize);

  // Instantiate convolution filter
  irtkConvolution_1D<VoxelType> *convolution = NULL;
  if (this->_input->HasBackgroundValue()) {
    convolution = new irtkConvolutionWithPadding_1D<VoxelType>(this->_input->GetBackgroundValueAsDouble());
  } else {
    convolution = new irtkConvolution_1D<VoxelType>();
  }

  // Create scalar function which corresponds to a 1D Gaussian function in X
  irtkScalarGaussianDx gaussianDx(this->_Sigma/xsize, 1, 1, 0, 0, 0);
  // Create filter kernel for 1D Gaussian function in X
  irtkGenericImage<irtkRealPixel> kernelX(2*round(4*this->_Sigma/xsize)+1, 1, 1);
  // Do conversion from  scalar function to filter kernel
  irtkScalarFunctionToImage<irtkRealPixel> gaussianSourceX;
  gaussianSourceX.SetInput (&gaussianDx);
  gaussianSourceX.SetOutput(&kernelX);
  gaussianSourceX.Run();

  // Do convolution
  convolution->SetInput ( this->_input);
  convolution->SetInput2(&kernelX);
  convolution->SetOutput(this->_output);
  convolution->irtkImageToImage<VoxelType>::Run();

  // Flip x and y axis of image
  this->_output->FlipXY(1);

  // Create scalar function which corresponds to a 1D Gaussian function in Y
  irtkScalarGaussian gaussianY(this->_Sigma/ysize, 1, 1, 0, 0, 0);

  // Create filter kernel for 1D Gaussian function in Y
  irtkGenericImage<irtkRealPixel> kernelY(2*round(4*this->_Sigma/ysize)+1, 1, 1);

  // Do conversion from  scalar function to filter kernel
  irtkScalarFunctionToImage<irtkRealPixel> gaussianSourceY;
  gaussianSourceY.SetInput (&gaussianY);
  gaussianSourceY.SetOutput(&kernelY);
  gaussianSourceY.Run();

  // Do convolution
  convolution->SetInput (this->_output);
  convolution->SetInput2(&kernelY);
  convolution->SetOutput(this->_output);
  convolution->SetNormalization(true);
  convolution->irtkImageToImage<VoxelType>::Run();

  // Flip x and z axis of image
  this->_output->FlipXZ(1);

  if (this->_output->GetX() != 1) {
    // Create scalar function which corresponds to a 1D Gaussian function in Z
    irtkScalarGaussian gaussianZ(this->_Sigma/zsize, 1, 1, 0, 0, 0);

    // Create filter kernel for 1D Gaussian function in Z
    irtkGenericImage<irtkRealPixel> kernelZ(2*round(4*this->_Sigma/zsize)+1, 1, 1);

    // Do conversion from  scalar function to filter kernel
    irtkScalarFunctionToImage<irtkRealPixel> gaussianSourceZ;
    gaussianSourceZ.SetInput (&gaussianZ);
    gaussianSourceZ.SetOutput(&kernelZ);
    gaussianSourceZ.Run();

    // Do convolution
    convolution->SetInput (this->_output);
    convolution->SetInput2(&kernelZ);
    convolution->SetOutput(this->_output);
    convolution->SetNormalization(true);
    convolution->irtkImageToImage<VoxelType>::Run();
  }

  // Flip image back, first x and z axis, then x and y axis
  this->_output->FlipXZ(1);
  this->_output->FlipXY(1);

  // Do the final cleaning up
  delete convolution;
  this->Finalize();

}

template <class VoxelType> void irtkConvolutionWithGaussianDerivative<VoxelType>::Iy()
{
  double xsize, ysize, zsize;

  // Do the initial set up
  this->Initialize();

  // Get voxel dimensions
  this->_input->GetPixelSize(&xsize, &ysize, &zsize);

  // Instantiate convolution filter
  irtkConvolution_1D<VoxelType> *convolution = NULL;
  if (this->_input->HasBackgroundValue()) {
    convolution = new irtkConvolutionWithPadding_1D<VoxelType>(this->_input->GetBackgroundValueAsDouble());
  } else {
    convolution = new irtkConvolution_1D<VoxelType>();
  }

  // Create scalar function which corresponds to a 1D Gaussian function in X
  irtkScalarGaussian gaussianX(this->_Sigma/xsize, 1, 1, 0, 0, 0);
  // Create filter kernel for 1D Gaussian function in X
  irtkGenericImage<irtkRealPixel> kernelX(2*round(4*this->_Sigma/xsize)+1, 1, 1);

  // Do conversion from  scalar function to filter kernel
  irtkScalarFunctionToImage<irtkRealPixel> gaussianSourceX;
  gaussianSourceX.SetInput (&gaussianX);
  gaussianSourceX.SetOutput(&kernelX);
  gaussianSourceX.Run();

  // Do convolution
  convolution->SetInput ( this->_input);
  convolution->SetInput2(&kernelX);
  convolution->SetOutput(this->_output);
  convolution->SetNormalization(true);
  convolution->irtkImageToImage<VoxelType>::Run();

  // Flip x and y axis of image
  this->_output->FlipXY(1);

  // Create scalar function which corresponds to a 1D Gaussian function in Y
  irtkScalarGaussianDx gaussianDy(this->_Sigma/ysize, 1, 1, 0, 0, 0);

  // Create filter kernel for 1D Gaussian function in Y
  irtkGenericImage<irtkRealPixel> kernelY(2*round(4*this->_Sigma/ysize)+1, 1, 1);

  // Do conversion from  scalar function to filter kernel
  irtkScalarFunctionToImage<irtkRealPixel> gaussianSourceY;
  gaussianSourceY.SetInput (&gaussianDy);
  gaussianSourceY.SetOutput(&kernelY);
  gaussianSourceY.Run();

  // Do convolution
  convolution->SetInput (this->_output);
  convolution->SetInput2(&kernelY);
  convolution->SetOutput(this->_output);
  convolution->irtkImageToImage<VoxelType>::Run();

  // Flip x and z axis of image
  this->_output->FlipXZ(1);

  if (this->_output->GetX() != 1) {
    // Create scalar function which corresponds to a 1D Gaussian function in Z
    irtkScalarGaussian gaussianZ(this->_Sigma/zsize, 1, 1, 0, 0, 0);

    // Create filter kernel for 1D Gaussian function in Z
    irtkGenericImage<irtkRealPixel> kernelZ(2*round(4*this->_Sigma/zsize)+1, 1, 1);

    // Do conversion from  scalar function to filter kernel
    irtkScalarFunctionToImage<irtkRealPixel> gaussianSourceZ;
    gaussianSourceZ.SetInput (&gaussianZ);
    gaussianSourceZ.SetOutput(&kernelZ);
    gaussianSourceZ.Run();

    // Do convolution
    convolution->SetInput (this->_output);
    convolution->SetInput2(&kernelZ);
    convolution->SetOutput(this->_output);
    convolution->SetNormalization(true);
    convolution->irtkImageToImage<VoxelType>::Run();
  }

  // Flip image back, first x and z axis, then x and y axis
  this->_output->FlipXZ(1);
  this->_output->FlipXY(1);

  // Do the final cleaning up
  delete convolution;
  this->Finalize();
}

template <class VoxelType> void irtkConvolutionWithGaussianDerivative<VoxelType>::Iz()
{
  double xsize, ysize, zsize;

  // Do the initial set up
  this->Initialize();

  // Get voxel dimensions
  this->_input->GetPixelSize(&xsize, &ysize, &zsize);

  // Instantiate convolution filter
  irtkConvolution_1D<VoxelType> *convolution = NULL;
  if (this->_input->HasBackgroundValue()) {
    convolution = new irtkConvolutionWithPadding_1D<VoxelType>(this->_input->GetBackgroundValueAsDouble());
  } else {
    convolution = new irtkConvolution_1D<VoxelType>();
  }

  // Create scalar function which corresponds to a 1D Gaussian function in X
  irtkScalarGaussian gaussianX(this->_Sigma/xsize, 1, 1, 0, 0, 0);
  // Create filter kernel for 1D Gaussian function in X
  irtkGenericImage<irtkRealPixel> kernelX(2*round(4*this->_Sigma/xsize)+1, 1, 1);

  // Do conversion from  scalar function to filter kernel
  irtkScalarFunctionToImage<irtkRealPixel> gaussianSourceX;
  gaussianSourceX.SetInput (&gaussianX);
  gaussianSourceX.SetOutput(&kernelX);
  gaussianSourceX.Run();

  // Do convolution
  convolution->SetInput ( this->_input);
  convolution->SetInput2(&kernelX);
  convolution->SetOutput(this->_output);
  convolution->SetNormalization(true);
  convolution->irtkImageToImage<VoxelType>::Run();

  // Flip x and y axis of image
  this->_output->FlipXY(1);

  // Create scalar function which corresponds to a 1D Gaussian function in Y
  irtkScalarGaussian gaussianY(this->_Sigma/ysize, 1, 1, 0, 0, 0);

  // Create filter kernel for 1D Gaussian function in Y
  irtkGenericImage<irtkRealPixel> kernelY(2*round(4*this->_Sigma/ysize)+1, 1, 1);

  // Do conversion from  scalar function to filter kernel
  irtkScalarFunctionToImage<irtkRealPixel> gaussianSourceY;
  gaussianSourceY.SetInput (&gaussianY);
  gaussianSourceY.SetOutput(&kernelY);
  gaussianSourceY.Run();

  // Do convolution
  convolution->SetInput (this->_output);
  convolution->SetInput2(&kernelY);
  convolution->SetOutput(this->_output);
  convolution->SetNormalization(true);
  convolution->irtkImageToImage<VoxelType>::Run();

  // Flip x and z axis of image
  this->_output->FlipXZ(1);

  if (this->_output->GetX() != 1) {
    // Create scalar function which corresponds to a 1D Gaussian function in Z
    irtkScalarGaussianDx gaussianDz(this->_Sigma/zsize, 1, 1, 0, 0, 0);

    // Create filter kernel for 1D Gaussian function in Z
    irtkGenericImage<irtkRealPixel> kernelZ(2*round(4*this->_Sigma/zsize)+1, 1, 1);

    // Do conversion from  scalar function to filter kernel
    irtkScalarFunctionToImage<irtkRealPixel> gaussianSourceZ;
    gaussianSourceZ.SetInput (&gaussianDz);
    gaussianSourceZ.SetOutput(&kernelZ);
    gaussianSourceZ.Run();

    // Do convolution
    convolution->SetInput (this->_output);
    convolution->SetInput2(&kernelZ);
    convolution->SetOutput(this->_output);
    convolution->irtkImageToImage<VoxelType>::Run();
  }

  // Flip image back, first x and z axis, then x and y axis
  this->_output->FlipXZ(1);
  this->_output->FlipXY(1);

  // Do the final cleaning up
  delete convolution;
  this->Finalize();
}


template class irtkConvolutionWithGaussianDerivative<irtkBytePixel>;
template class irtkConvolutionWithGaussianDerivative<irtkGreyPixel>;
template class irtkConvolutionWithGaussianDerivative<float>;
template class irtkConvolutionWithGaussianDerivative<double>;
