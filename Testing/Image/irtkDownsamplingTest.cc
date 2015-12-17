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

#include <gtest/gtest.h>

#include <irtkImage.h>
#include <irtkDownsampling.h>
#include <irtkScalarGaussianDx.h>
#include <irtkGaussianBlurring.h>
#include <irtkResampling.h>
#include <irtkScalarFunctionToImage.h>
#include <irtkConvolutionFunction.h>
using namespace irtkConvolutionFunction;

static const char *test_file   = NULL;
static const char *result_file = NULL;
static double      sigma       = 0.0;
static double      bgvalue     = 0.0;
static int         factor      = 2;

// ===========================================================================
// Tests
// ===========================================================================

// ---------------------------------------------------------------------------
TEST(irtkDownsampling, GaussianPyramid1)
{
  char buffer[128];
  const int                       levels = 4;
  irtkGenericImage<irtkRealPixel> pyramid[levels];
  pyramid[0].Read(test_file);
  pyramid[0].PutBackgroundValueAsDouble(bgvalue);
  for (int l = 1; l < levels; ++l) {
    irtkGaussianBlurring<irtkRealPixel> blur(.7 * pyramid[l-1].GetXSize(),
                                             .7 * pyramid[l-1].GetYSize(),
                                             .7 * pyramid[l-1].GetZSize());
    blur.SetInput (&pyramid[l-1]);
    blur.SetOutput(&pyramid[l]);
    blur.Run();
    irtkResampling<irtkRealPixel> resample(pyramid[l-1].GetX() / 2,
                                           pyramid[l-1].GetY() / 2,
                                           pyramid[l-1].GetZ() / 2);
    resample.SetInterpolator(irtkInterpolateImageFunction::New(Interpolation_Linear, &pyramid[l]));
    resample.SetInput (&pyramid[l]);
    resample.SetOutput(&pyramid[l]);
    resample.Run();
    snprintf(buffer, 128, "pyramid1_level_%d.nii.gz", l);
    delete resample.GetInterpolator();
    pyramid[l].Write(buffer);
  }
}

// ---------------------------------------------------------------------------
TEST(irtkDownsampling, GaussianPyramid2)
{
  char buffer[128];
  const int                       levels = 4;
  irtkGenericImage<irtkRealPixel> pyramid[levels];
  pyramid[0].Read(test_file);
  pyramid[0].PutBackgroundValueAsDouble(bgvalue, true);
  irtkScalarGaussian   gaussian (0.7, 1.0, 1.0, .0, .0, .0);
  irtkScalarGaussianDx gaussianD(0.7, 1.0, 1.0, .0, .0, .0);
  irtkGenericImage<irtkRealPixel> kernel (2*round(4.*0.7)+1, 1, 1);
  irtkGenericImage<irtkRealPixel> kernelD(2*round(4.*0.7)+1, 1, 1);
  irtkScalarFunctionToImage<irtkRealPixel> gaussianSource;
  gaussianSource.SetInput (&gaussian);
  gaussianSource.SetOutput(&kernel);
  gaussianSource.Run();
  gaussianSource.SetInput (&gaussianD);
  gaussianSource.SetOutput(&kernelD);
  gaussianSource.Run();
  irtkRealPixel sum = .0;
  for (int i = 0; i < kernel.GetX(); ++i) sum += kernel(i, 0, 0);
  for (int i = 0; i < kernel.GetX(); ++i) kernel(i, 0, 0) /= sum;
  for (int l = 1; l < levels; ++l) {
    irtkImageAttributes attr = pyramid[l-1].GetImageAttributes();
    irtkGenericImage<irtkRealPixel> tmp;
    tmp = pyramid[l] = pyramid[l-1];
    ConvolveExtendedForegroundInX<irtkRealPixel> convx(&pyramid[l-1], &kernel);
    ParallelForEachVoxel(attr, pyramid[l-1], pyramid[l], convx);
    ConvolveExtendedForegroundInY<irtkRealPixel> convy(&pyramid[l-1], &kernel);
    ParallelForEachVoxel(attr, pyramid[l], tmp, convy);
    ConvolveExtendedForegroundInZ<irtkRealPixel> convz(&pyramid[l-1], &kernel);
    ParallelForEachVoxel(attr, tmp, pyramid[l], convz);
    irtkResampling<irtkRealPixel> resample(pyramid[l-1].GetX() / 2,
                                           pyramid[l-1].GetY() / 2,
                                           pyramid[l-1].GetZ() / 2);
    resample.SetInterpolator(irtkInterpolateImageFunction::New(Interpolation_Linear, &pyramid[l]));
    resample.SetInput (&pyramid[l]);
    resample.SetOutput(&pyramid[l]);
    resample.Run();
    snprintf(buffer, 128, "pyramid2_level_%d.nii.gz", l);
    delete resample.GetInterpolator();
    pyramid[l].Write(buffer);
  }
}

/*
// ---------------------------------------------------------------------------
TEST(irtkDownsampling, Run)
{
  irtkGreyImage      image(test_file);
  irtkScalarGaussian gaussian(sigma, 1.0, 1.0, .0, .0, .0);
  image.PutBackgroundValueAsDouble(bgvalue, true);
  irtkDownsampling<irtkGreyPixel> downsampler(factor);
  double x = -0.5, y = -0.5, z = -0.5;
  image.ImageToWorld(x, y, z);
  cout << "Start of image region before: (" << x << ", " << y << ", " << z << ")" << endl;
  downsampler.Kernel   (&gaussian, round(4.0 * sigma));
  downsampler.SetInput (&image);
  downsampler.SetOutput(&image);
  downsampler.Run();
  x = -0.5, y = -0.5, z = -0.5;
  image.ImageToWorld(x, y, z);
  cout << "Start of image region after:  (" << x << ", " << y << ", " << z << ")" << endl;
  image.Write(result_file);
}
*/

// ===========================================================================
// Main
// ===========================================================================

// ---------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  if (argc < 3 || argc > 6) {
    cerr << "usage: " << argv[0] << " <input> <output> [factor] [sigma] [background]" << endl;
    exit(1);
  }
  test_file   = argv[1];
  result_file = argv[2];
  if (argc > 3) factor  = atof(argv[3]);
  if (argc > 4) sigma   = atof(argv[4]);
  if (argc > 5) bgvalue = atof(argv[5]);
  return RUN_ALL_TESTS();
}
