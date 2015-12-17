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
#include <irtkConvolutionFunction.h>
#include <irtkScalarGaussian.h>
#include <irtkScalarFunctionToImage.h>
using namespace irtkConvolutionFunction;

static const char *test_file   = NULL;
static const char *result_file = NULL;
static double      sigma       = 1.0;
static double      bgvalue     = 0.0;

// ===========================================================================
// Kernels
// ===========================================================================

irtkGenericImage<irtkRealPixel> Gaussian(double sigma, double dx)
{
  sigma = sigma / dx;
  irtkScalarGaussian                       gaussian(sigma, 1, 1, 0, 0, 0);
  irtkGenericImage<irtkRealPixel>          kernel(2 * int(round(4.0 * sigma)) + 1, 1, 1);
  irtkScalarFunctionToImage<irtkRealPixel> generator;
  generator.SetInput (&gaussian);
  generator.SetOutput(&kernel);
  generator.Run();
  irtkRealPixel sum = .0;
  for (int i = 0; i < kernel.GetX(); ++i) sum += kernel(i, 0, 0, 0);
  for (int i = 0; i < kernel.GetX(); ++i) kernel(i, 0, 0, 0) /= sum;
  return kernel;
}

// ===========================================================================
// Tests
// ===========================================================================

// ---------------------------------------------------------------------------
TEST(irtkConvolutionFunction, ConvolveMirroredForeground)
{
  irtkGenericImage<irtkRealPixel> image (test_file);
  irtkGenericImage<irtkRealPixel> output(image.GetImageAttributes());
  irtkGenericImage<irtkRealPixel> xkernel = Gaussian(sigma, image.GetXSize());
  irtkGenericImage<irtkRealPixel> ykernel = Gaussian(sigma, image.GetYSize());
  irtkGenericImage<irtkRealPixel> zkernel = Gaussian(sigma, image.GetZSize());
  image.PutBackgroundValueAsDouble(bgvalue, true);
  ConvolveMirroredForegroundInX<irtkRealPixel> xconv(&image, &xkernel);
  ParallelForEachVoxel(image.GetImageAttributes(), image, output, xconv);
  ConvolveMirroredForegroundInY<irtkRealPixel> yconv(&image, &ykernel);
  ParallelForEachVoxel(image.GetImageAttributes(), output, image, yconv);
  ConvolveMirroredForegroundInZ<irtkRealPixel> zconv(&image, &zkernel);
  ParallelForEachVoxel(image.GetImageAttributes(), image, output, zconv);
  output.Write(result_file);
}

// ===========================================================================
// Main
// ===========================================================================

// ---------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  if (argc < 3 || argc > 5) {
    cerr << "usage: " << argv[0] << " <input> <output> [sigma] [background]" << endl;
    exit(1);
  }
  test_file   = argv[1];
  result_file = argv[2];
  if (argc > 3) sigma   = atof(argv[3]);
  if (argc > 4) bgvalue = atof(argv[4]);
  return RUN_ALL_TESTS();
}
