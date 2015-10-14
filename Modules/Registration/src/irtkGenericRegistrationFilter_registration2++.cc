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

#include <irtkGenericRegistrationFilter.h>

#include <irtkResampling.h>
#include <irtkResamplingWithPadding.h>
#include <irtkGaussianBlurring.h>
#include <irtkGaussianBlurringWithPadding.h>


// =============================================================================
// Adjust number of histogram bins - copied from registration2++/src/irtkUtil.cc
// =============================================================================

// -----------------------------------------------------------------------------
static int irtkCalculateNumberOfBins(irtkGreyImage *image, int maxbin, int min, int max)
{
  int range = max - min + 1;
  int nbins = max - min + 1;
  int width = 1;

  // Calculate number of bins to use
  if (maxbin > 0) {
    while (int(ceil(range/(double)width)) > maxbin) ++width;
    nbins = int(ceil(range/(double)width));
    cout << "Using " << nbins << " out of " << maxbin << " bin(s) with width " << width << endl;
  } else {
    cout << "Using " << nbins << " bin(s) with width " << width << endl;
  }

  // Rescale intensities to the number of bins
  irtkGreyPixel *ptr = image->GetPointerToVoxels();
  for (int i = 0; i < image->GetNumberOfVoxels(); i++) {
    if (*ptr > 0) *ptr = int(*ptr/(double)width);
    ptr++;
  }

  return nbins;
}

// =============================================================================
// Implementation of InitializeInput as implemented by registration2++
// =============================================================================

// -----------------------------------------------------------------------------
// Initial implementation of irtkGenericRegistration::Initialize
//
// The only difference between the very first implementation, here referred to
// as v2.0 of the IRTK image registration, and a following modification with the
// objective to improve performance is the order of the blurring and resampling
// steps. This second initialization, where the images are first resampled and
// then blurred is referred to as v2.1 of the image registration.
void irtkGenericRegistrationFilter::InitializePyramid_v20_or_v21()
{
  string    msg; // Log message
  const int N = NumberOfImages();

  if (N != 2) {
    cerr << "irtkGenericRegistrationFilter::InitializeInput: Only implemented for"
         << " 2 input images (target and source) if Version set to 2.0 or 2.1" << endl;
    exit(1);
  }

  _Image.resize(_NumberOfLevels + 1);
  for (int l = 1; l < _NumberOfLevels; ++l) {

    // Copy/cast input images
    // Note: The registration2++ package was initially limited to grey values
    irtkGreyImage *image = Allocate<irtkGreyImage>(N);
    for (int n = 0; n < N; ++n) image[n] = *_Input[n];

    // Resample/blur images
    for (int n = 0; n < N; ++n) {
      const double                padding = _Padding[n];
      const double                sigma   = _Blurring  [_CurrentLevel][n];
      const irtkVector3D<double> &res     = _Resolution[_CurrentLevel][n];
      // Blur image before resampling (v2.0)
      if (version == irtkVersion(2, 0) && sigma > 0) {
        msg = string("Blurring image ") + ToString(n + 1, 3, '0') + " ... ";
        Broadcast(LogEvent, msg.c_str());
        irtkGaussianBlurringWithPadding<irtkGreyPixel> blurring(sigma, padding);
        blurring.SetInput (&image[n]);
        blurring.SetOutput(&image[n]);
        blurring.Run();
        Broadcast(LogEvent, "done\n");
      }
      // Resample image using linear interpolation
      double dx, dy, dz;
      image[n].GetPixelSize(&dx, &dy, &dz);
      if (!FinalLevel() || (fabs(res._x - dx) + fabs(res._y - dy) + fabs(res._z - dz)) > 1.0e-6) {
        msg = string("Resampling image ") + ToString(n + 1, 3, '0') + " ... ";
        Broadcast(LogEvent, msg.c_str());
        irtkResamplingWithPadding<irtkGreyPixel> resample(res._x, res._y, res._z, padding);
        resample.SetInput (&image[n]);
        resample.SetOutput(&image[n]);
        resample.Run();
        Broadcast(LogEvent, "done\n");
      }
      // Blur image after resampling (v2.1)
      if (version == irtkVersion(2, 1) && sigma > 0) {
        msg = string("Blurring image ") + ToString(n + 1, 3, '0') + " ... ";
        Broadcast(LogEvent, msg.c_str());
        irtkGaussianBlurringWithPadding<irtkGreyPixel> blurring(sigma, padding);
        blurring.SetInput (&image[n]);
        blurring.SetOutput(&image[n]);
        blurring.Run();
        cout << "done" << endl;
      }
    }

    // Find out the intensity range of the images
    // and apply background threshold (padding)
    double min[2] = {static_cast<double>(MAX_GREY), static_cast<double>(MAX_GREY)};
    double max[2] = {static_cast<double>(MIN_GREY), static_cast<double>(MIN_GREY)};

    for (int n = 0; n < N; ++n) {
      irtkGreyPixel *p    = image[n].GetPointerToVoxels();
      const int      nvox = image[n].GetNumberOfVoxels();
      for (int i = 0; i < nvox; ++i, ++p) {
        if (*p > _Padding[n]) {
          if (*p < min[n]) min[n] = *p;
          if (*p > max[n]) max[n] = *p;
        } else {
          *p = _Padding[n];
        }
      }
      if (max[n] == MIN_GREY && min[n] == MAX_GREY) {
        max[n] = _Padding[n];
        min[n] = _Padding[n];
      }
    }

    // Shift intensities such that min == 0
    irtkGreyPixel *p;
    int            nvox;

    p    = image[0].GetPointerToVoxels();
    nvox = image[0].GetNumberOfVoxels();
    if (max[0] - min[0] > MAX_GREY) {
      cerr << "irtkGenericRegistrationFilter::Initialize: Dynamic range of target is too large" << endl;
      delete[] image;
      exit(1);
    }
    for (int i = 0; i < nvox; ++i, ++p) {
      if (*p > _Padding[0]) *p -= min[0];
      else                  *p  = -1;
    }

    double shift = (_SimilarityMeasure == SSD ? min[0] : min[1]);
    if (max[1] - shift > MAX_GREY) {
      cerr << "irtkGenericRegistrationFilter::Initialize: Dynamic range of source is too large" << endl;
      delete[] image;
      exit(1);
    }
    p    = image[1].GetPointerToVoxels();
    nvox = image[1].GetNumberOfVoxels();
    for (int i = 0; i < nvox; ++i, ++p) {
      if (*p > _Padding[1]) *p -= shift;
      else                  *p  = -1;
    }

    // Compute offsets of background voxels (not used by FFD registration)
    for (int n = 0; n < N; ++n) {
      // Note: The remaining components of the registration3++ package assume
      //       that the background has a constant value. The offsets were
      //       only used to speed up the linear registration anyways.
      //irtkPadding(*image[n], ImagePadding(n));
      image[n].PutBackgroundValueAsDouble(-1.0);
    }

    // Adjust number of histogram bins and rescale intensities to [0, nbins-1]
    // Note: Actually the bins should be rescaled to [2, nbins-3] because the
    //       smoothing kernel for the Parzen window approximation is 5 bins wide.
    if (_SimilarityMeasure == NMI) {
      int nbins = 64;
      irtkParameterConstIterator it = Find(_Parameter[_CurrentLevel], "No. of bins");
      if (it != _Parameter[_CurrentLevel].end()) FromString(it->second.c_str(), nbins);
      else {
        it = Find(_Parameter[0], "No. of bins");
        if (it != _Parameter[0].end()) FromString(it->second.c_str(), nbins);
      }
      int tbins = irtkCalculateNumberOfBins(&image[0], nbins, min[0], max[0]);
      int sbins = irtkCalculateNumberOfBins(&image[1], nbins, min[1], max[1]);
      Insert(_Parameter[_CurrentLevel], "No. of target bins", ToString(tbins));
      Insert(_Parameter[_CurrentLevel], "No. of source bins", ToString(sbins));
    }

    // Copy grey values to internally used images
    _Image[l].resize(N);
    for (int n = 0; n < N; ++n) _Image[l][n] = image[n];
    Deallocate(image);

  }
}

// -----------------------------------------------------------------------------
void irtkGenericRegistrationFilter::InitializePyramid_v22()
{
  cerr << "irtkGenericRegistrationFilter::InitializePyramid_v22: Initialization of v2.2 currently not implemented" << endl;
  exit(1);
}
