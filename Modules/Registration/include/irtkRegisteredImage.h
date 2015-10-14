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

#ifndef _IRTKREGISTEREDIMAGE_H
#define _IRTKREGISTEREDIMAGE_H

#include <irtkTestProd.h>

#include <irtkCommon.h>
#include <irtkImage.h>
#include <irtkTransformation.h>


/**
 * Registered image such as fixed target image or transformed source image
 *
 * An instance of this (multi-channel) image class stores the intensities
 * of a registered image after alignment using the current transformation
 * estimate. If no (changing) transformation is set, the image is considered
 * to be fixed and is thus only updated once during the initialization.
 * If a transformation is set, however, upon every update of the input
 * transformation, the new moving image intensities are computed by applying
 * the current transformation and interpolating the input image intensities.
 * For the computation of the similarity measure gradient, the warped image
 * derivatives are stored along with the intensities in the consecutive image
 * channels of the registered image instance (i.e., t dimension). In case of
 * a gradient field similarity measure, these derivatives include not only the
 * 1st order, but also 2nd order derivatives needed for the similarity measure
 * gradient computation. An image similarity specifies which channels are
 * required for its evaluation by calling the Initialize function with a suitable
 * number of channels (i.e., 1, 4, or 10) and request their udpate using the
 * appropriate parameters of the Update function.
 *
 * The assignment of the transformed intensities and derivatives to the fourth
 * dimension (i.e. channels) of the registered image is as follows:
 *
 * - t=0: Transformed intensity
 * - t=1: Transformed 1st order derivative w.r.t x
 * - t=2: Transformed 1st order derivative w.r.t y
 * - t=3: Transformed 1st order derivative w.r.t z
 * - t=4: Transformed 2nd order derivative w.r.t xx
 * - t=5: Transformed 2nd order derivative w.r.t xy
 * - t=6: Transformed 2nd order derivative w.r.t xz
 * - t=7: Transformed 2nd order derivative w.r.t yy
 * - t=8: Transformed 2nd order derivative w.r.t yz
 * - t=9: Transformed 2nd order derivative w.r.t zz
 */
class irtkRegisteredImage : public irtkGenericImage<double>
{
  irtkObjectMacro(irtkRegisteredImage);

public:

  // Do not override other base class overloads
  using irtkGenericImage<double>::ImageToWorld;

  // ---------------------------------------------------------------------------
  // Types
public:

  /// Enumeration of registered image channel indices
  enum Channel { I = 0, Dx, Dy, Dz, Dxx, Dxy, Dxz, Dyx, Dyy, Dyz, Dzx, Dzy, Dzz };

  /// Type of untransformed input image
  typedef irtkGenericImage<double>   InputImageType;

  /// Type of untransformed gradient image
  typedef irtkGenericImage<double>   GradientImageType;

  /// Type of untransformed Hessian image
  typedef irtkGenericImage<double>   HessianImageType;

  /// Type of cached displacement fields
  typedef irtkGenericImage<double>   DisplacementImageType;

  // ---------------------------------------------------------------------------
  // Attributes

  /// Untransformed input image
  irtkPublicAggregateMacro(InputImageType, InputImage);

  /// Untransformed input gradient image
  irtkPublicComponentMacro(GradientImageType, InputGradient);

  /// Untransformed input Hessian image
  irtkPublicComponentMacro(HessianImageType, InputHessian);

  /// Current transformation estimate
  irtkPublicAggregateMacro(const irtkTransformation, Transformation);

  /// Interpolation mode
  irtkPublicAttributeMacro(irtkInterpolationMode, InterpolationMode);

  /// Extrapolation mode
  irtkPublicAttributeMacro(irtkExtrapolationMode, ExtrapolationMode);

  /// Pre-computed image to world coordinates
  irtkPublicAggregateMacro(irtkWorldCoordsImage, WorldCoordinates);

  /// Pre-computed image to world coordinates
  irtkPublicComponentMacro(irtkWorldCoordsImage, ImageToWorld);

  /// Externally pre-computed displacements to use
  irtkPublicAggregateMacro(DisplacementImageType, ExternalDisplacement);

  /// Pre-computed fixed displacements
  irtkComponentMacro(DisplacementImageType, FixedDisplacement);

  /// Pre-computed displacements
  irtkComponentMacro(DisplacementImageType, Displacement);

  /// Whether to use pre-computed world coordinates
  irtkAttributeMacro(bool, CacheWorldCoordinates);

  /// Whether to use pre-computed fixed displacements
  irtkAttributeMacro(bool, CacheFixedDisplacement);

  /// Whether to use pre-computed displacements
  irtkAttributeMacro(bool, CacheDisplacement);

  /// Whether self-update is enabled
  irtkPublicAttributeMacro(bool, SelfUpdate);

  /// Minimum foreground intensity of warped image or NaN
  irtkPublicAttributeMacro(double, MinIntensity);

  /// Maximum foreground intensity of warped image or NaN
  irtkPublicAttributeMacro(double, MaxIntensity);

  /// Standard deviation of Gaussian smoothing kernel applied before
  /// 1st order derivative computations in voxel units
  irtkPublicAttributeMacro(double, GradientSigma);

  /// Standard deviation of Gaussian smoothing kernel applied before
  /// 2nd order derivative computations in voxel units
  irtkPublicAttributeMacro(double, HessianSigma);

  /// Whether to precompute image derivatives
  /// true:  Compute derivatives of input image and transform these
  /// false: Use derivative of interpolation kernel to evaluate image derivative
  irtkPublicAttributeMacro(bool, PrecomputeDerivatives);

protected:

  /// Number of active levels
  int _NumberOfActiveLevels;

  /// Number of passive levels (fixed transformation)
  int _NumberOfPassiveLevels;

  /// Offsets of the different registered image channels
  int _Offset[13];

  /// (Pre-)compute gradient of input image
  /// \param[in] sigma Standard deviation of Gaussian smoothing filter in voxels.
  void ComputeInputGradient(double sigma);

  /// (Pre-)compute Hessian of input image
  /// \param[in] sigma Standard deviation of Gaussian smoothing filter in voxels.
  void ComputeInputHessian(double sigma);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  irtkRegisteredImage();

  /// Copy constructor
  irtkRegisteredImage(const irtkRegisteredImage &);

  /// Assignment operator
  irtkRegisteredImage &operator =(const irtkRegisteredImage &);

  /// Destructor
  ~irtkRegisteredImage();

  // ---------------------------------------------------------------------------
  // Channels

  /// Number of registered channels
  int NumberOfChannels() const;

  /// Number of voxels per registered channel
  int NumberOfVoxels() const;

  /// Offset of channel \p c with respect to the start of the image data
  int Offset(int) const;

  // ---------------------------------------------------------------------------
  // Initialization/Update

  /// Initialize image
  ///
  /// This function is called by the image similarity to set the attributes
  /// of the registered image, i.e., the attributes of the discrete image domain
  /// on which similarity is to be evaluated. It initializes the image to zero.
  /// The Update function must be called to initialize the output image
  /// intensities and derivatives before similarity can be evaluated.
  void Initialize(const irtkImageAttributes &, int = 0);

  /// Update image intensity and transformed derivatives
  ///
  /// This function only updates the first intensity channel within the image
  /// specified image region. If the image is not transformed or only warped by
  /// a fixed transformation, this function does nothing. Use \c force=true to
  /// initialize this image even if does not change over the course of the registration.
  ///
  /// \param[in] region    Image region to update.
  /// \param[in] intensity Request update of scalar intensities.
  /// \param[in] gradient  Request update of 1st order derivatives.
  /// \param[in] hessian   Request update of 2nd order derivatives.
  /// \param[in] force     Force update in any case.
  /// \param[in] disp      Pre-computed displacement field to use.
  ///                      If NULL, the default, the displacements corresponding
  ///                      to the set image transformation are used instead.
  void Update(const blocked_range3d<int> &region,
              bool intensity = true, bool gradient = false, bool hessian = false,
              bool force     = false);

  /// Update image intensity and transformed derivatives using custom displacements
  ///
  /// \param[in] region    Image region to update.
  /// \param[in] disp      Custom displacement field to use.
  /// \param[in] intensity Request update of scalar intensities.
  /// \param[in] gradient  Request update of 1st order derivatives.
  /// \param[in] hessian   Request update of 2nd order derivatives.
  void Update(const blocked_range3d<int> &region,
              const DisplacementImageType *disp,
              bool intensity = true, bool gradient = false, bool hessian = false);

  /// Update image intensity and transformed derivatives
  ///
  /// This function only updates the specified image channels if the self-update
  /// of this image is enabled and only if a (changing) transformation is set.
  /// If the image is not transformed or only warped by a fixed transformation,
  /// this function does nothing. Use \c force=true to initialize this image
  /// even if does not change over the course of the registration.
  ///
  /// \param[in] intensity Request update of scalar intensities.
  /// \param[in] gradient  Request update of 1st order derivatives.
  /// \param[in] hessian   Request update of 2nd order derivatives.
  /// \param[in] force     Force update in any case.
  void Update(bool intensity = true, bool gradient = false, bool hessian = false,
              bool force     = false);

private:

  template <class>                      void Update1(const blocked_range3d<int> &, bool, bool, bool);
  template <class, class, class, class> void Update2(const blocked_range3d<int> &, bool, bool, bool);
  template <class, class>               void Update3(const blocked_range3d<int> &, bool, bool, bool);

  FRIEND_TEST(irtkRegisteredImage, GlobalAndLocalTransformation);
};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline int irtkRegisteredImage::NumberOfChannels() const
{
  return _attr._t;
}

// -----------------------------------------------------------------------------
inline int irtkRegisteredImage::NumberOfVoxels() const
{
  return _attr._x * _attr._y * _attr._z;
}

// -----------------------------------------------------------------------------
inline int irtkRegisteredImage::Offset(int c) const
{
  return _Offset[c];
}


#endif
