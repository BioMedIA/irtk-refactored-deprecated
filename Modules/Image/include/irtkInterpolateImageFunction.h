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

#ifndef _IRTKINTERPOLATEIMAGEFUNCTION_H
#define _IRTKINTERPOLATEIMAGEFUNCTION_H

#include <irtkExtrapolateImageFunction.h>
#include <irtkVoxelFunction.h>


////////////////////////////////////////////////////////////////////////////////
// Abstract interpolation interface
////////////////////////////////////////////////////////////////////////////////

/**
 * Abstract base class for any general image interpolation function
 */

class irtkInterpolateImageFunction : public irtkImageFunction
{
  irtkAbstractMacro(irtkInterpolateImageFunction);

  // ---------------------------------------------------------------------------
  // Data members

  /// Number of dimensions to interpolate
  ///
  /// Determined either from input dimensions by default or set to a fixed
  /// number of dimension by specialized subclasses.
  irtkPublicAttributeMacro(int, NumberOfDimensions);

protected:

  /// Infinite discrete image obtained by extrapolation of finite input image.
  /// Unused by default, i.e., NULL which corresponds to extrapolation mode
  /// Extrapolation_None. If \c NULL, the interpolator has to deal with
  /// boundary cases itself either by only partially interpolating the
  /// available values or returning the _DefaultValue.
  irtkExtrapolateImageFunction *_InfiniteInput;

  /// Whether infinite discrete image was instantiated by this image function
  bool _InfiniteInputOwner;

  /// Domain of finite input image for which the interpolation is defined
  /// without requiring any extrapolation: [x1, x2]x[y1, y2]x[z1, z2]x[t1, t2]
  double _x1, _y1, _z1, _t1, _x2, _y2, _z2, _t2;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  irtkInterpolateImageFunction();

public:

  /// Destructor
  virtual ~irtkInterpolateImageFunction();

  /// Construct interpolator with default infinite extension of input image
  static irtkInterpolateImageFunction *New(irtkInterpolationMode, const irtkBaseImage * = NULL);

  /// Construct extrapolator which is compatible with this interpolator
  virtual irtkExtrapolateImageFunction *New(irtkExtrapolationMode,
                                            const irtkBaseImage * = NULL);

  /// Construct interpolator with specified infinite extension of input image
  ///
  /// The caller is required to set the input, initialize, and destroy the
  /// interpolator only, the extrapolator is initialized and destroyed by the
  /// interpolator unless the extrapolator has been replaced using the setter.
  static irtkInterpolateImageFunction *New(irtkInterpolationMode,
                                           irtkExtrapolationMode, const irtkBaseImage * = NULL);

  // ---------------------------------------------------------------------------
  // Initialization

  /// Set input image
  void Input(const irtkBaseImage *);

  /// Get input image
  const irtkBaseImage *Input() const;

  /// Get interpolation mode corresponding to this interpolator
  virtual irtkInterpolationMode InterpolationMode() const = 0;

  /// Get extrapolation mode used by this interpolator
  irtkExtrapolationMode ExtrapolationMode() const;

  /// Set extrapolate image function for evaluation outside of image domain
  virtual void Extrapolator(irtkExtrapolateImageFunction *, bool = false);

  /// Get extrapolate image function for evaluation outside of image domain
  /// or \c NULL if extrapolation mode is \c Extrapolation_None
  irtkExtrapolateImageFunction *Extrapolator();

  /// Get extrapolate image function for evaluation outside of image domain
  /// or \c NULL if extrapolation mode is \c Extrapolation_None
  const irtkExtrapolateImageFunction *Extrapolator() const;

  /// Initialize image function
  ///
  /// \note Override the virtual Initialize(bool) member function in subclasses,
  ///       but not this member function which is required by the abstract
  ///       image function interface (irtkImageFunction::Initialize).
  virtual void Initialize();

  /// Initialize image function
  ///
  /// \param[in] coeff Whether input image contains interpolation coefficients
  ///                  already. Otherwise, the interpolate image function will
  ///                  compute these coefficients from the input intensities.
  virtual void Initialize(bool coeff);

  // ---------------------------------------------------------------------------
  // Domain checks

  /// Convert world coordinates to image location (in pixels)
  void WorldToImage(double &, double &) const;

  /// Convert world coordinates to image location (in pixels)
  void WorldToImage(double &, double &, double &) const;

  /// Convert world coordinates to image location (in pixels)
  void WorldToImage(irtkPoint &) const;

  /// Returns the image domain for which this image interpolation function
  /// can be used without handling any form of boundary conditions
  void Inside(double &, double &,
              double &, double &) const;

  /// Returns the image domain for which this image interpolation function
  /// can be used without handling any form of boundary conditions
  void Inside(double &, double &, double &,
              double &, double &, double &) const;

  /// Returns the image domain for which this image interpolation function
  /// can be used without handling any form of boundary conditions
  void Inside(double &, double &, double &, double &,
              double &, double &, double &, double &) const;

  /// Returns interval of discrete image indices whose values are needed for
  /// interpolation of the image value at a given continuous coordinate
  virtual void BoundingInterval(double, int &, int &) const = 0;

  /// Returns discrete boundaries of local 2D image region needed for interpolation
  virtual void BoundingBox(double, double, int &, int &,
                                           int &, int &) const;

  /// Returns discrete boundaries of local 3D image region needed for interpolation
  virtual void BoundingBox(double, double, double, int &, int &, int &,
                                                   int &, int &, int &) const;

  /// Returns discrete boundaries of local 4D image region needed for interpolation
  virtual void BoundingBox(double, double, double, double,
                           int &, int &, int &, int &,
                           int &, int &, int &, int &) const;

  /// Check if the location (in pixels) is inside the domain for which this image
  /// interpolation can be used without handling any form of boundary condition
  bool IsInside(double, double) const;

  /// Check if the location (in pixels) is inside the domain for which this image
  /// interpolation can be used without handling any form of boundary condition
  bool IsInside(double, double, double) const;

  /// Check if the location (in pixels) is inside the domain for which this image
  /// interpolation can be used without handling any form of boundary condition
  bool IsInside(double, double, double, double) const;

  /// Check if the location (in pixels) is outside the domain for which this image
  /// interpolation can be used without handling any form of boundary condition
  bool IsOutside(double, double) const;

  /// Check if the location (in pixels) is outside the domain for which this image
  /// interpolation can be used without handling any form of boundary condition
  bool IsOutside(double, double, double) const;

  /// Check if the location (in pixels) is outside the domain for which this image
  /// interpolation can be used without handling any form of boundary condition
  bool IsOutside(double, double, double, double) const;

  /// Check if the location is fully inside the foreground of the image, i.e.,
  /// including all discrete image locations required for interpolation
  bool IsForeground(double, double) const;

  /// Check if the location is fully inside the foreground of the image, i.e.,
  /// including all discrete image locations required for interpolation
  bool IsForeground(double, double, double) const;

  /// Check if the location is fully inside the foreground of the image, i.e.,
  /// including all discrete image locations required for interpolation
  bool IsForeground(double, double, double, double) const;

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Evaluate scalar image without handling boundary conditions
  ///
  /// This version is faster than EvaluateOutside, but is only defined inside
  /// the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual double EvaluateInside(double, double, double = 0, double = 0) const = 0;

  /// Evaluate scalar image at an arbitrary location (in pixels)
  virtual double EvaluateOutside(double, double, double = 0, double = 0) const = 0;

  /// Evaluate scalar image at an arbitrary location (in pixels)
  ///
  /// If the location is inside the domain for which the filter can perform
  /// an interpolation without considering a particular boundary condition,
  /// the faster EvaluateInside method is called. Otherwise, the EvaluateOutside
  /// method which makes use of the extrapolation of the discrete image domain
  /// in order to interpolate also at boundary or outside locations is used.
  double Evaluate(double, double, double = 0, double = 0) const;

  /// Evaluate scalar image at an arbitrary location (in pixels)
  ///
  /// If the location is inside the domain for which the filter can perform
  /// an interpolation without considering a particular boundary condition,
  /// the faster EvaluateInside method is called. Otherwise, the EvaluateOutside
  /// method which makes use of the extrapolation of the discrete image domain
  /// in order to interpolate also at boundary or outside locations is used.
  ///
  /// \note This overloaded function corrects for the const-ness of the
  ///       corresponding virtual base class function irtkImageFunction::Evaluate.
  double Evaluate(double, double, double = 0, double = 0);

  /// Evaluate scalar image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  ///
  /// This version is faster than EvaluateWithPaddingOutside, but is only defined
  /// inside the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual double EvaluateWithPaddingInside(double, double, double = 0, double = 0) const = 0;

  /// Evaluate scalar image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  virtual double EvaluateWithPaddingOutside(double, double, double = 0, double = 0) const = 0;

  /// Evaluate scalar image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  ///
  /// If the location is inside the domain for which the filter can perform
  /// an interpolation without considering a particular boundary condition,
  /// the faster EvaluateWithPaddingInside method is called. Otherwise, the
  /// EvaluateWithPaddingOutside method which makes use of the extrapolation
  /// of the discrete image domain in order to interpolate also at boundary or
  /// outside locations is used.
  double EvaluateWithPadding(double, double, double = 0, double = 0) const;

  /// Evaluate multi-channel image without handling boundary conditions
  ///
  /// This version is faster than EvaluateOutside, but is only defined inside
  /// the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual void EvaluateInside(double *, double, double, double = 0, int = 1) const;

  /// Evaluate multi-channel image at an arbitrary location (in pixels)
  virtual void EvaluateOutside(double *, double, double, double = 0, int = 1) const;

  /// Evaluate multi-channel image at an arbitrary location (in pixels)
  ///
  /// If the location is inside the domain for which the filter can perform
  /// an interpolation without considering a particular boundary condition,
  /// the faster EvaluateInside method is called. Otherwise, the EvaluateOutside
  /// method which makes use of the extrapolation of the discrete image domain
  /// in order to interpolate also at boundary or outside locations is used.
  void Evaluate(double *, double, double, double = 0, int = 1) const;

  /// Evaluate multi-channel image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  ///
  /// This version is faster than EvaluateWithPaddingOutside, but is only defined
  /// inside the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual void EvaluateWithPaddingInside(double *, double, double, double = 0, int = 1) const;

  /// Evaluate multi-channel image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  virtual void EvaluateWithPaddingOutside(double *, double, double, double = 0, int = 1) const;

  /// Evaluate multi-channel image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  ///
  /// If the location is inside the domain for which the filter can perform
  /// an interpolation without considering a particular boundary condition,
  /// the faster EvaluateWithPaddingInside method is called. Otherwise, the
  /// EvaluateWithPaddingOutside method which makes use of the extrapolation
  /// of the discrete image domain in order to interpolate also at boundary or
  /// outside locations is used.
  void EvaluateWithPadding(double *, double, double, double = 0, int = 1) const;

  /// Evaluate vector image without handling boundary conditions
  ///
  /// This version is faster than EvaluateOutside, but is only defined inside
  /// the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual void EvaluateInside(irtkVector &, double, double, double = 0, double = 0) const = 0;

  /// Evaluate vector image at an arbitrary location (in pixels)
  virtual void EvaluateOutside(irtkVector &, double, double, double = 0, double = 0) const = 0;

  /// Evaluate vector image at an arbitrary location (in pixels)
  ///
  /// If the location is inside the domain for which the filter can perform
  /// an interpolation without considering a particular boundary condition,
  /// the faster EvaluateInside method is called. Otherwise, the EvaluateOutside
  /// method which makes use of the extrapolation of the discrete image domain
  /// in order to interpolate also at boundary or outside locations is used.
  void Evaluate(irtkVector &, double, double, double = 0, double = 0) const;

  /// Evaluate vector image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, a vector set to
  /// the _DefaultValue is returned.
  ///
  /// This version is faster than EvaluateWithPaddingOutside, but is only defined
  /// inside the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual void EvaluateWithPaddingInside(irtkVector &, double, double, double = 0, double = 0) const = 0;

  /// Evaluate vector image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, a vector set to
  /// the _DefaultValue is returned.
  virtual void EvaluateWithPaddingOutside(irtkVector &, double, double, double = 0, double = 0) const = 0;

  /// Evaluate vector image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, a vector set to
  /// the _DefaultValue is returned.
  ///
  /// If the location is inside the domain for which the filter can perform
  /// an interpolation without considering a particular boundary condition,
  /// the faster EvaluateWithPaddingInside method is called. Otherwise, the
  /// EvaluateWithPaddingOutside method which makes use of the extrapolation
  /// of the discrete image domain in order to interpolate also at boundary or
  /// outside locations is used.
  void EvaluateWithPadding(irtkVector &, double, double, double = 0, double = 0) const;

  /// Evaluate image function at all locations of the output image
  template <class TOutputImage>
  void Evaluate(TOutputImage &) const;

};

////////////////////////////////////////////////////////////////////////////////
// Generic interpolation interface
////////////////////////////////////////////////////////////////////////////////

/**
 * Abstract base class of generic interpolation functions
 *
 * This interpolation interface is templated over the input image type and
 * thus can access the image data using non-virtual getters which return
 * the image values with the appropriate voxel type. Therefore, it is more
 * efficient and type safe to use this interpolation interface whenever the
 * image type is known. Otherwise, use the abstract irtkInterpolateImageFunction
 * interface instead.
 *
 * \sa irtkInterpolateImageFunction
 */

template <class TImage>
class irtkGenericInterpolateImageFunction : public irtkInterpolateImageFunction
{
  irtkAbstractMacro(irtkGenericInterpolateImageFunction);

public:

  // ---------------------------------------------------------------------------
  // Types

  typedef TImage                                         ImageType;
  typedef typename ImageType::VoxelType                  VoxelType;
  typedef typename ImageType::RealType                   RealType;
  typedef typename voxel_info<RealType>::ScalarType      Real;
  typedef irtkGenericExtrapolateImageFunction<ImageType> ExtrapolatorType;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  irtkGenericInterpolateImageFunction();

public:

  /// Destructor
  virtual ~irtkGenericInterpolateImageFunction();

  /// Construct interpolator with default infinite extension of input image
  static irtkGenericInterpolateImageFunction *New(irtkInterpolationMode,
                                                  const TImage * = NULL);

  /// Construct extrapolator which is compatible with this interpolator
  virtual irtkExtrapolateImageFunction *New(irtkExtrapolationMode,
                                            const irtkBaseImage * = NULL);

  /// Construct interpolator with specified infinite extension of input image
  ///
  /// The caller is required to set the input, initialize, and destroy the
  /// interpolator only, the extrapolator is initialized and destroyed by the
  /// interpolator unless the extrapolator has been replaced using the setter.
  static irtkGenericInterpolateImageFunction *New(irtkInterpolationMode,
                                                  irtkExtrapolationMode,
                                                  const TImage * = NULL);

  // ---------------------------------------------------------------------------
  // Initialization

  /// Set input image
  virtual void Input(const irtkBaseImage *);

  /// Get input image
  const ImageType *Input() const;

  /// Set extrapolate image function for evaluation outside of image domain
  virtual void Extrapolator(irtkExtrapolateImageFunction *, bool = false);

  /// Get extrapolate image function for evaluation outside of image domain
  /// or \c NULL if extrapolation mode is \c Extrapolation_None
  ExtrapolatorType *Extrapolator();

  /// Get extrapolate image function for evaluation outside of image domain
  /// or \c NULL if extrapolation mode is \c Extrapolation_None
  const ExtrapolatorType *Extrapolator() const;

  /// Initialize image function
  ///
  /// \param[in] coeff Whether input image contains interpolation coefficients
  ///                  already. Otherwise, the interpolate image function will
  ///                  compute these coefficients from the input intensities.
  virtual void Initialize(bool coeff = false);

  // ---------------------------------------------------------------------------
  // Evaluation (generic type)

  /// Evaluate generic image without handling boundary conditions
  ///
  /// This version is faster than GetOutside, but is only defined inside
  /// the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual VoxelType GetInside(double, double, double = 0, double = 0) const = 0;

  /// Evaluate generic image at an arbitrary location (in pixels)
  virtual VoxelType GetOutside(double, double, double = 0, double = 0) const = 0;

  /// Evaluate generic image without handling boundary conditions
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  ///
  /// This version is faster than GetWithPaddingOutside, but is only defined
  /// inside the domain for which all image values required for interpolation
  /// are defined and thus require no extrapolation of the finite image.
  virtual VoxelType GetWithPaddingInside(double, double, double = 0, double = 0) const = 0;

  /// Evaluate generic image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  virtual VoxelType GetWithPaddingOutside(double, double, double = 0, double = 0) const = 0;

  /// Evaluate generic image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  VoxelType GetWithPadding(double, double, double = 0, double = 0) const;

  /// Evaluate generic image at an arbitrary location (in pixels)
  ///
  /// If the location is inside the domain for which the filter can perform
  /// an interpolation without considering a particular boundary condition,
  /// the faster GetInside method is called. Otherwise, the GetOutside method
  /// which makes use of the extrapolation of the discrete image domain in
  /// order to interpolate also at boundary or outside locations is used.
  VoxelType Get(double, double, double = 0, double = 0) const;

  /// Evaluate generic image at an arbitrary location (in pixels)
  VoxelType operator ()(double, double, double = 0, double = 0) const;

  // ---------------------------------------------------------------------------
  // Evaluation (general type)

  /// Evaluate scalar image without handling boundary conditions
  ///
  /// This version is faster than EvaluateOutside, but is only defined inside
  /// the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual double EvaluateInside(double, double, double = 0, double = 0) const;

  /// Evaluate scalar image at an arbitrary location (in pixels)
  virtual double EvaluateOutside(double, double, double = 0, double = 0) const;

  /// Evaluate scalar image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  ///
  /// This version is faster than EvaluateWithPaddingOutside, but is only defined
  /// inside the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual double EvaluateWithPaddingInside(double, double, double = 0, double = 0) const;

  /// Evaluate scalar image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  virtual double EvaluateWithPaddingOutside(double, double, double = 0, double = 0) const;

  /// Evaluate multi-channel image without handling boundary conditions
  ///
  /// This version is faster than EvaluateOutside, but is only defined inside
  /// the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual void EvaluateInside(double *, double, double, double = 0, int = 1) const;

  /// Evaluate multi-channel image at an arbitrary location (in pixels)
  virtual void EvaluateOutside(double *, double, double, double = 0, int = 1) const;

  /// Evaluate multi-channel image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  ///
  /// This version is faster than EvaluateWithPaddingOutside, but is only defined
  /// inside the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual void EvaluateWithPaddingInside(double *, double, double, double = 0, int = 1) const;

  /// Evaluate multi-channel image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  virtual void EvaluateWithPaddingOutside(double *, double, double, double = 0, int = 1) const;

  /// Evaluate vector image without handling boundary conditions
  ///
  /// This version is faster than EvaluateOutside, but is only defined inside
  /// the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual void EvaluateInside(irtkVector &, double, double, double = 0, double = 0) const;

  /// Evaluate vector image at an arbitrary location (in pixels)
  virtual void EvaluateOutside(irtkVector &, double, double, double = 0, double = 0) const;

  /// Evaluate vector image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, a vector set to
  /// the _DefaultValue is returned.
  ///
  /// This version is faster than EvaluateWithPaddingOutside, but is only defined
  /// inside the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual void EvaluateWithPaddingInside(irtkVector &, double, double, double = 0, double = 0) const;

  /// Evaluate vector image at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, a vector set to
  /// the _DefaultValue is returned.
  virtual void EvaluateWithPaddingOutside(irtkVector &, double, double, double = 0, double = 0) const;

};

////////////////////////////////////////////////////////////////////////////////
// Auxiliary macro for subclass implementation
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
#define irtkInterpolatorMacro(clsname, mode)                                   \
    irtkObjectMacro(clsname);                                                  \
  public:                                                                      \
    /** Get interpolation mode implemented by this interpolator */             \
    inline virtual irtkInterpolationMode InterpolationMode() const             \
    { return mode; }                                                           \
    /** Get interpolation mode implemented by this class */                    \
    inline static  irtkInterpolationMode InterpolationType()                   \
    { return mode; }                                                           \
  private:                                                                     \
    static void _irtkInterpolatorMacro_needs_trailing_semicolon()

// -----------------------------------------------------------------------------
#define irtkGenericInterpolatorTypes(superclsname)                             \
  public:                                                                      \
    typedef superclsname<TImage>                           Superclass;         \
    typedef typename Superclass::ImageType                 ImageType;          \
    typedef typename Superclass::VoxelType                 VoxelType;          \
    typedef typename Superclass::RealType                  RealType;           \
    typedef typename Superclass::Real                      Real;               \
    typedef typename Superclass::ExtrapolatorType          ExtrapolatorType;   \
  private:                                                                     \
    static void _irtkGenericInterpolatorTypedefsMacro_needs_trailing_semicolon()

// -----------------------------------------------------------------------------
#define irtkGenericInterpolatorMacro(clsname, mode)                            \
    irtkInterpolatorMacro(clsname, mode);                                      \
    irtkGenericInterpolatorTypes(irtkGenericInterpolateImageFunction);         \
  private:                                                                     \
     static void _irtkGenericInterpolatorMacro_needs_trailing_semicolon()

////////////////////////////////////////////////////////////////////////////////
// Inline definitions -- irtkInterpolateImageFunction
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
inline void irtkInterpolateImageFunction::Input(const irtkBaseImage *input)
{
  this->_input = const_cast<irtkBaseImage *>(input);
}

// -----------------------------------------------------------------------------
inline const irtkBaseImage *irtkInterpolateImageFunction::Input() const
{
  return this->_input;
}

// -----------------------------------------------------------------------------
inline void irtkInterpolateImageFunction::Initialize()
{
  this->Initialize(false);
}

// -----------------------------------------------------------------------------
inline void irtkInterpolateImageFunction
::Extrapolator(irtkExtrapolateImageFunction *input, bool owner)
{
  if (_InfiniteInputOwner) Delete(_InfiniteInput);
  _InfiniteInput      = input;
  _InfiniteInputOwner = _InfiniteInput && owner;
}

// -----------------------------------------------------------------------------
inline irtkExtrapolateImageFunction *irtkInterpolateImageFunction::Extrapolator()
{
  return _InfiniteInput;
}

// -----------------------------------------------------------------------------
inline const irtkExtrapolateImageFunction *irtkInterpolateImageFunction::Extrapolator() const
{
  return _InfiniteInput;
}

// -----------------------------------------------------------------------------
inline irtkExtrapolationMode irtkInterpolateImageFunction::ExtrapolationMode() const
{
  return (_InfiniteInput ? _InfiniteInput->ExtrapolationMode() : Extrapolation_None);
}

// =============================================================================
// Domain checks
// =============================================================================

// ----------------------------------------------------------------------------
inline void irtkInterpolateImageFunction::WorldToImage(double &x, double &y) const
{
  this->_input->WorldToImage(x, y);
}

// ----------------------------------------------------------------------------
inline void irtkInterpolateImageFunction::WorldToImage(double &x, double &y, double &z) const
{
  this->_input->WorldToImage(x, y, z);
}

// ----------------------------------------------------------------------------
inline void irtkInterpolateImageFunction::WorldToImage(irtkPoint &p) const
{
  this->_input->WorldToImage(p);
}

// -----------------------------------------------------------------------------
inline void irtkInterpolateImageFunction::Inside(double &x1, double &y1,
                                                 double &x2, double &y2) const
{
  x1 = _x1, y1 = _y1;
  x2 = _x2, y2 = _y2;
}

// -----------------------------------------------------------------------------
inline void irtkInterpolateImageFunction::Inside(double &x1, double &y1, double &z1,
                                                 double &x2, double &y2, double &z2) const
{
  x1 = _x1, y1 = _y1, z1 = _z1;
  x2 = _x2, y2 = _y2, z2 = _z2;
}

// -----------------------------------------------------------------------------
inline void irtkInterpolateImageFunction::Inside(double &x1, double &y1, double &z1, double &t1,
                                                 double &x2, double &y2, double &z2, double &t2) const
{
  x1 = _x1, y1 = _y1, z1 = _z1, t1 = _t1;
  x2 = _x2, y2 = _y2, z2 = _z2, t2 = _t2;
}

// -----------------------------------------------------------------------------
inline void irtkInterpolateImageFunction::BoundingBox(double x, double y,
                                                      int &i1, int &j1,
                                                      int &i2, int &j2) const
{
  this->BoundingInterval(x, i1, i2);
  this->BoundingInterval(y, j1, j2);
}

// -----------------------------------------------------------------------------
inline void irtkInterpolateImageFunction::BoundingBox(double x, double y, double z,
                                                      int &i1, int &j1, int &k1,
                                                      int &i2, int &j2, int &k2) const
{
  if (this->NumberOfDimensions() >= 3) {
    this->BoundingInterval(x, i1, i2);
    this->BoundingInterval(y, j1, j2);
    this->BoundingInterval(z, k1, k2);
  } else {
    this->BoundingInterval(x, i1, i2);
    this->BoundingInterval(y, j1, j2);
    k1 = k2 = round(z);
  }
}

// -----------------------------------------------------------------------------
inline void irtkInterpolateImageFunction::BoundingBox(double x, double y, double z, double t,
                                                      int &i1, int &j1, int &k1, int &l1,
                                                      int &i2, int &j2, int &k2, int &l2) const
{
  if (this->NumberOfDimensions() >= 4) {
    this->BoundingInterval(x, i1, i2);
    this->BoundingInterval(y, j1, j2);
    this->BoundingInterval(z, k1, k2);
    this->BoundingInterval(t, l1, l2);
  } else if (this->NumberOfDimensions() == 3) {
    this->BoundingInterval(x, i1, i2);
    this->BoundingInterval(y, j1, j2);
    this->BoundingInterval(z, k1, k2);
    l1 = l2 = round(t);
  } else {
    this->BoundingInterval(x, i1, i2);
    this->BoundingInterval(y, j1, j2);
    k1 = k2 = round(z);
    l1 = l2 = round(t);
  }
}

// -----------------------------------------------------------------------------
inline bool irtkInterpolateImageFunction::IsInside(double x, double y) const
{
  return (_x1 <= x && x <= _x2) && (_y1 <= y && y <= _y2);
}

// -----------------------------------------------------------------------------
inline bool irtkInterpolateImageFunction::IsInside(double x, double y, double z) const
{
  return (_x1 <= x && x <= _x2) && (_y1 <= y && y <= _y2) && (_z1 <= z && z <= _z2);
}

// -----------------------------------------------------------------------------
inline bool irtkInterpolateImageFunction::IsInside(double x, double y, double z, double t) const
{
  return (_x1 <= x && x <= _x2) && (_y1 <= y && y <= _y2) && (_z1 <= z && z <= _z2) && (_t1 <= t && t <= _t2);
}

// -----------------------------------------------------------------------------
inline bool irtkInterpolateImageFunction::IsOutside(double x, double y) const
{
  return !IsInside(x, y);
}

// -----------------------------------------------------------------------------
inline bool irtkInterpolateImageFunction::IsOutside(double x, double y, double z) const
{
  return !IsInside(x, y, z);
}

// -----------------------------------------------------------------------------
inline bool irtkInterpolateImageFunction::IsOutside(double x, double y, double z, double t) const
{
  return !IsInside(x, y, z, t);
}

// -----------------------------------------------------------------------------
inline bool irtkInterpolateImageFunction::IsForeground(double x, double y) const
{
  int i1, j1, i2, j2;
  BoundingBox(x, y,                             i1, j1, i2, j2);
  return Input()->IsBoundingBoxInsideForeground(i1, j1, i2, j2);
}

// -----------------------------------------------------------------------------
inline bool irtkInterpolateImageFunction::IsForeground(double x, double y, double z) const
{
  int i1, j1, k1, i2, j2, k2;
  BoundingBox(x, y, z,                          i1, j1, k1, i2, j2, k2);
  return Input()->IsBoundingBoxInsideForeground(i1, j1, k1, i2, j2, k2);
}

// -----------------------------------------------------------------------------
inline bool irtkInterpolateImageFunction::IsForeground(double x, double y, double z, double t) const
{
  int i1, j1, k1, l1, i2, j2, k2, l2;
  BoundingBox(x, y, z, t,                       i1, j1, k1, l1, i2, j2, k2, l2);
  return Input()->IsBoundingBoxInsideForeground(i1, j1, k1, l1, i2, j2, k2, l2);
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
inline double irtkInterpolateImageFunction::Evaluate(double x, double y, double z, double t) const
{
  if (IsInside(x, y, z, t)) return this->EvaluateInside (x, y, z, t);
  else                      return this->EvaluateOutside(x, y, z, t);
}

// -----------------------------------------------------------------------------
inline double irtkInterpolateImageFunction::Evaluate(double x, double y, double z, double t)
{
  return const_cast<const irtkInterpolateImageFunction *>(this)->Evaluate(x, y, z, t);
}

// -----------------------------------------------------------------------------
inline double irtkInterpolateImageFunction::EvaluateWithPadding(double x, double y, double z, double t) const
{
  if (IsInside(x, y, z, t)) return this->EvaluateWithPaddingInside (x, y, z, t);
  else                      return this->EvaluateWithPaddingOutside(x, y, z, t);
}

// -----------------------------------------------------------------------------
inline void irtkInterpolateImageFunction::EvaluateInside(double *v, double x, double y, double z, int vt) const
{
  for (int t = 0; t < Input()->T(); ++t, v += vt) {
    (*v) = this->EvaluateInside(x, y, z, static_cast<double>(t));
  }
}

// -----------------------------------------------------------------------------
inline void irtkInterpolateImageFunction::EvaluateOutside(double *v, double x, double y, double z, int vt) const
{
  for (int t = 0; t < Input()->T(); ++t, v += vt) {
    (*v) = this->EvaluateOutside(x, y, z, static_cast<double>(t));
  }
}

// -----------------------------------------------------------------------------
inline void irtkInterpolateImageFunction::Evaluate(double *v, double x, double y, double z, int vt) const
{
  if (IsInside(x, y, z)) this->EvaluateInside (v, x, y, z, vt);
  else                   this->EvaluateOutside(v, x, y, z, vt);
}

// -----------------------------------------------------------------------------
inline void irtkInterpolateImageFunction::EvaluateWithPaddingInside(double *v, double x, double y, double z, int vt) const
{
  for (int t = 0; t < Input()->T(); ++t, v += vt) {
    (*v) = this->EvaluateWithPaddingInside(x, y, z, static_cast<double>(t));
  }
}

// -----------------------------------------------------------------------------
inline void irtkInterpolateImageFunction::EvaluateWithPaddingOutside(double *v, double x, double y, double z, int vt) const
{
  for (int t = 0; t < Input()->T(); ++t, v += vt) {
    (*v) = this->EvaluateWithPaddingOutside(x, y, z, static_cast<double>(t));
  }
}

// -----------------------------------------------------------------------------
inline void irtkInterpolateImageFunction::EvaluateWithPadding(double *v, double x, double y, double z, int vt) const
{
  if (IsInside(x, y, z)) this->EvaluateWithPaddingInside (v, x, y, z, vt);
  else                   this->EvaluateWithPaddingOutside(v, x, y, z, vt);
}

// -----------------------------------------------------------------------------
inline void irtkInterpolateImageFunction::Evaluate(irtkVector &v, double x, double y, double z, double t) const
{
  if (IsInside(x, y, z, t)) this->EvaluateInside (v, x, y, z, t);
  else                      this->EvaluateOutside(v, x, y, z, t);
}

// -----------------------------------------------------------------------------
inline void irtkInterpolateImageFunction::EvaluateWithPadding(irtkVector &v, double x, double y, double z, double t) const
{
  if (IsInside(x, y, z, t)) this->EvaluateWithPaddingInside (v, x, y, z, t);
  else                      this->EvaluateWithPaddingOutside(v, x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TOutputImage>
inline void irtkInterpolateImageFunction::Evaluate(TOutputImage &output) const
{
  // Interpolate multi-channel image (or 3D+t vector image)
  if (output.GetTSize() == .0) {
    irtkUnaryVoxelFunction::InterpolateMultiChannelImage<irtkInterpolateImageFunction> eval(this, &output);
    ParallelForEachVoxel(output.Attributes(), output, eval);
  // Interpolate scalar image
  } else if (output.N() == 1) {
    irtkUnaryVoxelFunction::InterpolateScalarImage<irtkInterpolateImageFunction> eval(this, &output);
    ParallelForEachVoxel(output.Attributes(), output, eval);
  // Interpolate vector image
  } else {
    irtkUnaryVoxelFunction::InterpolateImage<irtkInterpolateImageFunction> eval(this, &output);
    ParallelForEachVoxel(output.Attributes(), output, eval);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Inline definitions -- irtkGenericInterpolateImageFunction
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Getters
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
inline const TImage *irtkGenericInterpolateImageFunction<TImage>::Input() const
{
  return reinterpret_cast<const TImage *>(this->_input);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename irtkGenericInterpolateImageFunction<TImage>::ExtrapolatorType *
irtkGenericInterpolateImageFunction<TImage>::Extrapolator()
{
  return reinterpret_cast<ExtrapolatorType *>(_InfiniteInput);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline const typename irtkGenericInterpolateImageFunction<TImage>::ExtrapolatorType *
irtkGenericInterpolateImageFunction<TImage>::Extrapolator() const
{
  return reinterpret_cast<const ExtrapolatorType *>(_InfiniteInput);
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
inline typename TImage::VoxelType
irtkGenericInterpolateImageFunction<TImage>
::Get(double x, double y, double z, double t) const
{
  if (IsInside(x, y, z, t)) return this->GetInside (x, y, z, t);
  else                      return this->GetOutside(x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename TImage::VoxelType
irtkGenericInterpolateImageFunction<TImage>
::GetWithPadding(double x, double y, double z, double t) const
{
  if (IsInside(x, y, z, t)) return this->GetWithPaddingInside (x, y, z, t);
  else                      return this->GetWithPaddingOutside(x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename TImage::VoxelType irtkGenericInterpolateImageFunction<TImage>
::operator ()(double x, double y, double z, double t) const
{
  return this->Get(x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline double irtkGenericInterpolateImageFunction<TImage>
::EvaluateInside(double x, double y, double z, double t) const
{
  return voxel_cast<double>(this->GetInside(x, y, z, t));
}

// -----------------------------------------------------------------------------
template <class TImage>
inline double irtkGenericInterpolateImageFunction<TImage>
::EvaluateOutside(double x, double y, double z, double t) const
{
  return voxel_cast<double>(this->GetOutside(x, y, z, t));
}

// -----------------------------------------------------------------------------
template <class TImage>
inline double irtkGenericInterpolateImageFunction<TImage>
::EvaluateWithPaddingInside(double x, double y, double z, double t) const
{
  return voxel_cast<double>(this->GetWithPaddingInside(x, y, z, t));
}

// -----------------------------------------------------------------------------
template <class TImage>
inline double irtkGenericInterpolateImageFunction<TImage>
::EvaluateWithPaddingOutside(double x, double y, double z, double t) const
{
  return voxel_cast<double>(this->GetWithPaddingOutside(x, y, z, t));
}

// -----------------------------------------------------------------------------
template <class TImage>
inline void irtkGenericInterpolateImageFunction<TImage>
::EvaluateInside(double *v, double x, double y, double z, int vt) const
{
  for (int t = 0; t < Input()->T(); ++t, v += vt) {
    (*v) = voxel_cast<double>(this->GetInside(x, y, z, static_cast<double>(t)));
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
inline void irtkGenericInterpolateImageFunction<TImage>
::EvaluateOutside(double *v, double x, double y, double z, int vt) const
{
  for (int t = 0; t < Input()->T(); ++t, v += vt) {
    (*v) = voxel_cast<double>(this->GetOutside(x, y, z, static_cast<double>(t)));
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
inline void irtkGenericInterpolateImageFunction<TImage>
::EvaluateWithPaddingInside(double *v, double x, double y, double z, int vt) const
{
  for (int t = 0; t < Input()->T(); ++t, v += vt) {
    (*v) = voxel_cast<double>(this->GetWithPaddingInside(x, y, z, static_cast<double>(t)));
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
inline void irtkGenericInterpolateImageFunction<TImage>
::EvaluateWithPaddingOutside(double *v, double x, double y, double z, int vt) const
{
  for (int t = 0; t < Input()->T(); ++t, v += vt) {
    (*v) = voxel_cast<double>(this->GetWithPaddingOutside(x, y, z, static_cast<double>(t)));
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
inline void irtkGenericInterpolateImageFunction<TImage>
::EvaluateInside(irtkVector &v, double x, double y, double z, double t) const
{
  v = voxel_cast<irtkVector>(this->GetInside(x, y, z, t));
}

// -----------------------------------------------------------------------------
template <class TImage>
inline void irtkGenericInterpolateImageFunction<TImage>
::EvaluateOutside(irtkVector &v, double x, double y, double z, double t) const
{
  v = voxel_cast<irtkVector>(this->GetOutside(x, y, z, t));
}

// -----------------------------------------------------------------------------
template <class TImage>
inline void irtkGenericInterpolateImageFunction<TImage>
::EvaluateWithPaddingInside(irtkVector &v, double x, double y, double z, double t) const
{
  v = voxel_cast<irtkVector>(this->GetWithPaddingInside(x, y, z, t));
}

// -----------------------------------------------------------------------------
template <class TImage>
inline void irtkGenericInterpolateImageFunction<TImage>
::EvaluateWithPaddingOutside(irtkVector &v, double x, double y, double z, double t) const
{
  v = voxel_cast<irtkVector>(this->GetWithPaddingOutside(x, y, z, t));
}

////////////////////////////////////////////////////////////////////////////////
// Interpolate image functions
////////////////////////////////////////////////////////////////////////////////

// ND
#include <irtkNearestNeighborInterpolateImageFunction.h>
#include <irtkLinearInterpolateImageFunction.h>
#include <irtkBSplineInterpolateImageFunction.h>
#include <irtkCubicBSplineInterpolateImageFunction.h>
#include <irtkFastCubicBSplineInterpolateImageFunction.h>
#include <irtkCSplineInterpolateImageFunction.h>
#include <irtkGaussianInterpolateImageFunction.h>
#include <irtkSincInterpolateImageFunction.h>

// 2D
#include <irtkLinearInterpolateImageFunction2D.h>
#include <irtkBSplineInterpolateImageFunction2D.h>
#include <irtkCubicBSplineInterpolateImageFunction2D.h>
#include <irtkFastCubicBSplineInterpolateImageFunction2D.h>
#include <irtkCSplineInterpolateImageFunction2D.h>
#include <irtkGaussianInterpolateImageFunction2D.h>
#include <irtkSincInterpolateImageFunction2D.h>

// 3D
#include <irtkLinearInterpolateImageFunction3D.h>
#include <irtkBSplineInterpolateImageFunction3D.h>
#include <irtkCubicBSplineInterpolateImageFunction3D.h>
#include <irtkFastCubicBSplineInterpolateImageFunction3D.h>
#include <irtkCSplineInterpolateImageFunction3D.h>
#include <irtkGaussianInterpolateImageFunction3D.h>
#include <irtkSincInterpolateImageFunction3D.h>

#include <irtkShapeBasedInterpolateImageFunction.h> // 3D scalar image only

// 4D
#include <irtkLinearInterpolateImageFunction4D.h>
#include <irtkBSplineInterpolateImageFunction4D.h>
#include <irtkCubicBSplineInterpolateImageFunction4D.h>
#include <irtkFastCubicBSplineInterpolateImageFunction4D.h>
#include <irtkCSplineInterpolateImageFunction4D.h>
#include <irtkGaussianInterpolateImageFunction4D.h>
#include <irtkSincInterpolateImageFunction4D.h>


#endif
