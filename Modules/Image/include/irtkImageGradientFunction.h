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

#ifndef _IRTKIMAGEGRADIENTFUNCTION_H
#define _IRTKIMAGEGRADIENTFUNCTION_H

#include <irtkImageFunction.h>


/**
 * Abstract base class of image gradient functions
 *
 * An image gradient function applies the derivative of the interpolation
 * kernel instead of the interpolation kernel itself to evaluate the image
 * gradient at an arbitrary image location. It interpolates the 1st order
 * image derivatives instead of the image intensity value. Note that it
 * cannot be applied to multi-channel or vector-valued images.
 */
class irtkImageGradientFunction
{
  irtkAbstractMacro(irtkImageGradientFunction);

  // ---------------------------------------------------------------------------
  // Types

public:

  /// Type of interpolated gradient vectors
  typedef irtkVector3D<double> GradientType;

private:

  // ---------------------------------------------------------------------------
  // Data members

  /// Whether to evaluate image gradient w.r.t. world coordinate system
  irtkPublicAttributeMacro(bool, WrtWorld);

  /// Default value to return
  irtkPublicAttributeMacro(double, DefaultValue);

  /// Number of dimensions
  ///
  /// Determined either from input dimensions by default or set to a fixed
  /// number of dimension by specialized subclasses.
  irtkPublicAttributeMacro(int, NumberOfDimensions);

protected:

  /// Input image for filter
  irtkBaseImage *_Input;

  /// Image resolution
  GradientType _VoxelSize;

  /// Image orientation matrix
  irtkMatrix _Orientation;

  /// Infinite discrete image obtained by extrapolation of finite input image.
  /// Unused by default, i.e., NULL which corresponds to extrapolation mode
  /// Extrapolation_None. If \c NULL, the interpolator has to deal with
  /// boundary cases itself either by only partially interpolating the
  /// available values or returning the _DefaultValue.
  irtkExtrapolateImageFunction *_InfiniteInput;

  /// Whether infinite discrete image was instantiated by this image function
  bool _InfiniteInputOwner;

  /// Domain of finite input image for which the image gradient is defined
  /// without requiring any extrapolation: [x1, x2]x[y1, y2]x[z1, z2]x[t1, t2]
  double _x1, _y1, _z1, _t1, _x2, _y2, _z2, _t2;

  // ---------------------------------------------------------------------------
  // Auxiliary functions

  /// Orient and scale image gradient by world to image matrix
  void ImageGradientToWorld(GradientType &) const;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  irtkImageGradientFunction();

public:

  /// Destructor
  virtual ~irtkImageGradientFunction();

  /// Construct image gradient function with default infinite extension of input image
  static irtkImageGradientFunction *New(irtkInterpolationMode, const irtkBaseImage * = NULL);

  /// Construct extrapolator which is compatible with this image gradient function
  virtual irtkExtrapolateImageFunction *New(irtkExtrapolationMode,
                                            const irtkBaseImage * = NULL);

  /// Construct image gradient function with specified infinite extension of input image
  ///
  /// The caller is required to set the input, initialize, and destroy the
  /// image gradient function only, the extrapolator is initialized and destroyed
  /// by the image function unless the extrapolator has been replaced using the setter.
  static irtkImageGradientFunction *New(irtkInterpolationMode,
                                        irtkExtrapolationMode,
                                        const irtkBaseImage * = NULL);

  // ---------------------------------------------------------------------------
  // Initialization

  /// Set input image
  virtual void Input(const irtkBaseImage *);

  /// Get input image
  const irtkBaseImage *Input() const;

  /// Get interpolation mode corresponding to this image gradient function
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
  /// \param[in] coeff Whether input image contains interpolation coefficients
  ///                  already. Otherwise, the interpolate image function will
  ///                  compute these coefficients from the input intensities.
  virtual void Initialize(bool coeff = false);

  // ---------------------------------------------------------------------------
  // Domain checks

  /// Get temporal origin of input image
  double GetTOrigin() const;

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
  // Evaluation (C++ style return value, always w.r.t. image coordinates)

  /// Evaluate image gradient without handling boundary conditions
  ///
  /// This version is faster than GetOutside, but is only defined inside
  /// the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  virtual GradientType GetInside(double, double, double = 0, double = 0) const = 0;

  /// Evaluate generic image gradient at an arbitrary location (in pixels)
  virtual GradientType GetOutside(double, double, double = 0, double = 0) const = 0;

  /// Evaluate image gradient at an arbitrary location (in pixels)
  ///
  /// If the location is inside the domain for which the filter can perform
  /// an interpolation without considering a particular boundary condition,
  /// the faster GetInside method is called. Otherwise, the GetOutside method
  /// which makes use of the extrapolation of the discrete image domain in
  /// order to interpolate also at boundary or outside locations is used.
  GradientType Get(double, double, double = 0, double = 0) const;

  /// Evaluate image gradient without handling boundary conditions
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  ///
  /// This version is faster than GetWithPaddingOutside, but is only defined
  /// inside the domain for which all image values required for interpolation
  /// are defined and thus require no extrapolation of the finite image.
  virtual GradientType GetWithPaddingInside(double, double, double = 0, double = 0) const = 0;

  /// Evaluate image gradient at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  virtual GradientType GetWithPaddingOutside(double, double, double = 0, double = 0) const = 0;

  /// Evaluate image gradient at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  GradientType GetWithPadding(double, double, double = 0, double = 0) const;

  /// Evaluate image gradient at an arbitrary location (in pixels)
  GradientType operator ()(double, double, double = 0, double = 0) const;

  // ---------------------------------------------------------------------------
  // Evaluation (C style return value)

  // Note: This interface is intentionally identical to the
  //       irtkInterpolateImageFunction Evaluate* functions.

  /// Evaluate image gradient without handling boundary conditions
  ///
  /// This version is faster than EvaluateOutside, but is only defined inside
  /// the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  ///
  /// If \c _WrtWorld is \c true, the returned gradient vector is in
  /// world units and oriented according to the image orientation.
  virtual void EvaluateInside(double *, double, double, double = 0, int = 1) const;

  /// Evaluate image gradient at an arbitrary location (in pixels)
  ///
  /// If \c _WrtWorld is \c true, the returned gradient vector is in
  /// world units and oriented according to the image orientation.
  virtual void EvaluateOutside(double *, double, double, double = 0, int = 1) const;

  /// Evaluate image gradient at an arbitrary location (in pixels)
  ///
  /// If the location is inside the domain for which the filter can perform
  /// an interpolation without considering a particular boundary condition,
  /// the faster EvaluateInside method is called. Otherwise, the EvaluateOutside
  /// method which makes use of the extrapolation of the discrete image domain
  /// in order to interpolate also at boundary or outside locations is used.
  ///
  /// If \c _WrtWorld is \c true, the returned gradient vector is in
  /// world units and oriented according to the image orientation.
  void Evaluate(double *, double, double, double = 0, int = 1) const;

  /// Evaluate image gradient at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  ///
  /// This version is faster than EvaluateWithPaddingOutside, but is only defined
  /// inside the domain for which all image values required for interpolation are
  /// defined and thus require no extrapolation of the finite image.
  ///
  /// If \c _WrtWorld is \c true, the returned gradient vector is in
  /// world units and oriented according to the image orientation.
  virtual void EvaluateWithPaddingInside(double *, double, double, double = 0, int = 1) const;

  /// Evaluate image gradient at an arbitrary location (in pixels)
  ///
  /// If the location is partially inside the foreground region of the image,
  /// only the foreground values are interpolated. Otherwise, the _DefaultValue
  /// is returned.
  ///
  /// If \c _WrtWorld is \c true, the returned gradient vector is in
  /// world units and oriented according to the image orientation.
  virtual void EvaluateWithPaddingOutside(double *, double, double, double = 0, int = 1) const;

  /// Evaluate image gradient at an arbitrary location (in pixels)
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
  ///
  /// If \c _WrtWorld is \c true, the returned gradient vector is in
  /// world units and oriented according to the image orientation.
  void EvaluateWithPadding(double *, double, double, double = 0, int = 1) const;

};


////////////////////////////////////////////////////////////////////////////////
// Generic image gradient function interface
////////////////////////////////////////////////////////////////////////////////

/**
 * Abstract base class of generic image gradient functions
 *
 * This image function interface is templated over the input image type and
 * thus can access the image data using non-virtual getters which return
 * the image values with the appropriate voxel type. Therefore, it is more
 * efficient and type safe to use this image function interface whenever the
 * image type is known. Otherwise, use the abstract irtkImageGradientFunction
 * interface instead.
 *
 * \sa irtkImageGradientFunction
 */

template <class TImage>
class irtkGenericImageGradientFunction : public irtkImageGradientFunction
{
  irtkAbstractMacro(irtkGenericImageGradientFunction);

public:

  // ---------------------------------------------------------------------------
  // Types

  typedef TImage                                         ImageType;
  typedef typename ImageType::VoxelType                  VoxelType;
  typedef GradientType::ComponentType                    Real;
  typedef irtkGenericExtrapolateImageFunction<ImageType> ExtrapolatorType;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Default constructor
  irtkGenericImageGradientFunction();

public:

  /// Destructor
  virtual ~irtkGenericImageGradientFunction();

  /// Construct interpolator with default infinite extension of input image
  static irtkGenericImageGradientFunction *New(irtkInterpolationMode,
                                               const TImage * = NULL);

  /// Construct extrapolator which is compatible with this interpolator
  virtual irtkExtrapolateImageFunction *New(irtkExtrapolationMode,
                                            const irtkBaseImage * = NULL);

  /// Construct interpolator with specified infinite extension of input image
  ///
  /// The caller is required to set the input, initialize, and destroy the
  /// interpolator only, the extrapolator is initialized and destroyed by the
  /// interpolator unless the extrapolator has been replaced using the setter.
  static irtkGenericImageGradientFunction *New(irtkInterpolationMode,
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

};

////////////////////////////////////////////////////////////////////////////////
// Auxiliary macro for subclass implementation
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
#define irtkGradientInterpolatorMacro(clsname, mode)                           \
    irtkObjectMacro(clsname);                                                  \
  public:                                                                      \
    /** Get interpolation mode implemented by this image gradient function */  \
    inline virtual irtkInterpolationMode InterpolationMode() const             \
    { return mode; }                                                           \
    /** Get interpolation mode implemented by this class */                    \
    inline static  irtkInterpolationMode InterpolationType()                   \
    { return mode; }                                                           \
  private:                                                                     \
    static void _irtkGradientInterpolatorMacro_needs_trailing_semicolon()

// -----------------------------------------------------------------------------
#define irtkGenericGradientInterpolatorTypes(superclsname)                     \
  public:                                                                      \
    typedef superclsname<TImage>                           Superclass;         \
    typedef typename Superclass::ImageType                 ImageType;          \
    typedef typename Superclass::VoxelType                 VoxelType;          \
    typedef typename Superclass::Real                      Real;               \
    typedef typename Superclass::ExtrapolatorType          ExtrapolatorType;   \
    typedef typename Superclass::GradientType              GradientType;       \
  private:                                                                     \
    static void _irtkGenericGradientInterpolatorTypedefsMacro_needs_trailing_semicolon()

// -----------------------------------------------------------------------------
#define irtkGenericGradientInterpolatorMacro(clsname, mode)                    \
    irtkGradientInterpolatorMacro(clsname, mode);                              \
    irtkGenericGradientInterpolatorTypes(irtkGenericImageGradientFunction);    \
  private:                                                                     \
     static void _irtkGenericGradientInterpolatorMacro_needs_trailing_semicolon()

////////////////////////////////////////////////////////////////////////////////
// Inline definitions -- irtkImageGradientFunction
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Auxiliary functions
// =============================================================================

// -----------------------------------------------------------------------------
inline void irtkImageGradientFunction::ImageGradientToWorld(GradientType &g) const
{
  const GradientType v = g / _VoxelSize;
  g._x = v._x * _Orientation(0, 0) + v._y * _Orientation(1, 0) + v._z * _Orientation(2, 0);
  g._y = v._x * _Orientation(0, 1) + v._y * _Orientation(1, 1) + v._z * _Orientation(2, 1);
  g._z = v._x * _Orientation(0, 2) + v._y * _Orientation(1, 2) + v._z * _Orientation(2, 2);
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
inline void irtkImageGradientFunction::Input(const irtkBaseImage *input)
{
  this->_Input = const_cast<irtkBaseImage *>(input);
}

// -----------------------------------------------------------------------------
inline const irtkBaseImage *irtkImageGradientFunction::Input() const
{
  return this->_Input;
}

// -----------------------------------------------------------------------------
inline void irtkImageGradientFunction
::Extrapolator(irtkExtrapolateImageFunction *input, bool owner)
{
  if (_InfiniteInputOwner) Delete(_InfiniteInput);
  _InfiniteInput      = input;
  _InfiniteInputOwner = _InfiniteInput && owner;
}

// -----------------------------------------------------------------------------
inline irtkExtrapolateImageFunction *irtkImageGradientFunction::Extrapolator()
{
  return _InfiniteInput;
}

// -----------------------------------------------------------------------------
inline const irtkExtrapolateImageFunction *irtkImageGradientFunction::Extrapolator() const
{
  return _InfiniteInput;
}

// -----------------------------------------------------------------------------
inline irtkExtrapolationMode irtkImageGradientFunction::ExtrapolationMode() const
{
  return (_InfiniteInput ? _InfiniteInput->ExtrapolationMode() : Extrapolation_None);
}

// =============================================================================
// Domain checks
// =============================================================================

// ----------------------------------------------------------------------------
inline double irtkImageGradientFunction::GetTOrigin() const
{
  return this->_Input->GetTOrigin();
}

// ----------------------------------------------------------------------------
inline void irtkImageGradientFunction::WorldToImage(double &x, double &y) const
{
  this->_Input->WorldToImage(x, y);
}

// ----------------------------------------------------------------------------
inline void irtkImageGradientFunction::WorldToImage(double &x, double &y, double &z) const
{
  this->_Input->WorldToImage(x, y, z);
}

// ----------------------------------------------------------------------------
inline void irtkImageGradientFunction::WorldToImage(irtkPoint &p) const
{
  this->_Input->WorldToImage(p);
}

// -----------------------------------------------------------------------------
inline void irtkImageGradientFunction::Inside(double &x1, double &y1,
                                              double &x2, double &y2) const
{
  x1 = _x1, y1 = _y1;
  x2 = _x2, y2 = _y2;
}

// -----------------------------------------------------------------------------
inline void irtkImageGradientFunction::Inside(double &x1, double &y1, double &z1,
                                              double &x2, double &y2, double &z2) const
{
  x1 = _x1, y1 = _y1, z1 = _z1;
  x2 = _x2, y2 = _y2, z2 = _z2;
}

// -----------------------------------------------------------------------------
inline void irtkImageGradientFunction::Inside(double &x1, double &y1, double &z1, double &t1,
                                              double &x2, double &y2, double &z2, double &t2) const
{
  x1 = _x1, y1 = _y1, z1 = _z1, t1 = _t1;
  x2 = _x2, y2 = _y2, z2 = _z2, t2 = _t2;
}

// -----------------------------------------------------------------------------
inline void irtkImageGradientFunction::BoundingBox(double x, double y,
                                                   int &i1, int &j1,
                                                   int &i2, int &j2) const
{
  this->BoundingInterval(x, i1, i2);
  this->BoundingInterval(y, j1, j2);
}

// -----------------------------------------------------------------------------
inline void irtkImageGradientFunction::BoundingBox(double x, double y, double z,
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
inline void irtkImageGradientFunction::BoundingBox(double x, double y, double z, double t,
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
inline bool irtkImageGradientFunction::IsInside(double x, double y) const
{
  return (_x1 <= x && x <= _x2) && (_y1 <= y && y <= _y2);
}

// -----------------------------------------------------------------------------
inline bool irtkImageGradientFunction::IsInside(double x, double y, double z) const
{
  return (_x1 <= x && x <= _x2) && (_y1 <= y && y <= _y2) && (_z1 <= z && z <= _z2);
}

// -----------------------------------------------------------------------------
inline bool irtkImageGradientFunction::IsInside(double x, double y, double z, double t) const
{
  return (_x1 <= x && x <= _x2) && (_y1 <= y && y <= _y2) && (_z1 <= z && z <= _z2) && (_t1 <= t && t <= _t2);
}

// -----------------------------------------------------------------------------
inline bool irtkImageGradientFunction::IsOutside(double x, double y) const
{
  return !IsInside(x, y);
}

// -----------------------------------------------------------------------------
inline bool irtkImageGradientFunction::IsOutside(double x, double y, double z) const
{
  return !IsInside(x, y, z);
}

// -----------------------------------------------------------------------------
inline bool irtkImageGradientFunction::IsOutside(double x, double y, double z, double t) const
{
  return !IsInside(x, y, z, t);
}

// -----------------------------------------------------------------------------
inline bool irtkImageGradientFunction::IsForeground(double x, double y) const
{
  int i1, j1, i2, j2;
  BoundingBox(x, y,                             i1, j1, i2, j2);
  return Input()->IsBoundingBoxInsideForeground(i1, j1, i2, j2);
}

// -----------------------------------------------------------------------------
inline bool irtkImageGradientFunction::IsForeground(double x, double y, double z) const
{
  int i1, j1, k1, i2, j2, k2;
  BoundingBox(x, y, z,                          i1, j1, k1, i2, j2, k2);
  return Input()->IsBoundingBoxInsideForeground(i1, j1, k1, i2, j2, k2);
}

// -----------------------------------------------------------------------------
inline bool irtkImageGradientFunction::IsForeground(double x, double y, double z, double t) const
{
  int i1, j1, k1, l1, i2, j2, k2, l2;
  BoundingBox(x, y, z, t,                       i1, j1, k1, l1, i2, j2, k2, l2);
  return Input()->IsBoundingBoxInsideForeground(i1, j1, k1, l1, i2, j2, k2, l2);
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
inline irtkImageGradientFunction::GradientType
irtkImageGradientFunction::Get(double x, double y, double z, double t) const
{
  if (IsInside(x, y, z, t)) return this->GetInside (x, y, z, t);
  else                      return this->GetOutside(x, y, z, t);
}

// -----------------------------------------------------------------------------
inline irtkImageGradientFunction::GradientType
irtkImageGradientFunction::GetWithPadding(double x, double y, double z, double t) const
{
  if (IsInside(x, y, z, t)) return this->GetWithPaddingInside (x, y, z, t);
  else                      return this->GetWithPaddingOutside(x, y, z, t);
}

// -----------------------------------------------------------------------------
inline irtkImageGradientFunction::GradientType
irtkImageGradientFunction::operator ()(double x, double y, double z, double t) const
{
  return this->Get(x, y, z, t);
}

// -----------------------------------------------------------------------------
inline void
irtkImageGradientFunction
::EvaluateInside(double *v, double x, double y, double z, int vt) const
{
  GradientType g = this->GetInside(x, y, z, .0);
  if (_WrtWorld) ImageGradientToWorld(g);
  for (int i = 0; i < GradientType::Rows(); ++i, v += vt) *v = g(i);
}

// -----------------------------------------------------------------------------
inline void
irtkImageGradientFunction
::EvaluateOutside(double *v, double x, double y, double z, int vt) const
{
  GradientType g = this->GetOutside(x, y, z, .0);
  if (_WrtWorld) ImageGradientToWorld(g);
  for (int i = 0; i < GradientType::Rows(); ++i, v += vt) *v = g(i);
}

// -----------------------------------------------------------------------------
inline void
irtkImageGradientFunction
::Evaluate(double *v, double x, double y, double z, int vt) const
{
  if (IsInside(x, y, z)) return this->EvaluateInside (v, x, y, z, vt);
  else                   return this->EvaluateOutside(v, x, y, z, vt);
}

// -----------------------------------------------------------------------------
inline void
irtkImageGradientFunction
::EvaluateWithPaddingInside(double *v, double x, double y, double z, int vt) const
{
  GradientType g = this->GetWithPaddingInside(x, y, z, .0);
  if (_WrtWorld) ImageGradientToWorld(g);
  for (int i = 0; i < GradientType::Rows(); ++i, v += vt) *v = g(i);
}

// -----------------------------------------------------------------------------
inline void
irtkImageGradientFunction
::EvaluateWithPaddingOutside(double *v, double x, double y, double z, int vt) const
{
  GradientType g = this->GetWithPaddingOutside(x, y, z, .0);
  if (_WrtWorld) ImageGradientToWorld(g);
  for (int i = 0; i < GradientType::Rows(); ++i, v += vt) *v = g(i);
}

// -----------------------------------------------------------------------------
inline void
irtkImageGradientFunction
::EvaluateWithPadding(double *v, double x, double y, double z, int vt) const
{
  if (IsInside(x, y, z)) return this->EvaluateWithPaddingInside (v, x, y, z, vt);
  else                   return this->EvaluateWithPaddingOutside(v, x, y, z, vt);
}

////////////////////////////////////////////////////////////////////////////////
// Inline definitions -- irtkGenericImageGradientFunction
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
template <class TImage>
inline const TImage *irtkGenericImageGradientFunction<TImage>::Input() const
{
  return reinterpret_cast<const TImage *>(this->_Input);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename irtkGenericImageGradientFunction<TImage>::ExtrapolatorType *
irtkGenericImageGradientFunction<TImage>::Extrapolator()
{
  return reinterpret_cast<ExtrapolatorType *>(_InfiniteInput);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline const typename irtkGenericImageGradientFunction<TImage>::ExtrapolatorType *
irtkGenericImageGradientFunction<TImage>::Extrapolator() const
{
  return reinterpret_cast<const ExtrapolatorType *>(_InfiniteInput);
}

////////////////////////////////////////////////////////////////////////////////
// Image gradient functions
////////////////////////////////////////////////////////////////////////////////

// ND
#include <irtkLinearImageGradientFunction.h>
#include <irtkFastLinearImageGradientFunction.h>

// 2D
#include <irtkLinearImageGradientFunction2D.h>
#include <irtkFastLinearImageGradientFunction2D.h>

// 3D
#include <irtkLinearImageGradientFunction3D.h>
#include <irtkFastLinearImageGradientFunction3D.h>


#endif
