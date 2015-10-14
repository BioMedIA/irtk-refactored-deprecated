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

#ifndef _IRTKREGISTRATIONFILTER_H
#define _IRTKREGISTRATIONFILTER_H

#include <irtkCommon.h>
#include <irtkTransformation.h>


////////////////////////////////////////////////////////////////////////////////
// Interface
////////////////////////////////////////////////////////////////////////////////

/**
 * Base class for registration filters
 *
 * \code
 * irtkMyRegistrationFilter registration;
 * irtkTransformation      *transformation = NULL;
 * registration.Input(...);
 * registration.Output(&transformation);
 * registration.Read(config_file);
 * registration.Run();
 * \endcode
 */
class irtkRegistrationFilter : public irtkObservable
{
  irtkAbstractMacro(irtkRegistrationFilter);

private:

  /// Pointer to output transformation
  irtkTransformation **_Output;

private:

  /// Copy constructor
  /// \note Intentionally not implemented
  irtkRegistrationFilter(const irtkRegistrationFilter &);

  /// Assignment operator
  /// \note Intentionally not implemented
  void operator =(const irtkRegistrationFilter &);

protected:

  /// Constructor
  irtkRegistrationFilter();

public:

  /// Destructor
  virtual ~irtkRegistrationFilter() = 0;

  /// Read registration parameters from input stream
  virtual bool Read(const char *, bool = false);

  /// Read registration parameters from input stream
  virtual bool Read(istream &, bool = false) = 0;

  /// Write registration parameters to file
  virtual void Write(const char *) const = 0;

  /// Runs the registration filter
  virtual void Run() = 0;

  /// Set pointer to output transformation
  void Output(irtkTransformation **);

  /// Get output transformation
  irtkTransformation *Output();

protected:

  /// Set current output transformation
  bool Output(irtkTransformation *);

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline irtkRegistrationFilter::irtkRegistrationFilter()
{
}

// -----------------------------------------------------------------------------
inline irtkRegistrationFilter::~irtkRegistrationFilter()
{
}

// -----------------------------------------------------------------------------
inline bool irtkRegistrationFilter::Read(const char *fname, bool echo)
{
  ifstream from(fname);
  if (!from) return false;
  bool ok = this->Read(from, echo);
  from.close();
  return ok;
}

// -----------------------------------------------------------------------------
inline void irtkRegistrationFilter::Output(irtkTransformation **output)
{
  _Output = output;
}

// -----------------------------------------------------------------------------
inline bool irtkRegistrationFilter::Output(irtkTransformation *output)
{
  if (_Output) {
    (*_Output) = output;
    return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
inline irtkTransformation *irtkRegistrationFilter::Output()
{
  return _Output ? *_Output : NULL;
}

////////////////////////////////////////////////////////////////////////////////
// Constants
////////////////////////////////////////////////////////////////////////////////

#ifdef MAX_NO_RESOLUTIONS
#  undef MAX_NO_RESOLUTIONS
#endif
const int MAX_NO_RESOLUTIONS = 10;

////////////////////////////////////////////////////////////////////////////////
// Transformation models
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
/// Enumeration of transformation models
///
/// A transformation model is implemented by one or more transformation
/// classes and futhermore may be used to determine additional registration
/// settings, such as hard and soft transformation constraints. Thus, the
/// transformation model differs in semantics from the transformation type.
///
/// \sa irtkTransformationType, ToTransformationType
enum irtkTransformationModel
{
  TM_Unknown,                       ///< Unknown/invalid transformation model
  // Add new enumeration values below
  TM_Rigid,                         ///< Linear transformation with up to  6 DoFs (rotate, translate)
  TM_Similarity,                    ///< Linear transformation with up to  7 DoFs (rotate, translate, global scale)
  TM_Affine,                        ///< Linear transformation with up to 12 DoFs (rotate, translate, scale, skew)
  TM_LinearFFD,                     ///< Displacement field with linear interpolation
  TM_BSplineFFD,                    ///< Displacement field with B-spline interpolation
  TM_BSplineStatFFD,                ///< Displacement field with B-spline interpolation using a statistical model
  TM_BSplineSVFFD,                  ///< Stationary velocity field with B-spline interpolation
  TM_BSplineTDFFD,                  ///< Non-stationary velocity field with B-spline interpolation
  // Add new enumeration values above
  TM_Last                           ///< Number of available transformation models + 1
};

// -----------------------------------------------------------------------------
inline string ToString(const irtkTransformationModel &m)
{
  switch (m) {
    case TM_Rigid:           return "Rigid";
    case TM_Similarity:      return "Similarity";
    case TM_Affine:          return "Affine";
    case TM_LinearFFD:       return "LinearFFD";
    case TM_BSplineFFD:      return "BSplineFFD";
    case TM_BSplineStatFFD:  return "BSplineStatFFD";
    case TM_BSplineSVFFD:    return "BSplineSVFFD";
    case TM_BSplineTDFFD:    return "BSplineTDFFD";
    default:                 return "Unknown";
  }
}

// -----------------------------------------------------------------------------
inline string ToPrettyString(const irtkTransformationModel &m)
{
  switch (m) {
    case TM_Rigid:           return "rigid transformation";
    case TM_Similarity:      return "similarity transformation";
    case TM_Affine:          return "affine transformation";
    case TM_LinearFFD:       return "non-parametric displacement field";
    case TM_BSplineFFD:      return "free-form deformation";
    case TM_BSplineStatFFD:  return "statistical free-form deformation";
    case TM_BSplineSVFFD:    return "parametric stationary velocity field transformation";
    case TM_BSplineTDFFD:    return "temporal diffeomorphic free-form deformation";
    default:                 return "unknown transformation";
  }
}

// -----------------------------------------------------------------------------
inline bool FromString(const char *str, irtkTransformationModel &m)
{
  m = TM_Unknown;
  // Alternative names of transformation models
  if      (strcmp(str, "FFD")   == 0) m = TM_BSplineFFD;
  else if (strcmp(str, "SVFFD") == 0) m = TM_BSplineSVFFD;
  else if (strcmp(str, "TDFFD") == 0) m = TM_BSplineTDFFD;
  // Default names of transformation models
  if (m == TM_Unknown) {
    m = static_cast<irtkTransformationModel>(TM_Last - 1);
    while (m != TM_Unknown) {
      if (ToString(m) == str) break;
      m = static_cast<irtkTransformationModel>(m - 1);
    }
  }
  return (m != TM_Unknown);
}

// -----------------------------------------------------------------------------
/// Whether a given transformation model is linear
inline bool IsLinear(irtkTransformationModel model)
{
  return model == TM_Rigid || model == TM_Similarity || model == TM_Affine;
}

// -----------------------------------------------------------------------------
/// Whether any given transformation model is linear
inline bool IsLinear(const vector<irtkTransformationModel> &model)
{
  for (size_t i = 0; i < model.size(); ++i) {
    if (IsLinear(model[i])) return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
/// Whether any given transformation model is non-linear
inline bool IsNonLinear(const vector<irtkTransformationModel> &model)
{
  for (size_t i = 0; i < model.size(); ++i) {
    if (!IsLinear(model[i])) return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
/// Whether a given transformation model is a FFD with a linear interpolation kernel
inline bool IsLinearFFD(irtkTransformationModel model)
{
  return model == TM_LinearFFD;
}

// -----------------------------------------------------------------------------
/// Whether any given transformation model is a FFD with a linear interpolation kernel
inline bool IsLinearFFD(const vector<irtkTransformationModel> &model)
{
  for (size_t i = 0; i < model.size(); ++i) {
    if (IsLinearFFD(model[i])) return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
/// Whether a given transformation model is diffeomorphic (velocity parameterization)
inline bool IsDiffeo(irtkTransformationModel model)
{
  return IsLinear(model) || model == TM_BSplineSVFFD || model == TM_BSplineTDFFD;
}

// -----------------------------------------------------------------------------
/// Whether any given transformation model is a diffeomorphic model
inline bool IsDiffeo(const vector<irtkTransformationModel> &model)
{
  for (size_t i = 0; i < model.size(); ++i) {
    if (IsDiffeo(model[i])) return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
/// Whether a given transformation model is 3D+t
inline bool IsSpatioTemporal(irtkTransformationModel model)
{
  return model == TM_BSplineTDFFD;
}

// -----------------------------------------------------------------------------
/// Whether any given transformation model is 3D+t
inline bool IsSpatioTemporal(const vector<irtkTransformationModel> &model)
{
  for (size_t i = 0; i < model.size(); ++i) {
    if (IsSpatioTemporal(model[i])) return true;
  }
  return false;
}

// -----------------------------------------------------------------------------
/// Get type of (default) transformation which implements a specific model
///
/// The mapping from transformation model to transformtion type is not
/// one-to-one. More then one transformation can be suitable for a
/// transformation model. This function defines the default type used
/// for each model. The base registration filter implementations make use
/// of it, but a specialized registration filter can choose another
/// transformation for a given model if desired.
///
/// For example, see irtkImageRegistrationFilter::TransformationType.
inline irtkTransformationType
ToTransformationType(irtkTransformationModel    model,
                     const irtkImageAttributes &domain)
{
  // 2D/3D
  if (domain._t == 1) {
    switch (model) {
      case TM_LinearFFD:     return IRTKTRANSFORMATION_LINEAR_FFD_3D_v3;
      case TM_BSplineFFD:    return IRTKTRANSFORMATION_BSPLINE_FFD_3D_v3;
      default: break;
    }
  // 4D
  } else {
    switch (model) {
      case TM_LinearFFD:     return IRTKTRANSFORMATION_LINEAR_FFD_4D_v2;
      case TM_BSplineFFD:    return IRTKTRANSFORMATION_BSPLINE_FFD_4D_v2;
      default: break;
    }
  }
  // nD
  switch (model) {
    case TM_Rigid:           return IRTKTRANSFORMATION_RIGID;
    case TM_Similarity:      return IRTKTRANSFORMATION_SIMILARITY;
    case TM_Affine:          return IRTKTRANSFORMATION_AFFINE;
    case TM_BSplineStatFFD:  return IRTKTRANSFORMATION_BSPLINE_FFD_STATISTICAL;
    case TM_BSplineSVFFD:    return IRTKTRANSFORMATION_BSPLINE_FFD_SV_v5;
    case TM_BSplineTDFFD:    return IRTKTRANSFORMATION_BSPLINE_FFD_TD_v3;
    default:                 return IRTKTRANSFORMATION_UNKNOWN;
  }
}

// -----------------------------------------------------------------------------
/// Enumeration of available multi-level transformation modes
enum irtkMFFDMode
{
  MFFD_Default,  ///< Choose suitable default multi-level transformation model
  MFFD_None,     ///< Use single transformation without additional global or local transformations
  MFFD_Sum,      ///< One transformation for each resolution level with additive composition
  MFFD_Fluid,    ///< One transformation for each resolution level with fluid composition
  MFFD_LogSum    ///< Additive multi-level stationary velocity field
};

// -----------------------------------------------------------------------------
inline string ToString(const irtkMFFDMode &m)
{
  switch (m) {
    case MFFD_Default:  return "Default";
    case MFFD_None:     return "None";
    case MFFD_Sum:      return "Sum";
    case MFFD_LogSum:   return "LogSum";
    case MFFD_Fluid:    return "Fluid";
    default:            return "Unknown";
  }
}

// -----------------------------------------------------------------------------
inline bool FromString(const char *str, irtkMFFDMode &m)
{
  if      (strcmp(str, "Default")  == 0) m = MFFD_Default;
  else if (strcmp(str, "None")     == 0) m = MFFD_None;
  else if (strcmp(str, "Sum")      == 0) m = MFFD_Sum;
  else if (strcmp(str, "LogSum")   == 0) m = MFFD_LogSum;
  else if (strcmp(str, "Fluid")    == 0) m = MFFD_Fluid;
  else return false;
  return true;
}


#endif
