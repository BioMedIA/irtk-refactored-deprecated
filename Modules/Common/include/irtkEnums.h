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

#ifndef _IRTKENUMS_H
#define _IRTKENUMS_H


// =============================================================================
// Orientation codes (same as NIFTI)
// =============================================================================

#define IRTK_L2R  1    // Left to Right
#define IRTK_R2L  2    // Right to Left
#define IRTK_P2A  3    // Posterior to Anterior
#define IRTK_A2P  4    // Anterior to Posterior
#define IRTK_I2S  5    // Inferior to Superior
#define IRTK_S2I  6    // Superior to Inferior

// =============================================================================
// Inter-/Extrapolation
// =============================================================================

// ----------------------------------------------------------------------------
/// Image interpolation modes
enum irtkInterpolationMode {
  Interpolation_NN,
  Interpolation_Linear,
  Interpolation_FastLinear,
  Interpolation_BSpline,
  Interpolation_CubicBSpline,
  Interpolation_FastCubicBSpline,
  Interpolation_CSpline,
  Interpolation_SBased,
  Interpolation_Sinc,
  Interpolation_Gaussian,
  Interpolation_NNWithPadding,
  Interpolation_LinearWithPadding,
  Interpolation_FastLinearWithPadding,
  Interpolation_BSplineWithPadding,
  Interpolation_CubicBSplineWithPadding,
  Interpolation_FastCubicBSplineWithPadding,
  Interpolation_CSplineWithPadding,
  Interpolation_SBasedWithPadding,
  Interpolation_SincWithPadding,
  Interpolation_GaussianWithPadding,
  // Add new enumeration values above
  Interpolation_Last
};

// ----------------------------------------------------------------------------
inline string ToString(const irtkInterpolationMode &m)
{
  switch(m) {
    case Interpolation_NN:                          return "NN";
    case Interpolation_Linear:                      return "Linear";
    case Interpolation_FastLinear:                  return "Fast linear";
    case Interpolation_BSpline:                     return "BSpline";
    case Interpolation_CSpline:                     return "CSpline";
    case Interpolation_CubicBSpline:                return "Cubic BSpline";
    case Interpolation_FastCubicBSpline:            return "Fast cubic BSpline";
    case Interpolation_SBased:                      return "SBased";
    case Interpolation_Sinc:                        return "Sinc";
    case Interpolation_Gaussian:                    return "Gaussian";
    case Interpolation_NNWithPadding:               return "NN with padding";
    case Interpolation_LinearWithPadding:           return "Linear with padding";
    case Interpolation_FastLinearWithPadding:       return "Fast linear with padding";
    case Interpolation_BSplineWithPadding:          return "BSpline with padding";
    case Interpolation_CubicBSplineWithPadding:     return "Cubic BSpline with padding";
    case Interpolation_FastCubicBSplineWithPadding: return "Fast cubic BSpline with padding";
    case Interpolation_CSplineWithPadding:          return "CSpline with padding";
    case Interpolation_SBasedWithPadding:           return "SBased with padding";
    case Interpolation_SincWithPadding:             return "Sinc with padding";
    case Interpolation_GaussianWithPadding:         return "Gaussian with padding";
    default:                                        return "Unknown";
  }
}

// ----------------------------------------------------------------------------
inline bool FromString(const char *str, irtkInterpolationMode &m)
{
  if      (strcmp(str, "NN") == 0) m = Interpolation_NN;
  else if (strcmp(str, "Linear") == 0) m = Interpolation_Linear;
  else if (strcmp(str, "Fast linear") == 0) m = Interpolation_FastLinear;
  else if (strcmp(str, "BSpline") == 0) m = Interpolation_BSpline;
  else if (strcmp(str, "CSpline") == 0) m = Interpolation_CSpline;
  else if (strcmp(str, "Cubic BSpline") == 0) m = Interpolation_CubicBSpline;
  else if (strcmp(str, "Fast cubic BSpline") == 0) m = Interpolation_FastCubicBSpline;
  else if (strcmp(str, "SBased") == 0) m = Interpolation_SBased;
  else if (strcmp(str, "Sinc") == 0) m = Interpolation_Sinc;
  else if (strcmp(str, "Gaussian") == 0) m = Interpolation_Gaussian;
  else if (strcmp(str, "NN with padding") == 0) m = Interpolation_NNWithPadding;
  else if (strcmp(str, "NN (with padding)") == 0) m = Interpolation_NNWithPadding;
  else if (strcmp(str, "Linear with padding") == 0) m = Interpolation_LinearWithPadding;
  else if (strcmp(str, "Linear (with padding)") == 0) m = Interpolation_LinearWithPadding;
  else if (strcmp(str, "Fast linear with padding") == 0) m = Interpolation_FastLinearWithPadding;
  else if (strcmp(str, "Fast linear (with padding)") == 0) m = Interpolation_FastLinearWithPadding;
  else if (strcmp(str, "BSpline with padding") == 0) m = Interpolation_BSplineWithPadding;
  else if (strcmp(str, "BSpline (with padding)") == 0) m = Interpolation_BSplineWithPadding;
  else if (strcmp(str, "Cubic BSpline with padding") == 0) m = Interpolation_CubicBSplineWithPadding;
  else if (strcmp(str, "Cubic BSpline (with padding)") == 0) m = Interpolation_CubicBSplineWithPadding;
  else if (strcmp(str, "Fast cubic BSpline with padding") == 0) m = Interpolation_FastCubicBSplineWithPadding;
  else if (strcmp(str, "Fast cubic BSpline (with padding)") == 0) m = Interpolation_FastCubicBSplineWithPadding;
  else if (strcmp(str, "CSpline with padding") == 0) m = Interpolation_CSplineWithPadding;
  else if (strcmp(str, "CSpline (with padding)") == 0) m = Interpolation_CSplineWithPadding;
  else if (strcmp(str, "SBased with padding") == 0) m = Interpolation_SBasedWithPadding;
  else if (strcmp(str, "SBased (with padding)") == 0) m = Interpolation_SBasedWithPadding;
  else if (strcmp(str, "Sinc with padding") == 0) m = Interpolation_SincWithPadding;
  else if (strcmp(str, "Sinc (with padding)") == 0) m = Interpolation_SincWithPadding;
  else if (strcmp(str, "Gaussian with padding") == 0) m = Interpolation_GaussianWithPadding;
  else if (strcmp(str, "Gaussian (with padding)") == 0) m = Interpolation_GaussianWithPadding;
  else return false;
  return true;
}

// ----------------------------------------------------------------------------
/// Get corresponding interpolation with padding
inline irtkInterpolationMode InterpolationWithPadding(irtkInterpolationMode m)
{
  switch(m) {
    case Interpolation_NN:               return Interpolation_NNWithPadding;
    case Interpolation_Linear:           return Interpolation_LinearWithPadding;
    case Interpolation_FastLinear:       return Interpolation_FastLinearWithPadding;
    case Interpolation_BSpline:          return Interpolation_BSplineWithPadding;
    case Interpolation_CubicBSpline:     return Interpolation_CubicBSplineWithPadding;
    case Interpolation_FastCubicBSpline: return Interpolation_FastCubicBSplineWithPadding;
    case Interpolation_CSpline:          return Interpolation_CSplineWithPadding;
    case Interpolation_SBased:           return Interpolation_SBasedWithPadding;
    case Interpolation_Sinc:             return Interpolation_SincWithPadding;
    case Interpolation_Gaussian:         return Interpolation_GaussianWithPadding;
    default:                             return m;
  }
}

// ----------------------------------------------------------------------------
/// Get corresponding interpolation without padding
inline irtkInterpolationMode InterpolationWithoutPadding(irtkInterpolationMode m)
{
  switch(m) {
    case Interpolation_NNWithPadding:               return Interpolation_NN;
    case Interpolation_LinearWithPadding:           return Interpolation_Linear;
    case Interpolation_FastLinearWithPadding:       return Interpolation_FastLinear;
    case Interpolation_BSplineWithPadding:          return Interpolation_BSpline;
    case Interpolation_CubicBSplineWithPadding:     return Interpolation_CubicBSpline;
    case Interpolation_FastCubicBSplineWithPadding: return Interpolation_FastCubicBSpline;
    case Interpolation_CSplineWithPadding:          return Interpolation_CSpline;
    case Interpolation_SBasedWithPadding:           return Interpolation_SBased;
    case Interpolation_SincWithPadding:             return Interpolation_Sinc;
    case Interpolation_GaussianWithPadding:         return Interpolation_Gaussian;
    default:                                        return m;
  }
}

// ----------------------------------------------------------------------------
/// Image extrapolation modes
enum irtkExtrapolationMode {
  Extrapolation_Default, // i.e., not specified
  Extrapolation_None,    // i.e., partial interpolation or default value
  Extrapolation_Const,
  Extrapolation_NN,
  Extrapolation_Repeat,  // i.e., periodic
  Extrapolation_Mirror,
  Extrapolation_ConstWithPeriodicTime,
  // Add new enumeration values above
  Extrapolation_Last
};

// ----------------------------------------------------------------------------
inline string ToString(const irtkExtrapolationMode &m)
{
  switch(m) {
    case Extrapolation_Default:               return "Default";
    case Extrapolation_None:                  return "None";
    case Extrapolation_Const:                 return "Const";
    case Extrapolation_NN:                    return "NN";
    case Extrapolation_Repeat:                return "Repeat";
    case Extrapolation_Mirror:                return "Mirror";
    case Extrapolation_ConstWithPeriodicTime: return "ConstWithPeriodicTime";
    default:                                  return "Unknown";
  }
}

// ----------------------------------------------------------------------------
/// Get corresponding extrapolation with periodic time
inline irtkExtrapolationMode ExtrapolationWithPeriodicTime(irtkExtrapolationMode m)
{
  switch(m) {
    case Extrapolation_ConstWithPeriodicTime:
    case Extrapolation_Const:  return Extrapolation_ConstWithPeriodicTime;
    case Extrapolation_Repeat: return Extrapolation_Repeat;
    default:                   return Extrapolation_None;
  }
}

// ----------------------------------------------------------------------------
/// Get corresponding extrapolation without periodic time
inline irtkExtrapolationMode ExtrapolationWithoutPeriodicTime(irtkExtrapolationMode m)
{
  switch(m) {
    case Extrapolation_Const:
    case Extrapolation_ConstWithPeriodicTime: return Extrapolation_Const;
    // Note: Extrapolation_Repeat remains periodic in all dimensions
    default: return m;
  }
}

// ----------------------------------------------------------------------------
inline bool FromString(const char *str, irtkExtrapolationMode &m)
{
  if      (strcmp(str, "Default")               == 0) m = Extrapolation_Default;
  else if (strcmp(str, "None")                  == 0) m = Extrapolation_None;
  else if (strcmp(str, "Const")                 == 0) m = Extrapolation_Const;
  else if (strcmp(str, "NN")                    == 0) m = Extrapolation_NN;
  else if (strcmp(str, "Repeat")                == 0) m = Extrapolation_Repeat;
  else if (strcmp(str, "Mirror")                == 0) m = Extrapolation_Mirror;
  else if (strcmp(str, "ConstWithPeriodicTime") == 0) m = Extrapolation_ConstWithPeriodicTime;
  else return false;
  return true;
}

// =============================================================================
// Similarity measures
// =============================================================================

// -----------------------------------------------------------------------------
/// Enumeration of available similarity measures
enum irtkSimilarityMeasure
{
  JE,     ///< Joint entropy
  CC,     ///< Cross-correlation
  MI,     ///< Mutual information
  NMI,    ///< Normalized mutual information
  SSD,    ///< Sum of squared differences
  CR_XY,  ///< Correlation ratio
  CR_YX,  ///< Correlation ratio
  LC,
  K,
  ML,
  NGF_COS,  ///< Cosine of normalzed gradient field
  NCC    ///< Normalized/Local cross-correlation
};

// -----------------------------------------------------------------------------
inline string ToString(const irtkSimilarityMeasure &m)
{
  switch (m) {
    case JE:       return "JE";
    case CC:       return "CC";
    case MI:       return "MI";
    case NMI:      return "NMI";
    case SSD:      return "SSD";
    case CR_XY:    return "CR_XY";
    case CR_YX:    return "CR_YX";
    case LC:       return "LC";
    case K:        return "K";
    case ML:       return "ML";
    case NGF_COS:  return "NGF_COS";
    case NCC:      return "LNCC";
    default:       return "Unknown";
  }
}

// -----------------------------------------------------------------------------
inline bool FromString(const char *str, irtkSimilarityMeasure &m)
{
  if      (strcmp(str, "JE")       == 0) m = JE;
  else if (strcmp(str, "CC")       == 0) m = CC;
  else if (strcmp(str, "MI")       == 0) m = MI;
  else if (strcmp(str, "NMI")      == 0) m = NMI;
  else if (strcmp(str, "SSD")      == 0) m = SSD;
  else if (strcmp(str, "CR_XY")    == 0) m = CR_XY;
  else if (strcmp(str, "CR_YX")    == 0) m = CR_YX;
  else if (strcmp(str, "LC")       == 0) m = LC;
  else if (strcmp(str, "K")        == 0) m = K;
  else if (strcmp(str, "ML")       == 0) m = ML;
  else if (strcmp(str, "NGF_COS")  == 0) m = NGF_COS;
  else if (strcmp(str, "NCC")      == 0) m = NCC;
  else if (strcmp(str, "LCC")      == 0) m = NCC;
  else if (strcmp(str, "LNCC")     == 0) m = NCC;
  else return false;
  return true;
}

// =============================================================================
// Point set (dis-)similarity measures
// =============================================================================

// -----------------------------------------------------------------------------
/// Enumeration of available polydata distance measures
enum irtkPointSetDistanceMeasure
{
  PDM_Unknown,                   ///< Unknown/invalid point set distance measure
  // Add new enumeration values below
  PDM_FRE,                       ///< Fiducial registration error (FRE) measure
  PDM_CorrespondenceDistance,    ///< Point correspondence distance measure
  PDM_CurrentsDistance,          ///< Distance measure based on currents representation
  PDM_VarifoldDistance,          ///< Distance measure based on varifold representation
  // Add new enumeration values above
  PDM_Last                       ///< Number of available point set distance measures + 1
};

// -----------------------------------------------------------------------------
inline string ToString(const irtkPointSetDistanceMeasure &m)
{
  switch (m) {
    case PDM_FRE:                    return "FRE";
    case PDM_CorrespondenceDistance: return "PCD";
    case PDM_CurrentsDistance:       return "CurrentsDistance";
    case PDM_VarifoldDistance:       return "VarifoldDistance";
    default:                         return "Unknown";
  }
}

// -----------------------------------------------------------------------------
inline bool FromString(const char *str, irtkPointSetDistanceMeasure &m)
{
  m = PDM_Unknown;
  // Alternative names for point set distance measures
  if (strcmp(str, "Fiducial Registration Error") == 0 ||
      strcmp(str, "Fiducial registration error") == 0 ||
      strcmp(str, "Fiducial Error")              == 0 ||
      strcmp(str, "Fiducial error")              == 0 ||
      strcmp(str, "Landmark Registration Error") == 0 ||
      strcmp(str, "Landmark registration error") == 0 ||
      strcmp(str, "Landmark Error")              == 0 ||
      strcmp(str, "Landmark error")              == 0) {
    m = PDM_FRE;
  } else if (strcmp(str, "Point Correspondence Distance") == 0 ||
             strcmp(str, "Point correspondence distance") == 0 ||
             strcmp(str, "Correspondence Distance")       == 0 ||
             strcmp(str, "Correspondence distance")       == 0) {
    m = PDM_CorrespondenceDistance;
  } else if (strcmp(str, "Currents distance") == 0 ||
             strcmp(str, "Currents Distance") == 0) {
    m = PDM_CurrentsDistance;
  } else if (strcmp(str, "Varifold distance") == 0 ||
             strcmp(str, "Varifold Distance") == 0) {
    m = PDM_VarifoldDistance;
  }
  // Default names of point set distance measures
  if (m == PDM_Unknown) {
    m = static_cast<irtkPointSetDistanceMeasure>(PDM_Last - 1);
    while (m != PDM_Unknown) {
      if (ToString(m) == str) break;
      m = static_cast<irtkPointSetDistanceMeasure>(m - 1);
    }
  }
  return (m != PDM_Unknown);
}

// =============================================================================
// External simplicial complex forces (resp., "distance" of point set to image)
// =============================================================================

// -----------------------------------------------------------------------------
enum irtkPointSetForceMeasure
{
  PFM_Unknown,                 ///< Unknown/invalid external force
  // Add new enumeration values below
  PFM_BalloonForce,            ///< Balloon/inflation force
  PFM_EdgeForce,               ///< Image edge force
  // Add new enumeration values above
  PFM_Last                     ///< Number of available external forces + 1
};

// -----------------------------------------------------------------------------
inline string ToString(const irtkPointSetForceMeasure &m)
{
  switch (m) {
    case PFM_BalloonForce: return "BalloonForce";
    case PFM_EdgeForce:    return "EdgeForce";
    default:               return "Unknown";
  }
}

// -----------------------------------------------------------------------------
inline bool FromString(const char *str, irtkPointSetForceMeasure &m)
{
  m = PFM_Unknown;
  // Alternative names for constraint measures
  if      (strcmp(str, "InflationForce") == 0) m = PFM_BalloonForce;
  else if (strcmp(str, "ImageEdgeForce") == 0) m = PFM_EdgeForce;
  // Default names of constraint measures
  if (m == PFM_Unknown) {
    m = static_cast<irtkPointSetForceMeasure>(PFM_Last - 1);
    while (m != PFM_Unknown) {
      if (ToString(m) == str) break;
      m = static_cast<irtkPointSetForceMeasure>(m - 1);
    }
  }
  return (m != PFM_Unknown);
}

// =============================================================================
// Internal simplicial complex forces
// =============================================================================

// -----------------------------------------------------------------------------
enum irtkPointSetConstraintMeasure
{
  PCM_Unknown,                 ///< Unknown/invalid regularizer
  // Add new enumeration values below
  PCM_Stretching,              ///< Penalize changes of edge lengths (e.g., metric distortion)
  PCM_Curvature,               ///< Minimize curvature of point set surface
  PCM_NonSelfIntersection,     ///< Repels too close non-neighboring elements
  PCM_Inflation,               ///< Inflate point set surface
  // Add new enumeration values above
  PCM_Last                     ///< Number of available regularizers + 1
};

// -----------------------------------------------------------------------------
inline string ToString(const irtkPointSetConstraintMeasure &m)
{
  switch (m) {
    case PCM_Stretching:          return "Stretching";
    case PCM_Curvature:           return "Curvature";
    case PCM_NonSelfIntersection: return "NSI";
    case PCM_Inflation:           return "Inflation";
    default:                      return "Unknown";
  }
}

// -----------------------------------------------------------------------------
inline bool FromString(const char *str, irtkPointSetConstraintMeasure &m)
{
  m = PCM_Unknown;
  // Alternative names for constraint measures
  if      (strcmp(str, "EdgeLength")          == 0) m = PCM_Stretching;
  else if (strcmp(str, "Distortion")          == 0) m = PCM_Stretching;
  else if (strcmp(str, "MetricDistortion")    == 0) m = PCM_Stretching;
  else if (strcmp(str, "Bending")             == 0) m = PCM_Curvature;
  else if (strcmp(str, "SurfaceBending")      == 0) m = PCM_Curvature;
  else if (strcmp(str, "SurfaceCurvature")    == 0) m = PCM_Curvature;
  else if (strcmp(str, "SurfaceInflation")    == 0) m = PCM_Inflation;
  else if (strcmp(str, "NonSelfIntersection") == 0) m = PCM_NonSelfIntersection;
  // Default names of constraint measures
  if (m == PCM_Unknown) {
    m = static_cast<irtkPointSetConstraintMeasure>(PCM_Last - 1);
    while (m != PCM_Unknown) {
      if (ToString(m) == str) break;
      m = static_cast<irtkPointSetConstraintMeasure>(m - 1);
    }
  }
  return (m != PCM_Unknown);
}

// =============================================================================
// Transformation regularization
// =============================================================================

// -----------------------------------------------------------------------------
enum irtkConstraintMeasure
{
  CM_Unknown,                 ///< Unknown/invalid regularizer
  // Add new enumeration values below
  CM_Smoothness,              ///< Default smoothness constraint
  CM_VolumePreservation,      ///< Default volume preservation constraint
  CM_TopologyPreservation,    ///< Default topology preservation constraint
  CM_Sparsity,                ///< Default sparsity constraint
  CM_BendingEnergy,           ///< Thin-plate spline bending energy
  CM_L0Norm,                  ///< Sparsity constraint based on l0-norm
  CM_L1Norm,                  ///< Sparsity constraint based on l1-norm
  CM_L2Norm,                  ///< Sparsity constraint based on l2-norm
  CM_SqLogDetJac,             ///< Squared logarithm of the Jacobian determinant
  CM_MinDetJac,               ///< Constrain minimum Jacobian determinant
  // Add new enumeration values above
  CM_Last                     ///< Number of available regularizers + 1
};

// -----------------------------------------------------------------------------
inline string ToString(const irtkConstraintMeasure &m)
{
  switch (m) {
    case CM_BendingEnergy:        return "BE";
    case CM_Smoothness:           return "S";
    case CM_VolumePreservation:   return "VP";
    case CM_TopologyPreservation: return "TP";
    case CM_Sparsity:             return "Sparsity";
    case CM_L0Norm:               return "L0";
    case CM_L1Norm:               return "L1";
    case CM_L2Norm:               return "L2";
    case CM_SqLogDetJac:          return "SqLogDetJac";
    case CM_MinDetJac:            return "MinDetJac";
    default:                      return "Unknown";
  }
}

// -----------------------------------------------------------------------------
inline bool FromString(const char *str, irtkConstraintMeasure &m)
{
  m = CM_Unknown;
  // Alternative names for constraint measures
  if      (strcmp(str, "JAC")    == 0) m = CM_SqLogDetJac;
  else if (strcmp(str, "MinJac") == 0) m = CM_MinDetJac;
  // Default names of constraint measures
  if (m == CM_Unknown) {
    m = static_cast<irtkConstraintMeasure>(CM_Last - 1);
    while (m != CM_Unknown) {
      if (ToString(m) == str) break;
      m = static_cast<irtkConstraintMeasure>(m - 1);
    }
  }
  return (m != CM_Unknown);
}

// =============================================================================
// Optimization methods
// =============================================================================

// -----------------------------------------------------------------------------
/// Enumeration of available states for individual DoFs
enum DOFStatus { _Active  = 0, Active  = _Active,
                 _Passive = 1, Passive = _Passive,
                 _Unknown = 2 };

// -----------------------------------------------------------------------------
/// Enumeration of available optimization methods
enum irtkOptimizationMethod
{
  DownhillDescent,
  GradientDescent,
  GradientDescentConstrained,
  SteepestGradientDescent,
  ConjugateGradientDescent,
  ClosedForm,
  LBFGS,
  LineSearch,
  EulerMethod,            ///< Explicit Euler method for deformable surface models
  EulerMethodWithMomentum ///< Explicit Euler method with momentum for deformable surface models
};

// -----------------------------------------------------------------------------
inline string ToString(const irtkOptimizationMethod &m)
{
  switch (m) {
    case DownhillDescent:            return "DownhillDescent";
    case GradientDescent:            return "GradientDescent";
    case GradientDescentConstrained: return "GradientDescentConstrained";
    case SteepestGradientDescent:    return "SteepestGradientDescent";
    case ConjugateGradientDescent:   return "ConjugateGradientDescent";
    case ClosedForm:                 return "ClosedForm";
    case LBFGS:                      return "LBFGS";
    case LineSearch:                 return "LineSearch";
    case EulerMethod:                return "EulerMethod";
    case EulerMethodWithMomentum:    return "EulerMethodWithMomentum";
    default:                         return "Unknown";
  }
}

// -----------------------------------------------------------------------------
inline bool FromString(const char *str, irtkOptimizationMethod &m)
{
  if      (strcmp(str, "DownhillDescent")            == 0 || strcmp(str, "Downhill") == 0) m = DownhillDescent;
  else if (strcmp(str, "GradientDescent")            == 0) m = GradientDescent;
  else if (strcmp(str, "ConstrainedGradientDescent") == 0 ||
           strcmp(str, "GradientDescentConstrained") == 0) m = GradientDescentConstrained;
  else if (strcmp(str, "SteepestGradientDescent")    == 0 || strcmp(str, "SteepestGradient")  == 0 || strcmp(str, "SGD") == 0) m = SteepestGradientDescent;
  else if (strcmp(str, "ConjugateGradientDescent")   == 0 || strcmp(str, "ConjugateGradient") == 0 || strcmp(str, "CGD") == 0) m = ConjugateGradientDescent;
  else if (strcmp(str, "ClosedForm")                 == 0) m = ClosedForm;
  else if (strcmp(str, "LBFGS")                      == 0) m = LBFGS;
  else if (strcmp(str, "LineSearch")                 == 0) m = LineSearch;
  else if (strcmp(str, "EulerMethod")                == 0) m = EulerMethod;
  else if (strcmp(str, "EulerMethodWithMomentum")    == 0) m = EulerMethodWithMomentum;
  else return false;
  return true;
}

// ============================================================================
// FFD integration methods
// ============================================================================

/// Enumeration of implemented numerical integration methods
enum FFDIntegrationMethod
{
  FFDIM_UNKNOWN = 0,                                                   // keep at begin
  FFDIM_RKE1,   FFDIM_RKE2,   FFDIM_RKH2,  FFDIM_RK4,                  // explicit RK
  FFDIM_RKEH12, FFDIM_RKBS23, FFDIM_RKF45, FFDIM_RKDP45, FFDIM_RKCK45, // embedded RK
  FFDIM_SS, FFDIM_FastSS,                                              // scaling and squaring
  FFDIM_NUM                                                            // keep at end
};

// ---------------------------------------------------------------------------
inline const char *ToString(FFDIntegrationMethod m)
{
  switch (m) {
    case FFDIM_SS:     return "SS";
    case FFDIM_FastSS: return "FastSS";
    case FFDIM_RKE1:   return "RKE1";
    case FFDIM_RKE2:   return "RKE2";
    case FFDIM_RKH2:   return "RKH2";
    case FFDIM_RK4:    return "RK4";
    case FFDIM_RKEH12: return "RKEH12";
    case FFDIM_RKBS23: return "RKBS23";
    case FFDIM_RKF45:  return "RKF45";
    case FFDIM_RKDP45: return "RKDP45";
    case FFDIM_RKCK45: return "RKCK45";
    default:           return "Unknown";
  }
}

// ---------------------------------------------------------------------------
inline bool FromString(const char *str, FFDIntegrationMethod &m)
{
  // Convert string to all uppercase
  char *STR = new char[strlen(str) + 3]; // extra space for alias replacement
  strcpy(STR, str);
  for (char *p = STR; *p != '\0'; p++) *p = toupper(*p);
  // Translate aliases of methods to actual identifiers
  if      (strcmp(STR, "SCALINGANDSQUARING")        == 0) strcpy(STR, "SS");
  else if (strcmp(STR, "SCALING AND SQUARING")      == 0) strcpy(STR, "SS");
  else if (strcmp(STR, "FASTSCALINGANDSQUARING")    == 0) strcpy(STR, "FastSS");
  else if (strcmp(STR, "FAST SCALING AND SQUARING") == 0) strcpy(STR, "FastSS");
  else if (strcmp(STR, "EULER")                     == 0) strcpy(STR, "RKE1");
  else if (strcmp(STR, "FORWARDEULER")              == 0) strcpy(STR, "RKE1");
  else if (strcmp(STR, "FORWARD EULER")             == 0) strcpy(STR, "RKE1");
  else if (strcmp(STR, "FWDEULER")                  == 0) strcpy(STR, "RKE1");
  else if (strcmp(STR, "MODIFIEDEULER")             == 0) strcpy(STR, "RKE2");
  else if (strcmp(STR, "MODIFIED EULER")            == 0) strcpy(STR, "RKE2");
  else if (strcmp(STR, "MODEULER")                  == 0) strcpy(STR, "RKE2");
  else if (strcmp(STR, "HEUN")                      == 0) strcpy(STR, "RKH2");
  else if (strcmp(STR, "IMPROVEDEULER")             == 0) strcpy(STR, "RKH2");
  else if (strcmp(STR, "IMPROVED EULER")            == 0) strcpy(STR, "RKH2");
  else if (strcmp(STR, "IMPEULER")                  == 0) strcpy(STR, "RKH2");
  else if (strcmp(STR, "RK1")                       == 0) strcpy(STR, "RKE1");
  else if (strcmp(STR, "RK2")                       == 0) strcpy(STR, "RKH2");
  else if (strcmp(STR, "RK12")                      == 0) strcpy(STR, "RKEH12");
  else if (strcmp(STR, "RK23")                      == 0) strcpy(STR, "RKBS23");
  else if (strcmp(STR, "RK45")                      == 0) strcpy(STR, "RKCK45");
  // Convert string to enumeration value
  m = FFDIM_NUM;
  const char *s;
  do {
    m = static_cast<FFDIntegrationMethod>(static_cast<int>(m) - 1);
    if (m == FFDIM_UNKNOWN) break;
    s = ToString(m);
  } while (strcmp(s, str) != 0 && strcmp(s, STR) != 0);
  // Clean up
  delete[] STR;
  // Return whether conversion was successful
  return (m != FFDIM_UNKNOWN);
}


#endif
