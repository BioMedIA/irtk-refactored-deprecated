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

#ifndef _IRTKBALLOONFORCE_H
#define _IRTKBALLOONFORCE_H

#include <irtkSurfaceForce.h>

#include <vtkSmartPointer.h>
class vtkDataArray;


/**
 * Balloon/inflation force
 *
 * The balloon force inflates the surface while inside the object and deflates
 * it when outside the object assuming that the image intensities within the
 * object are within a certain range. It is therefore also referred to as inflation
 * force and is mainly used while the surface is not in the vicinity of an edge.
 *
 *   McInerney & Terzopoulos, Topology adaptive deformable surfaces for medical
 *   image volume segmentation. IEEE Transactions on Medical Imaging, 18(10),
 *   840–850. doi:10.1109/42.811261 (1999)
 *
 * This balloon force has also been used in
 *
 *   Park et al., A non-self-intersecting adaptive deformable surface for complex
 *   boundary extraction from volumetric images. Computer & Graphics, 25, 421–440,
 *   (2001)
 */
class irtkBalloonForce : public irtkSurfaceForce
{
  irtkObjectMacro(irtkBalloonForce);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Lower global intensity threshold for inside region
  irtkPublicAttributeMacro(double, LowerIntensity);

  /// Upper global intensity threshold for inside regin
  irtkPublicAttributeMacro(double, UpperIntensity);

  /// Multiplier of intensity standard deviation
  irtkPublicAttributeMacro(double, SigmaFactor);

  /// Multiplier of standard deviation of foreground intensities
  irtkPublicAttributeMacro(double, ForegroundSigmaFactor);

  /// Multiplier of standard deviation of background intensities
  irtkPublicAttributeMacro(double, BackgroundSigmaFactor);

  /// Box window radius within which to compute local image intensity statistics
  irtkPublicAttributeMacro(double, Radius);

  /// Positive multiplicative factor in (0, 1] used to progressively reduce
  /// the magnitude of the balloon force whenever its sign changes
  irtkPublicAttributeMacro(double, DampingFactor);

  /// Minimum magnitude threshold, once the magnitude falls below this threshold,
  /// it is set to zero such that the node is no longer affected by this force
  irtkPublicAttributeMacro(double, MagnitudeThreshold);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
protected:

  /// Copy attributes of this class from another instance
  void Copy(const irtkBalloonForce &);

public:

  /// Constructor
  irtkBalloonForce(const char * = "", double = 1.0);

  /// Copy constructor
  irtkBalloonForce(const irtkBalloonForce &);

  /// Assignment operator
  irtkBalloonForce &operator =(const irtkBalloonForce &);

  /// Destructor
  virtual ~irtkBalloonForce();

  // ---------------------------------------------------------------------------
  // Configuration

  // Import other overloads
  using irtkExternalForce::Parameter;

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  /// Get parameter name/value pairs
  virtual irtkParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Initialization
protected:

  /// Compute local intensity thresholds and/or background/foreground statistics
  void ComputeLocalIntensityAttributes(bool, bool);

public:

  /// Initialize external force once input and parameters have been set
  virtual void Initialize();

  /// Update moving input points and internal state of force term
  virtual void Update(bool = true);

  // ---------------------------------------------------------------------------
  // Evaluation
protected:

  /// Evaluate external force
  virtual void EvaluateGradient(double *, double, double);

};


#endif
