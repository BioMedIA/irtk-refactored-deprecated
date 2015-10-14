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

#ifndef _IRTKEULERMETHOD_H
#define _IRTKEULERMETHOD_H

#include <irtkLocalOptimizer.h>

class irtkDeformableSurfaceModel;


/**
 * Minimizes deformable surface model using Euler's method
 *
 * This method has for example been used in the following paper on which this
 * implementation is also based on:
 *
 *   McInerney et al. (1999), Topology adaptive deformable surfaces for medical
 *   image volume segmentation. IEEE Transactions on Medical Imaging, 18(10),
 *   840–850. doi:10.1109/42.811261
 */
class irtkEulerMethod : public irtkLocalOptimizer
{
  irtkOptimizerMacro(irtkEulerMethod, EulerMethod);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Deformable surface model
  irtkReadOnlyAggregateMacro(irtkDeformableSurfaceModel, Model);

  /// Maximum number of integration steps
  irtkPublicAttributeMacro(int, NumberOfSteps);

  /// Damping coefficient
  irtkPublicAttributeMacro(double, DampingFactor);

  /// Length of each integration step (\delta t)
  irtkPublicAttributeMacro(double, StepLength);

  /// Current force acting on each node
  irtkComponentMacro(double, Force);

  /// Current displacement of each node
  irtkComponentMacro(double, Displacement);

private:

  /// Size of allocated vectors, may be larger than actual number of model DoFs!
  int _NumberOfDOFs;

  /// Copy attributes of this class from another instance
  void Copy(const irtkEulerMethod &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
public:

  /// Constructor
  irtkEulerMethod(irtkObjectiveFunction * = NULL);

  /// Copy constructor
  irtkEulerMethod(const irtkEulerMethod &);

  /// Assignment operator
  irtkEulerMethod &operator =(const irtkEulerMethod &);

  /// Destructor
  virtual ~irtkEulerMethod();

  // ---------------------------------------------------------------------------
  // Parameters
  using irtkLocalOptimizer::Parameter;

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  /// Get parameters as key/value as string map
  virtual irtkParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Execution

  /// Initialize optimization
  ///
  /// This member funtion is implicitly called by Run.
  virtual void Initialize();

  /// Integrate deformable surface model
  virtual double Run();

  /// Compute node velocity given its last velocity and the sum of node forces
  virtual void ComputeVelocity(double *, const double *) const;

protected:

  /// Update node velocities
  virtual void UpdateVelocity();

  /// Finalize optimization
  virtual void Finalize();

};


#endif
