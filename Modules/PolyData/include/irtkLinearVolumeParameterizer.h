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

#ifndef _IRTKLINEARVOLUMEPARAMETERIZER_H
#define _IRTKLINEARVOLUMEPARAMETERIZER_H

#include <irtkVolumeParameterizer.h>

class irtkMatrix3x3;


namespace irtk { namespace polydata {


/**
 * Base class of volume re-parameterization filters based on discrete linear operators
 *
 * The reparameterization is the solution of a system of linear equations.
 */
class irtkLinearVolumeParameterizer : public irtkVolumeParameterizer
{
  irtkAbstractMacro(irtkLinearVolumeParameterizer);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Maximum number of linear solver iterations
  irtkPublicAttributeMacro(int, NumberOfIterations);

  /// Linear solver tolerance
  irtkPublicAttributeMacro(double, Tolerance);

  /// Relaxation factor
  irtkPublicAttributeMacro(double, RelaxationFactor);

  /// Original ID of n-th interior point
  irtkReadOnlyAttributeMacro(vector<int>, InteriorPointId);

  /// Variable/linear equation offset of n-th interior point
  irtkReadOnlyAttributeMacro(vector<int>, InteriorPointPos);

  /// Copy attributes of this class from another instance
  void Copy(const irtkLinearVolumeParameterizer &);

public:

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  irtkLinearVolumeParameterizer();

  /// Copy constructor
  irtkLinearVolumeParameterizer(const irtkLinearVolumeParameterizer &);

  /// Assignment operator
  irtkLinearVolumeParameterizer &operator =(const irtkLinearVolumeParameterizer &);

  /// Destructor
  virtual ~irtkLinearVolumeParameterizer();

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Initialize filter after input and parameters are set
  virtual void Initialize();

  /// Parameterize interior points
  virtual void Parameterize();

  /// Solve linear system with operator weights computed using the passed object
  void Solve(const irtkLinearVolumeParameterizer *);

  // ---------------------------------------------------------------------------
  // Auxiliary functions

public:

  /// Calculate operator weight for given tetrahadron
  ///
  /// \param[in] cellId ID of tetrahedron.
  /// \param[in] v0     First  vertex/point of tetrahedron.
  /// \param[in] v1     Second vertex/point of tetrahedron.
  /// \param[in] v2     Third  vertex/point of tetrahedron.
  /// \param[in] v3     Fourth vertex/point of tetrahedron.
  /// \param[in] volume Volume of tetrahedron.
  ///
  /// \return Operator weight contribution of tetrahedron.
  virtual irtkMatrix3x3 GetWeight(vtkIdType cellId,
                                  const double v0[3],
                                  const double v1[3],
                                  const double v2[3],
                                  const double v3[3],
                                  double       volume) const = 0;

};


} } // namespace irtk::polydata

#endif
