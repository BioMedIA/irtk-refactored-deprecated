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

#ifndef _IRTKASCONFORMALASPOSSIBLEVOLUMEPARAMETERIZER_H
#define _IRTKASCONFORMALASPOSSIBLEVOLUMEPARAMETERIZER_H

#include <irtkLinearVolumeParameterizer.h>

#include <irtkMatrix3x3.h>


namespace irtk { namespace polydata {


/**
 * As-conformal-as-possible (ACAP) discrete volumetric map
 *
 * Paillé & Poulin (2012), As-conformal-as-possible discrete volumetric mapping,
 * Computers and Graphics (Pergamon), 36(5), 427–433. doi:10.1016/j.cag.2012.03.014
 */
class irtkAsConformalAsPossibleVolumeParameterizer : public irtkLinearVolumeParameterizer
{
  irtkObjectMacro(irtkAsConformalAsPossibleVolumeParameterizer);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Uniform weight of scale and angle conformality
  irtkPublicAttributeMacro(double, UniformWeight);

  /// Local orientation of tetrahedron (rotation matrix)
  irtkAttributeMacro(vector<irtkMatrix3x3>, Orientation);

  /// Copy attributes of this class from another instance
  void Copy(const irtkAsConformalAsPossibleVolumeParameterizer &);

public:

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Default constructor
  irtkAsConformalAsPossibleVolumeParameterizer();

  /// Copy constructor
  irtkAsConformalAsPossibleVolumeParameterizer(const irtkAsConformalAsPossibleVolumeParameterizer &);

  /// Assignment operator
  irtkAsConformalAsPossibleVolumeParameterizer &operator =(const irtkAsConformalAsPossibleVolumeParameterizer &);

  /// Destructor
  virtual ~irtkAsConformalAsPossibleVolumeParameterizer();

  /// Type of volumetric map produced by this parameterizer
  virtual MapType TypeOfMap() const { return ACAPMap; }

  // ---------------------------------------------------------------------------
  // Exection

protected:

  /// Initialize filter after input and parameters are set
  void Initialize();

  /// Finalize filter execution
  void Finalize();

  // ---------------------------------------------------------------------------
  // Auxiliary functions

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
                                  double       volume) const;

};


} } // namespace irtk::polydata

#endif
