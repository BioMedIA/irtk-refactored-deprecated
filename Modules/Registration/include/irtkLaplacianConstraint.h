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

#ifndef _IRTKLAPLACIANCONSTRAINT_H
#define _IRTKLAPLACIANCONSTRAINT_H

#include <irtkSurfaceConstraint.h>
#include <irtkEdgeTable.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkDataArray.h>


/**
 * Internal surface bending/tensile force based on discrete Laplace operator
 *
 * This internal deformable surface force term is based on the umbrella
 * operator which approximates the discrete Laplace operator (cf. the
 * "tensile force" used by McInerney & Terzopoulos).
 *
 * McInerney & Terzopoulos, Topology adaptive deformable surfaces for medical
 * image volume segmentation. IEEE Transactions on Medical Imaging, 18(10),
 * 840–850, doi:10.1109/42.811261 (1999)
 *
 * \sa irtkCurvatureConstraint
 */
class irtkLaplacianConstraint : public irtkSurfaceConstraint
{
  irtkObjectMacro(irtkLaplacianConstraint);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Edge table
  irtkAttributeMacro(irtk::polydata::irtkEdgeTable, EdgeTable);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  irtkLaplacianConstraint(const char * = "", double = 1.0);

  /// Copy constructor
  irtkLaplacianConstraint(const irtkLaplacianConstraint &);

  /// Assignment operator
  irtkLaplacianConstraint &operator =(const irtkLaplacianConstraint &);

  /// Destructor
  virtual ~irtkLaplacianConstraint();

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Initialize internal force term
  virtual void Initialize();

  /// Reinitialize internal force term after change of input topology
  virtual void Reinitialize();

protected:

  /// Common (re-)initialization code of this class only (non-virtual function!)
  void Init();

  /// Evaluate energy of internal force term
  /// \return Infinity if only the internal force (i.e., gradient) is provided.
  virtual double Evaluate();

  /// Evaluate internal force w.r.t. transformation parameters or surface nodes
  virtual void EvaluateGradient(double *, double, double);

};


#endif
