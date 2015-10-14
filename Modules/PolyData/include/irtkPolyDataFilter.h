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

#ifndef _IRTKPOLYDATAFILTER_H
#define _IRTKPOLYDATAFILTER_H

#include <irtkObject.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>


namespace irtk { namespace polydata {


/**
 * Base class for filters which process polygonal surface meshes
 */
class irtkPolyDataFilter : public irtkObject
{
  irtkAbstractMacro(irtkPolyDataFilter);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Input surface mesh
  irtkPublicAttributeMacro(vtkSmartPointer<vtkPolyData>, Input);

  /// Output surface mesh (NULL if filter does not produce polygonal output)
  irtkReadOnlyAttributeMacro(vtkSmartPointer<vtkPolyData>, Output);

  /// Whether to output floating point results in double precision
  /// If \c true, output data arrays are of type vtkDoubleArray.
  /// If \c false, output data arrays are of type vtkFloatArray.
  irtkPublicAttributeMacro(bool, DoublePrecision);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Copy attributes of this filter from another instance
  void Copy(const irtkPolyDataFilter &);

protected:

  /// Default constructor
  irtkPolyDataFilter();

  /// Copy constructor
  irtkPolyDataFilter(const irtkPolyDataFilter &);

  /// Assignment operator
  irtkPolyDataFilter &operator =(const irtkPolyDataFilter &);

  /// Destructor
  virtual ~irtkPolyDataFilter();

  // ---------------------------------------------------------------------------
  // Execution

public:

  /// Run filter
  ///
  /// \note Most filters assume that irtkPolyData::BuildLinks of the input
  ///       surface mesh has been invoked before execution of the filter.
  ///       See the documentation/implementation of each individual filter.
  ///       Generally assume that this is required unless documented otherwise.
  virtual void Run();

protected:

  /// Initialize filter after input and parameters are set
  virtual void Initialize();

  /// Execute filter
  virtual void Execute() = 0;

  /// Finalize filter execution
  virtual void Finalize();

  // ---------------------------------------------------------------------------
  // Alternative VTK-like API

public:

  /// Enable/disable double precision output
  irtkOnOffMacro(DoublePrecision);

  /// Set input surface mesh
  void SetInputData(vtkPolyData *);

  /// Set input surface mesh
  void SetInput(vtkPolyData *);

  /// Run filter
  void Update();

  /// Get output surface mesh
  vtkPolyData *GetOutput();

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Alternative VTK-like API
// =============================================================================

// -----------------------------------------------------------------------------
inline void irtkPolyDataFilter::SetInputData(vtkPolyData *poly)
{
  this->Input(poly);
}

// -----------------------------------------------------------------------------
inline void irtkPolyDataFilter::SetInput(vtkPolyData *poly)
{
  this->Input(poly);
}

// -----------------------------------------------------------------------------
inline void irtkPolyDataFilter::Update()
{
  this->Run();
}

// -----------------------------------------------------------------------------
inline vtkPolyData *irtkPolyDataFilter::GetOutput()
{
  return this->Output();
}


} } // namespace irtk::polydata

#endif
