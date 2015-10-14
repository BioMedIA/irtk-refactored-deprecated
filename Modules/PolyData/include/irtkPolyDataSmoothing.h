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

#ifndef IRTKPOLYDATASMOOTHING_H_
#define IRTKPOLYDATASMOOTHING_H_

#include <irtkPolyDataFilter.h>

class vtkDataArray;


namespace irtk { namespace polydata {

class irtkEdgeTable;


/**
 * Smooth scalars and/or points of triangulated surface mesh
 */
class irtkPolyDataSmoothing : public irtkPolyDataFilter
{
  irtkObjectMacro(irtkPolyDataSmoothing);

  // ---------------------------------------------------------------------------
  // Types

public:

  /// Enumeration of smoothing kernel functions
  enum WeightFunction
  {
    Combinatorial,      ///< Uniform node weights
    InverseDistance,    ///< Inverse node distance
    Gaussian,           ///< Gaussian node weights
    AnisotropicGaussian ///< Anisotropic Gaussian node weights
  };

  /// List of point data array names
  typedef vector<string> ArrayNames;

  /// Vector of point data arrays to be smoothed
  typedef vector<vtkSmartPointer<vtkDataArray> > DataArrays;

  // ---------------------------------------------------------------------------
  // Attributes

protected:

  /// Edge table of input mesh
  irtkPublicAggregateMacro(const irtkEdgeTable, EdgeTable);

  /// Whether edge table is to be destroyed by this filter
  irtkPublicAttributeMacro(bool, EdgeTableOwner);

  /// Number of smoothing iterations
  irtkPublicAttributeMacro(int, NumberOfIterations);

  /// Relaxation factor
  irtkPublicAttributeMacro(double, Lambda);

  /// Smoothing kernel parameter
  ///
  /// In case of a Gaussian smoothing kernel, if the sigma value is negative,
  /// the standard deviation of the Gaussian kernel is set to the average
  /// edge length times the absolute value of \c _Sigma. If _Sigma is zero,
  /// the average edge length is used as standard deviation. The _Sigma attribute
  /// is set to the actual used standard deviation after the filter execution.
  ///
  /// In case of the inverse distance weighting, the _Sigma value is added to
  /// the edge length before computing the inverse value.
  irtkPublicAttributeMacro(double, Sigma);

  /// Smoothing kernel parameter in direction of maximum curvature
  ///
  /// \note The direction of maximum principle curvature is orthogonal to the
  ///       direction in which the surface is most bended! It is the direction
  ///       with the most variance, i.e., along ridges, not orthogonal to these.
  ///
  /// This parameter is only used by the anisotropic Gaussian kernel.
  /// By default, i.e., when _Sigma2 = 0, the standard deviation along the
  /// direction of maximum change is half the standard deviation along the
  /// direction of minimum change. Hence, the surface points or data values are
  /// smoothed less in the direction of maximum change (i.e., maximum curvature).
  /// If the sigma value is negative, the standard deviation is set to the
  /// average edge length times the absolute value of \c _Sigma2.
  ///
  /// When an array of local geometry tensors is used instead of the direction
  /// of minimum and/or maximum change, the default is to use an isotropic
  /// Gaussian kernel in the local coordinate system defined by the tensor.
  /// In this case the axes of the local coordinate system are expected to be
  /// scaled anisotropically as in case of the curvature tensor, for example.
  irtkPublicAttributeMacro(double, MaximumDirectionSigma);

  /// Smoothing kernel function
  irtkPublicAttributeMacro(WeightFunction, Weighting);

  /// Name of input point data array with local geometry tensor used for anisotropic smoothing
  ///
  /// For example, the local curvature tensors can be computed using
  /// irtkPolyDataCurvature and used for anisotropic Gaussian smoothing.
  /// The input point data array is then named irtkPolyDataCurvature::TENSOR.
  irtkPublicAttributeMacro(string, GeometryTensorName);

  /// Name of input point data array with direction along which to smooth less
  /// \note This array is only used if no _GeometryTensorName is specified.
  irtkPublicAttributeMacro(string, MinimumDirectionName);

  /// Name of input point data array with direction along which to smooth more
  /// \note This array is only used if no _GeometryTensorName is specified.
  irtkPublicAttributeMacro(string, MaximumDirectionName);

  /// Whether to average values of adjacent nodes only or to
  /// also include the node's values themselves in the average
  ///
  /// \note In case of an InverseDistance node weighting, the values of the
  ///       node itself are only included in the average if _Sigma > .0.
  irtkPublicAttributeMacro(bool, AdjacentValuesOnly);

  /// Whether to smooth the node positions, i.e., input geometry
  irtkPublicAttributeMacro(bool, SmoothPoints);

  /// Names of input point data arrays to be smoothed
  irtkPublicAttributeMacro(ArrayNames, SmoothArrays);

  /// Input point data arrays to be smoothed
  irtkAttributeMacro(DataArrays, InputArrays);

  /// Output point data arrays
  irtkAttributeMacro(DataArrays, OutputArrays);

  /// Copy attributes of this class from another instance
  void Copy(const irtkPolyDataSmoothing &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  irtkPolyDataSmoothing();

  /// Copy constructor
  irtkPolyDataSmoothing(const irtkPolyDataSmoothing &);

  /// Assignment operator
  irtkPolyDataSmoothing &operator =(const irtkPolyDataSmoothing &);

  /// Destructor
  virtual ~irtkPolyDataSmoothing();

  /// Add named point data array to list of arrays to be smoothed
  void SmoothArray(const char *);

  // ---------------------------------------------------------------------------
  // Execution

protected:

  /// Initialize filter after input and parameters are set
  virtual void Initialize();

  /// Execute filter
  virtual void Execute();

  /// Finalize filter execution
  virtual void Finalize();

  // ---------------------------------------------------------------------------
  // Alternative VTK-like interface

public:

  /// Set number of smoothing iterations
  irtkSetMacro(NumberOfIterations, int);

  /// Get number of smoothing iterations
  irtkGetMacro(NumberOfIterations, int);

  /// Set relaxation factor
  irtkSetMacro(Lambda, double);

  /// Get relaxation factor
  irtkGetMacro(Lambda, double);

  /// Set smoothing kernel standard deviation
  irtkSetMacro(Sigma, double);

  /// Get smoothing kernel standard deviation
  irtkGetMacro(Sigma, double);

  /// Set smoothing kernel standard deviation in direction of maximum curvature
  irtkSetMacro(MaximumDirectionSigma, double);

  /// Get smoothing kernel standard deviation in direction of maximum curvature
  irtkGetMacro(MaximumDirectionSigma, double);

  /// Enable/disable averaging of adjacent node values only
  irtkOnOffMacro(AdjacentValuesOnly);

  /// Enable/disable smoothing of node positions
  irtkOnOffMacro(SmoothPoints);

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline void irtkPolyDataSmoothing::SmoothArray(const char *name)
{
  _SmoothArrays.push_back(name);
}


} } // namespace irtk::polydata

#endif
