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

#ifndef _IRTKPOINTCORRESPONDENCE_H
#define _IRTKPOINTCORRESPONDENCE_H

#include <vtkSmartPointer.h>
#include <vtkPointSet.h>

#include <irtkObservable.h>
#include <irtkPointLocator.h>
#include <irtkRegisteredPointSet.h>


/**
 * Locator of corresponding points of point clouds, curves, or surfaces
 *
 * An instance of a subclass derived from irtkPointCorrespondence
 * finds the point of the source data set which corresponds to the specified
 * point in the target data set, i.e., a solution of the linear assignment
 * problem. In case of curves or surfaces, the corresponding point need not be
 * an explicit point of the source data set, but can also be a point on the
 * line segment or surface element (cell), respectively.
 */

class irtkPointCorrespondence : public irtkObservable
{
  irtkAbstractMacro(irtkPointCorrespondence);

  // ---------------------------------------------------------------------------
  // Types
public:

  /// Enumeration of available point correspondence maps
  enum TypeId
  {
    Unknown,            ///< Unknown/invalid/default correspondence
    FiducialMatch,      ///< Match points by their respective indices
    ClosestPoint,       ///< Find closest point in source data set
    ClosestPointLabel,  ///< Find closest point in source data set with identical label
    ClosestCell,        ///< Find closest cell in source data set
    SpectralMatch,      ///< Closest point interpolated using spectral similarity
    RobustClosestPoint, ///< Find closest points within proximity in both directions
    RobustPointMatch    ///< Robust point matching (RPM) correspondence
  };

  /// Enumeration of correspondence mapping direction
  enum Direction
  {
    TargetToSource,
    SourceToTarget
  };

  /// List of point features to use for point matching
  typedef irtkPointLocator::FeatureList FeatureList;

  /// Info structure of point feature to use for point matching
  typedef irtkPointLocator::FeatureInfo FeatureInfo;

  // ---------------------------------------------------------------------------
  // Static helpers for subclass implementation
public:

  /// Get number of points
  ///
  /// \param[in] dataset Dataset.
  /// \param[in] sample  Indices of points in \p dataset to consider only or NULL for all.
  ///
  /// \returns Number of (sample) points.
  static int GetNumberOfPoints(vtkPointSet *dataset, const vector<int> *sample = NULL);

  /// Get number of points
  ///
  /// \param[in] dataset Dataset.
  /// \param[in] sample  Indices of points in \p dataset to consider only or NULL for all.
  ///
  /// \returns Number of (sample) points.
  static int GetNumberOfPoints(const irtkRegisteredPointSet *dataset, const vector<int> *sample = NULL);

  /// Get size of feature vectors
  ///
  /// \param[in] dataset Dataset.
  /// \param[in] feature Indices/names and weights/rescaling parameters of features.
  static int GetPointDimension(vtkPointSet *dataset, const FeatureList *feature);

  /// Get size of feature vectors
  ///
  /// \param[in] dataset Dataset.
  /// \param[in] feature Indices/names and weights/rescaling parameters of features.
  static int GetPointDimension(const irtkRegisteredPointSet *dataset, const FeatureList *feature);

  /// Get index of point
  ///
  /// \param[in] dataset Dataset.
  /// \param[in] sample  Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in] index   Index of point in \p dataset or \p sample (if not NULL).
  ///
  /// \returns Index of feature vector/point in \p dataset.
  static int GetPointIndex(vtkPointSet *dataset, const vector<int> *sample, int index);

  /// Get index of point
  ///
  /// \param[in] dataset Dataset.
  /// \param[in] sample  Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in] index   Index of point in \p dataset or \p sample (if not NULL).
  ///
  /// \returns Index of feature vector/point in \p dataset.
  static int GetPointIndex(const irtkRegisteredPointSet *dataset, const vector<int> *sample, int index);

  /// Get spatial (sample) point
  ///
  /// \param[out] point   Feature vector.
  /// \param[in]  dataset Dataset.
  /// \param[in]  sample  Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in]  index   Index of point in \p dataset or \p sample (if not NULL).
  static void GetPoint(irtkPoint &point, vtkPointSet *dataset, const vector<int> *sample, int index);

  /// Get spatial (sample) point
  ///
  /// \param[out] point   Feature vector.
  /// \param[in]  dataset Dataset.
  /// \param[in]  sample  Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in]  index   Index of point in \p dataset or \p sample (if not NULL).
  static void GetPoint(irtkPoint &point, const irtkRegisteredPointSet *dataset, const vector<int> *sample, int index);

  /// Get feature vector
  ///
  /// \param[out] point   Feature vector.
  /// \param[in]  dataset Dataset.
  /// \param[in]  sample  Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in]  index   Index of point in \p dataset or \p sample (if not NULL).
  /// \param[in]  feature Indices/names and weights/rescaling parameters of features.
  static void GetPoint(double *point, vtkPointSet *dataset, const vector<int> *sample,
                       int index, const FeatureList *feature = NULL);

  /// Get feature vector
  ///
  /// \param[out] point   Feature vector.
  /// \param[in]  dataset Dataset.
  /// \param[in]  sample  Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in]  index   Index of point in \p dataset or \p sample (if not NULL).
  /// \param[in]  feature Indices/names and weights/rescaling parameters of features.
  static void GetPoint(double *point, const irtkRegisteredPointSet *dataset,
                       const vector<int> *sample, int index, const FeatureList *feature = NULL);

  /// Calculate squared Euclidean distance between feature vectors
  ///
  /// \param[in] a First  feature vector.
  /// \param[in] b Second feature vector.
  /// \param[in] d Dimension of feature vectors.
  ///
  /// \returns Squared Euclidean distance, i.e., dot product of a and b.
  static double Distance2BetweenPoints(const double *a, const double *b, int d = 3);

  /// Get subset of sample points
  ///
  /// \param[in] dataset Dataset.
  /// \param[in] sample  Indices of points in \p dataset to extract or NULL for all.
  ///
  /// \returns Set of \p sample points only.
  static vtkSmartPointer<vtkPoints> GetPoints(vtkPointSet *dataset, const vector<int> *sample);

  /// Get subset of sample points
  ///
  /// \param[in] dataset Registered dataset.
  /// \param[in] sample  Indices of points in \p dataset to extract or NULL for all.
  ///
  /// \returns Set of \p sample points only.
  static vtkSmartPointer<vtkPoints> GetPoints(const irtkRegisteredPointSet *dataset, const vector<int> *sample);

  /// Get subset of sample points
  ///
  /// \param[in] dataset Dataset.
  /// \param[in] sample  Indices of points in \p dataset to extract or NULL for all.
  ///
  /// \returns Dataset consisting of \p sample points only.
  static vtkSmartPointer<vtkPointSet> GetPointSet(vtkPointSet *dataset, const vector<int> *sample);

  /// Get subset of sample points
  ///
  /// \param[in] dataset Registered dataset.
  /// \param[in] sample  Indices of points in \p dataset to extract or NULL for all.
  ///
  /// \returns Dataset consisting of \p sample points only.
  static vtkSmartPointer<vtkPointSet> GetPointSet(const irtkRegisteredPointSet *dataset, const vector<int> *sample);

  /// Get index of point data array using case insensitive name
  static int GetPointDataIndexByCaseInsensitiveName(vtkPointData *, string);

  // ---------------------------------------------------------------------------
  // Attributes
protected:

  /// Number of point (samples) in target data set
  irtkAttributeMacro(int, M);

  /// Number of point (samples) in source data set
  irtkAttributeMacro(int, N);

  /// Dimension of feature vectors
  irtkReadOnlyAttributeMacro(int, NumberOfFeatures);

  /// (Transformed) target data set
  irtkPublicAggregateMacro(const irtkRegisteredPointSet, Target);

  /// Indices of target point samples or NULL if all points are considered
  irtkPublicAggregateMacro(const vector<int>, TargetSample);

  /// Indices and rescaling parameters of target point features
  irtkPublicAttributeMacro(FeatureList, TargetFeatures);

  /// (Transformed) source data set
  irtkPublicAggregateMacro(const irtkRegisteredPointSet, Source);

  /// Indices of source point samples or NULL if all points are considered
  irtkPublicAggregateMacro(const vector<int>, SourceSample);

  /// Indices and rescaling parameters of source point features
  irtkPublicAttributeMacro(FeatureList, SourceFeatures);

  /// Eigenvalues corresponding to spectral coordinates of target
  irtkAttributeMacro(irtkVector, TargetEigenvalues);

  /// Eigenvalues corresponding to spectral coordinates of source
  irtkAttributeMacro(irtkVector, SourceEigenvalues);

  /// Number of eigenmodes used for spectral features if used
  irtkPublicAttributeMacro(int, DimensionOfSpectralPoints);

  /// Whether to use diffeomorphic spectral matching for spectral decomposition
  irtkPublicAttributeMacro(bool, DiffeomorphicSpectralDecomposition);

  /// Whether to update spectral coordinates of moving dataset
  irtkPublicAttributeMacro(bool, UpdateSpectralPoints);

  /// Whether target to source correspondences are needed (i.e., GetSourcePoint)
  irtkPublicAttributeMacro(bool, FromTargetToSource);

  /// Whether source to target correspondences are needed (i.e., GetTargetPoint)
  irtkPublicAttributeMacro(bool, FromSourceToTarget);

  /// Default direction (i.e., GetPoint)
  irtkPublicAttributeMacro(Direction, DefaultDirection);

  // ---------------------------------------------------------------------------
  // Construction/Destruction
protected:

  /// Constructor
  irtkPointCorrespondence(const irtkRegisteredPointSet * = NULL,
                          const irtkRegisteredPointSet * = NULL);

  /// Copy constructor
  irtkPointCorrespondence(const irtkPointCorrespondence &);

  /// Fill in missing point feature info
  ///
  /// This function is called by Initialize after the input data has
  /// been set. The AddFeature function cannot inspect the input data
  /// to determine the index of a named feature because the input may
  /// not be set yet when AddFeature is called.
  void CompleteFeatureInfo(const irtkRegisteredPointSet *, FeatureList &);

public:

  /// Construct new instance of specified type
  static irtkPointCorrespondence *New(TypeId);

  /// Construct new instance of specified type
  static irtkPointCorrespondence *New(const char *);

  /// Copy construct a new instance
  virtual irtkPointCorrespondence *NewInstance() const = 0;

  /// Destructor
  virtual ~irtkPointCorrespondence();

  /// Type enumeration value
  virtual TypeId Type() const = 0;

  // ---------------------------------------------------------------------------
  // Parameters

  /// Set parameter value from string
  virtual bool Set(const char *, const char *);

  // Do not hide other overload
  using irtkObservable::Parameter;

  /// Get parameter key/value as string map
  virtual irtkParameterList Parameter() const;

  /// Add point data array to use as feature in point matching
  ///
  /// This function is used to specify the point data to use as features in
  /// the point matching. By default, if no features are specified, the
  /// spatial coordinates of the points are used. When a feature is added
  /// to the lists _TargetFeatures and _SourceFeatures, only the named features
  /// are used. The point matching is then done between n-dimensional vectors
  /// which are the concatenation of the selected features, where each feature
  /// may be linearly rescaled to normalize them or adjust their weight in the
  /// matching function. If the spatial coordinates are to be used as features
  /// along with other point data features, the feature "spatial coordinates"
  /// or "spatial point" must be added explicitly as well.
  ///
  /// The feature value used is \p weight * (\p slope * value + \p intercept).
  ///
  /// \param[in] name      Case insensitive name of point data array in both
  ///                      target and source dataset.
  /// \param[in] weight    Weight of feature.
  /// \param[in] slope     Slope of feature rescaling function.
  /// \param[in] intercept Intercept of feature rescaling function, e.g., to
  ///                      normalize feature to zero mean.
  ///
  /// \note Whether a feature name is valid is determined by FinalizeFeatureLists
  ///       called by Initialize. This cannot be done by AddFeature because the
  ///       input poly data objects may not be set yet.
  bool AddFeature(const char *name, double weight = 1.0, double slope = 1.0, double intercept = .0);

  /// Remove point data array from list of features to use for point matching
  ///
  /// \param[in] name Name of point data array in both target and source dataset.
  void RemoveFeature(const char *name);

  // ---------------------------------------------------------------------------
  // Correspondences

  /// Initialize correspondence map
  virtual void Initialize();

  /// Reinitialize correspondence map after change of input topology
  virtual void Reinitialize();

protected:

  /// Common (re-)initialization steps of this class
  /// \note Must be a non-virtual function!
  void Init();

public:

  /// Update correspondence map after change of input data sets
  virtual void Update();

  /// Update correspondence map after convergence
  /// \returns Whether correspondences have changed.
  virtual bool Upgrade();

  /// Get untransformed target point corresponding to i-th source (sample) point
  virtual bool GetInputTargetPoint(int, irtkPoint &) const = 0;

  /// Get untransformed source point corresponding to i-th target (sample) point
  virtual bool GetInputSourcePoint(int, irtkPoint &) const = 0;

  /// Get (transformed) target point corresponding to i-th source (sample) point
  virtual bool GetTargetPoint(int, irtkPoint &) const = 0;

  /// Get (transformed) source point corresponding to i-th target (sample) point
  virtual bool GetSourcePoint(int, irtkPoint &) const = 0;

  /// Get untransformed output point corresponding to i-th input (sample) point
  bool GetInputPoint(int, irtkPoint &) const;

  /// Get (transformed) output point corresponding to i-th input (sample) point
  bool GetPoint(int, irtkPoint &) const;

  /// Get index of target point corresponding to i-th source (sample) point
  ///
  /// \returns Index of corresponding target point and -1 if point is an outlier
  ///          or its corresponding point is not a target vertex position.
  virtual int GetTargetIndex(int) const;

  /// Get index of source point corresponding to i-th target (sample) point
  ///
  /// \returns Index of corresponding source point and -1 if point is an outlier
  ///          or its corresponding point is not a source vertex position.
  virtual int GetSourceIndex(int) const;

  /// \returns Index of corresponding point and -1 if point is an outlier
  ///          or its corresponding point is not a vertex position.
  int GetIndex(int) const;

  // ---------------------------------------------------------------------------
  // Debugging

  /// Write input of data fidelity term
  virtual void WriteDataSets(const char *, const char *, bool = true) const;

  /// Write first three spectral coordinates as polydata points
  virtual void WriteSpectralPoints(const char *, vtkPointSet *) const;

};


// -----------------------------------------------------------------------------
template <> bool FromString(const char *str, irtkPointCorrespondence::TypeId &);
template <> string ToString(const irtkPointCorrespondence::TypeId &, int, char, bool);

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Static helpers for subclass implementation
// =============================================================================

// -----------------------------------------------------------------------------
inline int irtkPointCorrespondence
::GetNumberOfPoints(vtkPointSet *dataset, const vector<int> *sample)
{
  return irtkPointLocator::GetNumberOfPoints(dataset, sample);
}

// -----------------------------------------------------------------------------
inline int irtkPointCorrespondence
::GetNumberOfPoints(const irtkRegisteredPointSet *dataset, const vector<int> *sample)
{
  return GetNumberOfPoints(dataset->PointSet(), sample);
}

// -----------------------------------------------------------------------------
inline int irtkPointCorrespondence
::GetPointDimension(vtkPointSet *dataset, const FeatureList *features)
{
  return irtkPointLocator::GetPointDimension(dataset, features);
}

// -----------------------------------------------------------------------------
inline int irtkPointCorrespondence
::GetPointDimension(const irtkRegisteredPointSet *dataset, const FeatureList *features)
{
  return GetPointDimension(dataset->PointSet(), features);
}

// -----------------------------------------------------------------------------
inline int irtkPointCorrespondence
::GetPointIndex(vtkPointSet *dataset, const vector<int> *sample, int index)
{
  return irtkPointLocator::GetPointIndex(dataset, sample, index);
}

// -----------------------------------------------------------------------------
inline int irtkPointCorrespondence
::GetPointIndex(const irtkRegisteredPointSet *dataset, const vector<int> *sample, int index)
{
  return GetPointIndex(dataset->PointSet(), sample, index);
}

// -----------------------------------------------------------------------------
inline void irtkPointCorrespondence
::GetPoint(irtkPoint &point, vtkPointSet *dataset, const vector<int> *sample, int index)
{
  return irtkPointLocator::GetPoint(point, dataset, sample, index);
}

// -----------------------------------------------------------------------------
inline void irtkPointCorrespondence
::GetPoint(irtkPoint &point, const irtkRegisteredPointSet *dataset, const vector<int> *sample, int index)
{
  return GetPoint(point, dataset->PointSet(), sample, index);
}

// -----------------------------------------------------------------------------
inline void irtkPointCorrespondence
::GetPoint(double *point, vtkPointSet *dataset, const vector<int> *sample, int index, const FeatureList *features)
{
  return irtkPointLocator::GetPoint(point, dataset, sample, index, features);
}

// -----------------------------------------------------------------------------
inline void irtkPointCorrespondence
::GetPoint(double *point, const irtkRegisteredPointSet *dataset,
           const vector<int> *sample, int index, const FeatureList *features)
{
  return GetPoint(point, dataset->PointSet(), sample, index, features);
}

// -----------------------------------------------------------------------------
inline double irtkPointCorrespondence
::Distance2BetweenPoints(const double *a, const double *b, int dim)
{
  return irtkPointLocator::Distance2BetweenPoints(a, b, dim);
}

// -----------------------------------------------------------------------------
inline vtkSmartPointer<vtkPoints> irtkPointCorrespondence
::GetPoints(vtkPointSet *dataset, const vector<int> *sample)
{
  vtkSmartPointer<vtkPoints> points;
  if (sample && sample->size() > 0) {
    points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(sample->size());
    double p[3];
    for (size_t i = 0; i < sample->size(); ++i) {
      dataset->GetPoint(i, p);
      points ->SetPoint(i, p);
    }
  } else {
    points = dataset->GetPoints();
  }
  return points;
}

// -----------------------------------------------------------------------------
inline vtkSmartPointer<vtkPoints> irtkPointCorrespondence
::GetPoints(const irtkRegisteredPointSet *dataset, const vector<int> *sample)
{
  return GetPoints(dataset->PointSet(), sample);
}

// -----------------------------------------------------------------------------
inline vtkSmartPointer<vtkPointSet> irtkPointCorrespondence
::GetPointSet(vtkPointSet *dataset, const vector<int> *sample)
{
  vtkSmartPointer<vtkPointSet> output;
  if (sample && sample->size() > 0) {
    output = dataset->NewInstance();
    output->SetPoints(GetPoints(dataset, sample));
  } else {
    output = dataset;
  }
  return output;
}

// -----------------------------------------------------------------------------
inline vtkSmartPointer<vtkPointSet> irtkPointCorrespondence
::GetPointSet(const irtkRegisteredPointSet *dataset, const vector<int> *sample)
{
  return GetPointSet(dataset->PointSet(), sample);
}

// =============================================================================
// Correspondences
// =============================================================================

// -----------------------------------------------------------------------------
inline bool irtkPointCorrespondence::GetInputPoint(int i, irtkPoint &p) const
{
  if (_DefaultDirection == TargetToSource) return this->GetInputSourcePoint(i, p);
  else                                     return this->GetInputTargetPoint(i, p);
}

// -----------------------------------------------------------------------------
inline bool irtkPointCorrespondence::GetPoint(int i, irtkPoint &p) const
{
  if (_DefaultDirection == TargetToSource) return this->GetSourcePoint(i, p);
  else                                     return this->GetTargetPoint(i, p);
}

// -----------------------------------------------------------------------------
inline int irtkPointCorrespondence::GetIndex(int i) const
{
  if (_DefaultDirection == TargetToSource) return this->GetSourceIndex(i);
  else                                     return this->GetTargetIndex(i);
}

// =============================================================================
// Auxiliary functions for subclass implementation
// =============================================================================

namespace irtkPointCorrespondenceUtils {


// -----------------------------------------------------------------------------
/// Draw spatially stratified sample points from given dataset
///
/// \param[in]  dataset    Dataset from which to draw point samples.
/// \param[out] indices    Indices of drawn points.
/// \param[in]  maxnum     Maximum number of samples.
/// \param[in]  maxdist    Approximate maximum distance between point samples.
/// \param[in]  stratified Sample spatially stratified points.
void SamplePoints(vtkPointSet *dataset, vector<int> &indices, int maxnum, double maxdist = .0, bool stratified = true);

// -----------------------------------------------------------------------------
/// Draw spatially stratified sample points from given dataset
///
/// \param[in]  dataset    Dataset from which to draw point samples.
/// \param[out] indices    Indices of drawn points.
/// \param[in]  maxnum     Maximum number of samples.
/// \param[in]  maxdist    Approximate maximum distance between point samples.
/// \param[in]  stratified Sample spatially stratified points.
inline void SamplePoints(const irtkRegisteredPointSet *dataset, vector<int> &indices,
                         int maxnum, double maxdist = .0, bool stratified = true)
{
  SamplePoints(dataset->InputPointSet(), indices, maxnum, maxdist, stratified);
}


} // namespace irtkPointCorrespondenceUtils


#endif
