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

#ifndef _IRTKPOINTLOCATOR_H
#define _IRTKPOINTLOCATOR_H

#include <vtkSmartPointer.h>
#include <vtkPointSet.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkOctreePointLocator.h>

#include <irtkObservable.h>
#include <irtkPoint.h>

#ifdef HAS_FLANN
#  include <flann/flann.hpp>
#endif


/**
 * Point search structure for establishing point correspondences
 *
 * This point locator implementation is specialized for use by
 * irtkPointCorrespondence subclass implementations. It allows the search of
 * nearest neighbors within the n-dimensional feature space spanned by the
 * feature vectors used to establish point correspondences.
 *
 * The implementation uses either VTK for point search in three dimensions or
 * FLANN for higher dimensional feature spaces if available. Alternatively,
 * an ITK Kd tree is used if the library is available. As last resort, a brute
 * force search without actual Kd tree search structure is performed.
 *
 * \todo Consider actual implementation of own thread-safe N-dimensional Kd tree?
 */

class irtkPointLocator : public irtkObject
{
  irtkObjectMacro(irtkPointLocator);

  // ---------------------------------------------------------------------------
  // Feature vectors
public:

  /// Type storing name and/or index of point feature together with
  /// an optional linear rescaling function
  struct FeatureInfo
  {
    string _Name;      ///< Name of feature/point data array
    int    _Index;     ///< Index of feature/point data array
    double _Weight;    ///< Weight of feature/point data
    double _Slope;     ///< Rescaling slope of feature/point data
    double _Intercept; ///< Rescaling intercept of feature/point data

    /// Constructor
    FeatureInfo(int index = -2, double weight = 1.0, double slope = 1.0, double intercept = .0)
    :
      _Name(""), _Index(index), _Weight(weight), _Slope(slope), _Intercept(intercept)
    {}

    /// Constructor
    FeatureInfo(const char *name, double weight = 1.0, double slope = 1.0, double intercept = .0)
    :
      _Name(name), _Index(-2), _Weight(weight), _Slope(slope), _Intercept(intercept)
    {}

    /// Constructor
    FeatureInfo(const string &name, double weight = 1.0, double slope = 1.0, double intercept = .0)
    :
      _Name(name), _Index(-2), _Weight(weight), _Slope(slope), _Intercept(intercept)
    {}

    /// Constructor
    FeatureInfo(const char *name, int index, double weight = 1.0, double slope = 1.0, double intercept = .0)
    :
      _Name(name), _Index(index), _Weight(weight), _Slope(slope), _Intercept(intercept)
    {}

    /// Constructor
    FeatureInfo(const string &name, int index, double weight = 1.0, double slope = 1.0, double intercept = .0)
    :
      _Name(name), _Index(index), _Weight(weight), _Slope(slope), _Intercept(intercept)
    {}
  };

  /// List of point features to use for nearest neighbor search
  typedef vector<FeatureInfo> FeatureList;

  /// Get point data array
  ///
  /// \param[in] dataset Dataset.
  /// \param[in] feature Index/name and weight/rescaling parameters.
  ///
  /// \returns Pointer to point data array of dataset or NULL.
  static vtkDataArray *GetDataArray(vtkPointSet *dataset, const FeatureInfo &feature);

  /// Get number of points
  ///
  /// \param[in] dataset Dataset.
  /// \param[in] sample  Indices of points in \p dataset to consider only or NULL for all.
  ///
  /// \returns Number of (sample) points.
  static int GetNumberOfPoints(vtkPointSet *dataset, const vector<int> *sample = NULL);

  /// Get size of feature vectors
  ///
  /// \param[in] dataset Dataset.
  /// \param[in] feature Indices/names and weights/rescaling parameters of features.
  static int GetPointDimension(vtkPointSet *dataset, const FeatureList *feature);

  /// Get index of feature vector
  ///
  /// \param[in] dataset Dataset.
  /// \param[in] sample  Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in] index   Index of point in \p dataset or \p sample (if not NULL).
  ///
  /// \returns Index of feature vector/point in \p dataset.
  static int GetPointIndex(vtkPointSet *dataset, const vector<int> *sample, int index);

  /// Get spatial (sample) point
  ///
  /// \param[out] point   Spatial point.
  /// \param[in]  dataset Dataset.
  /// \param[in]  sample  Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in]  index   Index of point in \p dataset or \p sample (if not NULL).
  static void GetPoint(irtkPoint &point, vtkPointSet *dataset, const vector<int> *sample, int index);

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
  /// \param[in]  index   Index of point in \p dataset.
  /// \param[in]  feature Indices/names and weights/rescaling parameters of features.
  static void GetPoint(double *point, vtkPointSet *dataset, int index, const FeatureList *feature = NULL);

  /// Calculate squared Euclidean distance between feature vectors
  ///
  /// \param[in] a First  feature vector.
  /// \param[in] b Second feature vector.
  /// \param[in] d Dimension of feature vectors.
  ///
  /// \returns Squared Euclidean distance, i.e., dot product of a and b.
  static double Distance2BetweenPoints(const double *a, const double *b, int d = 3);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Dataset for which search structure is build
  irtkPublicAggregateMacro(vtkPointSet, DataSet);

  /// Indices of points to consider only or NULL
  irtkPublicAggregateMacro(const vector<int>, Sample);

  /// Indices/names and rescaling parameters of point data arrays
  irtkPublicAttributeMacro(FeatureList, Features);

  /// Number of points in search structure
  irtkReadOnlyAttributeMacro(int, NumberOfPoints);

  /// Dimension of feature vectors/points
  irtkReadOnlyAttributeMacro(int, PointDimension);

  /// VTK point locator used for three-dimensional feature spaces
  vtkSmartPointer<vtkOctreePointLocator> _VtkLocator;

#ifdef HAS_FLANN

  /// Datatype used for internal FLANN index
  typedef float FlannType;

  /// FLANN Kd tree used for n-dimensional feature spaces
  flann::Index<flann::L2_Simple<FlannType> > _FlannTree;

  /// Create FLANN matrix containing feature vectors of points of given dataset
  flann::Matrix<FlannType> FlannMatrix(vtkPointSet *, const vector<int> * = NULL,
                                       const FeatureList * = NULL);

  /// Create FLANN matrix from single feature vector
  flann::Matrix<FlannType> FlannMatrix(double *);

#elif HAS_ITK

  /// TODO: ITK Kd tree used for n-dimensional feature spaces
  //itk::Statistics::KdTree::Pointer _ItkTree;

#endif

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Constructor
  irtkPointLocator();

  /// Initialize point locator
  void Initialize();

private:

  /// Copy constructor
  /// \note Intentionally not implemented
  irtkPointLocator(const irtkPointLocator &);

  /// Assignment operator
  /// \note Intentionally not implemented
  void operator =(const irtkPointLocator &);

public:

  /// Construct new point locator for search in given dataset
  ///
  /// \param[in] dataset Dataset in which points are searched.
  /// \param[in] sample  Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in] feature Indices and weights of point data in \p dataset to use.
  static irtkPointLocator *New(vtkPointSet       *dataset,
                               const vector<int> *sample  = NULL,
                               const FeatureList *feature = NULL);

  /// Destructor
  virtual ~irtkPointLocator();

  // ---------------------------------------------------------------------------
  // Closest point

  /// Find closest point
  ///
  /// \param[in]  point Query point.
  /// \param[out] dist2 Squared Euclidean distance of closest point.
  ///
  /// \returns Index of point/sample in dataset closest to the query \p point.
  int FindClosestPoint(double *point, double *dist2 = NULL);

  /// Find closest point
  ///
  /// \param[in]  dataset  Dataset whose nearest neighbor is queried.
  /// \param[in]  sample   Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in]  index    Index of point in \p dataset or \p sample (if not NULL).
  /// \param[in]  features Indices and weights of point data in \p dataset to use.
  /// \param[out] dist2    Squared Euclidean distance of closest point.
  ///
  /// \returns Index of point/sample in dataset closest to the i-th point in \p dataset.
  int FindClosestPoint(vtkPointSet       *dataset,
                       const vector<int> *sample,
                       int                index,
                       const FeatureList *features,
                       double            *dist2 = NULL);

  /// Find closest point
  ///
  /// \param[in]  dataset  Dataset whose nearest neighbor is queried.
  /// \param[in]  index    Index of point in \p dataset or \p sample (if not NULL).
  /// \param[in]  features Indices and weights of point data in \p dataset to use.
  /// \param[out] dist2    Squared Euclidean distance of closest point.
  ///
  /// \returns Index of point/sample in dataset closest to the i-th point in \p dataset.
  int FindClosestPoint(vtkPointSet       *dataset,
                       int                index,
                       const FeatureList *features,
                       double            *dist2 = NULL);

  /// Find closest point
  ///
  /// \param[in]  dataset Dataset whose nearest neighbor is queried.
  /// \param[in]  sample  Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in]  index   Index of point in \p dataset or \p sample (if not NULL).
  /// \param[out] dist2   Squared Euclidean distance of closest point.
  ///
  /// \returns Index of point/sample in dataset closest to the i-th point in \p dataset.
  int FindClosestPoint(vtkPointSet       *dataset,
                       const vector<int> *sample,
                       int                index,
                       double            *dist2 = NULL);

  /// Find closest point
  ///
  /// \param[in]  dataset Dataset whose nearest neighbor is queried.
  /// \param[in]  index   Index of point in \p dataset.
  /// \param[out] dist2   Squared Euclidean distance of closest point.
  ///
  /// \returns Index of point/sample in dataset closest to the i-th point in \p dataset.
  int FindClosestPoint(vtkPointSet *dataset,
                       int          index,
                       double      *dist2 = NULL);

  /// Find closest point
  ///
  /// \param[in]  dataset  Dataset whose nearest neighbor is queried.
  /// \param[in]  sample   Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in]  features Indices and weights of point data in \p dataset to use.
  /// \param[out] dist2    Squared Euclidean distance of closest point.
  ///
  /// \returns Indices of points/samples closest to each point in \p dataset.
  vector<int> FindClosestPoint(vtkPointSet       *dataset,
                               const vector<int> *sample,
                               const FeatureList *features,
                               vector<double>    *dist2 = NULL);

  /// Find closest point
  ///
  /// \param[in]  dataset  Dataset whose nearest neighbor is queried.
  /// \param[in]  features Indices and weights of point data in \p dataset to use.
  /// \param[out] dist2    Squared Euclidean distance of closest point.
  ///
  /// \returns Indices of points/samples closest to each point in \p dataset.
  vector<int> FindClosestPoint(vtkPointSet       *dataset,
                               const FeatureList *features,
                               vector<double>    *dist2 = NULL);

  /// Find closest point
  ///
  /// \param[in]  dataset Dataset whose nearest neighbor is queried.
  /// \param[in]  sample  Indices of points in \p dataset to consider only or NULL for all.
  /// \param[out] dist2   Squared Euclidean distance of closest point.
  ///
  /// \returns Indices of points/samples closest to each point in \p dataset.
  vector<int> FindClosestPoint(vtkPointSet       *dataset,
                               const vector<int> *sample,
                               vector<double>    *dist2 = NULL);

  /// Find closest point
  ///
  /// \param[in]  dataset Dataset whose nearest neighbor is queried.
  /// \param[out] dist2   Squared Euclidean distance of closest point.
  ///
  /// \returns Indices of points/samples closest to each point in \p dataset.
  vector<int> FindClosestPoint(vtkPointSet    *dataset,
                               vector<double> *dist2 = NULL);

  /// Find closest point
  ///
  /// \param[in]  dataset1  Dataset whose nearest neighbor is queried.
  /// \param[in]  sample1   Indices of points in \p dataset1 to consider only or NULL for all.
  /// \param[in]  features1 Indices and weights of point data in \p dataset1 to use.
  /// \param[in]  dataset2  Dataset in which to perform nearest neighbor search.
  /// \param[in]  sample2   Indices of points in \p dataset2 to consider only or NULL for all.
  /// \param[in]  features2 Indices and weights of point data in \p dataset2 to use.
  /// \param[out] dist2     Squared Euclidean distance of closest point.
  ///
  /// \returns Indices of points/samples in \p dataset2 which are closest to each point in \p dataset1.
  static vector<int> FindClosestPoint(vtkPointSet       *dataset1,
                                      const vector<int> *sample1,
                                      const FeatureList *features1,
                                      vtkPointSet       *dataset2,
                                      const vector<int> *sample2,
                                      const FeatureList *features2,
                                      vector<double>    *dist2 = NULL);

  // ---------------------------------------------------------------------------
  // Nearest neighbors (kNN)

  /// Find k nearest neighbors
  ///
  /// \param[in]  k       Number of nearest neighbors to find.
  /// \param[in]  dataset Dataset whose nearest neighbors are queried.
  /// \param[in]  sample  Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in]  index   Index of point in \p dataset or \p sample (if not NULL).
  /// \param[in]  feature Indices and weights of point data in \p dataset to use.
  /// \param[out] dist2   Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns Indices of points in dataset closests to the i-th point in \p d.
  vector<int> FindClosestNPoints(int k, double *point, vector<double> *dist2 = NULL);

  /// Find k nearest neighbors
  ///
  /// \param[in]  k        Number of nearest neighbors to find.
  /// \param[in]  dataset  Dataset whose nearest neighbors are queried.
  /// \param[in]  sample   Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in]  index    Index of point in \p dataset or \p sample (if not NULL).
  /// \param[in]  features Indices and weights of point data in \p dataset to use.
  /// \param[out] dist2    Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns Indices of points in dataset closests to the i-th point in \p d.
  vector<int> FindClosestNPoints(int k, vtkPointSet       *dataset,
                                        const vector<int> *sample,
                                        int                index,
                                        const FeatureList *features,
                                 vector<double> *dist2 = NULL);

  /// Find k nearest neighbors
  ///
  /// \param[in]  k        Number of nearest neighbors to find.
  /// \param[in]  dataset  Dataset whose nearest neighbors are queried.
  /// \param[in]  sample   Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in]  index    Index of point in \p dataset or \p sample (if not NULL).
  /// \param[out] dist2    Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns Indices of points in dataset closests to the i-th point in \p d.
  vector<int> FindClosestNPoints(int k, vtkPointSet       *dataset,
                                        const vector<int> *sample,
                                        int                index,
                                 vector<double> *dist2 = NULL);

  /// Find k nearest neighbors
  ///
  /// \param[in]  k        Number of nearest neighbors to find.
  /// \param[in]  dataset  Dataset whose nearest neighbors are queried.
  /// \param[in]  index    Index of point in \p dataset.
  /// \param[in]  features Indices and weights of point data in \p dataset to use.
  /// \param[out] dist2    Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns Indices of points in dataset closests to the i-th point in \p d.
  vector<int> FindClosestNPoints(int k, vtkPointSet       *dataset,
                                        int                index,
                                        const FeatureList *features,
                                 vector<double> *dist2 = NULL);

  /// Find k nearest neighbors
  ///
  /// \param[in]  k       Number of nearest neighbors to find.
  /// \param[in]  dataset Dataset whose nearest neighbors are queried.
  /// \param[in]  index   Index of point in \p dataset.
  /// \param[out] dist2   Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns Indices of points in dataset closests to the i-th point in \p d.
  vector<int> FindClosestNPoints(int k, vtkPointSet *dataset, int index,
                                 vector<double> *dist2 = NULL);

  /// Find k nearest neighbors
  ///
  /// \param[in]  k        Number of nearest neighbors to find.
  /// \param[in]  dataset  Dataset whose nearest neighbors are queried.
  /// \param[in]  sample   Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in]  features Indices and weights of point data in \p dataset to use.
  /// \param[out] dist2    Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns For each point in \p dataset, indices of closest \p k points.
  vector<vector<int> > FindClosestNPoints(int k, vtkPointSet       *dataset,
                                                 const vector<int> *sample,
                                                 const FeatureList *features,
                                          vector<vector<double> > *dist2 = NULL);

  /// Find k nearest neighbors
  ///
  /// \param[in]  k        Number of nearest neighbors to find.
  /// \param[in]  dataset  Dataset whose nearest neighbors are queried.
  /// \param[in]  sample   Indices of points in \p dataset to consider only or NULL for all.
  /// \param[out] dist2    Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns For each point in \p dataset, indices of closest \p k points.
  vector<vector<int> > FindClosestNPoints(int k, vtkPointSet       *dataset,
                                                 const vector<int> *sample,
                                          vector<vector<double> > *dist2 = NULL);

  /// Find k nearest neighbors
  ///
  /// \param[in]  k        Number of nearest neighbors to find.
  /// \param[in]  dataset  Dataset whose nearest neighbors are queried.
  /// \param[in]  index    Index of point in \p dataset.
  /// \param[in]  features Indices and weights of point data in \p dataset to use.
  /// \param[out] dist2    Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns For each point in \p dataset, indices of closest \p k points.
  vector<vector<int> > FindClosestNPoints(int k, vtkPointSet       *dataset,
                                                 const FeatureList *features,
                                          vector<vector<double> > *dist2 = NULL);

  /// Find k nearest neighbors
  ///
  /// \param[in]  k       Number of nearest neighbors to find.
  /// \param[in]  dataset Dataset whose nearest neighbors are queried.
  /// \param[out] dist2   Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns For each point in \p dataset, indices of closest \p k points.
  vector<vector<int> > FindClosestNPoints(int k, vtkPointSet *dataset,
                                          vector<vector<double> > *dist2 = NULL);

  /// Find k nearest neighbors
  ///
  /// \param[in]  k         Number of nearest neighbors to find.
  /// \param[in]  dataset1  Dataset whose nearest neighbor is queried.
  /// \param[in]  sample1   Indices of points in \p dataset1 to consider only or NULL for all.
  /// \param[in]  features1 Indices and weights of point data in \p dataset1 to use.
  /// \param[in]  dataset2  Dataset in which to perform nearest neighbor search.
  /// \param[in]  sample2   Indices of points in \p dataset2 to consider only or NULL for all.
  /// \param[in]  features2 Indices and weights of point data in \p dataset2 to use.
  /// \param[out] dist2     Squared Euclidean distance of closest point.
  ///
  /// \returns For each point in \p dataset1, indices of closest \p k points in \p dataset2.
  static vector<vector<int> > FindClosestNPoints(int                      k,
                                                 vtkPointSet             *dataset1,
                                                 const vector<int>       *sample1,
                                                 const FeatureList       *features1,
                                                 vtkPointSet             *dataset2,
                                                 const vector<int>       *sample2,
                                                 const FeatureList       *features2,
                                                 vector<vector<double> > *dist2 = NULL);

  // ---------------------------------------------------------------------------
  // Radius search

  /// Find neighbors within search radius
  ///
  /// \param[in]  radius  Search radius in N-D point feature space.
  /// \param[in]  dataset Dataset whose nearest neighbors are queried.
  /// \param[in]  sample  Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in]  index   Index of point in \p dataset or \p sample (if not NULL).
  /// \param[in]  feature Indices and weights of point data in \p dataset to use.
  /// \param[out] dist2   Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns Indices of points in dataset closests to the i-th point in \p d.
  vector<int> FindPointsWithinRadius(double radius, double *point, vector<double> *dist2 = NULL);

  /// Find neighbors within search radius
  ///
  /// \param[in]  radius   Search radius in N-D point feature space.
  /// \param[in]  dataset  Dataset whose nearest neighbors are queried.
  /// \param[in]  sample   Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in]  index    Index of point in \p dataset or \p sample (if not NULL).
  /// \param[in]  features Indices and weights of point data in \p dataset to use.
  /// \param[out] dist2    Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns Indices of points in dataset closests to the i-th point in \p d.
  vector<int> FindPointsWithinRadius(double radius, vtkPointSet       *dataset,
                                                    const vector<int> *sample,
                                                    int                index,
                                                    const FeatureList *features,
                                     vector<double> *dist2 = NULL);

  /// Find neighbors within search radius
  ///
  /// \param[in]  radius  Search radius in N-D point feature space.
  /// \param[in]  dataset Dataset whose nearest neighbors are queried.
  /// \param[in]  sample  Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in]  index   Index of point in \p dataset or \p sample (if not NULL).
  /// \param[out] dist2   Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns Indices of points in dataset closests to the i-th point in \p d.
  vector<int> FindPointsWithinRadius(double radius, vtkPointSet       *dataset,
                                                    const vector<int> *sample,
                                                    int                index,
                                     vector<double> *dist2 = NULL);

  /// Find neighbors within search radius
  ///
  /// \param[in]  radius   Search radius in N-D point feature space.
  /// \param[in]  dataset  Dataset whose nearest neighbors are queried.
  /// \param[in]  index    Index of point in \p dataset.
  /// \param[in]  features Indices and weights of point data in \p dataset to use.
  /// \param[out] dist2    Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns Indices of points in dataset closests to the i-th point in \p d.
  vector<int> FindPointsWithinRadius(double radius, vtkPointSet       *dataset,
                                                    int                index,
                                                    const FeatureList *features,
                                     vector<double> *dist2 = NULL);

  /// Find neighbors within search radius
  ///
  /// \param[in]  radius  Search radius in N-D point feature space.
  /// \param[in]  dataset Dataset whose nearest neighbors are queried.
  /// \param[in]  index   Index of point in \p dataset.
  /// \param[out] dist2   Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns Indices of points in dataset closests to the i-th point in \p d.
  vector<int> FindPointsWithinRadius(double radius, vtkPointSet *dataset, int index,
                                     vector<double> *dist2 = NULL);

  /// Find neighbors within search radius
  ///
  /// \param[in]  radius   Search radius in N-D point feature space.
  /// \param[in]  dataset  Dataset whose nearest neighbors are queried.
  /// \param[in]  sample   Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in]  features Indices and weights of point data in \p dataset to use.
  /// \param[out] dist2    Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns For each point in \p dataset, indices of points within search \p radius.
  vector<vector<int> > FindPointsWithinRadius(double radius, vtkPointSet       *dataset,
                                                             const vector<int> *sample,
                                                             const FeatureList *features,
                                              vector<vector<double> > *dist2 = NULL);

  /// Find neighbors within search radius
  ///
  /// \param[in]  radius  Search radius in N-D point feature space.
  /// \param[in]  dataset Dataset whose nearest neighbors are queried.
  /// \param[in]  sample  Indices of points in \p dataset to consider only or NULL for all.
  /// \param[out] dist2   Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns For each point in \p dataset, indices of points within search \p radius.
  vector<vector<int> > FindPointsWithinRadius(double radius, vtkPointSet       *dataset,
                                                             const vector<int> *sample,
                                              vector<vector<double> > *dist2 = NULL);

  /// Find neighbors within search radius
  ///
  /// \param[in]  radius   Search radius in N-D point feature space.
  /// \param[in]  dataset  Dataset whose nearest neighbors are queried.
  /// \param[in]  features Indices and weights of point data in \p dataset to use.
  /// \param[out] dist2    Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns For each point in \p dataset, indices of points within search \p radius.
  vector<vector<int> > FindPointsWithinRadius(double radius, vtkPointSet       *dataset,
                                                             const FeatureList *features,
                                              vector<vector<double> > *dist2 = NULL);

  /// Find neighbors within search radius
  ///
  /// \param[in]  radius  Search radius in N-D point feature space.
  /// \param[in]  dataset Dataset whose nearest neighbors are queried.
  /// \param[out] dist2   Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns For each point in \p dataset, indices of points within search \p radius.
  vector<vector<int> > FindPointsWithinRadius(double radius, vtkPointSet *dataset,
                                              vector<vector<double> > *dist2 = NULL);

  /// Find neighbors within search radius
  ///
  /// \param[in]  radius    Search radius in N-D point feature space.
  /// \param[in]  dataset1  Dataset whose nearest neighbor is queried.
  /// \param[in]  sample1   Indices of points in \p dataset1 to consider only or NULL for all.
  /// \param[in]  features1 Indices and weights of point data in \p dataset1 to use.
  /// \param[in]  dataset2  Dataset in which to perform nearest neighbor search.
  /// \param[in]  sample2   Indices of points in \p dataset2 to consider only or NULL for all.
  /// \param[in]  features2 Indices and weights of point data in \p dataset2 to use.
  /// \param[out] dist2     Squared Euclidean distance of closest point.
  ///
  /// \returns For each point in \p dataset1, indices of points in \p dataset2 within search \p radius.
  static vector<vector<int> > FindPointsWithinRadius(double                   radius,
                                                     vtkPointSet             *dataset1,
                                                     const vector<int>       *sample1,
                                                     const FeatureList       *features1,
                                                     vtkPointSet             *dataset2,
                                                     const vector<int>       *sample2,
                                                     const FeatureList       *features2,
                                                     vector<vector<double> > *dist2 = NULL);

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Static helpers
// =============================================================================

// -----------------------------------------------------------------------------
inline vtkDataArray *irtkPointLocator::GetDataArray(vtkPointSet *dataset, const FeatureInfo &feature)
{
  vtkPointData * const pd = dataset->GetPointData();
  if (feature._Index >= 0) return pd->GetArray(feature._Index);
  return pd->GetArray(feature._Name.c_str());
}

// -----------------------------------------------------------------------------
inline int irtkPointLocator::GetNumberOfPoints(vtkPointSet *dataset, const vector<int> *sample)
{
  return (sample && sample->size() > 0 ? sample->size() : dataset->GetNumberOfPoints());
}

// -----------------------------------------------------------------------------
inline int irtkPointLocator::GetPointDimension(vtkPointSet *dataset, const FeatureList *features)
{
  int dim;
  if (features && features->size() > 0) {
    dim = 0;
    vtkDataArray *feature_array;
    for (FeatureList::const_iterator feature = features->begin(); feature != features->end(); ++feature) {
      if (feature->_Weight != .0) {
        if (feature->_Index == -1) {
          dim += 3;
        } else {
          feature_array = GetDataArray(dataset, *feature);
          dim += feature_array->GetNumberOfComponents();
        }
      }
    }
  } else {
    dim = 3;
  }
  return dim;
}

// -----------------------------------------------------------------------------
inline int irtkPointLocator::GetPointIndex(vtkPointSet *dataset, const vector<int> *sample, int index)
{
  return (sample && sample->size() > 0 ? (*sample)[index] : index);
}

// -----------------------------------------------------------------------------
inline void irtkPointLocator
::GetPoint(irtkPoint &point, vtkPointSet *dataset, const vector<int> *sample, int index)
{
  double pt[3];
  const int i = GetPointIndex(dataset, sample, index);
  dataset->GetPoint(i, pt);
  point._x = pt[0], point._y = pt[1], point._z = pt[2];
}

// -----------------------------------------------------------------------------
inline void irtkPointLocator
::GetPoint(double *point, vtkPointSet *dataset, const vector<int> *sample, int index, const FeatureList *features)
{
  const int i = GetPointIndex(dataset, sample, index);
  if (features && features->size() > 0) {
    vtkDataArray *feature_array;
    for (FeatureList::const_iterator feature = features->begin(); feature != features->end(); ++feature) {
      if (feature->_Weight != .0) {
        if (feature->_Index == -1) {
          dataset->GetPoint(i, point);
          for (int d = 0; d < 3; ++d, ++point) {
            (*point) = feature->_Weight * (feature->_Slope * (*point) + feature->_Intercept);
          }
        } else {
          feature_array = GetDataArray(dataset, *feature);
          feature_array->GetTuple(i, point);
          for (int d = 0; d < feature_array->GetNumberOfComponents(); ++d, ++point) {
            (*point) = feature->_Weight * (feature->_Slope * (*point) + feature->_Intercept);
          }
        }
      }
    }
  } else {
    dataset->GetPoint(i, point);
  }
}

// -----------------------------------------------------------------------------
inline void irtkPointLocator
::GetPoint(double *point, vtkPointSet *dataset, int index, const FeatureList *features)
{
  GetPoint(point, dataset, NULL, index, features);
}

// -----------------------------------------------------------------------------
inline double irtkPointLocator
::Distance2BetweenPoints(const double *a, const double *b, int d)
{
  double dx, dist2 = .0;
  for (int i = 0; i < d; ++i) {
    dx = b[i] - a[i];
    dist2 += dx * dx;
  }
  return dist2;
}

// =============================================================================
// Closest point
// =============================================================================

// -----------------------------------------------------------------------------
inline int irtkPointLocator
::FindClosestPoint(vtkPointSet *dataset, const vector<int> *sample, int index,
                   const FeatureList *features, double *dist2)
{
  double *point = new double[_PointDimension];
  GetPoint(point, dataset, sample, index, features);
  int i = this->FindClosestPoint(point, dist2);
  delete[] point;
  return i;
}

// -----------------------------------------------------------------------------
inline int irtkPointLocator
::FindClosestPoint(vtkPointSet *dataset, int index,
                   const FeatureList *features, double *dist2)
{
  return FindClosestPoint(dataset, NULL, index, features, dist2);
}

// -----------------------------------------------------------------------------
inline int irtkPointLocator
::FindClosestPoint(vtkPointSet *dataset, const vector<int> *sample,
                   int index, double *dist2)
{
  return FindClosestPoint(dataset, sample, index, NULL, dist2);
}

// -----------------------------------------------------------------------------
inline int irtkPointLocator
::FindClosestPoint(vtkPointSet *dataset, int index, double *dist2)
{
  return FindClosestPoint(dataset, NULL, index, dist2);
}

// -----------------------------------------------------------------------------
inline vector<int> irtkPointLocator
::FindClosestPoint(vtkPointSet *dataset, const FeatureList *features, vector<double> *dist2)
{
  return FindClosestPoint(dataset, NULL, features, dist2);
}

// -----------------------------------------------------------------------------
inline vector<int> irtkPointLocator
::FindClosestPoint(vtkPointSet *dataset, const vector<int> *sample, vector<double> *dist2)
{
  return FindClosestPoint(dataset, sample, NULL, dist2);
}

// -----------------------------------------------------------------------------
inline vector<int> irtkPointLocator
::FindClosestPoint(vtkPointSet *dataset, vector<double> *dist2)
{
  return FindClosestPoint(dataset, NULL, NULL, dist2);
}

// -----------------------------------------------------------------------------
inline vector<int> irtkPointLocator
::FindClosestPoint(vtkPointSet *dataset1, const vector<int> *sample1, const FeatureList *features1,
                   vtkPointSet *dataset2, const vector<int> *sample2, const FeatureList *features2,
                   vector<double> *dist2)
{
  auto_ptr<irtkPointLocator> locator(irtkPointLocator::New(dataset2, sample2, features2));
  return locator->FindClosestPoint(dataset1, sample1, features1, dist2);
}

// =============================================================================
// Nearest neighbors (kNN)
// =============================================================================

// -----------------------------------------------------------------------------
inline vector<int> irtkPointLocator
::FindClosestNPoints(int k, vtkPointSet *dataset, const vector<int> *sample, int index,
                     const FeatureList *features, vector<double> *dist2)
{
  double *point = new double[_PointDimension];
  GetPoint(point, dataset, sample, index, features);
  vector<int> indices = this->FindClosestNPoints(k, point, dist2);
  delete[] point;
  return indices;
}

// -----------------------------------------------------------------------------
inline vector<int> irtkPointLocator
::FindClosestNPoints(int k, vtkPointSet *dataset, int index,
                     const FeatureList *features, vector<double> *dist2)
{
  return FindClosestNPoints(k, dataset, NULL, index, features, dist2);
}

// -----------------------------------------------------------------------------
inline vector<int> irtkPointLocator
::FindClosestNPoints(int k, vtkPointSet *dataset, const vector<int> *sample,
                     int index, vector<double> *dist2)
{
  return FindClosestNPoints(k, dataset, sample, index, NULL, dist2);
}

// -----------------------------------------------------------------------------
inline vector<int> irtkPointLocator
::FindClosestNPoints(int k, vtkPointSet *dataset, int index, vector<double> *dist2)
{
  return FindClosestNPoints(k, dataset, NULL, index, dist2);
}

// -----------------------------------------------------------------------------
inline vector<vector<int> > irtkPointLocator
::FindClosestNPoints(int k, vtkPointSet *dataset, const FeatureList *features, vector<vector<double> > *dist2)
{
  return FindClosestNPoints(k, dataset, NULL, features, dist2);
}

// -----------------------------------------------------------------------------
inline vector<vector<int> > irtkPointLocator
::FindClosestNPoints(int k, vtkPointSet *dataset, const vector<int> *sample, vector<vector<double> > *dist2)
{
  return FindClosestNPoints(k, dataset, sample, NULL, dist2);
}

// -----------------------------------------------------------------------------
inline vector<vector<int> > irtkPointLocator
::FindClosestNPoints(int k, vtkPointSet *dataset, vector<vector<double> > *dist2)
{
  return FindClosestNPoints(k, dataset, NULL, NULL, dist2);
}

// -----------------------------------------------------------------------------
inline vector<vector<int> > irtkPointLocator
::FindClosestNPoints(int k, vtkPointSet *dataset1, const vector<int> *sample1, const FeatureList *features1,
                            vtkPointSet *dataset2, const vector<int> *sample2, const FeatureList *features2,
                     vector<vector<double> > *dist2)
{
  auto_ptr<irtkPointLocator> locator(irtkPointLocator::New(dataset2, sample2, features2));
  return locator->FindClosestNPoints(k, dataset1, sample1, features1, dist2);
}

// =============================================================================
// Radius search
// =============================================================================

// -----------------------------------------------------------------------------
inline vector<int> irtkPointLocator
::FindPointsWithinRadius(double radius, vtkPointSet *dataset, const vector<int> *sample,
                         int index, const FeatureList *features, vector<double> *dist2)
{
  double *point = new double[_PointDimension];
  GetPoint(point, dataset, sample, index, features);
  vector<int> indices = this->FindPointsWithinRadius(radius, point, dist2);
  delete[] point;
  return indices;
}

// -----------------------------------------------------------------------------
inline vector<int> irtkPointLocator
::FindPointsWithinRadius(double radius, vtkPointSet *dataset, int index,
                         const FeatureList *features, vector<double> *dist2)
{
  return FindPointsWithinRadius(radius, dataset, NULL, index, features, dist2);
}

// -----------------------------------------------------------------------------
inline vector<int> irtkPointLocator
::FindPointsWithinRadius(double radius, vtkPointSet *dataset, const vector<int> *sample,
                         int index, vector<double> *dist2)
{
  return FindPointsWithinRadius(radius, dataset, sample, index, NULL, dist2);
}

// -----------------------------------------------------------------------------
inline vector<int> irtkPointLocator
::FindPointsWithinRadius(double radius, vtkPointSet *dataset, int index, vector<double> *dist2)
{
  return FindPointsWithinRadius(radius, dataset, NULL, index, dist2);
}

// -----------------------------------------------------------------------------
inline vector<vector<int> > irtkPointLocator
::FindPointsWithinRadius(double radius, vtkPointSet *dataset,
                         const FeatureList *features, vector<vector<double> > *dist2)
{
  return FindPointsWithinRadius(radius, dataset, NULL, features, dist2);
}

// -----------------------------------------------------------------------------
inline vector<vector<int> > irtkPointLocator
::FindPointsWithinRadius(double radius, vtkPointSet *dataset, const vector<int> *sample, vector<vector<double> > *dist2)
{
  return FindPointsWithinRadius(radius, dataset, sample, NULL, dist2);
}

// -----------------------------------------------------------------------------
inline vector<vector<int> > irtkPointLocator
::FindPointsWithinRadius(double radius, vtkPointSet *dataset, vector<vector<double> > *dist2)
{
  return FindPointsWithinRadius(radius, dataset, NULL, NULL, dist2);
}

// -----------------------------------------------------------------------------
inline vector<vector<int> > irtkPointLocator
::FindPointsWithinRadius(double radius, vtkPointSet *dataset1, const vector<int> *sample1, const FeatureList *features1,
                                        vtkPointSet *dataset2, const vector<int> *sample2, const FeatureList *features2,
                         vector<vector<double> > *dist2)
{
  auto_ptr<irtkPointLocator> locator(irtkPointLocator::New(dataset2, sample2, features2));
  return locator->FindPointsWithinRadius(radius, dataset1, sample1, features1, dist2);
}


#endif
