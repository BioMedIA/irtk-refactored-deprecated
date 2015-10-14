// =============================================================================
// Project: Image Registration Toolkit (IRTK)
// Package: PolyData
//
// Copyright (c) 2015 Imperial College London
// Copyright (c) 2015 Andreas Schuh
// =============================================================================

#include <irtkCommon.h>


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
/// Print help screen
void PrintHelp(const char *name)
{
  cout << "usage: " << name << " <target> <source> [options]" << endl;
  cout << endl;
  cout << "Reports label statistics given two vtkPolyData files and a text file" << endl;
  cout << "listing the indices of corresponding points. The labels have to be" << endl;
  cout << "stored as named point data array in both datasets." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  target   File name of target dataset." << endl;
  cout << "  source   File name of source dataset." << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << "  -corr <file>          Name of text file listing the point-wise correspondences separated by spaces. (default: NN)" << endl;
  cout << "  -array       <name>   Name of array storing labels. (default: labels)" << endl;
  cout << "  -targetarray <name>   Name of target array storing labels. (default: labels)" << endl;
  cout << "  -sourcearray <name>   Name of source array storing labels. (default: labels)" << endl;
  cout << "  -sep <str>            Separating string between printed numbers. (default: ,)" << endl;
}

// =============================================================================
// Includes
// =============================================================================

#include <vtkCell.h>
#include <vtkTriangle.h>
#include <vtkTriangleFilter.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkKdTree.h>
#include <vtkKdTreePointLocator.h>

#ifdef HAS_FLANN
#  include <flann/flann.hpp>
#endif

// =============================================================================
// Types
// =============================================================================

#ifdef HAS_FLANN
typedef float                                      FlannType;
typedef flann::Matrix<FlannType>                   FlannMatrix;
typedef flann::Index<flann::L2_Simple<FlannType> > FlannIndex;
#endif

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
/// Read correspondences from text file
vector<vtkIdType>
ReadCorrespondences(const char *fname, vtkPolyData *target, vtkPolyData *source, vtkIdType n)
{
  ifstream ifs(fname);
  if (!ifs.is_open()) {
    cerr << "Error: Failed to open correspondence input file: " << fname << endl;
    exit(1);
  }
  string    line;
  vtkIdType t, s;
  vector<vtkIdType> corr(n, -1);
  while (getline(ifs, line)) {
    istringstream is(line);
    if (!(is >> t >> s)) {
      cerr << "Error: Failed to read correspondence map from file: " << fname << endl;
      exit(1);
    }
    if (t < 0 || t >= n) {
      cerr << "Error: Invalid target index (" << t << ") in correspondence map: " << fname << endl;
      exit(1);
    }
    if (s < 0 || s >= n) {
      cerr << "Error: Invalid source index (" << t << ") in correspondence map: " << fname << endl;
      exit(1);
    }
    corr[t] = s;
  }
  ifs.close();
  return corr;
}

// -----------------------------------------------------------------------------
/// Get indices of data arrays using case insensitive name
vector<vtkDataArray *>
GetArraysByCaseInsensitiveName(vtkDataSetAttributes *data, const vector<string> &name)
{
  vector<vtkDataArray *> array(name.size(), NULL);
  for (size_t i = 0; i < name.size(); ++i) {
    array[i] = GetArrayByCaseInsensitiveName(data, name[i].c_str());
  }
  return array;
}

// -----------------------------------------------------------------------------
/// Get pointers to point data arrays used as features to establish correspondences
///
/// @returns Pointers to point data arrays of point set. If an array is \c NULL,
///          the feature corresponds to the spatial point coordinates instead.
vector<vtkDataArray *> GetPointFeatureArrays(vtkPolyData *polydata, const vector<string> &name)
{
  vector<vtkDataArray *> array = GetArraysByCaseInsensitiveName(polydata->GetPointData(), name);
  for (size_t i = 0; i < array.size(); ++i) {
    if (array[i] == NULL) {
      if (name[i] != "spatial coordinates" &&
          name[i] != "spatial points" && name[i] != "points") {
        cerr << "Missing target point feature array named " << name[i] << endl;
        exit(1);
      }
    }
  }
  return array;
}

// -----------------------------------------------------------------------------
/// Get dimension of feature points used to establish nearest neighbor correspondences
inline int GetFeaturePointDimension(const vector<vtkDataArray *> &feature,
                                    const vector<double>         &weight,
                                    const vector<int>            &rows)
{
  int dim = 0;
  for (size_t i = 0; i < feature.size(); ++i) {
    if (!fequal(weight[i], .0)) {
      if (feature[i]) {
        if (0 <= rows[i] && rows[i] < feature[i]->GetNumberOfComponents()) dim += rows[i];
        else dim += feature[i]->GetNumberOfComponents();
      } else dim += 3; // spatial coordinates
    }
  }
  return dim;
}

// -----------------------------------------------------------------------------
/// Get weighted coordinates of feature point used to establish correspondences
inline void GetFeaturePoint(vtkIdType id, vtkPoints *points,
                            const vector<vtkDataArray *> &feature,
                            const vector<double>         &weight,
                            const vector<int>            &rows,
                            double                       *p)
{
  for (size_t i = 0; i < feature.size(); ++i) {
    if (!fequal(weight[i], .0)) {
      if (feature[i]) {
        if (rows[i] < 0 || rows[i] >= feature[i]->GetNumberOfComponents()) {
          feature[i]->GetTuple(id, p);
          for (int j = 0; j < feature[i]->GetNumberOfComponents(); ++j) {
            p[j] *= weight[i];
          }
          p += feature[i]->GetNumberOfComponents();
        } else {
          for (int j = 0; j < rows[i]; ++j) {
            p[j] = weight[i] * feature[i]->GetComponent(id, j);
          }
          p += rows[i];
        }
      } else {
        points->GetPoint(id, p);
        p[0] *= weight[i], p[1] *= weight[i], p[2] *= weight[i];
        p += 3;
      }
    }
  }
}

// -----------------------------------------------------------------------------
/// Get cell centers
vtkSmartPointer<vtkPoints> GetCellCenters(vtkPolyData *polydata)
{
  double  *weights = new double[polydata->GetMaxCellSize()];
  double   pcoords[3], p[3];
  int      subId;
  vtkCell *cell;

  vtkSmartPointer<vtkPoints> centers = vtkSmartPointer<vtkPoints>::New();
  centers->SetNumberOfPoints(polydata->GetNumberOfCells());
  for (vtkIdType cellId = 0; cellId < polydata->GetNumberOfCells(); ++cellId) {
    cell  = polydata->GetCell(cellId);
    subId = cell->GetParametricCenter(pcoords);
    cell->EvaluateLocation(subId, pcoords, p, weights);
    centers->SetPoint(cellId, p);
  }

  delete[] weights;
  return centers;
}

#ifdef HAS_FLANN
// -----------------------------------------------------------------------------
/// Get new FLANN matrix with feature points at data set points
FlannMatrix FlannFeaturesAtDataSetPoints(vtkPolyData                  *polydata,
                                         const vector<vtkDataArray *> &feature,
                                         const vector<double>         &weight,
                                         const vector<int>            &rows)
{
  const size_t num = polydata->GetNumberOfPoints();
  const size_t dim = GetFeaturePointDimension(feature, weight, rows);
  double *point = new double[dim];
  FlannMatrix m(new FlannType[num * dim], num, dim);
  for (vtkIdType ptId = 0; ptId < polydata->GetNumberOfPoints(); ++ptId) {
    GetFeaturePoint(ptId, polydata->GetPoints(), feature, weight, rows, point);
    FlannType *row = m[ptId];
    for (size_t d = 0; d < dim; ++d) {
      row[d] = static_cast<FlannType>(point[d]);
    }
  }
  delete[] point;
  return m;
}

// -----------------------------------------------------------------------------
/// Get new FLANN matrix with feature points at cell centers
FlannMatrix FlannFeaturesAtCellCenters(vtkPolyData                  *polydata,
                                       const vector<vtkDataArray *> &feature,
                                       const vector<double>         &weight,
                                       const vector<int>            &rows)
{
  vtkCell   *cell;
  vtkIdList *ptIds;
  double    *ptWeights = new double[polydata->GetMaxCellSize()];
  double     pcoords[3], p[3];
  int        subId;

  const int num = polydata->GetNumberOfCells();
  const int dim = GetFeaturePointDimension(feature, weight, rows);
  double *point = new double[dim];
  FlannMatrix m(new FlannType[num * dim], num, dim);
  memset(m.ptr(), 0, num * dim * sizeof(FlannType));

  for (vtkIdType cellId = 0; cellId < polydata->GetNumberOfCells(); ++cellId) {
    cell  = polydata->GetCell(cellId);
    ptIds = cell->GetPointIds();
    subId = cell->GetParametricCenter(pcoords);
    cell->EvaluateLocation(subId, pcoords, p, ptWeights);
    for (vtkIdType i = 0; i < ptIds->GetNumberOfIds(); ++i) {
      GetFeaturePoint(ptIds->GetId(i), polydata->GetPoints(), feature, weight, rows, point);
      FlannType *row = m[cellId];
      for (int d = 0; d < dim; ++d) {
        row[d] += static_cast<FlannType>(ptWeights[i] * point[d]);
      }
    }
  }

  delete[] point;
  delete[] ptWeights;
  return m;
}
#endif

// -----------------------------------------------------------------------------
/// Find source points closest to target points
vector<vtkIdType>
FindClosestPoints(vtkPolyData *target, vtkPolyData *source,
                  const vector<string> &feature,
                  const vector<double> &weight,
                  const vector<int>    &rows)
{
  vector<vtkIdType> corr(target->GetNumberOfPoints());
  if (feature.empty()) {
    vtkSmartPointer<vtkKdTreePointLocator> kdtree;
    kdtree = vtkSmartPointer<vtkKdTreePointLocator>::New();
    kdtree->SetDataSet(source);
    kdtree->BuildLocator();
    for (vtkIdType i = 0; i < target->GetNumberOfPoints(); ++i) {
      corr[i] = kdtree->FindClosestPoint(target->GetPoint(i));
    }
  } else {
    vector<vtkDataArray *> target_feature = GetPointFeatureArrays(target, feature);
    vector<vtkDataArray *> source_feature = GetPointFeatureArrays(source, feature);
#ifdef HAS_FLANN
    FlannIndex  index(flann::KDTreeSingleIndexParams(10));
    FlannMatrix points = FlannFeaturesAtDataSetPoints(source, source_feature, weight, rows);
    index.buildIndex(points);
    delete[] points.ptr();
    FlannMatrix queries = FlannFeaturesAtDataSetPoints(target, target_feature, weight, rows);
    vector<vector<int> >       indices;
    vector<vector<FlannType> > dists;
    index.knnSearch(queries, indices, dists, 1,
                    flann::SearchParams(flann::FLANN_CHECKS_UNLIMITED));
    for (size_t i = 0; i < indices.size(); ++i) {
      corr[i] = static_cast<vtkIdType>(indices[i][0]);
    }
    delete[] queries.ptr();
#else
    cerr << "The -feature option for establishing point correspondences is only implemented using FLANN" << endl;
    exit(1);
#endif
  }
  return corr;
}

// -----------------------------------------------------------------------------
/// Find source cells whose centers are closest to the target cell centers
vector<vtkIdType> FindClosestCells(vtkPolyData *target, vtkPolyData *source,
                                   const vector<string> &feature,
                                   const vector<double> &weight,
                                   const vector<int>    &rows)
{
  vector<vtkIdType> corr(target->GetNumberOfCells(), -1);

  if (feature.empty()) {
    vtkSmartPointer<vtkPoints> target_centers = GetCellCenters(target);
    vtkSmartPointer<vtkPoints> source_centers = GetCellCenters(source);
    vtkSmartPointer<vtkKdTree> locator = vtkSmartPointer<vtkKdTree>::New();
    locator->BuildLocatorFromPoints(source_centers);
    double dist2;
    for (vtkIdType i = 0; i < target_centers->GetNumberOfPoints(); ++i) {
      corr[i] = locator->FindClosestPoint(target_centers->GetPoint(i), dist2);
    }
  } else {
    vector<vtkDataArray *> target_feature = GetPointFeatureArrays(target, feature);
    vector<vtkDataArray *> source_feature = GetPointFeatureArrays(source, feature);
#ifdef HAS_FLANN
    FlannIndex  index(flann::KDTreeSingleIndexParams(10));
    FlannMatrix points = FlannFeaturesAtCellCenters(source, source_feature, weight, rows);
    index.buildIndex(points);
    delete[] points.ptr();
    FlannMatrix queries = FlannFeaturesAtCellCenters(target, target_feature, weight, rows);
    vector<vector<int> >       indices;
    vector<vector<FlannType> > dists;
    index.knnSearch(queries, indices, dists, 1,
                    flann::SearchParams(flann::FLANN_CHECKS_UNLIMITED));
    for (size_t i = 0; i < indices.size(); ++i) {
      corr[i] = static_cast<vtkIdType>(indices[i][0]);
    }
    delete[] queries.ptr();
#else
    cerr << "The -feature option for establishing cell correspondences is only implemented using FLANN" << endl;
    exit(1);
#endif
  }

  return corr;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  // Positional arguments
  EXPECTS_POSARGS(2);
  const char *target_name = POSARG(1);
  const char *source_name = POSARG(2);

  // Optional arguments
  const char *corrin_name    = NULL;
  const char *target_pdarray = NULL;
  const char *target_cdarray = NULL;
  const char *source_pdarray = NULL;
  const char *source_cdarray = NULL;
  const char *sep            = ",";

  vector<string> feature_name;
  vector<double> feature_weight;
  vector<int>    feature_rows;

  for (ALL_OPTIONS) {
    if      (OPTION("-corrin")) corrin_name = ARGUMENT;
    else if (OPTION("-pointdata") || OPTION("-pd")) {
      target_pdarray = source_pdarray = ARGUMENT;
    }
    else if (OPTION("-target-pointdata") || OPTION("-Tpd")) target_pdarray = ARGUMENT;
    else if (OPTION("-source-pointdata") || OPTION("-Spd")) source_pdarray = ARGUMENT;
    else if (OPTION("-celldata") || OPTION("-cd")) {
      target_cdarray = source_cdarray = ARGUMENT;
    }
    else if (OPTION("-target-celldata") || OPTION("-Tcd")) target_cdarray = ARGUMENT;
    else if (OPTION("-source-celldata") || OPTION("-Scd")) source_cdarray = ARGUMENT;
    else if (OPTION("-feature") || OPTION("-f")) {
      feature_name.push_back(ARGUMENT);
      if (HAS_ARGUMENT) feature_weight.push_back(atof(ARGUMENT));
      else              feature_weight.push_back(1.0);
      if (HAS_ARGUMENT) feature_rows.push_back(atoi(ARGUMENT));
      else              feature_rows.push_back(-1);
    }
    else if (OPTION("-sep")) sep = ARGUMENT;
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // Check that either point data only or cell data only is used
  if ((target_pdarray && target_cdarray) || (source_pdarray && source_cdarray)) {
    cerr << "Can only use either point data or cell data labels, not both!" << endl;
    exit(1);
  }
  if ((target_pdarray && source_cdarray) || (target_cdarray && source_pdarray)) {
    cerr << "Can only use either point data or cell data labels, not a combination of these!" << endl;
    exit(1);
  }

  // Set default array names
  const char *default_array = "labels";
  if      ( target_pdarray && !source_pdarray) source_pdarray = default_array;
  else if ( target_cdarray && !source_cdarray) source_cdarray = default_array;
  else if (!target_pdarray &&  source_pdarray) target_pdarray = default_array;
  else if (!target_cdarray &&  source_cdarray) target_cdarray = default_array;

  // If neither was specified, use point data array
  if (!target_pdarray && !target_cdarray) {
    target_pdarray = source_pdarray = default_array;
  }

  // Read input files
  if (verbose > 1) cout << "Reading target from " << target_name << endl;
  vtkSmartPointer<vtkPolyData> target = ReadPolyData(target_name);
  if (verbose > 1) cout << "Reading source from " << source_name << endl;
  vtkSmartPointer<vtkPolyData> source = ReadPolyData(source_name);

  // Triangulate surfaces if cell overlap is evaluated
  if (target_cdarray) {
    target = Triangulate(target);
    source = Triangulate(source);
  }

  // Establish correspondences
  vector<vtkIdType> corr;
  if (target_pdarray) {
    if (corrin_name) {
      if (verbose > 1) cout << "Reading point correspondences from " << corrin_name << endl;
      corr = ReadCorrespondences(corrin_name, target, source, target->GetNumberOfPoints());
    } else {
      if (verbose > 1) cout << "Using closest points as correspondences" << endl;
      corr = FindClosestPoints(target, source, feature_name, feature_weight, feature_rows);
    }
  } else {
    if (corrin_name) {
      if (verbose > 1) cout << "Reading cell correspondences from " << corrin_name << endl;
      corr = ReadCorrespondences(corrin_name, target, source, target->GetNumberOfCells());
    } else {
      if (verbose > 1) cout << "Using closest cell centers as correspondences" << endl;
      corr = FindClosestCells(target, source, feature_name, feature_weight, feature_rows);
    }
  }

  // Get data arrays storing labels
  vtkDataArray *target_labels = NULL;
  vtkDataArray *source_labels = NULL;
  if (target_pdarray) {
    target_labels = GetArrayByCaseInsensitiveName(target->GetPointData(), target_pdarray);
    source_labels = GetArrayByCaseInsensitiveName(source->GetPointData(), source_pdarray);
  } else {
    target_labels = GetArrayByCaseInsensitiveName(target->GetCellData(), target_cdarray);
    source_labels = GetArrayByCaseInsensitiveName(source->GetCellData(), source_cdarray);
  }
  if (!target_labels) {
    if (target_pdarray) {
      cerr << "Error: Missing target point data array named " << target_pdarray << " (case insensitive)" << endl;
    } else {
      cerr << "Error: Missing target cell data array named " << target_cdarray << " (case insensitive)" << endl;
    }
    exit(1);
  }
  if (!source_labels) {
    if (source_pdarray) {
      cerr << "Error: Missing source data array named " << source_pdarray << " (case insensitive)" << endl;
    } else {
      cerr << "Error: Missing source data array named " << source_cdarray << " (case insensitive)" << endl;
    }
    exit(1);
  }

  // Determine maximum possible label
  double range[2];
  target_labels->GetRange(range);
  int max_label = static_cast<int>(range[1]);
  source_labels->GetRange(range);
  if (range[1] > static_cast<double>(max_label)) {
    max_label = static_cast<int>(range[1]);
  }

  // Count label occurrences
  vector<size_t> target_hist(max_label + 1, 0);
  vector<size_t> source_hist(max_label + 1, 0);
  vector<size_t> joint_hist (max_label + 1, 0);
  int t, s;
  for (size_t i = 0; i < corr.size(); ++i) {
    const vtkIdType &j = corr[i];
    t = static_cast<int>(round(target_labels->GetTuple1(i)));
    s = static_cast<int>(round(source_labels->GetTuple1(j)));
    ++target_hist[t];
    ++source_hist[s];
    if (t == s) ++joint_hist[t];
  }

  // Print label statistics
  if (verbose > 1) cout << endl;
  double overlap, dice, total_overlap = .0, total_dice = .0;
  int    num_labels = 0;

  for (size_t l = 0; l <= max_label; ++l) {
    // Skip unused labels
    if (target_hist[l] == 0 && source_hist[l] == 0) continue;
    ++num_labels;
    // Calculate overlap
    overlap = double(joint_hist[l]) / double(target_hist[l] + source_hist[l] - joint_hist[l]);
    dice    = 2.0 * double(joint_hist[l]) / double(target_hist[l] + source_hist[l]);
    total_overlap += overlap;
    total_dice    += dice;
    // Print label overlap
    if (verbose) {
      cout << l << sep << target_hist[l] << sep << source_hist[l] << sep << joint_hist[l]
                << sep << overlap << sep << dice << endl;
    }
  }

  // Print summary
  if (num_labels > 0) {
    double mean_overlap = total_overlap / num_labels;
    double mean_dice    = total_dice    / num_labels;
    if (verbose) cout << endl;
    cout << "Mean overlap ratio: " << mean_overlap << endl;
    cout << "Mean dice overlap:  " << mean_dice << endl;
  } else {
    cout << "Zero overlap between datasets" << endl;
  }

  return 0;
}
