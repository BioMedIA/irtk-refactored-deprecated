// =============================================================================
// Project: Image Registration Toolkit (IRTK)
// Package: PolyData
//
// Copyright (c) 2015 Imperial College London
// Copyright (c) 2015 Andreas Schuh
// =============================================================================

#include <irtkCommon.h>
#include <memory>

// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
/// Print help screen
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "usage: " << name << " <input> <output> [options]" << endl;
  cout << "       " << name << " <input1> <input2> <output1> <output2> [options]" << endl;
  cout << endl;
  cout << "Performs a spectral analysis of the general graph Laplacian matrix" << endl;
  cout << "computed from the given input surface mesh(es). If a -target surface" << endl;
  cout << "is specified, the eigenvectors are reordered and their sign flipped" << endl;
  cout << "to match the eigenvectors of the target surface mesh." << endl;
  cout << endl;
  cout << "If two input and output surface meshes are specified, a spectral analysis" << endl;
  cout << "of the joint graph Laplacian is performed after an initial surface match." << endl;
  cout << "The type of initial point correspondences is defined by the -corr option." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  input   File name of input  surface (vtkPolyData)." << endl;
  cout << "  output  File name of output surface (vtkPolyData)." << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << "  -target <file>            File name of dataset (vtkPolyData) whose spectral components" << endl;
  cout << "                            should be used to adjust the sign and order of the eigenmodes." << endl;
  cout << "  -corr <string>            Type of point correspondences to use for initial correspondences:" << endl;
  cout << "                            \"closest point\" (\"CP\"), \"spectral match\" (\"SM\"). (default: SM)" << endl;
  cout << "  -corrpar <name>=<value>   Name/value pair of correspondence parameter." << endl;
  cout << "  -corrin <file>            Text file naming for each point in input1 the index of the corresponding" << endl;
  cout << "                            point in the second dataset input2, one index per line. (default: -corr)" << endl;
  cout << "  -points                   Replace output vertex coordinates by first three spectral coordinates." << endl;
  PrintCommonOptions(cout);
  cout << endl;
}

// =============================================================================
// Includes
// =============================================================================

#include <vtkMath.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>

#include <irtkRegisteredSurface.h>
#include <irtkPointCorrespondence.h>
#include <irtkFuzzyCorrespondence.h>
#include <irtkSpectralDecomposition.h>

#ifdef HAVE_MATLAB
#include <mclmcrrt.h>
#endif

using namespace irtkSpectralDecomposition;

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
/// Write point set with surface eigenmodes
void Write(const char *fname, vtkPointSet *pointset, vtkPolyData *surface, bool ascii = false, bool as_points = false)
{
  vtkSmartPointer<vtkPointSet> output = pointset;

  vtkDataArray *modes = surface->GetPointData()->GetArray("joint_eigenmodes");
  if (!modes)   modes = surface->GetPointData()->GetArray("eigenmodes");
  if (!modes) {
    cerr << "Output surface mesh has no eigenmodes point data!" << endl;
    exit(1);
  }

  if (as_points) {

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(surface->GetNumberOfPoints());

    double *p = new double[max(3, int(modes->GetNumberOfComponents()))];
    for (int c = 0; c < 3; ++c) p[c] = .0;
    for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i) {
      modes ->GetTuple(i, p);
      points->SetPoint(i, p);
    }
    delete[] p;

    output = vtkSmartPointer<vtkPolyData>::New();
    output->ShallowCopy(surface);
    output->SetPoints(points);

  } else if (pointset != surface) {

    vtkSmartPointer<vtkDataArray> output_modes = modes->NewInstance();
    output_modes->SetName(modes->GetName());
    output_modes->SetNumberOfComponents(modes->GetNumberOfComponents());
    output_modes->SetNumberOfTuples(pointset->GetNumberOfPoints());

    for (vtkIdType ptId = 0; ptId < pointset->GetNumberOfPoints(); ++ptId) {
      for (int j = 0; j < modes->GetNumberOfComponents(); ++j) {
        output_modes->SetComponent(ptId, j, .0);
      }
    }

    vtkIdType origPtId;
    vtkDataArray *origPtIds = surface->GetPointData()->GetArray("vtkOriginalPointIds");
    for (vtkIdType ptId = 0; ptId < surface->GetNumberOfPoints(); ++ptId) {
      origPtId = static_cast<vtkIdType>(origPtIds->GetComponent(ptId, 0));
      for (int j = 0; j < modes->GetNumberOfComponents(); ++j) {
        output_modes->SetComponent(origPtId, j, modes->GetComponent(ptId, j));
      }
    }

    output->GetPointData()->AddArray(output_modes);
  }

  WritePointSet(fname, output, true, ascii);
}

// -----------------------------------------------------------------------------
/// Compute eigenmodes of joint intra- and inter-mesh graph Laplacian matrix
int ComputeJointEigenmodes(vtkPolyData *target, const vector<int> &target_sample,
                           vtkPolyData *source, const vector<int> &source_sample,
                           const irtkPointCorrespondence *corr, int k,
                           irtkMatrix &modes, irtkVector &freq)
{
  // Set intra-mesh affinity weights
  const int m = target->GetNumberOfPoints();
  const int n = source->GetNumberOfPoints();
  SparseMatrix::Entries *cols = new SparseMatrix::Entries[m + n];
  AdjacencyMatrix(cols, SparseMatrix::CCS, 0, 0, target);
  AdjacencyMatrix(cols, SparseMatrix::CCS, m, m, source);
  // Set inter-mesh affinity weights
  const irtkFuzzyCorrespondence *fuzzy;
  if ((fuzzy = dynamic_cast<const irtkFuzzyCorrespondence *>(corr))) {
    ConnectivityMatrix(cols, SparseMatrix::CCS, 0, m, target, &target_sample,
                       source, &source_sample, &fuzzy->Weight(), false);
    ConnectivityMatrix(cols, SparseMatrix::CCS, m, 0, source, &source_sample,
                       target, &target_sample, &fuzzy->Weight(), true );
  } else {
    int i, j;
    double w, p1[3], p2[3];
    for (size_t s = 0; s < target_sample.size(); ++s) {
      i = target_sample[s];
      j = corr->GetSourceIndex(i);
      if (j < 0) {
        cerr << "Error: Point correspondence must be fuzzy or map point indices" << endl;
        exit(1);
      }
      target->GetPoint(i, p1);
      source->GetPoint(j, p2);
      w = 1.0 / (sqrt(vtkMath::Distance2BetweenPoints(p1, p2)) + EPSILON);
      cols[i  ].push_back(make_pair(j+m, w));
      cols[j+m].push_back(make_pair(i,   w));
    }
    for (size_t s = 0; s < source_sample.size(); ++s) {
      j = source_sample[s];
      i = corr->GetTargetIndex(j);
      if (i < 0) {
        cerr << "Error: Point correspondence must be fuzzy or map point indices" << endl;
        exit(1);
      }
      target->GetPoint(i, p1);
      source->GetPoint(j, p2);
      w = 1.0 / (sqrt(vtkMath::Distance2BetweenPoints(p1, p2)) + EPSILON);
      cols[j+m].push_back(make_pair(i,   w));
      cols[i  ].push_back(make_pair(j+m, w));
    }
  }
  // Compute graph Laplacian of joint connectivity graph
  SparseMatrix L(SparseMatrix::CCS);
  L.Initialize(m + n, m + n, cols);
  delete[] cols;
  NormalizedLaplacian(L, L);
  // Compute eigenmodes of joint graph Laplacian
  return ComputeEigenmodes(L, k+1, modes, freq);
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  REQUIRES_POSARGS(2);

  const char *input_name [2] = {NULL};
  const char *output_name[2] = {NULL};
  int ndatasets = 0;

  switch (NUM_POSARGS) {
    case 2:
      input_name [0] = POSARG(1);
      output_name[0] = POSARG(2);
      ndatasets = 1;
      break;
    case 4:
      input_name [0] = POSARG(1);
      input_name [1] = POSARG(2);
      output_name[0] = POSARG(3);
      output_name[1] = POSARG(4);
      ndatasets = 2;
      break;
    default:
      PrintHelp(EXECNAME);
      cerr << "Error: Invalid number of positional arguments (" << NUM_POSARGS << ")" << endl;
      exit(1);
  };

  // Optional arguments
  const char *target_name = NULL;
  const char *dofin_name  = NULL;
  int         k           = 5;
  bool        ascii       = false;
  bool        as_points   = false;

  irtkPointCorrespondence::TypeId ctype = irtkPointCorrespondence::ClosestPoint;
  irtkParameterList               cparam;
  vector<string>                  cfeature_name;
  vector<double>                  cfeature_weight;

  for (ALL_OPTIONS) {
    if      (OPTION("-target")) target_name  = ARGUMENT;
    else if (OPTION("-dofin"))  dofin_name   = ARGUMENT;
    else if (OPTION("-k") || OPTION("-dim")) k = atoi(ARGUMENT);
    else if (OPTION("-corr")) {
      const char *arg = ARGUMENT;
      if (!FromString(arg, ctype)) {
        cerr << "Invalid -corr argument: " << arg << endl;
        exit(1);
      }
    }
    else if (OPTION("-corrin")) {
      ctype = irtkPointCorrespondence::FiducialMatch;
      Insert(cparam, "Correspondence map", ARGUMENT);
    }
    else if (OPTION("-corrpar")) {
      Insert(cparam, ARGUMENT, ARGUMENT);
    }
    else if (OPTION("-feature")) {
      cfeature_name.push_back(ARGUMENT);
      if (HAS_ARGUMENT) cfeature_weight.push_back(atof(ARGUMENT));
      else              cfeature_weight.push_back(1.0);
    }
    else if (OPTION("-points")) as_points = true;
    else if (OPTION("-ascii")  || OPTION("-nobinary")) ascii = true;
    else if (OPTION("-binary") || OPTION("-noascii") ) ascii = false;
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // Default correspondence features and normalization of names to lowercase
  if (cfeature_name.empty()) {
    cfeature_name  .push_back("eigenmodes");
    cfeature_weight.push_back(1.0);
  } else {
    for (size_t i = 0; i < cfeature_name.size(); ++i) {
      transform(cfeature_name[i].begin(), cfeature_name[i].end(),
                cfeature_name[i].begin(), ::tolower);
    }
  }

  // Ensure number of spectral components is positive
  if (k <= 0) {
    cerr << "Error: Argument for option -k/-dim must be positive" << endl;
    exit(1);
  }
  if (as_points && k < 3) {
    cerr << "Error: Argument for option -k/-dim must be >= 3 when -points should be written" << endl;
    exit(1);
  }

  // Read input datasets
  vtkSmartPointer<vtkPointSet> dataset[2];
  vtkSmartPointer<vtkPolyData> surface[2];
  for (int i = 0; i < ndatasets; ++i) {
    if (verbose > 1) cout << "Reading point set " << (i+1) << " from " << input_name[i] << endl;
    dataset[i] = ReadPointSet(input_name[i]);
    if (dataset[i] == NULL || dataset[i]->GetNumberOfPoints() == 0) {
      cerr << "Error: Failed to read dataset " << (i+1) << " or dataset is empty!" << endl;
      exit(1);
    }
    surface[i] = DataSetSurface(dataset[i], true);
  }

  // Transform first surface by specified transformation
  if (dofin_name) {
    std::unique_ptr<irtkTransformation> dof(irtkTransformation::New(dofin_name));
    irtkRegisteredSurface transformed;
    transformed.InputSurface(surface[0]);
    transformed.Transformation(dof.get());
    transformed.Initialize();
    transformed.Update();
    surface[0] = transformed.Surface();
  }

  // Initialize MCR (for debug output of matrices as MAT files)
#ifdef HAVE_MATLAB
  if (debug) mclmcrInitialize();
#endif

  // Skip individual spectral analysis if joint analysis requested based on
  // pre-defined correspondences, non-spectral matches, or known eigenmodes
  bool individual_analysis[2] = {false};
  bool spectral_match         =  false;
  if (ndatasets == 1) {
    spectral_match = !as_points; // -points option can be applied to single input dataset
  } else if (ndatasets > 1) {
    switch (ctype) {
      case irtkPointCorrespondence::FiducialMatch:
        spectral_match = false;
        break;
      case irtkPointCorrespondence::SpectralMatch:
        spectral_match = true;
        break;
      default:
        spectral_match = false;
        for (size_t i = 0; i < cfeature_name.size(); ++i) {
          if (cfeature_name[i].find("spectral" ) != string::npos ||
              cfeature_name[i].find("eigenmode") != string::npos) {
            spectral_match = true;
            break;
          }
        }
    }
  }
  if (spectral_match) {
    for (int i = 0; i < ndatasets; ++i) {
      individual_analysis[i] = (!surface[i]->GetPointData()->HasArray("eigenmodes") &&
                                !surface[i]->GetPointData()->HasArray("joint_eigenmodes"));
    }
  }

  // Individual spectral analysis of input datasets
  irtkMatrix modes[2];
  irtkVector freq [2];
  for (int i = 0; i < ndatasets; ++i) {
    if (individual_analysis[i]) {
      if (ComputeEigenmodes(surface[i], k, modes[i], freq[i]) < k) {
        cerr << "Error: Failed to find " << k << " eigenmodes of dataset " << (i+1) << endl;
        exit(1);
      }
      if (verbose) {
        cout << "Frequencies";
        if (ndatasets > 1) cout << " of dataset " << (i+1);
        cout << ": ";
        for (int c = 0; c < k; ++c) {
          if (c > 0) cout << "  ";
          cout << scientific << freq[i](c);
        }
        cout << endl;
      }
    }
  }

  // Joint spectral analysis
  if (ndatasets == 2) {

    // Adjust sign and order of individually computed eigenmodes
    if (spectral_match) {
      irtkVector cost = MatchEigenmodes(surface[0]->GetPoints(), modes[0], freq[0],
                                        surface[1]->GetPoints(), modes[1], freq[1]);

      // TODO: Transform modes[1] using coherent point drift algorithm

      // Add individual eigenmodes to output point data
      for (int i = 0; i < ndatasets; ++i) {
        SetEigenmodes(surface[i], modes[i], "eigenmodes");
        modes[i].Clear();
        freq [i].Clear();
      }
    }

    // Obtain correspondence map for initial links
    irtkRegisteredSurface target, source;
    target.InputSurface(surface[0]);
    source.InputSurface(surface[1]);
    target.Initialize();
    source.Initialize();
    target.Update();
    source.Update();

    std::unique_ptr<irtkPointCorrespondence> cmap(irtkPointCorrespondence::New(ctype));
    cmap->FromTargetToSource(true);
    cmap->FromSourceToTarget(true);
    cmap->Parameter(cparam);
    cmap->Target(&target);
    cmap->Source(&source);
    for (size_t i = 0; i < cfeature_name.size(); ++i) {
      cmap->AddFeature(cfeature_name[i].c_str(), cfeature_weight[i]);
    }
    cmap->Initialize();
    cmap->Update();

    // Draw samples for which to add connectivity links
    vector<int> target_sample;
    vector<int> source_sample;
    using irtkPointCorrespondenceUtils::SamplePoints;
    SamplePoints(&target, target_sample, target.NumberOfPoints() / 10);
    SamplePoints(&source, source_sample, source.NumberOfPoints() / 10);

    // Perform spectral analysis of joint graph Laplacian
    irtkMatrix modes;
    irtkVector freq;
    if (ComputeJointEigenmodes(target.Surface(), target_sample,
                               source.Surface(), source_sample,
                               cmap.get(), k, modes, freq) < k) {
      cerr << "Error: Failed to find " << k << " eigenmodes of joint graph Laplacian" << endl;
      exit(1);
    }

    // Weight eigenmodes
    irtkVector w(k);
    double wsum = .0;
    for (int i = 0; i < k; ++i) {
      wsum += 1.0 / sqrt(freq(i));
    }
    for (int i = 0; i < k; ++i) {
      modes.ScaleCol(i, (1.0 / sqrt(freq(i))) / wsum);
    }

    // Add eigenmodes to output point data
    SetEigenmodes(surface[0], modes, 0,                       k, "joint_eigenmodes");
    SetEigenmodes(surface[1], modes, target.NumberOfPoints(), k, "joint_eigenmodes");

  } else {

    // Add individual eigenmodes to output point data
    for (int i = 0; i < ndatasets; ++i) {
      if (individual_analysis[i]) {
        SetEigenmodes(surface[i], modes[i], "eigenmodes");
        modes[i].Clear();
        freq [i].Clear();
      }
    }

  }

  // Match spectral components of input datasets with those of target dataset
  if (target_name) {
    if (verbose > 1) cout << "Reading target from " << target_name << endl;
    vtkSmartPointer<vtkPointSet> target = ReadPointSet(target_name);
    if (target == NULL || target->GetNumberOfPoints() == 0) {
      cerr << "Warning: Failed to read target dataset or dataset has no points!" << endl;
      cerr << "         Skipping matching of spectral components therefore." << endl;
    } else {
      for (int i = 0; i < ndatasets; ++i) {
        // TODO: Match eigenmodes
      }
    }
  }

  // Write datasets together with their spectral components
  for (int i = 0; i < ndatasets; ++i) {
    Write(output_name[i], dataset[i], surface[i], ascii, as_points);
    if (verbose > 1) cout << "Wrote result dataset to " << output_name[i] << endl;
  }

  return 0;
}
