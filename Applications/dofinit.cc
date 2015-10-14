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

#include <irtkTransformation.h>
#include <irtkPointSamples.h>

#ifdef HAS_VTK
#include <irtkPolyDataUtils.h>
using namespace irtk::polydata;
#endif


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <dofout> [output option] [options]" << endl;
  cout << endl;
  cout << "This tool either creates a new affine transformation with the given parameters" << endl;
  cout << "or approximates an affine transformation or non-rigid deformation given" << endl;
  cout << "a set of corresponding landmarks or point displacements (cf. -disp option)." << endl;
  cout << endl;
  cout << "The output transformation is by default an affine transformation." << endl;
  cout << "One of the output options below can be used to request a different" << endl;
  cout << "type of output transformation. This option should be given right" << endl;
  cout << "after the output transformation name as it may be altered using" << endl;
  cout << "one of the other output related options (e.g. -noshearing)." << endl;
  cout << endl;
  cout << "An affine transformation is the composition of a shearing followed by" << endl;
  cout << "a scaling, a rotation, and a translation. The homogeneous transformation" << endl;
  cout << "matrix, A, is given by the matrix product:" << endl;
  cout << endl;
  cout << "  A = Translate * Rotate * Scale * Shear, where Rotate = (Rz Ry Rx)^T" << endl;
  cout << endl;
  cout << "Output options:" << endl;
  cout << "  -rigid         Output rigid transformation." << endl;
  cout << "  -affine        Output affine transformation." << endl;
  cout << "  -affine-mffd   Output affine transformation as global component of MFFD." << endl;
  cout << "  -mffd          Output multi-level free-form deformation (MFFD)." << endl;
  cout << "  -ffd           Output free-form deformation (FFD)." << endl;
  cout << "  -svffd         Output diffeomorphic free-form deformation (SV FFD)." << endl;
  cout << endl;
  cout << "  -[no]translations   Allow/disallow translations." << endl;
  cout << "  -[no]rotations      Allow/disallow rotations." << endl;
  cout << "  -[no]scaling        Allow/disallow scaling." << endl;
  cout << "  -[no]shearing       Allow/disallow shearing." << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -dofin <dof>     Read affine parameters from input transformation." << endl;
  cout << "  -orientation <x1> <x2> <x3>  <y1> <y2> <y3>  <z1> <z2> <z3>" << endl;
  cout << "                   Set upper 3x3 (rotation) matrix entries" << endl;
  cout << "                   (e.g., image orientation from image info output)." << endl;
  cout << endl;
  cout << "  -tx    <tx>      Translation along x axis in mm.             (default: 0)" << endl;
  cout << "  -ty    <ty>      Translation along y axis in mm.             (default: 0)" << endl;
  cout << "  -tz    <tz>      Translation along z axis in mm.             (default: 0)" << endl;
  cout << "  -rx    <rx>      Rotation around x axis in degrees.          (default: 0)" << endl;
  cout << "  -ry    <ry>      Rotation around x axis in degrees.          (default: 0)" << endl;
  cout << "  -rz    <rz>      Rotation around x axis in degrees.          (default: 0)" << endl;
  cout << "  -sx    <sx>      Scaling of x axis in percentage.            (default: 100)" << endl;
  cout << "  -sy    <sy>      Scaling of y axis in percentage.            (default: 100)" << endl;
  cout << "  -sz    <sz>      Scaling of z axis in percentage.            (default: 100)" << endl;
  cout << "  -s     <s>       Isotropic scaling in percentage.            (default: 100)" << endl;
  cout << "  -sxy   <sxy>     Skew angle between x and y axis in degrees. (default: 0)" << endl;
  cout << "  -sxz   <sxz>     Skew angle between x and z axis in degrees. (default: 0)" << endl;
  cout << "  -syz   <syz>     Skew angle between y and z axis in degrees. (default: 0)" << endl;
  cout << endl;
  cout << "  -target <image>  Target image. Required when output is a FFD. (default: none)" << endl;
  cout << "  -source <image>  Source image. Used to initialize translation. (default: none)" << endl;
  cout << endl;
  cout << "  -disp <pset1> <pset2> [<cor12>]   Create transformation which minimizes the mean" << endl;
  cout << "                                    squared distance of (fiducial) points." << endl;
  cout << "                                    Option can be used multiple times." << endl;
  cout << "  -disp <pset>                      Create transformation which approximates the" << endl;
  cout << "                                    displacement VECTORS of the given VTK point set." << endl;
  cout << "                                    Option can be used multiple times." << endl;
  cout << "  -approximate <dofin>              Create transformation which approximates" << endl;
  cout << "                                    another given transformation." << endl;
  cout << "                                    Option can be used multiple times." << endl;
  cout << endl;
  cout << "Arguments for creating a free-form deformation (FFD):" << endl;
  cout << "  -dx <dx>         Control point spacing of FFD in x. (default: target voxel size)" << endl;
  cout << "  -dy <dy>         Control point spacing of FFD in y. (default: target voxel size)" << endl;
  cout << "  -dz <dz>         Control point spacing of FFD in z. (default: target voxel size)" << endl;
  cout << "  -ds <ds>         Control point spacing of FFD.      (default: target voxel size)" << endl;
  cout << "  -subdiv <n>      Number of subdivisions to use for approximating a non-rigid deformation" << endl;
  cout << "                   from input -displacements. (default: 4)" << endl;
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliary functions
// =============================================================================

// -----------------------------------------------------------------------------
/// Read correspondences from text file
vector<int> ReadCorrespondences(const char *fname, int n1, int n2)
{
  ifstream ifs(fname);
  if (!ifs.is_open()) {
    cerr << "Error: Failed to open correspondence input file: " << fname << endl;
    exit(1);
  }
  string line;
  int    t, s;
  vector<int> corr(n1, -1);
  while (getline(ifs, line)) {
    istringstream is(line);
    if (!(is >> t >> s)) {
      cerr << "Error: Failed to read correspondence map from file: " << fname << endl;
      exit(1);
    }
    if (t < 0 || t >= n1) {
      cerr << "Error: Invalid target index (" << t << ") in correspondence map: " << fname << endl;
      exit(1);
    }
    if (s < 0 || s >= n2) {
      cerr << "Error: Invalid source index (" << s << ") in correspondence map: " << fname << endl;
      exit(1);
    }
    corr[t] = s;
  }
  for (t = 0; t < n1; ++t) {
    if (corr[t] == -1) {
      cerr << "Error: Missing correspondence for target index " << t << endl;
      exit(1);
    }
  }
  ifs.close();
  return corr;
}

// -----------------------------------------------------------------------------
irtkFreeFormTransformation3D *
BSplineFFD(const irtkAffineTransformation &dof, const irtkImageAttributes &attr, double sx, double sy, double sz)
{
  irtkBSplineFreeFormTransformation3D *ffd = new irtkBSplineFreeFormTransformation3D(attr, sx, sy, sz);

  const int N = ffd->NumberOfCPs();
  double *dx = new double[N];
  double *dy = new double[N];
  double *dz = new double[N];

  int n = 0;
  for (int k = 0; k < ffd->GetZ(); ++k) {
    for (int j = 0; j < ffd->GetY(); ++j) {
      for (int i = 0; i < ffd->GetX(); ++i) {
        dx[n] = i, dy[n] = j, dz[n] = k;
        ffd->LatticeToWorld(dx[n], dy[n], dz[n]);
        dof .Displacement  (dx[n], dy[n], dz[n]);
        ++n;
      }
    }
  }

  ffd->Interpolate(dx, dy, dz);

  delete[] dx;
  delete[] dy;
  delete[] dz;

  return ffd;
}

// -----------------------------------------------------------------------------
irtkFreeFormTransformation3D *
BSplineSVFFD(const irtkAffineTransformation &dof, irtkImageAttributes attr, double sx, double sy, double sz)
{
  irtkBSplineFreeFormTransformationSV *ffd = new irtkBSplineFreeFormTransformationSV();

  ffd->NumberOfSteps(64);
  ffd->MaxScaledVelocity(.0);
  ffd->Initialize(attr, sx, sy, sz, &dof);

  if (verbose) {
    double error = .0, x, y, z;
    irtkGenericImage<double> disp(attr, 3);
    ffd->Displacement(disp);
    for (int k = 0; k < disp.GetZ(); ++k) {
      for (int j = 0; j < disp.GetY(); ++j) {
        for (int i = 0; i < disp.GetX(); ++i) {
          x = i, y = j, z = k;
          disp.ImageToWorld(x, y, z);
          dof .Displacement(x, y, z);
          x -= disp(i, j, k, 0);
          y -= disp(i, j, k, 1);
          z -= disp(i, j, k, 2);
          error += x*x + y*y + z*z;
        }
      }
    }
    error = sqrt(error / disp.GetX() * disp.GetY() * disp.GetZ());
    cout << "RMS error of B-spline SV FFD approximation = " << error << "\n" << endl;
  }

  return ffd;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  EXPECTS_POSARGS(1);

  // Number of point samples used for approximation of rigid/affine input transformation
  const int number_of_point_samples = 20;

  // Instantiate affine transformation
  irtkAffineTransformation dof;

  struct PointSetPair
  {
    const char *target_name; // Target fiducial points
    const char *source_name; // Source fiducial points
    const char *corrin_name; // Correspondences between target and source points
    PointSetPair() : target_name(NULL), source_name(NULL), corrin_name(NULL) {}
  };

  // Parse command-line arguments
  const char *dofout_name   = POSARG(1);
  const char *target_name   = NULL;  // Target image/FFD lattice read from NIfTI header
  const char *source_name   = NULL;  // Source image
  double      sx            = .0;    // Control point spacing in x = voxel size by default
  double      sy            = .0;    // Control point spacing in y = voxel size by default
  double      sz            = .0;    // Control point spacing in z = voxel size by default
  bool        multilevel    = false; // Create multi-level FFD
  bool        nonrigid      = false; // Create non-rigid transformation
  bool        diffeomorphic = false; // B-spline FFD or affine transformation by default
  int         nsubdiv       = 4;     // No. of subdivisions to approximate non-linear deformation
  vector<PointSetPair> psets;

  for (ALL_OPTIONS) {
    // Output options
    if (OPTION("-rigid")) {
      multilevel = false, nonrigid = false, diffeomorphic = true;
      dof.AllowTranslations(true);
      dof.AllowRotations(true);
      dof.AllowScaling(false);
      dof.AllowShearing(false);
    }
    else if (OPTION("-affine") || OPTION("-affine-mffd")) {
      multilevel    = (strstr(OPTNAME, "mffd") != NULL);
      nonrigid      = false;
      diffeomorphic = true;
      dof.AllowTranslations(true);
      dof.AllowRotations(true);
      dof.AllowScaling(true);
      dof.AllowShearing(true);
    }
    else if (OPTION("-mffd") || OPTION("-ffd") || OPTION("-svffd")) {
      multilevel    = (strstr(OPTNAME, "mffd") != NULL);
      nonrigid      = true;
      diffeomorphic = (strstr(OPTNAME, "sv") != NULL);
    }
    // Enable/disable affine transformation parameters
    // Note: Considered only for approximation of input displacements!
    else if (OPTION("-translations"  )) dof.AllowTranslations(true);
    else if (OPTION("-notranslations")) dof.AllowTranslations(false);
    else if (OPTION("-rotations"     )) dof.AllowRotations(true);
    else if (OPTION("-norotations"   )) dof.AllowRotations(false);
    else if (OPTION("-scaling"       )) dof.AllowScaling(true);
    else if (OPTION("-noscaling"     )) dof.AllowScaling(false);
    else if (OPTION("-shearing"      )) dof.AllowShearing(true);
    else if (OPTION("-noshearing"    )) dof.AllowShearing(false);
    // Initialize upper 3x3 (rotation) matrix
    else if (OPTION("-orientation")) {
      irtkMatrix m(4, 4);
      m(3, 3) = 1.0;
      const char *arg;
      for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
          arg = ARGUMENT;
          if (arg[0] == ',' && arg[1] == '\0') arg = ARGUMENT;
          m(r, c) = atof(arg);
        }
      }
      dof.PutMatrix(m);
    }
    // Read affine parameters from input transformation
    else if (OPTION("-dofin")) {
      irtkTransformation *dofin = irtkTransformation::New(ARGUMENT);
      irtkRigidTransformation      *dof6  = NULL;
      irtkSimilarityTransformation *dof7  = NULL;
      irtkAffineTransformation     *dof12 = NULL;
      irtkMultiLevelTransformation *mffd  = NULL;
      // Attention: Order of casts or if clauses below matters!
      (dof12 = dynamic_cast<irtkAffineTransformation     *>(dofin)) ||
      (dof7  = dynamic_cast<irtkSimilarityTransformation *>(dofin)) ||
      (dof6  = dynamic_cast<irtkRigidTransformation      *>(dofin)) ||
      (mffd  = dynamic_cast<irtkMultiLevelTransformation *>(dofin));
      if (mffd) dof12 = mffd->GetGlobalTransformation();
      if (dof12) {
        dof.PutTranslationX(dof12->GetTranslationX());
        dof.PutTranslationY(dof12->GetTranslationY());
        dof.PutTranslationZ(dof12->GetTranslationZ());
        dof.PutRotationX   (dof12->GetRotationX());
        dof.PutRotationY   (dof12->GetRotationY());
        dof.PutRotationZ   (dof12->GetRotationZ());
        dof.PutScaleX      (dof12->GetScaleX());
        dof.PutScaleY      (dof12->GetScaleY());
        dof.PutScaleZ      (dof12->GetScaleZ());
        dof.PutShearXY     (dof12->GetShearXY());
        dof.PutShearXZ     (dof12->GetShearXZ());
        dof.PutShearYZ     (dof12->GetShearYZ());
      } else if (dof7) {
        dof.PutTranslationX(dof7->GetTranslationX());
        dof.PutTranslationY(dof7->GetTranslationY());
        dof.PutTranslationZ(dof7->GetTranslationZ());
        dof.PutRotationX   (dof7->GetRotationX());
        dof.PutRotationY   (dof7->GetRotationY());
        dof.PutRotationZ   (dof7->GetRotationZ());
        dof.PutScale       (dof7->GetScale());
      } else if (dof6) {
        dof.PutTranslationX(dof6->GetTranslationX());
        dof.PutTranslationY(dof6->GetTranslationY());
        dof.PutTranslationZ(dof6->GetTranslationZ());
        dof.PutRotationX   (dof6->GetRotationX());
        dof.PutRotationY   (dof6->GetRotationY());
        dof.PutRotationZ   (dof6->GetRotationZ());
      } else {
        delete dofin;
        cerr << "Only a rigid, similarity, affine, or MFFD transformation can be used as -dofin" << endl;
        exit(1);
      }
      delete dofin;
    }
    // Affine transformation parameters
    else if (OPTION("-tx"))  dof.PutTranslationX(atof(ARGUMENT));
    else if (OPTION("-ty"))  dof.PutTranslationY(atof(ARGUMENT));
    else if (OPTION("-tz"))  dof.PutTranslationZ(atof(ARGUMENT));
    else if (OPTION("-rx"))  dof.PutRotationX   (atof(ARGUMENT));
    else if (OPTION("-ry"))  dof.PutRotationY   (atof(ARGUMENT));
    else if (OPTION("-rz"))  dof.PutRotationZ   (atof(ARGUMENT));
    else if (OPTION("-sx"))  dof.PutScaleX      (atof(ARGUMENT));
    else if (OPTION("-sy"))  dof.PutScaleY      (atof(ARGUMENT));
    else if (OPTION("-sz"))  dof.PutScaleZ      (atof(ARGUMENT));
    else if (OPTION("-s"))   dof.PutScale       (atof(ARGUMENT));
    else if (OPTION("-sxy")) dof.PutShearXY     (atof(ARGUMENT));
    else if (OPTION("-sxz")) dof.PutShearXZ     (atof(ARGUMENT));
    else if (OPTION("-syz")) dof.PutShearYZ     (atof(ARGUMENT));
    // Input displacements (unstructured)
    else if (OPTION("-approximate") || OPTION("-approx")) {
      PointSetPair p;
      p.target_name = ARGUMENT;
      psets.push_back(p);
    }
    else if (OPTION("-displacements") || OPTION("-disp")) {
      PointSetPair p;
      p.target_name = ARGUMENT;
      if (HAS_ARGUMENT) {
        p.source_name = ARGUMENT;
        if (HAS_ARGUMENT) p.corrin_name = ARGUMENT;
      }
      psets.push_back(p);
    }
    // FFD lattice attributes
    else if (OPTION("-dx"))  sx            = atof(ARGUMENT);
    else if (OPTION("-dy"))  sy            = atof(ARGUMENT);
    else if (OPTION("-dz"))  sz            = atof(ARGUMENT);
    else if (OPTION("-ds"))  sx = sy = sz  = atof(ARGUMENT);
    else if (OPTION("-subdiv")) nsubdiv = atoi(ARGUMENT);
    // Target/source image
    else if (OPTION("-target")) target_name = ARGUMENT;
    else if (OPTION("-source")) source_name = ARGUMENT;
    // Common or unknown option
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (nsubdiv <  0) nsubdiv = 4;
  if (nsubdiv == 0) nsubdiv = 1;

  if (nonrigid && !target_name) {
    cerr << "Target image required when creating a non-rigid FFD!" << endl;
    exit(1);
  }

  // Read attributes of target image
  irtkImageAttributes attr;
  if (target_name) {
    irtkGreyImage target(target_name);
    attr = target.Attributes();
  }

  // ---------------------------------------------------------------------------
  // Usage 1) Create affine (FFD) transformation from parameter values read
  //          from command-line or input transformation file
  if (psets.empty()) {

    // Initialize affine transformation from target and source image attributes
    if (source_name && dof.IsIdentity()) {
      irtkGreyImage source(source_name);
      const irtkImageAttributes &sattr = source.Attributes();
      dof.PutTranslationX(sattr._xorigin - attr._xorigin);
      dof.PutTranslationY(sattr._yorigin - attr._yorigin);
      dof.PutTranslationZ(sattr._zorigin - attr._zorigin);
    }

    // Reset passive parameters
    irtkAffineTransformation copy(dof);
    dof.CopyFrom(&copy);

    // Write affine transformation as global transformation of MFFD
    if (multilevel) {

      irtkMultiLevelFreeFormTransformation mffd(dof);
      if (verbose) mffd.Print();
      mffd.Write(dofout_name);

    // Write affine transformation as free-form deformation
    } else if (nonrigid) {

      auto_ptr<irtkFreeFormTransformation3D> ffd;
      if (diffeomorphic) ffd.reset(BSplineSVFFD(dof, attr, sx, sy, sz));
      else               ffd.reset(BSplineFFD  (dof, attr, sx, sy, sz));
      if (verbose) {
        dof.Print();
        ffd->Print();
      }
      ffd->Write(dofout_name);

    // Write affine transformation
    } else {

      if (verbose) dof.Print();
      dof.Write(dofout_name);

    }

  // ---------------------------------------------------------------------------
  // Usage 2) Create transformation which approximates given point displacements
  } else {

    // -------------------------------------------------------------------------
    // Combine input point sets into single unstructured displacement map

    irtkPointSet target_points, source_points;
    irtkHomogeneousTransformation *lin;
    irtkFreeFormTransformation    *ffd;
    irtkMultiLevelTransformation  *mffd;

    int n = 0;
    for (size_t i = 0; i < psets.size(); ++i) {
      if (psets[i].source_name == NULL && IsTransformation(psets[i].target_name)) {
        irtkTransformation *dof = irtkTransformation::New(psets[i].target_name);
        lin  = dynamic_cast<irtkHomogeneousTransformation *>(dof);
        ffd  = dynamic_cast<irtkFreeFormTransformation    *>(dof);
        mffd = dynamic_cast<irtkMultiLevelTransformation  *>(dof);
        if (mffd && mffd->NumberOfLevels() > 0) ffd = mffd->GetLocalTransformation(-1);
        if (lin) {
          n += number_of_point_samples;
        } else if (ffd) {
          n += ffd->NumberOfActiveCPs();
        } else {
          cerr << "Unsupported transformation file: " << psets[i].target_name << endl;
          delete dof;
          exit(1);
        }
        delete dof;
      } else {
        target_points.Read(psets[i].target_name);
        n += target_points.Size();
      }
    }

    // Allocate memory for unstructured displacement map
    double *x  = Allocate<double>(n);
    double *y  = Allocate<double>(n);
    double *z  = Allocate<double>(n);
    double *dx = Allocate<double>(n);
    double *dy = Allocate<double>(n);
    double *dz = Allocate<double>(n);

    int j = 0;
    for (size_t i = 0; i < psets.size(); ++i) {

      // Add displacements of corresponding landmark points
      if (psets[i].source_name) {

        target_points.Read(psets[i].target_name);
        source_points.Read(psets[i].source_name);

        // Calculate unstructured displacement map
        if (psets[i].corrin_name) {
          vector<int> corr = ReadCorrespondences(psets[i].corrin_name,
                                                 target_points.Size(),
                                                 source_points.Size());
          for (int t = 0; t < target_points.Size(); ++t, ++j) {
            const int &s = corr[t];
            x [j] = target_points(t)._x;
            y [j] = target_points(t)._y;
            z [j] = target_points(t)._z;
            dx[j] = source_points(s)._x - x[j];
            dy[j] = source_points(s)._y - y[j];
            dz[j] = source_points(s)._z - z[j];
          }
        } else {
          if (source_points.Size() < target_points.Size()) {
            cerr << "Second point set has different number of points! An explicit correspondence map is required." << endl;
            Deallocate(x);
            Deallocate(y);
            Deallocate(z);
            Deallocate(dx);
            Deallocate(dy);
            Deallocate(dz);
            exit(1);
          }
          for (int t = 0; t < target_points.Size(); ++t, ++j) {
            x [j] = target_points(t)._x;
            y [j] = target_points(t)._y;
            z [j] = target_points(t)._z;
            dx[j] = source_points(t)._x - x[j];
            dy[j] = source_points(t)._y - y[j];
            dz[j] = source_points(t)._z - z[j];
          }
        }

      // Add point displacements from input transformation file
      } else if (IsTransformation(psets[i].target_name)) {

        irtkTransformation *dof = irtkTransformation::New(psets[i].target_name);
        lin  = dynamic_cast<irtkHomogeneousTransformation *>(dof);
        ffd  = dynamic_cast<irtkFreeFormTransformation    *>(dof);
        mffd = dynamic_cast<irtkMultiLevelTransformation  *>(dof);
        if (mffd && mffd->NumberOfLevels() > 0) ffd = mffd->GetLocalTransformation(-1);
        if (lin) {

          irtkPointSamples samples(number_of_point_samples);
          samples.SampleGaussian(.0, 100.0);

          int j1 = j;
          for (int i = 0; i < samples.Size(); ++i, ++j) {
            x[j] = samples(i)._x;
            y[j] = samples(i)._y;
            z[j] = samples(i)._z;
          }

          lin->Transform(samples);

          j = j1;
          for (int i = 0; i < samples.Size(); ++i, ++j) {
            dx[j] = samples(i)._x - x[j];
            dy[j] = samples(i)._y - y[j];
            dz[j] = samples(i)._z - z[j];
          }

        } else if (ffd) {

          irtkGenericImage<double> disp(ffd->Attributes(), 3);
          ffd->Displacement(disp);

          for (int c = 0; c < disp.Z(); ++c)
          for (int b = 0; b < disp.Y(); ++b)
          for (int a = 0; a < disp.X(); ++a, ++j) {
            x[j] = a, y[j] = b, z[j] = c;
            disp.ImageToWorld(x[j], y[j], z[j]);
            dx[j] = disp(a, b, c, 0);
            dy[j] = disp(a, b, c, 1);
            dz[j] = disp(a, b, c, 2);
          }

        }
        delete dof;

      // Add point displacements read from VTK point set file
      } else {

#ifdef HAS_VTK
        vtkSmartPointer<vtkPointSet> pset = ReadPointSet(psets[i].target_name);
        bool pset_ok = true;
        if (pset->GetNumberOfPoints() == 0) {
          cerr << "Failed to read point set " << psets[i].target_name << " or it contains no points" << endl;
          pset_ok = false;
        }
        vtkDataArray *disp = pset->GetPointData()->GetVectors();
        if (pset_ok && disp == NULL) {
          cerr << "Point set " << psets[i].target_name << " has no displacement vectors!" << endl;
          pset_ok = false;
        }
        if (pset_ok && disp->GetNumberOfComponents() != 3) {
          cerr << "Point set " << psets[i].target_name << " displacement vectors must have dimension 3!" << endl;
          pset_ok = false;
        }
        if (!pset_ok) {
          Deallocate(x);
          Deallocate(y);
          Deallocate(z);
          Deallocate(dx);
          Deallocate(dy);
          Deallocate(dz);
          exit(1);
        }
        double p[3], v[3];
        for (vtkIdType ptId = 0; ptId < pset->GetNumberOfPoints(); ++ptId, ++j) {
          pset->GetPoint(ptId, p);
          disp->GetTuple(ptId, v);
          x [j] = p[0];
          y [j] = p[1];
          z [j] = p[2];
          dx[j] = v[0];
          dy[j] = v[1];
          dz[j] = v[2];
        }
#else
        cerr << EXECNAME << ": Must be compiled with VTK when displacements should be read from VTK point set file!" << endl;
        exit(1);
#endif
      }
    }

    // -------------------------------------------------------------------------
    // 2a) Find non-rigid (M)FFD which approximates the input displacements
    if (nonrigid) {

      sx *= pow(2.0, nsubdiv-1);
      sy *= pow(2.0, nsubdiv-1);
      sz *= pow(2.0, nsubdiv-1);

      double *rx = Allocate<double>(n);
      double *ry = Allocate<double>(n);
      double *rz = Allocate<double>(n);

      if (multilevel) {

        memcpy(rx, dx, n * sizeof(double));
        memcpy(ry, dy, n * sizeof(double));
        memcpy(rz, dz, n * sizeof(double));
        dof.ApproximateAsNew(x, y, z, rx, ry, rz, n);

        irtkMultiLevelFreeFormTransformation mffd(dof);
        irtkFreeFormTransformation *ffd;
        if (diffeomorphic) ffd = new irtkBSplineFreeFormTransformationSV(attr, sx, sy, sz);
        else               ffd = new irtkBSplineFreeFormTransformation3D(attr, sx, sy, sz);
        mffd.PushLocalTransformation(ffd);

        double rms;
        for (int l = 0; l < nsubdiv; ++l) {
          if (l > 0) ffd->Subdivide();
          memcpy(rx, dx, n * sizeof(double));
          memcpy(ry, dy, n * sizeof(double));
          memcpy(rz, dz, n * sizeof(double));
          rms = ffd->Approximate(x, y, z, rx, ry, rz, n);
        }

        if (verbose) {
          cout << "RMS error of approximation is " << rms << "\n" << endl;
          mffd.Print();
        }
        mffd.Write(dofout_name);

      } else {

        auto_ptr<irtkFreeFormTransformation3D> ffd;
        if (diffeomorphic) ffd.reset(new irtkBSplineFreeFormTransformationSV(attr, sx, sy, sz));
        else               ffd.reset(new irtkBSplineFreeFormTransformation3D(attr, sx, sy, sz));

        double rms;
        for (int l = 0; l < nsubdiv; ++l) {
          if (l > 0) ffd->Subdivide();
          memcpy(rx, dx, n * sizeof(double));
          memcpy(ry, dy, n * sizeof(double));
          memcpy(rz, dz, n * sizeof(double));
          rms = ffd->Approximate(x, y, z, rx, ry, rz, n);
        }

        if (verbose) {
          cout << "RMS error of approximation is " << rms << "\n" << endl;
          ffd->Print();
        }
        ffd->Write(dofout_name);

      }

      Deallocate(rx);
      Deallocate(ry);
      Deallocate(rz);

    // -------------------------------------------------------------------------
    // 2b) Find affine transformation which approximates the input displacements
    } else {

      // Minimize mean squared error of approximation
      double rms = dof.ApproximateAsNew(x, y, z, dx, dy, dz, n);
      if (verbose) cout << "RMS error of approximation is " << rms << "\n" << endl;

      // Write transformation
      if (multilevel) {
        irtkMultiLevelFreeFormTransformation mffd(dof);
        mffd.Print();
        mffd.Write(dofout_name);
      } else {
        if (verbose) dof.Print();
        dof.Write(dofout_name);
      }

    }

    // -------------------------------------------------------------------------
    // Free displacement map
    Deallocate(x);
    Deallocate(y);
    Deallocate(z);
    Deallocate(dx);
    Deallocate(dy);
    Deallocate(dz);
  }

  return 0;
}
