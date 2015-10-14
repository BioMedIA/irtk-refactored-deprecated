// =============================================================================
// Project: Image Registration Toolkit (IRTK)
// Package: PolyData
//
// Copyright (c) 2015 Imperial College London
// Copyright (c) 2015 Andreas Schuh
// =============================================================================

#include <irtkCommon.h>

// -----------------------------------------------------------------------------
// Configuration

// Choose self intersection checks to perform
// Best results are obtained with the triangle/triangle intersection checks.
#define SELFINTERSECTION_CHECK_NORMAL    0 // line segment of old and new vertex position
#define SEFLINTERSECTION_CHECK_EDGES     0 // edges formed by modified triangles
#define SELFINTERSECTION_CHECK_TRIANGLES 1 // exhaustive triangle/triangle intersection checks
#define SELFINTERSECTION_CHECK_TOUCHING  1 // avoid triangles getting too close

// -----------------------------------------------------------------------------
/// Print help screen
void PrintHelp(const char *name)
{
  cout << "usage: " << name << " <seg|psi> [options]" << endl;
  cout << endl;
  cout << "Extracts the cortical surface from a given tissue segmentation or" << endl;
  cout << "precomputed psi image." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  seg   Brain tissue segmentation with 0: Outside, 1: CSF, 2: cGM, 3: WM, 4: BG, 5: Vent, 6: Cerebellum+Brainstem, 7: dGM" << endl;
  cout << "  psi   Precomputed psi image (see -psi option). Useful when extracting" << endl;
  cout << "        multiple isosurfaces from the same tissue segmentation." << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -isseg               Whether input image is a label image. (default: yes)" << endl;
  cout << "  -ispsi               Whether input image is a psi image. (default: no)" << endl;
  cout << "  -psi  <file>         Name of psi image output file. (default: none)" << endl;
  cout << "  -psi0 <float>        Value of target isosurface. (default: 9000)" << endl;
  cout << "  -psi1 <float>        Initial isosurface value. (default: 2000)" << endl;
  cout << "  -init <file|hull>    Name of initial surface mesh. (default: none)" << endl;
  cout << "  -labels <file>       Text file with input labels corresponding to each tissue class," << endl;
  cout << "                       where line number corresponds to tissue label minus 1, i.e." << endl;
  cout << "                       the comma-separated input labels on the third line are associted" << endl;
  cout << "                       with cortical grey matter (cGM = 2). (default: input labels)" << endl;
  cout << "  -mesh <file>         Name of surface mesh output file. (default: none)" << endl;
  cout << "  -mask <file>         Mask image defining region of brain segmentation." << endl;
  cout << "  -rh/-lh              Whether to extract cortex of right (mask=1) or left (mask=2) hemisphere." << endl;
  cout << "  -lambda <float>      Laplacian smoothing parameter applied after each iteration. (default: 0)" << endl;
  cout << "  -minactive <float>   Minimum number of active vertices in %. (default: 10%)" << endl;
  cout << "  -maxn <n>            Maximum number of iterations. (default: 500)" << endl;
  cout << "  -normal              Move vertices along normal direction. (default)" << endl;
  cout << "  -gradient            Move vertices along gradient direction." << endl;
}

// -----------------------------------------------------------------------------
// Includes
#include <irtkConvolution.h>
#include <irtkGaussianBlurring.h>
#include <irtkResampling.h>
#include <irtkPolyDataSmoothing.h>
#include <irtkNearestNeighborInterpolateImageFunction.hxx>
#include <irtkFastCubicBSplineInterpolateImageFunction.hxx>
#include <irtkLinearInterpolateImageFunction.hxx>
#include <irtkGradientImageFilter.h>
#include <irtkRemesher.h>

#include <vtkMath.h>
#include <vtkPlane.h>
#include <vtkTriangle.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkCellArray.h>
#include <vtkCellLocator.h>
#include <vtkKdTreePointLocator.h>
#include <vtkMarchingCubes.h>
#include <vtkTriangleFilter.h>
#include <vtkModifiedBSPTree.h>
#include <vtkIntersectionPolyDataFilter.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkDecimatePro.h>
#include <vtkPolyDataNormals.h>
#include <vtkHull.h>
#include <vtkMaskPoints.h>
#include <vtkDelaunay3D.h>
#include <vtkDataSetSurfaceFilter.h>

// -----------------------------------------------------------------------------
// Typedefs

/// Type of psi image interpolator
typedef irtkGenericFastCubicBSplineInterpolateImageFunction<irtkRealImage> ImageInterpolator;

/// Type of psi image gradient interpolator
typedef irtkGenericFastCubicBSplineInterpolateImageFunction<irtkRealImage> GradientInterpolator;

// -----------------------------------------------------------------------------
/// Brain tissue segmentation labels (cf. NeoSeg)
enum Label
{
  Outside = 0, ///< Outside brain extraction
  CSF     = 1, ///< Cerebrospinal fluid
  cGM     = 2, ///< Cortical grey matter
  WM      = 3, ///< White matter
  BG      = 4, ///< Background
  Vent    = 5, ///< Lateral ventricles
  CB      = 6, ///< Cerebellum + brainstem
  dGM     = 7, ///< Deep grey matter
  NUM_LABELS   ///< Max. number of segments
};

// -----------------------------------------------------------------------------
/// Labels of brain hemispheres in hemispheres mask image (cf. NeoSeg cortmrf)
enum Hemisphere
{
  UH, ///< Unspecified hemishpere
  RH, ///< Right brain hemisphere
  LH  ///< Left  brain hemisphere
};

// -----------------------------------------------------------------------------
/// Type of mesh collision
enum TypeOfCollision
{
  NO_COLLISION,
  SELF_INTERSECTION,
  TOUCHING
};

// -----------------------------------------------------------------------------
/// Write surface mesh to disk
void Write(vtkPolyData *mesh, const char *fname, bool ascii = false)
{
  vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
  SetVTKInput(writer, mesh);
  writer->SetFileName(fname);
  if (ascii) writer->SetFileTypeToASCII();
  else       writer->SetFileTypeToBinary();
  writer->Write();
}

// -----------------------------------------------------------------------------
/// Resample image to isotropic voxel size
void MakeIsotropic(irtkByteImage &seg)
{
  const double ds = min(min(seg.GetXSize(), seg.GetYSize()), seg.GetZSize());
  irtkGenericNearestNeighborInterpolateImageFunction<irtkByteImage> interp;
  irtkResampling<irtkBytePixel> resampling(ds, ds, ds);
  resampling.SetInput(&seg);
  resampling.SetOutput(&seg);
  resampling.SetInterpolator(&interp);
  resampling.Run();
}

// -----------------------------------------------------------------------------
/// Resample mask or label image on specified image lattice
void Resample(irtkByteImage &image, const irtkImageAttributes &attr)
{
  irtkByteImage output(attr);
  irtkGenericNearestNeighborInterpolateImageFunction<irtkByteImage> interp;
  interp.Input(&image);
  interp.Initialize();
  irtkBytePixel *out = output.Data();
  double x, y, z;
  for (int k = 0; k < attr._z; ++k)
  for (int j = 0; j < attr._y; ++j)
  for (int i = 0; i < attr._x; ++i) {
    x = i, y = j, z = k;
    attr.LatticeToWorld(x, y, z);
    image.WorldToImage (x, y, z);
    (*out++) = interp.Get(x, y, z);
  }
  image = output;
}

// -----------------------------------------------------------------------------
/// Set value of voxels within specified segments to given padding value
void MaskSegments(irtkRealImage &image, const irtkByteImage &seg,
                  set<irtkBytePixel> labels, irtkRealPixel padding, bool inv = false)
{
  irtkRealPixel       *value = image.Data();
  const irtkBytePixel *label = seg.Data();
  for (int i = 0; i < image.NumberOfSpatialVoxels(); ++i, ++value, ++label) {
    if (labels.find(*label) == labels.end()) {
      if ( inv) *value = padding;
    } else {
      if (!inv) *value = padding;
    }
  }
}

// -----------------------------------------------------------------------------
/// Resample image on specified image lattice
void Resample(irtkRealImage &image, const irtkImageAttributes &attr)
{
  irtkRealImage output(attr);
  irtkGenericLinearInterpolateImageFunction<irtkRealImage> interp;
  interp.Input(&image);
  interp.Initialize();
  irtkRealPixel *out = output.Data();
  double x, y, z;
  for (int k = 0; k < attr._z; ++k)
  for (int j = 0; j < attr._y; ++j)
  for (int i = 0; i < attr._x; ++i) {
    x = i, y = j, z = k;
    attr.LatticeToWorld(x, y, z);
    image.WorldToImage (x, y, z);
    (*out++) = interp.Get(x, y, z);
  }
  image = output;
}

// -----------------------------------------------------------------------------
/// Voxel function for parallel application of image Laplace operator
struct LaplaceFunction : public irtkVoxelReduction
{
  irtkRealImage *_Psi;
  irtkBytePixel  _Label;
  double         _Delta;
  int            _Changed;

  LaplaceFunction(irtkRealImage *psi = NULL, irtkBytePixel label = cGM)
  :
    _Psi(psi), _Label(label), _Delta(.0), _Changed(0)
  {}

  void split(LaplaceFunction)
  {
    _Delta   = .0;
    _Changed =  0;
  }

  void join(LaplaceFunction &rhs)
  {
    _Delta   += rhs._Delta;
    _Changed += rhs._Changed;
  }

  void operator ()(int i, int j, int k, int, const irtkBytePixel *label, irtkRealPixel *psi)
  {
    if (*label != _Label) return;

    irtkRealPixel value = 0;
    int           num   = 0;

    if (i > 0           && _Psi->Get(i-1, j, k) != -1) { value += _Psi->Get(i-1, j, k); ++num; }
    if (i < _Psi->X()-1 && _Psi->Get(i+1, j, k) != -1) { value += _Psi->Get(i+1, j, k); ++num; }
    if (j > 0           && _Psi->Get(i, j-1, k) != -1) { value += _Psi->Get(i, j-1, k); ++num; }
    if (j < _Psi->Y()-1 && _Psi->Get(i, j+1, k) != -1) { value += _Psi->Get(i, j+1, k); ++num; }
    if (k > 0           && _Psi->Get(i, j, k-1) != -1) { value += _Psi->Get(i, j, k-1); ++num; }
    if (k < _Psi->Z()-1 && _Psi->Get(i, j, k+1) != -1) { value += _Psi->Get(i, j, k+1); ++num; }

    if (num > 0) {
      value /= num;
      irtkRealPixel preval = _Psi->Get(i, j, k);
      if (preval < 0) preval = 0;
      _Delta += abs(preval - value) / 10000;
      *psi = value;
      ++_Changed;
    }
  }
};

// -----------------------------------------------------------------------------
/// Applies Laplace image operator within specified tissue segment
double Laplace(const irtkByteImage &seg, irtkRealImage &psi, irtkBytePixel label)
{
  LaplaceFunction laplace(new irtkRealImage(psi), label);
  ParallelForEachVoxel(seg.Attributes(), seg, psi, laplace);
  delete laplace._Psi;
  if (laplace._Changed) laplace._Delta /= laplace._Changed;
  return laplace._Delta;
}

// -----------------------------------------------------------------------------
/// Iteratively applies Laplace operator within cortical grey matter
irtkRealImage ComputeLaplacian(const irtkByteImage &seg)
{
  IRTK_START_TIMING();

  // Initialize psi image
  irtkRealImage psi(seg.Attributes());
  const irtkBytePixel *segval = seg.Data();
  irtkRealPixel       *psival = psi.Data();

  for (int k = 0; k < seg.Z(); ++k)
  for (int j = 0; j < seg.Y(); ++j)
  for (int i = 0; i < seg.X(); ++i) {
    switch (*segval) {
      case Outside: case BG: case CSF: case CB: *psival = 10000; break;
      case cGM:                                 *psival =    -1; break;
      default:                                  *psival =     0; break;
    }
    ++segval, ++psival;
  }

  // Iteratively apply Laplace operator
  double precon      = Laplace(seg, psi, cGM);
  double convergence = 1.0;
  if (verbose) cout << "  0 Laplace convergence = " << precon << endl;
  for (int iter = 0; convergence > .0001 && iter < 2000; ++iter) {
    double con  = Laplace(seg, psi, cGM);
    convergence = 1.0 - con / precon;
    precon      = con;
    if (verbose) {
      cout << setw(3)<< right << (iter+1) << left
           << " Laplace convergence = " << convergence << " (delta = " << con << ")" << endl;
    }
    if (debug && iter % 20 == 0) {
      char fname[64];
      sprintf(fname, "polydatacortex_psi_%03d.nii.gz", iter+1);
      psi.Write(fname);
    }
  }

  // Finalize psi image
  psival = psi.Data();
  for (int k = 0; k < seg.Z(); ++k)
  for (int j = 0; j < seg.Y(); ++j)
  for (int i = 0; i < seg.X(); ++i) {
    if (*psival == -1) *psival = .0;
    ++psival;
  }

  IRTK_DEBUG_TIMING(1, "computing Laplacian");
  return psi;
}

// -----------------------------------------------------------------------------
/// Voxel function for parallel convolution of image with Laplacian kernel
struct LaplacianFunction : public irtkVoxelFunction
{
  const irtkRealImage *_Image;
  LaplacianFunction(const irtkRealImage *image = NULL)
  :
    _Image(image)
  {}

  void operator ()(int i, int j, int k, int, irtkRealPixel *out) const
  {
    const irtkRealImage &f = *_Image;
    if (0 < i && i < f.X()-1 && 0 < j && j < f.Y()-1 && 0 < k && k < f.Z()-1) {
      *out = .25 * (  f(i-1, j, k) + f(i+1, j, k)
                    + f(i, j-1, k) + f(i, j+1, k)
                    + f(i, j, k-1) + f(i, j, k+1) - 6.0 * f(i, j, k));
    }
  }
};

// -----------------------------------------------------------------------------
/// Convolve intensity image by the finite difference Laplacian kernel
irtkRealImage ComputeLaplacian(const irtkRealImage &image, double sigma = .0)
{
  irtkRealImage *input = const_cast<irtkRealImage *>(&image);
  irtkRealImage output(image);
  if (sigma > .0) {
    input = new irtkRealImage;
    irtkGaussianBlurring<irtkRealPixel> blur(sigma);
    blur.SetInput (const_cast<irtkRealImage *>(&image));
    blur.SetOutput(input);
    blur.Run();
  }
  ParallelForEachVoxel(LaplacianFunction(input), output.Attributes(), output);
  if (input != &image) delete input;
  return output;
}

// -----------------------------------------------------------------------------
/// http://en.wikipedia.org/wiki/Edge_detection#Differential_edge_detection
struct DifferentialEdgeFunction : public irtkVoxelFunction
{
  const irtkRealImage *_Image;
  DifferentialEdgeFunction(const irtkRealImage *image = NULL)
  :
    _Image(image)
  {}

  void operator ()(int i, int j, int k, int, irtkRealPixel *out) const
  {
    const irtkRealImage &f = *_Image;
    if (0 < i && i < f.X()-1 && 0 < j && j < f.Y()-1 && 0 < k && k < f.Z()-1) {
      // First order derivatives
      double fx = .5 * (f(i+1, j, k) - f(i-1, j, k));
      double fy = .5 * (f(i, j+1, k) - f(i, j-1, k));
      double fz = .5 * (f(i, j, k+1) - f(i, j, k-1));
      // Second order derivatives
      double fijk = f(i, j, k);
      double fxx = .25 * (f(i+1, j,   k  ) - fijk             - fijk             + f(i-1, j,   k  ));
      double fxy = .25 * (f(i+1, j+1, k  ) - f(i-1, j+1, k  ) - f(i+1, j-1, k  ) + f(i-1, j-1, k  ));
      double fxz = .25 * (f(i+1, j,   k+1) - f(i-1, j,   k+1) - f(i+1, j,   k-1) + f(i-1, j,   k-1));
      double fyy = .25 * (f(i,   j+1, k  ) - fijk             - fijk             + f(i,   j-1, k  ));
      double fyz = .25 * (f(i,   j+1, k+1) - f(i,   j-1, k+1) - f(i,   j+1, k-1) + f(i,   j-1, k-1));
      double fzz = .25 * (f(i,   j,   k+1) - fijk             - fijk             + f(i,   j,   k-1));
      // Second order derivative in the direction of the image gradient
      *out =        fx * fx * fxx + fy * fy * fyy + fz * fz * fzz
           + 2.0 * (fx * fy * fxy + fx * fz * fxz + fy * fz * fyz);
      if (fabs(*out) > 1000.) *out = copysign(1000.0, *out);
    }
  }
};

// -----------------------------------------------------------------------------
/// http://en.wikipedia.org/wiki/Edge_detection#Differential_edge_detection
irtkRealImage ComputeDifferentialEdgeFunction(const irtkRealImage &image, double sigma = .0)
{
  irtkRealImage *input = const_cast<irtkRealImage *>(&image);
  irtkRealImage output(image);
  if (sigma > .0) {
    irtkGaussianBlurring<irtkRealPixel> blur(sigma);
    blur.SetInput (input);
    blur.SetOutput(&output);
    blur.Run();
    input = new irtkRealImage(output);
  }
  ParallelForEachVoxel(DifferentialEdgeFunction(input), output.Attributes(), output);
  if (input != &image) delete input;
  return output;
}

// -----------------------------------------------------------------------------
/// Get convex hull of input points
vtkSmartPointer<vtkPolyData> ConvexHull(vtkSmartPointer<vtkPolyData> input)
{
  // Spatially stratify points to prevent "Unable to factor linear system"
  // warning of vtkDelaunay3D filter due to numerical imprecisions
  vtkSmartPointer<vtkMaskPoints> stratify;
  stratify = vtkSmartPointer<vtkMaskPoints>::New();
  stratify->RandomModeOn();
  stratify->SetRandomModeType(2);
  stratify->SetMaximumNumberOfPoints(.75 * input->GetNumberOfPoints());
  SetVTKInput(stratify, input);

  // Get convex hull of largest component
  vtkSmartPointer<vtkHull> hull;
  hull = vtkSmartPointer<vtkHull>::New();
  hull->AddRecursiveSpherePlanes(3);
  SetVTKConnection(hull, stratify);

  // Compute Delaunay triangulation
  vtkSmartPointer<vtkDelaunay3D> delaunay;
  delaunay = vtkSmartPointer<vtkDelaunay3D>::New();
  SetVTKConnection(delaunay, hull);

  // Construct surface mesh
  vtkSmartPointer<vtkDataSetSurfaceFilter> mesher;
  mesher = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
  SetVTKConnection(mesher, delaunay);

  mesher->Update();
  return mesher->GetOutput();
}

// -----------------------------------------------------------------------------
/// Extract initial isosurface from psi image
vtkSmartPointer<vtkPolyData>
InitialSurface(const irtkRealImage &psi, double psi1, bool hull = false)
{
  IRTK_START_TIMING();
  vtkSmartPointer<vtkPolyData> mesh;

  // Use image with default attributes because ImageToVTK ignores the orientation
  irtkRealImage mypsi(irtkImageAttributes(psi.X(), psi.Y(), psi.Z()),
                      const_cast<irtkRealPixel *>(psi.Data()));

  // Convert image to structured points
  vtkSmartPointer<vtkStructuredPoints> input;
  input = vtkSmartPointer<vtkStructuredPoints>::New();
  mypsi.ImageToVTK(input);

  // Extract surface using marching cubes
  vtkSmartPointer<vtkMarchingCubes> mcubes;
  mcubes = vtkSmartPointer<vtkMarchingCubes>::New();
  mcubes->SetValue(0, psi1);
  mcubes->ComputeNormalsOff();
  mcubes->ComputeGradientsOff();
  SetVTKInput(mcubes, input);

  // Extract largest component
  vtkSmartPointer<vtkPolyDataConnectivityFilter> lcc;
  lcc = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
  lcc->SetExtractionModeToLargestRegion();
  SetVTKConnection(lcc, mcubes);

  // Clean surface mesh
  vtkSmartPointer<vtkCleanPolyData> cleaner;
  cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
  cleaner->ConvertLinesToPointsOn();
  cleaner->ConvertPolysToLinesOn();
  cleaner->ConvertStripsToPolysOn();
  cleaner->PointMergingOn();
  SetVTKConnection(cleaner, lcc);

  // Triangulate surface mesh
  vtkSmartPointer<vtkTriangleFilter> triangulator;
  triangulator = vtkSmartPointer<vtkTriangleFilter>::New();
  triangulator->PassVertsOff();
  triangulator->PassLinesOff();
  SetVTKConnection(triangulator, cleaner);

  // Get cleaned and triangulated largest component
  triangulator->Update();
  mesh = triangulator->GetOutput();
  if (mesh == NULL || mesh->GetNumberOfPoints() == 0 || mesh->GetPolys()->GetNumberOfCells() == 0) {
    cerr << "Error: Failed to extract initial surface mesh (psi=" << psi1 << ")" << endl;
    exit(1);
  }

  // Optionally, compute convex hull of surface points
  if (hull) mesh = ConvexHull(mesh);

  // Smooth initial surface mesh
  vtkSmartPointer<vtkSmoothPolyDataFilter> smoother;
  smoother = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
  smoother->SetNumberOfIterations(100);
  SetVTKInput(smoother, mesh);

  // Get initial surface mesh
  smoother->Update();
  mesh = smoother->GetOutput();

  // Convert point coordinates to original world coordinates
  double w[3];
  for (vtkIdType id = 0; id < mesh->GetNumberOfPoints(); ++id) {
    mesh->GetPoint(id, w);
    mypsi.WorldToImage(w[0], w[1], w[2]);
    psi  .ImageToWorld(w[0], w[1], w[2]);
    mesh->GetPoints()->SetPoint(id, w);
  }

  IRTK_DEBUG_TIMING(1, "initializing surface");
  return mesh;
}

// -----------------------------------------------------------------------------
/// Compute direction of propagation
irtkRealImage PropagationDirection(vtkPolyData *init, const irtkRealImage &psi, irtkRealPixel psi0, bool outside = false)
{
  IRTK_START_TIMING();
  irtkRealImage grad;

  // Compute normalized gradient of input image
  int type = irtkGradientImageFilter<irtkRealPixel>::NORMALISED_GRADIENT_VECTOR;
  irtkGradientImageFilter<irtkRealPixel> gradient(type);
  gradient.SetInput(const_cast<irtkRealImage *>(&psi));
  gradient.SetOutput(&grad);
  gradient.UseVoxelSize  (false);
  gradient.UseOrientation(false);
  gradient.Run();

  if (outside) {

    // Convert image to structured points
    irtkRealImage mypsi(irtkImageAttributes(psi.X(), psi.Y(), psi.Z()),
                        const_cast<irtkRealPixel *>(psi.Data()));
    vtkSmartPointer<vtkStructuredPoints> input;
    input = vtkSmartPointer<vtkStructuredPoints>::New();
    mypsi.ImageToVTK(input);

    // Extract surface using marching cubes
    vtkSmartPointer<vtkMarchingCubes> mcubes;
    mcubes = vtkSmartPointer<vtkMarchingCubes>::New();
    mcubes->SetValue(0, psi0);
    mcubes->ComputeNormalsOff();
    mcubes->ComputeGradientsOff();
    SetVTKInput(mcubes, input);

    mcubes->Update();

    // TODO Convert extracted surface to binary image instead and use
    //      much more efficient Euclidean distance transform (EDT).
    vtkSmartPointer<vtkPolyData> mesh = mcubes->GetOutput();
    mesh->BuildLinks();
    mesh->BuildCells();
    mcubes = NULL;

    // Convert points to actual world coordinates
    double w[3];
    for (vtkIdType id = 0; id < mesh->GetNumberOfPoints(); ++id) {
      mesh->GetPoint(id, w);
      mypsi.WorldToImage(w[0], w[1], w[2]);
      psi  .ImageToWorld(w[0], w[1], w[2]);
      mesh->GetPoints()->SetPoint(id, w);
    }

    // Use normalized distance to closest point on isosurface if gradient is zero
    vtkSmartPointer<vtkCellLocator> locator;
    locator = vtkSmartPointer<vtkCellLocator>::New();
    locator->SetDataSet(mesh);
    locator->BuildLocator();

    irtkRealPixel *gx = grad.Data(0, 0, 0, 0);
    irtkRealPixel *gy = grad.Data(0, 0, 0, 1);
    irtkRealPixel *gz = grad.Data(0, 0, 0, 2);

    double bbox[6];
    init->GetBounds(bbox);

    vtkIdType cellId;
    int       subId;
    double p1[3], p2[3], gm, dist2;
    for (int k = 0; k < psi.Z(); ++k)
    for (int j = 0; j < psi.Y(); ++j)
    for (int i = 0; i < psi.X(); ++i) {
      if (*gx == 0 && *gy == 0 && *gz == 0) {
        p1[0] = i, p1[1] = j, p1[2] = k;
        psi.ImageToWorld(p1[0], p1[1], p1[2]);
        if (bbox[0] <= p1[0] && p1[0] <= bbox[1] &&
            bbox[2] <= p1[1] && p1[1] <= bbox[3] &&
            bbox[4] <= p1[2] && p1[2] <= bbox[5]) {
          locator->FindClosestPoint(p1, p2, cellId, subId, dist2);
          *gx = p2[0] - p1[0];
          *gy = p2[1] - p1[1];
          *gz = p2[2] - p1[2];
           gm = sqrt((*gx) * (*gx) + (*gy) * (*gy) + (*gz) * (*gz));
          if (gm) *gx /= gm, *gy /= gm, *gz /= gm;
        }
      }
      ++gx, ++gy, ++gz;
    }

  }

  IRTK_DEBUG_TIMING(1, "computing propagation direction");
  return grad;
}

// -----------------------------------------------------------------------------
/// Get IDs of adjacent points, ensure BuildLinks and BuildCells has been called!
void GetAdjacentPoints(vtkPolyData *poly, vtkIdType ptId, set<vtkIdType> &otherPtIds)
{
  otherPtIds.clear();
  vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
  poly->GetPointCells(ptId, cellIds);
  vtkIdType cellId, npts, *ptIds = NULL;
  for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {
    cellId = cellIds->GetId(i);
    poly->GetCellPoints(cellId, npts, ptIds);
    for (vtkIdType j = 0; j < npts; ++j) otherPtIds.insert(ptIds[j]);
  }
  otherPtIds.erase(ptId);
}

// -----------------------------------------------------------------------------
/// Get IDs of triangles adjacent to a point
void GetPointTriangles(vtkPolyData *poly, vtkIdType ptId, vector<irtkVector3D<vtkIdType> > &triPtIds)
{
  vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
  poly->GetPointCells(ptId, cellIds);
  vtkIdType cellId, npts, *ptIds = NULL;
  triPtIds.resize(cellIds->GetNumberOfIds());
  for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {
    cellId = cellIds->GetId(i);
    poly->GetCellPoints(cellId, npts, ptIds);
    triPtIds[i] = irtkVector3D<vtkIdType>(ptIds[0], ptIds[1], ptIds[2]);
  }
}

// -----------------------------------------------------------------------------
/// Get IDs of triangles close to a point
///
/// \attention The vtkAbstractCellLocator::FindCellsWithinBounds function
///            is currently not implemented by any of the subclasses! Therefore,
///            this function is at the moment useless and will raise an error.
void GetNearbyTriangles(vtkPolyData *poly, vtkCellLocator *cellLoc, vtkIdType ptId,
                        vector<irtkVector3D<vtkIdType> > &triPtIds, double radius)
{
  double bbox[6];
  poly->GetPoint(ptId, bbox);
  bbox[5] = bbox[2] + radius;
  bbox[4] = bbox[2] - radius;
  bbox[3] = bbox[1] + radius;
  bbox[2] = bbox[1] - radius;
  bbox[1] = bbox[0] + radius;
  bbox[0] = bbox[0] - radius;

  vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
  cellLoc->FindCellsWithinBounds(bbox, cellIds);

  vtkIdType cellId, npts, *ptIds = NULL;
  triPtIds.resize(cellIds->GetNumberOfIds());
  for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {
    cellId = cellIds->GetId(i);
    poly->GetCellPoints(cellId, npts, ptIds);
    triPtIds[i] = irtkVector3D<vtkIdType>(ptIds[0], ptIds[1], ptIds[2]);
  }
}

// -----------------------------------------------------------------------------
/// Move points without checking for collisions
struct MovePoints
{
  GradientInterpolator *_Direction;
  ImageInterpolator    *_Psi;
  double                _Psi0;
  vtkPolyData          *_Mesh;
  double                _MinStep;

  void operator ()(const blocked_range<vtkIdType> &ids) const
  {
    vtkPoints    * const points = _Mesh->GetPoints();
    vtkDataArray * const normals = _Mesh->GetPointData()->GetNormals();
    vtkDataArray * const step   = _Mesh->GetPointData()->GetArray("step");

    double h, w1[3], w2[3], v[3], n[3], psival1, psival2;
    for (vtkIdType id = ids.begin(); id != ids.end(); ++id) {

      // Get current step length
      step->GetTuple(id, &h);

      // Skip vertices that cannot be moved further
      if (h == .0) continue;
      if (fabs(h) < _MinStep) {
        step->SetTuple1(id, .0);
        continue;
      }

      // Get current vertex position
      points->GetPoint(id, w1);
      memcpy(v, w1, 3 * sizeof(double));
      _Psi->WorldToImage(v[0], v[1], v[2]);

      // Keep position if close enough to isosurface
      psival1 = _Psi->Get(v[0], v[1], v[2]);
      if (fabs(_Psi0 - psival1) < 10.0) {
        step->SetTuple1(id, .0);
        continue;
      }

      // Move vertex outwards (psi0 > value) or inwards (psi0 < value)
      h = copysign(h, _Psi0 - psival1);
      if (_Direction) {
        _Direction->Evaluate(n, v[0], v[1], v[2]);
        if (n[0] == .0 && n[1] == .0 && n[2] == .0) {
          normals->GetTuple(id, n);
        }
      } else {
        normals->GetTuple(id, n);
      }
      for (int d = 0; d < 3; ++d) w2[d] = w1[d] + h * n[d];
      memcpy(v, w2, 3 * sizeof(double));
      _Psi->WorldToImage (v[0], v[1], v[2]);
      psival2 = _Psi->Get(v[0], v[1], v[2]);

      // Ensure step is not exceeding the target psi value
      while ((_Psi0 - psival1) * (_Psi0 - psival2) < .0) {
        h *= .9;
        if (fabs(h) < _MinStep) break;
        for (int d = 0; d < 3; ++d) w2[d] = w1[d] + h * n[d];
        memcpy(v, w2, 3 * sizeof(double));
        _Psi->WorldToImage (v[0], v[1], v[2]);
        psival2 = _Psi->Get(v[0], v[1], v[2]);
      }
      if (fabs(h) < _MinStep) {
        step->SetTuple1(id, .0);
        continue;
      }

      // Set new vertex position
      points->SetPoint (id, w2);
      step  ->SetTuple1(id, h );
    }
  }
};

// -----------------------------------------------------------------------------
/// Check for collisions such as self-intersections and faces too close
struct FindCollisions
{
  #if SELFINTERSECTION_CHECK_NORMALS || SELFINTERSECTION_CHECK_EDGES
  vtkAbstractCellLocator *_CellLocator;
  #endif

  #if SELFINTERSECTION_CHECK_TRIANGLES
  const vector<set<irtkVector3D<vtkIdType> > > *_Triangles;
  const vector<set<irtkVector3D<vtkIdType> > > *_NearbyTriangles;
  #endif

  GradientInterpolator *_Direction;
  double                _Psi0;
  vtkPolyData          *_Mesh;
  double                _MinStep;
  double                _MinGap;
  double                _MinWidth;

  void operator ()(const blocked_range<vtkIdType> &ids) const
  {
    vtkPoints    * const points  = _Mesh->GetPoints();
    #if SELFINTERSECTION_CHECK_NORMAL || SELFINTERSECTION_CHECK_EDGES
    vtkDataArray * const normals = _Mesh->GetPointData()->GetNormals();
    #endif
    vtkDataArray * const step    = _Mesh->GetPointData()->GetArray("step");
    vtkDataArray * const coll    = _Mesh->GetPointData()->GetArray("coll");

    double h, p1[3], p2[3];

    #if SELFINTERSECTION_CHECK_NORMAL
    double w[3], n[3], v[3];
    #endif // SELFINTERSECTION_CHECK_NORMAL

    #if SELFINTERSECTION_CHECK_EDGES
    double t, x[3], pcoords[3];
    vtkIdType cellId;
    int subId;
    #endif // SELFINTERSECTION_CHECK_EDGES

    #if SELFINTERSECTION_CHECK_TRIANGLES
    const vector<set<irtkVector3D<vtkIdType> > > &triangles        = *_Triangles;
    const vector<set<irtkVector3D<vtkIdType> > > &nearby_triangles = *_NearbyTriangles;
    set<irtkVector3D<vtkIdType> >::const_iterator triangle,        triangle_end;
    set<irtkVector3D<vtkIdType> >::const_iterator nearby_triangle, nearby_triangle_end;
    bool   intersects, have_shared_edge, have_shared_vert;
    double tri1[3][3], tri2[3][3], n1[3], n2[3];
    int    tri12[3], coplanar;
    #endif // SELFINTERSECTION_CHECK_TRIANGLES

    for (vtkIdType id = ids.begin(); id != ids.end(); ++id) {

      // Get current step length
      step->GetTuple(id, &h);

      // Quick check for self intersection (inaccurate!)
      #if SELFINTERSECTION_CHECK_NORMAL || SELFINTERSECTION_CHECK_EDGES
      points->GetPoint(id, w);
      if (_Direction) {
        memcpy(v, w, 3 * sizeof(double));
        _Psi      ->WorldToImage(v[0], v[1], v[2]);
        _Direction->Evaluate(n,  v[0], v[1], v[2]);
        if (n[0] == .0 && n[1] == .0 && n[2] == .0) {
          normals->GetTuple(id, n);
        }
      } else {
        normals->GetTuple(id, n);
      }
      p2[0] = w[0] - (_MinStep + 1e-12) * n[0]; // near new position
      p2[1] = w[1] - (_MinStep + 1e-12) * n[1];
      p2[2] = w[2] - (_MinStep + 1e-12) * n[2];
      #endif // SELFINTERSECTION_CHECK_NORMAL || SELFINTERSECTION_CHECK_EDGES

      #if SELFINTERSECTION_CHECK_NORMAL
      // Find intersection of line segment along vertex displacement
      p1[0] = w[0] - (_MinStep + h) * n[0]; // near old position
      p1[1] = w[1] - (_MinStep + h) * n[1];
      p1[2] = w[2] - (_MinStep + h) * n[2];
      if (_CellLocator->IntersectWithLine(p1, p2, _MinStep, t, x, pcoords, subId, cellId) != 0) {
        step->SetTuple1(id, .5 * h);
        coll->SetTuple1(id, SELF_INTERSECTION);
        continue;
      }
      #endif // SELFINTERSECTION_CHECK_NORMAL

      #if SELFINTERSECTION_CHECK_EDGES
      GetAdjacentPoints(_Mesh, id, neighbors);
      for (neighbor = neighbors.begin(); neighbor != neighbors.end(); ++neighbor) {
        vertices->GetPoint(*neighbor, p1);
        p1[0] += (_MinStep + 1e-12) * (p2[0] - p1[0]);
        p1[1] += (_MinStep + 1e-12) * (p2[1] - p1[1]);
        p1[2] += (_MinStep + 1e-12) * (p2[2] - p1[2]);
        if (_CellLocator->IntersectWithLine(p1, p2, _MinStep, t, x, pcoords, subId, cellId) != 0) {
          break;
        }
      }
      if (neighbor != neighbors.end()) {
        step->SetTuple1(id, .5 * h);
        coll->SetTuple1(id, SELF_INTERSECTION);
        continue;
      }
      #endif // SELFINTERSECTION_CHECK_EDGES

      // More thorough check for any self intersection
      #if SELFINTERSECTION_CHECK_TRIANGLES
      //GetPointTriangles(_Mesh, id, triangles); // custom structure more efficient
      triangle     = triangles[id].begin();
      triangle_end = triangles[id].end();
      intersects   = false;
      while (!intersects && triangle != triangle_end) {
        points->GetPoint(triangle->_x, tri1[0]);
        points->GetPoint(triangle->_y, tri1[1]);
        points->GetPoint(triangle->_z, tri1[2]);
        vtkTriangle::ComputeNormal(tri1[0], tri1[1], tri1[2], n1);
        //GetNearbyTriangles(_Mesh, _CellLocator, id, nearby_triangles, _Radius);
        nearby_triangle     = nearby_triangles[id].begin();
        nearby_triangle_end = nearby_triangles[id].end();
        while (nearby_triangle != nearby_triangle_end) {
          points->GetPoint(nearby_triangle->_x, tri2[0]);
          points->GetPoint(nearby_triangle->_y, tri2[1]);
          points->GetPoint(nearby_triangle->_z, tri2[2]);
          vtkTriangle::ComputeNormal(tri2[0], tri2[1], tri2[2], n2);
          // Get corresponding indices of shared vertices
          tri12[0] = -1, tri12[1] = -1, tri12[2] = -1;
          if (triangle->_x == nearby_triangle->_x) tri12[0] = 0;
          if (triangle->_x == nearby_triangle->_y) tri12[0] = 1;
          if (triangle->_x == nearby_triangle->_z) tri12[0] = 2;
          if (triangle->_y == nearby_triangle->_x) tri12[1] = 0;
          if (triangle->_y == nearby_triangle->_y) tri12[1] = 1;
          if (triangle->_y == nearby_triangle->_z) tri12[1] = 2;
          if (triangle->_z == nearby_triangle->_x) tri12[2] = 0;
          if (triangle->_z == nearby_triangle->_y) tri12[2] = 1;
          if (triangle->_z == nearby_triangle->_z) tri12[2] = 2;
          // Determine whether triangles share one or two vertices
          have_shared_edge = false;
          have_shared_vert = false;
          for (int i = 0; i < 3; ++i) {
            if (tri12[i] != -1) {
              have_shared_vert = true;
              for (int j = i+1; j < 3; ++j) {
                if (tri12[j] != -1) {
                  have_shared_edge = true;
                  break;
                }
              }
              if (have_shared_edge) break;
            }
          }
          // Compute intersection of both triangles
          if (!have_shared_edge &&
              vtkIntersectionPolyDataFilter::TriangleTriangleIntersection(
                  tri1[0], tri1[1], tri1[2],
                  tri2[0], tri2[1], tri2[2], coplanar, p1, p2)) {
            // Ignore valid intersection of single shared vertex
            if (!have_shared_vert || fabs(p1[0] - p2[0]) > 1e-6
                                  || fabs(p1[1] - p2[1]) > 1e-6
                                  || fabs(p1[2] - p2[2]) > 1e-6) {
              intersects = true;
              break;
            }
          }
          coll->SetTuple1(id, SELF_INTERSECTION);
          // Check if triangles are nearly touching
          #if SELFINTERSECTION_CHECK_TOUCHING
          if (!have_shared_edge && !have_shared_vert) {
            const double mingap = (vtkMath::Dot(n1, n2) < .0 ? _MinWidth : _MinGap);
            for (int i = 0; i < 3; ++i) {
              if (vtkPlane::DistanceToPlane(tri1[i], n2, tri2[0]) < mingap &&
                  vtkTriangle::PointInTriangle(tri1[i], tri2[0], tri2[1], tri2[2], 1e-6)) {
                intersects = true;
                break;
              }
            }
            if (!intersects) {
              for (int i = 0; i < 3; ++i) {
                if (vtkPlane::DistanceToPlane(tri2[i], n1, tri1[0]) < mingap &&
                    vtkTriangle::PointInTriangle(tri2[i], tri1[0], tri1[1], tri1[2], 1e-6)) {
                  intersects = true;
                  break;
                }
              }
            }
            if (intersects) {
              coll->SetTuple1(id, TOUCHING);
              break;
            }
          }
          #endif
          ++nearby_triangle;
        }
        ++triangle;
      }
      if (intersects) {
        step->SetTuple1(id, .5 * h);
        continue;
      }
      #endif // SELFINTERSECTION_CHECK_TRIANGLES

      // No collision detected
      coll->SetTuple1(id, NO_COLLISION);
    }
  }
};

// -----------------------------------------------------------------------------
/// Smooth step lengths by averaging step lengths of adjacent vertices
struct SmoothSteps
{
  vtkPolyData  *_Mesh;
  vtkDataArray *_InputStep;
  vtkDataArray *_OutputStep;
  double        _MinStep;

  void operator ()(const blocked_range<vtkIdType> &ids) const
  {
    set<vtkIdType>                 neighbors;
    set<vtkIdType>::const_iterator neighbor;
    double h, nh;

    for (vtkIdType id = ids.begin(); id != ids.end(); ++id) {
      _InputStep->GetTuple(id, &h);
      if (h != .0) {
        GetAdjacentPoints(_Mesh, id, neighbors);
        for (neighbor = neighbors.begin(); neighbor != neighbors.end(); ++neighbor) {
          _InputStep->GetTuple(*neighbor, &nh);
          h += nh;
        }
        h /= (neighbors.size() + 1);
        if (fabs(h) < _MinStep) h = copysign(_MinStep, h);
      }
      _OutputStep->SetTuple1(id, h);
    }
  }
};

// -----------------------------------------------------------------------------
/// Update positions of vertices
struct UpdatePoints
{
  GradientInterpolator *_Direction;
  const irtkRealImage  *_Psi;
  vtkPolyData          *_Mesh;
  double                _MaxTotalDistance;
  vtkDataArray         *_PreStep;
  double                _SumStep, _MinStep, _MaxStep;
  double                _SumDist, _MinDist, _MaxDist;
  int                   _NumActive, _NumAccept;

  UpdatePoints()
  :
    _SumStep(.0),
    _MinStep(+ numeric_limits<double>::infinity()),
    _MaxStep(- numeric_limits<double>::infinity()),
    _SumDist(.0),
    _MinDist(+ numeric_limits<double>::infinity()),
    _MaxDist(- numeric_limits<double>::infinity()),
    _NumActive(0), _NumAccept(0)
  {}

  UpdatePoints(UpdatePoints &lhs, split)
  :
    _Direction(lhs._Direction),
    _Psi(lhs._Psi),
    _Mesh(lhs._Mesh),
    _MaxTotalDistance(lhs._MaxTotalDistance),
    _PreStep(lhs._PreStep),
    _SumStep(.0),
    _MinStep(+ numeric_limits<double>::infinity()),
    _MaxStep(- numeric_limits<double>::infinity()),
    _SumDist(.0),
    _MinDist(+ numeric_limits<double>::infinity()),
    _MaxDist(- numeric_limits<double>::infinity()),
    _NumActive(0),
    _NumAccept(0)
  {}

  void join(UpdatePoints &rhs)
  {
    _SumStep   += rhs._SumStep;
    _SumDist   += rhs._SumDist;
    _MinStep    = min(_MinStep, rhs._MinStep);
    _MaxStep    = max(_MaxStep, rhs._MaxStep);
    _MinDist    = min(_MinDist, rhs._MinDist);
    _MaxDist    = max(_MaxDist, rhs._MaxDist);
    _NumActive += rhs._NumActive;
    _NumAccept += rhs._NumAccept;
  }

  void operator()(const blocked_range<vtkIdType> &ids)
  {
    vtkPoints    * const points  = _Mesh->GetPoints();
    vtkDataArray * const normals = _Mesh->GetPointData()->GetNormals();
    vtkDataArray * const step    = _Mesh->GetPointData()->GetArray("step");
    vtkDataArray * const dist    = _Mesh->GetPointData()->GetArray("dist");
    vtkDataArray * const coll    = _Mesh->GetPointData()->GetArray("coll");

    double w[3], v[3], n[3], h0, h, d, c, total_dist;

    for (vtkIdType id = ids.begin(); id != ids.end(); ++id) {

      // Skip passive vertices and those that collided with the mesh
      step->GetTuple(id, &h);
      if (h == .0) continue;

      // Reject step if maximum allowed displacement exceeded
      d = sqrt(pow(h * n[0], 2) + pow(h * n[1], 2) + pow(h * n[2], 2));
      dist->GetTuple(id, &total_dist);
      total_dist += d;
      if (_MaxTotalDistance > .0 && total_dist > _MaxTotalDistance) {
        step->SetTuple1(id, .9 * h);
        continue;
      }
      dist->SetTuple1(id, total_dist);

      ++_NumActive;

      // Check if collision was detected
      coll->GetTuple(id, &c);
      if (static_cast<int>(c) == NO_COLLISION) {
        ++_NumAccept;
      } else {
        h = .0; // revert previous step
      }

      // Adjust vertex position
      points->GetPoint(id, w);
      if (_Direction) {
        memcpy(v, w, 3 * sizeof(double));
        _Psi      ->WorldToImage(v[0], v[1], v[2]);
        _Direction->Evaluate(n,  v[0], v[1], v[2]);
        if (n[0] == .0 && n[1] == .0 && n[2] == .0) {
          normals->GetTuple(id, n);
        }
      } else {
        normals->GetTuple(id, n);
      }
      _PreStep->GetTuple(id, &h0);
      w[0] += (h - h0) * n[0];
      w[1] += (h - h0) * n[1];
      w[2] += (h - h0) * n[2];
      points->SetPoint(id, w);

      // Do not update ranges if previous step was rejected
      if (h == .0) continue;

      // Update range of steps performed
      h = fabs(h);
      if (h < _MinStep) _MinStep = h;
      if (h > _MaxStep) _MaxStep = h;
      _SumStep += h;

      // Update range of displacements
      if (d > _MaxDist) _MaxDist = d;
      if (d < _MinDist) _MinDist = d;
      _SumDist += d;
    }
  }
};

// -----------------------------------------------------------------------------
/// Extract isosurface from specified image
vtkSmartPointer<vtkPolyData>
RefineSurface(vtkSmartPointer<vtkPolyData> mesh,
              const irtkRealImage &psi, const irtkRealImage &grad,
              double psi0, double maxh, double minw, double maxd, int maxit,
              double minactive, double edgelength, double smoothing,
              bool recalcnormals)
{
  const int _verbose = verbose;
  IRTK_START_TIMING();

  // Minimum voxel size (should be isotropic though)
  const double ds = min(min(psi.GetXSize(), psi.GetYSize()), psi.GetZSize());
  // Minimum step length
  const double minh = 0.25 * (maxh /= ds); // [voxel]
  #if SELFINTERSECTION_CHECK_TRIANGLES
  #if SELFINTERSECTION_CHECK_TOUCHING
  // Minimum triangle distance
  const double dtol = 0.1 /* [mm] */ / ds; // [voxel]
  #endif
  // Nearby triangles search radius
  //
  // Auxiliary data structures will be updated whenever necessary during
  // the propagation of the surface. Choosing a larger radius results in
  // a long time for initializing the data structures and a large memory
  // requirement, while a smaller radius requires often updates of the
  // auxiliary data structures. As usual, it's a trade off...
  //
  // FIXME: The memory occupied by the auxiliary data structures grows
  //        massively with an increase of the search radius. Instead of
  //        only a few 100MB, it can easily go up to >10GB! This also
  //        depends on the chosen average edgelength.
  const double radius = 25 * maxh; // [voxel]
  #endif // SELFINTERSECTION_CHECK_TRIANGLES

  irtkRemesher remesher;
  if (edgelength > .0) {
    IRTK_RESET_TIMING();
    // Ensure we are dealing with triangles only
    vtkSmartPointer<vtkTriangleFilter> triangulator;
    triangulator = vtkSmartPointer<vtkTriangleFilter>::New();
    triangulator->PassVertsOff();
    triangulator->PassLinesOff();
    SetVTKInput(triangulator, mesh);
    triangulator->Update();
    mesh = triangulator->GetOutput();
    // Iteratively remesh input surface
    verbose = 0;
    remesher.MaxEdgeLength(2.0 * edgelength / ds);
    remesher.MinEdgeLength(0.8 * edgelength / ds);
    int num_edges_processed = 0;
    for (int iter = 0; iter < 20; ++iter) {
      mesh = remesher.Remesh(mesh, &num_edges_processed);
      if (num_edges_processed == 0) break;
    }
    verbose = _verbose;
    IRTK_DEBUG_TIMING(1, "initial remeshing");
  }

  // Custom nearby triangles structure used because none of the
  // vtkAbstractCellLocator subclasses implements FindCellsWithinBounds
  // Moreover, first looking for the nearby vertices and inversely
  // build the nearby triangles set should be more efficient anyway.
  #if SELFINTERSECTION_CHECK_TRIANGLES
    vector<set<irtkVector3D<vtkIdType> > >  triangles;
    vector<set<vtkIdType> >                 nearby_vertices;
    vector<set<irtkVector3D<vtkIdType> > >  nearby_triangles;
    set<vtkIdType>::const_iterator          vertex, vertex_end;
  #endif // SELFINTERSECTION_CHECK_TRIANGLES

  // Initialize BSP tree for self-intersection check
  #if SELFINTERSECTION_CHECK_NORMAL || SELFINTERSECTION_CHECK_EDGES
    vtkSmartPointer<vtkModifiedBSPTree> bsptree;
    bsptree = vtkSmartPointer<vtkModifiedBSPTree>::New();
    bsptree->SetMaxLevel(100);
    bsptree->SetNumberOfCellsPerNode(1);
    bsptree->LazyEvaluationOff();
    bsptree->CacheCellBoundsOff();
    bsptree->UseExistingSearchStructureOff();
    bsptree->AutomaticOn();
    bsptree->SetDataSet(mesh);
  #endif // SELFINTERSECTION_CHECK_NORMAL || SELFINTERSECTION_CHECK_EDGES

  // Initialize interpolators
  ImageInterpolator f;
  f.Input(&psi);
  f.Initialize();

  // Initialize normal calculation filter
  vtkSmartPointer<vtkPolyDataNormals> calc_normals;
  calc_normals = vtkSmartPointer<vtkPolyDataNormals>::New();
  calc_normals->SplittingOff();
  calc_normals->ConsistencyOn();
  calc_normals->AutoOrientNormalsOn();

  GradientInterpolator *df = NULL;
  if (!grad.IsEmpty()) {
    df = new GradientInterpolator();
    df->Input(&grad);
    df->Initialize();
  }

  // Initialize auxiliary point data
  vtkIdType id, npts, *pts = NULL;
  vtkSmartPointer<vtkDataArray> prestep, newstep, step, dist, coll;

  prestep = vtkSmartPointer<vtkFloatArray>::New();
  newstep = vtkSmartPointer<vtkFloatArray>::New();
  step    = vtkSmartPointer<vtkFloatArray>::New();
  dist    = vtkSmartPointer<vtkFloatArray>::New();
  coll    = vtkSmartPointer<vtkIntArray  >::New();

  step->SetName("step");
  dist->SetName("dist");
  coll->SetName("coll");

  step->SetNumberOfComponents(1);
  dist->SetNumberOfComponents(1);
  coll->SetNumberOfComponents(1);

  step->SetNumberOfTuples(mesh->GetNumberOfPoints());
  dist->SetNumberOfTuples(mesh->GetNumberOfPoints());
  coll->SetNumberOfTuples(mesh->GetNumberOfPoints());

  for (vtkIdType id = 0; id < mesh->GetNumberOfPoints(); ++id) {
    step->SetTuple1(id, maxh);
    dist->SetTuple1(id,   .0);
    coll->SetTuple1(id,    0);
  }

  mesh->GetPointData()->AddArray(step);
  mesh->GetPointData()->AddArray(dist);
  mesh->GetPointData()->AddArray(coll);

  // Iteratively update position of vertices
  IRTK_RESET_TIMING();
  for (int niter = 1; niter <= maxit; ++niter) {

    IRTK_START_TIMING();
    if (verbose && debug_time > 1) cout << endl;

    // Remesh surface
    if (edgelength > .0) {
      IRTK_RESET_TIMING();
      verbose = 0;
      mesh    = remesher.Remesh(mesh);
      verbose = _verbose;
      IRTK_DEBUG_TIMING(2, "remeshing the surface");
    }

    // Smooth surface mesh
    if (smoothing > .0) {
      IRTK_RESET_TIMING();
      irtkPolyDataSmoothing smoother;
      smoother.Input(mesh);
      smoother.NumberOfIterations(1);
      smoother.Lambda(smoothing);
      smoother.Sigma(.0);
      smoother.Run();
      mesh = smoother.Output();
      IRTK_DEBUG_TIMING(2, "smoothing the surface");
      // Reset step lengths
      // TODO Figure out if there is a way to both relax the mesh and at the
      //      same time control the individual step length of each vertex.
      step = mesh->GetPointData()->GetArray("step");
      if (step) {
        for (id = 0; id < mesh->GetNumberOfPoints(); ++id) {
          step->SetTuple1(id, maxh);
        }
      }
    }

    // Recalculate normals
    if (!df && recalcnormals) {
      SetVTKInput(calc_normals, mesh);
      IRTK_RESET_TIMING();
      calc_normals->Update();
      IRTK_DEBUG_TIMING(2, "computing surface normals");
      mesh = calc_normals->GetOutput();
    }

    // Rebuild mesh links
    mesh->BuildCells();
    mesh->BuildLinks();

    // Get auxiliary point data
    step = mesh->GetPointData()->GetArray("step");
    dist = mesh->GetPointData()->GetArray("dist");
    coll = mesh->GetPointData()->GetArray("coll");

    if (!step) {
      cerr << "Error: Expected filters to preserve point data 'step'" << endl;
      exit(1);
    }
    if (!dist) {
      cerr << "Error: Expected filters to preserve point data 'dist'" << endl;
      exit(1);
    }
    if (!mesh->GetPointData()->HasArray("coll")) {
      coll = vtkSmartPointer<vtkIntArray>::New();
      coll->SetName("coll");
      coll->SetNumberOfComponents(1);
      coll->SetNumberOfTuples(mesh->GetNumberOfPoints());
      mesh->GetPointData()->AddArray(coll);
    }

    // Update auxiliary data structures
    #if SELFINTERSECTION_CHECK_TRIANGLES
      // Reset data structures
      for (size_t id = 0; id < triangles.size(); ++id) {
        triangles       [id].clear();
        nearby_vertices [id].clear();
        nearby_triangles[id].clear();
      }
      triangles       .resize(mesh->GetNumberOfPoints());
      nearby_vertices .resize(mesh->GetNumberOfPoints());
      nearby_triangles.resize(mesh->GetNumberOfPoints());
      // Find vertices in vicinity of each vertex first. Using this information,
      // find for each vertex the triangles within its vicinity which do not
      // contain this point as vertex. To prevent self intersections, each
      // triangle adjacent to a vertex is tested for intersection with any
      // other triangle in the vicinity of the vertex.
      IRTK_RESET_TIMING();
      vtkSmartPointer<vtkIdList>             ptIds;
      vtkSmartPointer<vtkKdTreePointLocator> kdtree;
      ptIds  = vtkSmartPointer<vtkIdList>::New();
      kdtree = vtkSmartPointer<vtkKdTreePointLocator>::New();
      kdtree->SetDataSet(mesh);
      kdtree->BuildLocator();
      double p[3];
      for (id = 0; id < mesh->GetNumberOfPoints(); ++id) {
        mesh->GetPoint(id, p);
        kdtree->FindPointsWithinRadius(radius, p, ptIds);
        for (vtkIdType j = 0; j < ptIds->GetNumberOfIds(); ++j) {
          nearby_vertices[id].insert(ptIds->GetId(j));
        }
      }
      IRTK_DEBUG_TIMING(2, "finding nearby vertices");
      IRTK_RESET_TIMING();
      vtkCellArray * const cells = mesh->GetPolys();
      cells->InitTraversal();
      irtkVector3D<vtkIdType> triPts;
      while (cells->GetNextCell(npts, pts)) {
        triPts._x = pts[0], triPts._y = pts[1], triPts._z = pts[2];
        for (int i = 0; i < 3; ++i) {
          id = pts[i];
          triangles[id].insert(triPts);
          vertex_end = nearby_vertices[id].end();
          for (vertex = nearby_vertices[id].begin(); vertex != vertex_end; ++vertex) {
            if (*vertex != pts[0] && *vertex != pts[1] && *vertex != pts[2]) {
              nearby_triangles[*vertex].insert(triPts);
            }
          }
        }
      }
      IRTK_DEBUG_TIMING(2, "finding nearby triangles");
    #endif // SELFINTERSECTION_CHECK_TRIANGLES

    if (debug > 0) {
      // Frequency of debug output depends on debug level (-debug <n>)
      //   1: Every 10th iteration
      //   2: Every  9th iteration
      // >10: Every iteration
      const int nskip = 10 - min(debug, 10) + 1;
      if (niter == 1 || niter == maxit || niter % nskip == 0) {
        char fname[64];
        sprintf(fname, "polydatacortex_%03d.vtk", niter);
        Write(mesh, fname);
      }
    }

    // Range of point IDs
    const vtkIdType npoints = mesh->GetNumberOfPoints();
    blocked_range<vtkIdType> ids(0, npoints);

    // Update point positions, possibly introducing self intersections
    IRTK_RESET_TIMING();
    MovePoints move;
    move._Direction = df;
    move._Psi     = &f;
    move._Psi0    = psi0;
    move._Mesh    = mesh;
    move._MinStep = minh;
    parallel_for(ids, move);
    prestep->DeepCopy(step);
    IRTK_DEBUG_TIMING(2, "moving all vertices");

    // Find any collisions such as self intersections
    IRTK_RESET_TIMING();
    FindCollisions check;
    check._Mesh     = mesh;
    check._MinStep  = minh;
    check._MinGap   = dtol;
    check._MinWidth = minw;
    #if SELFINTERSECTION_CHECK_NORMAL || SELFINTERSECTION_CHECK_EDGES
      bsptree->BuildLocator();
      check._Locator = bsptree;
    #endif
    #if SELFINTERSECTION_CHECK_TRIANGLES
      check._Triangles       = &triangles;
      check._NearbyTriangles = &nearby_triangles;
    #endif
    parallel_for(ids, check);
    IRTK_DEBUG_TIMING(2, "collision detection");

    // Smooth step lengths such that neighboring vertices move at similar speed
    IRTK_RESET_TIMING();
    SmoothSteps smooth;
    smooth._Mesh       = mesh;
    smooth._InputStep  = step;
    smooth._OutputStep = newstep;
    smooth._MinStep    = minh;

    newstep->SetNumberOfTuples(mesh->GetNumberOfPoints());
    for (int iter = 0; iter < 1; ++iter) {
      parallel_for(ids, smooth);
      step->DeepCopy(newstep);
    }
    IRTK_DEBUG_TIMING(2, "adjusting the step lengths");

    // Update position of vertices, reverting collisions
    IRTK_RESET_TIMING();
    UpdatePoints update;
    update._Direction        = df;
    update._Psi              = &psi;
    update._Mesh             = mesh;
    update._PreStep          = prestep;
    update._MaxTotalDistance = maxd;
    parallel_reduce(ids, update);
    IRTK_DEBUG_TIMING(2, "correcting vertex positions");

    // Update statistics
    const int &nactive = update._NumActive;
    const int &naccept = update._NumAccept;
    if (update._MinStep > update._MaxStep) update._MinStep = update._MaxStep = .0;
    if (update._MinDist > update._MaxDist) update._MinDist = update._MaxDist = .0;
    const double &minstep = update._MinStep;
    const double  avgdist = (nactive ? update._SumDist / nactive : .0);
    const double &mindist = update._MinDist;
    const double &maxdist = update._MaxDist;

    // Report progress
    if (verbose) {
      int nintersect = 0, ntouching = 0;
      for (id = 0; id < mesh->GetNumberOfPoints(); ++id) {
        switch (static_cast<int>(coll->GetTuple1(id))) {
          case SELF_INTERSECTION: ++nintersect; break;
          case TOUCHING:          ++ntouching;  break;
        }
      }
      if (debug_time > 1) cout << endl;
      cout << setw(3) << niter
           << "  #points: " << setw(6) << npoints
           << ", #active: " << setw(6) << nactive << " ("
           << setw(6) << fixed << setprecision(2) << (100.0 * nactive / npoints)
           << "%), #accepts: " << setw(6) << naccept << " ("
           << setw(6) << fixed << setprecision(2)
           << (nactive ? 100.0 * naccept / nactive : .0)
           << "%), #collisions: " << setw(6) << nintersect << " ("
           << setw(6) << fixed << setprecision(2)
           << (nactive ? 100.0 * nintersect / nactive : .0)
           << "%), #repulsions: " << setw(6) << ntouching  << " ("
           << setw(6) << fixed << setprecision(2)
           << (nactive ? 100.0 * ntouching / nactive : .0)
           << "%), delta: ";
      cout << setw(7) << fixed << setprecision(5) << avgdist << " mm";
      if (verbose > 1) {
        cout << " in [" << setw(7) << fixed << setprecision(5) << mindist << " "
                        << setw(7) << fixed << setprecision(5) << maxdist << "]";
      }
      cout << endl;
    }

    // Stop surface propagation if only few active vertices remain
    if (nactive <= round(minactive * npoints) || minstep < minh) break;
  }
  if (verbose && debug_time) cout << endl;
  IRTK_DEBUG_TIMING(1, "converging towards final surface");

  // Smooth final surface mesh
  vtkSmartPointer<vtkSmoothPolyDataFilter> smoother;
  smoother = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
  smoother->SetNumberOfIterations(20);
  SetVTKInput(smoother, mesh);

  // Compute normals of final surface mesh
  SetVTKConnection(calc_normals, smoother);
  calc_normals->Update();
  mesh = calc_normals->GetOutput();

  // Remove auxiliary point data
  mesh->GetPointData()->RemoveArray("step");
  mesh->GetPointData()->RemoveArray("dist");
  mesh->GetPointData()->RemoveArray("coll");

  return mesh;
}

// -----------------------------------------------------------------------------
/// Set psi value outside masked region to 10000
void MaskRegion(irtkRealImage &psi, const irtkByteImage &mask)
{
  irtkGenericNearestNeighborExtrapolateImageFunction<irtkByteImage> interp;
  interp.Input(&mask);
  interp.Initialize();

  double x, y, z;
  for (int k = 0; k < psi.Z(); ++k)
  for (int j = 0; j < psi.Y(); ++j)
  for (int i = 0; i < psi.X(); ++i) {
    x = i, y = j, z = k;
    psi.ImageToWorld(x, y, z);
    if (interp.Get(x, y, z) <= 0) psi(i, j, k) = 10000;
  }
}

// -----------------------------------------------------------------------------
/// Naive cutting plane calculation assuming symmetric brain hemispheres mask
double CuttingPlane(const irtkByteImage &mask, Hemisphere hemisphere, double n[3])
{
  irtkBytePixel h;
  irtkPoint center[2];
  int       numvox[2] = {0};
  double p[3];
  for (int k = 0; k < mask.Z(); ++k)
  for (int j = 0; j < mask.Y(); ++j)
  for (int i = 0; i < mask.X(); ++i) {
    h = mask(i, j, k) - 1;
    if (h < 0 || h > 2) continue;
    p[0] = i, p[1] = j, p[2] = k;
    mask.ImageToWorld(p[0], p[1], p[2]);
    center[h] += p;
    numvox[h] += 1;
  }
  if (numvox[0] == 0 || numvox[1] == 0) {
    cerr << "Error: Expected mask to contain label 1 for right hemisphere and 2 for left hemisphere" << endl;
    exit(1);
  }
  center[0] /= numvox[0];
  center[1] /= numvox[1];
  p[0] = (center[0]._x + center[1]._x) / 2.0;
  p[1] = (center[0]._y + center[1]._y) / 2.0;
  p[2] = (center[0]._z + center[1]._z) / 2.0;
  n[0] =  center[0]._x - center[1]._x;
  n[1] =  center[0]._y - center[1]._y;
  n[2] =  center[0]._z - center[1]._z;
  double m = sqrt(vtkMath::Dot(n, n));
  n[0] /= m, n[1] /= m, n[2] /= m;
  if (hemisphere == LH) n[0] = -n[0], n[1] = -n[1], n[2] = -n[2];
  return -vtkMath::Dot(p, n);
}

// -----------------------------------------------------------------------------
/// Set (subcortical) psi value outside specified brain hemisphere to 10000
void MaskHemisphere(irtkRealImage &psi, const irtkByteImage &mask, Hemisphere hemisphere)
{
  // Initialize mask interpolator
  irtkGenericNearestNeighborExtrapolateImageFunction<irtkByteImage> interp;
  interp.Input(&mask);
  interp.Initialize();

  // Get cutting plane of specified brain hemisphere
  double n[3];
  double b = CuttingPlane(mask, hemisphere, n);

  if (verbose) {
    cout << "\nMasking " << (hemisphere == RH ? "left" : "right") << " hemisphere:"
            " cutting plane b=" << b << ", n=[" << n[0] << " " << n[1] << " " << n[2] << "]\n" << endl;
  }

  // Mask values outside specified brain hemisphere
  double wx, wy, wz, vx, vy, vz;
  irtkBytePixel maskval;
  for (int k = 0; k < psi.Z(); ++k)
  for (int j = 0; j < psi.Y(); ++j)
  for (int i = 0; i < psi.X(); ++i) {
    wx = i, wy = j, wz = k;
    psi.ImageToWorld(wx, wy, wz);
    vx = wx, vy = wy, vz = wz;
    mask.WorldToImage(vx, vy, vz);
    maskval = interp.Get(vx, vy, vz);
    if (maskval == 0) {
      if ((wx * n[0] + wy * n[1] + wz * n[2] + b) < .0) {
        psi(i, j, k) = 10000;
      }
    } else if (maskval != hemisphere) {
      psi(i, j, k) = 10000;
    }
  }
}

// -----------------------------------------------------------------------------
/// Read map from input segmentation to tissue segmentation label
void ReadLabelToTissueMap(const char *csv_name, map<irtkBytePixel, irtkBytePixel> &tissue)
{
  ifstream in(csv_name);
  string line, token, first, last;
  int l = 0, label, from, to;
  while (getline(in, line)) {
    istringstream ss(line);
    while(getline(ss, token, ',')) {
      size_t pos = token.find("-");
      if (pos != string::npos) {
        first = token.substr(0, pos);
        last  = token.substr(pos+1);
        if (FromString(first.c_str(), from) && FromString(last.c_str(), to)) {
          for (label = from; label <= to; ++label) tissue[label] = l;
        } else {
          cerr << "Failed to convert token " << token << " to range of integers" << endl;
          exit(1);
        }
      }
      if (FromString(token.c_str(), label)) {
        tissue[label] = l;
      } else {
        cerr << "Failed to convert token " << token << " to integer" << endl;
        exit(1);
      }
    }
    ++l;
  }
}

// -----------------------------------------------------------------------------
/// Map input segmentation labels to their respective tissue segmentation label
void MakeTissueSegmentation(irtkByteImage &seg, const map<irtkBytePixel, irtkBytePixel> &tissue)
{
  map<irtkBytePixel, irtkBytePixel>::const_iterator entry;
  const int nvox = seg.NumberOfVoxels();
  irtkBytePixel *l = seg.Data();
  for (int i = 0; i < nvox; ++i, ++l) {
    entry = tissue.find(*l);
    *l = (entry == tissue.end()) ? BG : entry->second;
  }
}

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  REQUIRES_POSARGS(1);

  // Optional arguments
  bool        isseg          = true;  // Input is segmentation image
  const char *mesh_name      = NULL;  // Name of output mesh file
  const char *init_name      = NULL;  // Name of initial surface file
  const char *image_name     = NULL;  // Greyscale image (brain MR image)
  const char *labels_name    = NULL;  // CSV file mapping input label to tissue
  const char *mask_name      = NULL;  // Region/hemispheres mask
  const char *psi_name       = NULL;  // Name of output psi image file
  Hemisphere  hemisphere     = UH;    // Which hemishere/region?
  double      psi0           = 9000;  // Target  psi value
  double      psi1           = 2000;  // Initial psi value
  double      maxh           = 0.05;  // Initial step size in mm
  double      minw           = -1.0;  // Minimum width of mesh
  double      maxd           =  .0;   // Maximum displacement in mm
  int         maxn           = 500;   // Maximum no. of iterations
  double      minactive      =  10;   // Minimum number of active points in %
  double      smoothing      =  .0;   // Laplacian smoothing lambda value
  bool        usenormals     = true;  // Whether to propate in normal direction
  bool        recalcnormals  = true;  // Recompute normals after each step
  bool        ascii          = false; // Write vtkPolyData in ASCII format
  double      edgelength;             // Average edge length in mm
  edgelength = numeric_limits<double>::quiet_NaN();

  for (ALL_OPTIONS) {
    if      (OPTION("-mesh"))  mesh_name  = ARGUMENT;
    else if (OPTION("-init"))  init_name  = ARGUMENT;
    else if (OPTION("-mask"))  mask_name  = ARGUMENT;
    else if (OPTION("-image")) image_name = ARGUMENT;
    else if (OPTION("-labels")) labels_name = ARGUMENT;
    else if (OPTION("-rh"  ))  hemisphere = RH;
    else if (OPTION("-lh"  ))  hemisphere = LH;
    else if (OPTION("-isseg")) isseg = true;
    else if (OPTION("-ispsi")) isseg = false;
    else if (OPTION("-psi"))   psi_name   = ARGUMENT;
    else if (OPTION("-psi0"))  psi0       = atof(ARGUMENT);
    else if (OPTION("-psi1"))  psi1       = atof(ARGUMENT);
    else if (OPTION("-maxh"))  maxh       = atof(ARGUMENT);
    else if (OPTION("-minw"))  minw       = atof(ARGUMENT);
    else if (OPTION("-maxd"))  maxd       = atof(ARGUMENT);
    else if (OPTION("-maxn"))  maxn       = atoi(ARGUMENT);
    else if (OPTION("-minactive"))  minactive  = atof(ARGUMENT);
    else if (OPTION("-edgelength")) edgelength = atof(ARGUMENT);
    else if (OPTION("-lambda") || OPTION("-smoothing")) smoothing = atof(ARGUMENT);
    else if (OPTION("-normal"  ))        usenormals    = true;
    else if (OPTION("-gradient"))        usenormals    = false;
    else if (OPTION("-recalcnormals"  )) recalcnormals = true;
    else if (OPTION("-norecalcnormals")) recalcnormals = false;
    else if (OPTION("-ascii" ) || OPTION("-nobinary")) ascii = true;
    else if (OPTION("-binary") || OPTION("-noascii" )) ascii = false;
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (maxd <= .0 && !init_name) maxd = 5.0;
  if      (minactive > 100.0) minactive = 100.0;
  else if (minactive <    .0) minactive =    .0;

  // Read tissue segmentation labels and compute psi image or read
  // precomputed input psi image directly
  irtkRealImage psi;
  if (isseg) {
    if (verbose) {
      cout << "Reading segmentation labels from file " << POSARG(1) << endl;
    }
    irtkByteImage seg(POSARG(1));
    if (labels_name) {
      if (verbose) cout << "Reading tissue labels from file " << labels_name << endl;
      map<irtkBytePixel, irtkBytePixel> label2tissue;
      ReadLabelToTissueMap(labels_name, label2tissue);
      MakeTissueSegmentation(seg, label2tissue);
      if (debug) seg.Write("polydatacortex_labels.nii.gz");
    }
    MakeIsotropic(seg);
    psi = ComputeLaplacian(seg);
  } else {
    if (verbose) {
      cout << "Reading precomputed psi image from file " << POSARG(1) << endl;
    }
    psi.Read(POSARG(1));
  }

  // Set minimum distance between folds
  if (minw < .0) minw = psi.GetXSize();

  // Write psi image if requested
  if (psi_name) {
    psi.Write(psi_name);
    if (!mesh_name) return 0;
  }

  // Extract specified region / brain hemisphere
  if (mask_name) {
    irtkByteImage mask(mask_name);
    if (hemisphere == UH) MaskRegion    (psi, mask);
    else                  MaskHemisphere(psi, mask, hemisphere);
  }

  // Initialize surface
  vtkSmartPointer<vtkPolyData> mesh;

  if (init_name) {
    if (strcmp(init_name, "convex")      == 0 ||
        strcmp(init_name, "hull")        == 0 ||
        strcmp(init_name, "convex hull") == 0 ||
        strcmp(init_name, "chull")       == 0) {
      if (IsNaN(edgelength)) edgelength = .5; // [mm]
      mesh = InitialSurface(psi, psi1, true);
    } else {
      vtkSmartPointer<vtkPolyDataReader> reader;
      reader = vtkSmartPointer<vtkPolyDataReader>::New();
      reader->SetFileName(init_name);
      reader->Update();
      mesh = reader->GetOutput();
    }
  } else {
    mesh = InitialSurface(psi, psi1);
  }
  if (debug) Write(mesh, "polydatacortex_init.vtk");

  // By default, do not remesh (if no smoothing is used) as it can alter topology
  // FIXME Fix the irtkRemesher filter that it preserves topology! In particular
  //       the collapsing of edges causes undesired changes in topology.
  if (IsNaN(edgelength)) edgelength = .0;

  // Surface mesh refinement
  if (image_name) {

    // Compute Laplacian of intensity image
    irtkRealImage image(image_name);
    Resample(image, psi.Attributes());
    image = ComputeLaplacian(image, psi.GetXSize());
    if (debug) image.Write("polydatacortex_laplacian_of_image.nii.gz");

    // Refine surface mesh based on Laplacian of intensity image
    irtkRealImage grad; // no gradient used
    mesh = RefineSurface(mesh, image, grad, .0, maxh, minw, maxd, maxn,
                         minactive / 100.0, edgelength, smoothing, recalcnormals);

  } else {

    // Compute gradient of psi image
    irtkRealImage grad;
    if (!usenormals) {
      grad = PropagationDirection(mesh, psi, psi0);
      if (debug) grad.Write("polydatacortex_gradient.nii.gz");
    }

    // Refine surface mesh based on psi image computed from segmentation
    mesh = RefineSurface(mesh, psi, grad, psi0, maxh, minw, maxd, maxn,
                         minactive / 100.0, edgelength, smoothing, recalcnormals);

  }

  // Write surface mesh
  if (mesh_name) Write(mesh, mesh_name, ascii);

  return 0;
}
