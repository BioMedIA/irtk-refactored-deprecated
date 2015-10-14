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
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "usage: " << name << " <in> <out> [options]" << endl;
  cout << endl;
  cout << "Calculate surface curvature." << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << "  -H [<name>]               Mean curvature." << endl;
  cout << "  -K [<name>]               Gauss curvature." << endl;
  cout << "  -C [<name>]               Curvedness." << endl;
  cout << "  -k1 [<name>]              Minimum curvature." << endl;
  cout << "  -k2 [<name>]              Maximum curvature." << endl;
  cout << "  -k1k2 [<name>] [<name>]   Principal curvatures." << endl;
  cout << "  -e1 [<name>]              Direction of minimum curvature." << endl;
  cout << "  -e2 [<name>]              Direction of maximum curvature." << endl;
  cout << "  -normalize                Normalize curvature using volume of convex hull." << endl;
  cout << "  -vtk                      Use vtkCurvatures when possible." << endl;
  cout << "  -robust                   Do not use vtkCurvatures. Instead, estimate the curvature" << endl;
  cout << "                            tensor field and decompose it to obtain principle curvatures. (default)" << endl;
  cout << endl;
  cout << "  -smooth [<niter>] [<sigma>] [<sigma2>]" << endl;
  cout << "      Smooth scalar curvature measures using a Gaussian smoothing kernel." << endl;
  cout << "      If sigma2 is specified, an anisotropic kernel with standard deviation" << endl;
  cout << "      sigma along the direction of minimum curvature, and sigma2 in the" << endl;
  cout << "      direction of maximum curvature is used. If the value of sigma2 is \"tensor\"" << endl;
  cout << "      instead of a numeric value, the isotropic Gaussian kernel is oriented" << endl;
  cout << "      and scaled along each local geometry axis using the curvature tensor." << endl;
  PrintCommonOptions(cout);
  cout << endl;
}

// =============================================================================
// Includes
// =============================================================================
#ifdef HAS_VTK

#include <irtkPolyDataUtils.h>
#include <irtkPolyDataCurvature.h>
#include <irtkPolyDataSmoothing.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>

using namespace irtk::polydata;

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  EXPECTS_POSARGS(2);

  const char *input_name  = POSARG(1);
  const char *output_name = POSARG(2);

  const char *kmin_name       = irtkPolyDataCurvature::MINIMUM;
  const char *kmax_name       = irtkPolyDataCurvature::MAXIMUM;
  const char *gauss_name      = irtkPolyDataCurvature::GAUSS;
  const char *mean_name       = irtkPolyDataCurvature::MEAN;
  const char *curvedness_name = irtkPolyDataCurvature::CURVEDNESS;
  const char *e1_name         = irtkPolyDataCurvature::MINIMUM_DIRECTION;
  const char *e2_name         = irtkPolyDataCurvature::MAXIMUM_DIRECTION;
  const char *tensor_name     = irtkPolyDataCurvature::TENSOR;
  const char *inverse_name    = irtkPolyDataCurvature::INVERSE_TENSOR;

  bool   use_vtkCurvatures    = false;
  int    type                 = 0;
  int    tensor_averaging     = 3;
  bool   normalize            = false;
  int    smooth_iterations    = 0;
  double smooth_sigma         = .0;
  double smooth_sigma2        = .0;
  bool   smooth_along_tensor  = false;
  bool   smooth_anisotropic   = false;

  for (ALL_OPTIONS) {
    if (OPTION("-k1")) {
      type |= irtkPolyDataCurvature::Minimum;
      if (HAS_ARGUMENT) kmin_name = ARGUMENT;
    }
    else if (OPTION("-k2")) {
      type |= irtkPolyDataCurvature::Maximum;
      if (HAS_ARGUMENT) kmax_name = ARGUMENT;
    }
    else if (OPTION("-k1k2")) {
      type |= (irtkPolyDataCurvature::Minimum | irtkPolyDataCurvature::Maximum);
      if (HAS_ARGUMENT) {
        kmin_name = ARGUMENT;
        kmax_name = ARGUMENT;
      }
    }
    else if (OPTION("-H")) {
      type |= irtkPolyDataCurvature::Mean;
      if (HAS_ARGUMENT) mean_name = ARGUMENT;
    }
    else if (OPTION("-K")) {
      type |= irtkPolyDataCurvature::Gauss;
      if (HAS_ARGUMENT) gauss_name = ARGUMENT;
    }
    else if (OPTION("-C")) {
      type |= irtkPolyDataCurvature::Curvedness;
      if (HAS_ARGUMENT) curvedness_name = ARGUMENT;
    }
    else if (OPTION("-n") || OPTION("-normal")) {
      type |= irtkPolyDataCurvature::Normal;
    }
    else if (OPTION("-e1")) {
      type |= irtkPolyDataCurvature::MinimumDirection;
      if (HAS_ARGUMENT) e1_name = ARGUMENT;
    }
    else if (OPTION("-e2")) {
      type |= irtkPolyDataCurvature::MaximumDirection;
      if (HAS_ARGUMENT) e2_name = ARGUMENT;
    }
    else if (OPTION("-tensor")) {
      type |= irtkPolyDataCurvature::Tensor;
      if (HAS_ARGUMENT) tensor_name = ARGUMENT;
    }
    else if (OPTION("-inverse-tensor")) {
      type |= irtkPolyDataCurvature::InverseTensor;
      if (HAS_ARGUMENT) inverse_name = ARGUMENT;
    }
    else if (OPTION("-tensor-averaging")) {
      tensor_averaging = atoi(ARGUMENT);
    }
    else if (OPTION("-normalize")) normalize = true;
    else if (OPTION("-vtk"))    use_vtkCurvatures = true;
    else if (OPTION("-robust")) use_vtkCurvatures = false;
    else if (OPTION("-smooth")) {
      smooth_iterations = 1;
      if (HAS_ARGUMENT) smooth_iterations = atoi(ARGUMENT);
      if (HAS_ARGUMENT) smooth_sigma      = atof(ARGUMENT);
      if (HAS_ARGUMENT) {
        const char *arg = ARGUMENT;
        if (strcmp(arg, "tensor") == 0) {
          smooth_along_tensor = true;
        } else {
          smooth_along_tensor = false;
          smooth_sigma2       = atof(arg);
        }
        smooth_anisotropic = true;
      }
    }
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (type == 0) type = irtkPolyDataCurvature::Scalars;

  int curvature_type = type;
  if (smooth_along_tensor) {
    curvature_type |= irtkPolyDataCurvature::Tensor;
  } else if (smooth_sigma2 != .0) {
    curvature_type |= irtkPolyDataCurvature::MinimumDirection;
    curvature_type |= irtkPolyDataCurvature::MaximumDirection;
  }

  // Read input surface and triangulate it
  vtkSmartPointer<vtkPolyData> surface = Triangulate(ReadPolyData(input_name));
  surface->BuildLinks();

  // Compute curvature
  irtkPolyDataCurvature curvature;
  curvature.Input(surface);
  curvature.CurvatureType(curvature_type);
  curvature.VtkCurvatures(use_vtkCurvatures);
  curvature.TensorAveraging(tensor_averaging);
  curvature.Normalize(normalize);
  curvature.Run();
  surface = curvature.Output();

  // Get output arrays
  vtkPointData *pd         = surface->GetPointData();
  vtkDataArray *kmin       = pd->GetArray(irtkPolyDataCurvature::MINIMUM);
  vtkDataArray *kmax       = pd->GetArray(irtkPolyDataCurvature::MAXIMUM);
  vtkDataArray *mean       = pd->GetArray(irtkPolyDataCurvature::MEAN);
  vtkDataArray *gauss      = pd->GetArray(irtkPolyDataCurvature::GAUSS);
  vtkDataArray *curvedness = pd->GetArray(irtkPolyDataCurvature::CURVEDNESS);
  vtkDataArray *e1         = pd->GetArray(irtkPolyDataCurvature::MINIMUM_DIRECTION);
  vtkDataArray *e2         = pd->GetArray(irtkPolyDataCurvature::MAXIMUM_DIRECTION);
  vtkDataArray *tensor     = pd->GetArray(irtkPolyDataCurvature::TENSOR);
  vtkDataArray *inverse    = pd->GetArray(irtkPolyDataCurvature::INVERSE_TENSOR);

  // Rename output arrays
  if (kmin)       kmin      ->SetName(kmin_name);
  if (kmax)       kmax      ->SetName(kmax_name);
  if (mean)       mean      ->SetName(mean_name);
  if (gauss)      gauss     ->SetName(gauss_name);
  if (curvedness) curvedness->SetName(curvedness_name);
  if (e1)         e1        ->SetName(e1_name);
  if (e2)         e2        ->SetName(e2_name);
  if (tensor)     tensor    ->SetName(tensor_name);
  if (inverse)    inverse   ->SetName(inverse_name);

  // Smooth curvature measures
  if (smooth_iterations) {
    irtkPolyDataSmoothing smoother;
    smoother.Input(surface);
    smoother.SmoothPointsOff();
    if (kmin)       smoother.SmoothArray(kmin_name);
    if (kmax)       smoother.SmoothArray(kmax_name);
    if (mean)       smoother.SmoothArray(mean_name);
    if (gauss)      smoother.SmoothArray(gauss_name);
    if (curvedness) smoother.SmoothArray(curvedness_name);
    if (e1)         smoother.SmoothArray(e1_name);
    if (e2)         smoother.SmoothArray(e2_name);
    smoother.NumberOfIterations(smooth_iterations);
    smoother.Sigma(-smooth_sigma); // negative: multiple of avg. edge length
    if (smooth_anisotropic) {
      smoother.Weighting(irtkPolyDataSmoothing::AnisotropicGaussian);
      if (smooth_along_tensor) {
        smoother.GeometryTensorName(tensor_name);
      } else {
        smoother.MinimumDirectionName(e2_name);
        smoother.MaximumDirectionName(e1_name);
      }
      smoother.MaximumDirectionSigma(-smooth_sigma2);
    } else {
      smoother.Weighting(irtkPolyDataSmoothing::Gaussian);
    }
    smoother.Run();
    vtkPointData *pd = smoother.Output()->GetPointData();
    if (kmin) kmin->DeepCopy(pd->GetArray(kmin_name));
    if (kmax) kmax->DeepCopy(pd->GetArray(kmax_name));
    if (mean) mean->DeepCopy(pd->GetArray(mean_name));
    if (gauss) gauss->DeepCopy(pd->GetArray(gauss_name));
    if (curvedness) curvedness->DeepCopy(pd->GetArray(curvedness_name));
    if (e1) e1->DeepCopy(pd->GetArray(e1_name));
    if (e2) e2->DeepCopy(pd->GetArray(e2_name));
  }

  // Remove not requested output arrays which were used for anisotropic smoothing
  if ((type & irtkPolyDataCurvature::Tensor) == 0) {
    surface->GetPointData()->RemoveArray(tensor_name);
  }
  if ((type & irtkPolyDataCurvature::MinimumDirection) == 0) {
    surface->GetPointData()->RemoveArray(e1_name);
  }
  if ((type & irtkPolyDataCurvature::MaximumDirection) == 0) {
    surface->GetPointData()->RemoveArray(e2_name);
  }

  // Write output surface
  WritePolyData(output_name, surface);
  return 0;
}

#else // HAS_VTK

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  EXPECTS_POSARGS(2); // handle help options
  cerr << EXECNAME << ": Needs to be compiled with VTK" << endl;
  exit(1);
}

#endif // HAS_VTK
