// =============================================================================
// Project: Image Registration Toolkit (IRTK)
// Package: PolyData
//
// Copyright (c) 2015 Imperial College London
// Copyright (c) 2015 Andreas Schuh
// Copyright (c) Paul Aljabar
// =============================================================================

#ifdef HAS_VTK

#include <irtkCommon.h>
#include <irtkPolyDataSmoothing.h>
#include <irtkPolyDataCurvature.h>
#include <irtkPolyDataUtils.h>

#include <vtkDataSetAttributes.h>

using namespace irtk::polydata;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "usage: " << name << " <input> <output> [iterations] [relaxation] [options]" << endl;
  cout << endl;
  cout << "Smooths the node positions and/or scalar data of a surface mesh." << endl;
  cout << endl;
  cout << "When node positions are smoothed, iteratively move mesh nodes towards the" << endl;
  cout << "centroid of the adjacent nodes. The relaxation factor (RF) determines how" << endl;
  cout << "far each node moves towards the local centroid (0<= RF <= 1). The new position" << endl;
  cout << "of a node is a weighted combination of its previous position and the neighbours' centroid." << endl;
  cout << "RF = 1 moves a node all the way to the centroid while RF = 0 keeps a node at its" << endl;
  cout << "previous position." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  input                File name of input surface mesh." << endl;
  cout << "  output               File name of output surface mesh." << endl;
  cout << "  iterations           No. of smoothing iterations. (default: 1)" << endl;
  cout << "  relaxation           Relaxation factor (0 <= RF <= 1). (default: 1)" << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -points              Smooth node positions. (default if no -scalars specified)" << endl;
  cout << "  -scalars [<name>]    Name of scalar point data array to smooth. If no <name> given," << endl;
  cout << "                       all non-attribute scalar arrays and the SCALARS if any are smoothed." << endl;
  cout << "  -exclnode            Exclude node itself from smoothing, only average adjacent values." << endl;
  cout << "  -areaweighted        Area weighted Laplacian relaxation of surface node positions. (default)" << endl;
  cout << "  -combinatorial       Combinatorial weighting of adjacent nodes." << endl;
  cout << "  -inverse-distance [<sigma>]" << endl;
  cout << "                       Inverse distance weighted smoothing kernel." << endl;
  cout << "  -gaussian [<sigma>]  Gaussian smoothing kernel." << endl;
  cout << "  -anisotropic [<sigma>] [<tensor_array>]" << endl;
  cout << "                       Anisotropic Gaussian smoothing kernel given input tensor field point data array." << endl;
  cout << "                       The default tensor field used is the curvature tensor output of polydatacurvature." << endl;
  cout << "  -anisotropic <sigma> <sigma2> [<e2_array>]" << endl;
  cout << "                       Anisotropic Gaussian smoothing kernel given direction of maximum curvature." << endl;
  cout << "  -anisotropic <sigma> <sigma2> <e1_array> <e2_array>" << endl;
  cout << "                       Anisotropic Gaussian smoothing kernel given directions of minimum and maximum curvature." << endl;
  cout << "  -track [<name>]      Track the signed distances traversed by points (positive is outwards)." << endl;
  cout << "                       The optional arguments specifies the name of output output point data array." << endl;
  cout << "                       (default: \"smoothingDists\" if option is given without <name>)." << endl;
  cout << "  -threshold <value>   A value indicating the smoothness as the L_2 norm of H^2 over the surface." << endl;
  cout << "                       See Tosun, MedIA, 2004. Iterations stop if the norm for the surface drops" << endl;
  cout << "                       below the given value. Only used for -areaweighted smoothing." << endl;
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Tosun, MedIA, 2004
// =============================================================================

#include <vtkSmartPointer.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkTriangle.h>
#include <vtkIdList.h>
#include <vtkPolyDataNormals.h>
#include <vtkCurvatures.h>

// -----------------------------------------------------------------------------
double hSquareRobustMean(vtkPolyData* surface)
{
  const vtkIdType noOfPoints = surface->GetNumberOfPoints();

  const double lo =  5.0;
  const double hi = 95.0;

  // Compute mean curvature
  vtkCurvatures *curve = vtkCurvatures::New();
  SetVTKInput(curve, surface);
  curve->SetCurvatureTypeToMean();
  curve->Update();

  vtkPolyData  *output  = curve->GetOutput();
  vtkDataArray *scalars = output->GetPointData()->GetScalars("Mean_Curvature");

  // Sort mean curvature values
  vector<double> data(noOfPoints);
  for (vtkIdType i = 0; i < noOfPoints; ++i) data[i] = scalars->GetTuple1(i);
  sort(data.begin(), data.end());

  // Get 5th and 95th percentile
  double minVal = data[int(round(lo * (noOfPoints - 1) / 100.0))];
  double maxVal = data[int(round(hi * (noOfPoints - 1) / 100.0))];

  // Compute robust squared mean within percentile range
  double sum = .0, val;
  int    num = 0;

  for (vtkIdType i = 0; i < noOfPoints; ++i) {
    val = scalars->GetTuple1(i);
    if (minVal <= val && val <= maxVal) {
      sum += val * val;
      ++num;
    }
  }

  return sum / num;
}

// -----------------------------------------------------------------------------
double surfaceArea(vtkPolyData* surface)
{
  vtkCellArray* facets = surface->GetPolys();
  vtkTriangle* facet = vtkTriangle::New();

  double A = 0.0;
  double v0[3], v1[3], v2[3];

  vtkIdType f, *vert=0;
  facets->InitTraversal();
  while (facets->GetNextCell(f,vert)){

    surface->GetPoint(vert[0],v0);
    surface->GetPoint(vert[1],v1);
    surface->GetPoint(vert[2],v2);

    A += double(facet->TriangleArea(v0,v1,v2));
  }

  return A;
}

// -----------------------------------------------------------------------------
void getCoG(vtkPolyData *input, double*cog)
{
  int i, noOfPoints;
  double cofgx, cofgy, cofgz;
  double point[3];

  cofgx = cofgy = cofgz = 0.0;

  noOfPoints = input->GetNumberOfPoints();

  for (i = 0; i < noOfPoints; i++){
    input->GetPoint (i, point);
    cofgx += point[0];
    cofgy += point[1];
    cofgz += point[2];
  }

  if (noOfPoints > 0){
    cofgx /= noOfPoints;
    cofgy /= noOfPoints;
    cofgz /= noOfPoints;
  }

  cog[0] = cofgx;
  cog[1] = cofgy;
  cog[2] = cofgz;
}

// -----------------------------------------------------------------------------
void shiftAndScalePolyData(vtkPolyData* input, double *shift, double factor)
{
  int i, j, noOfPoints;
  noOfPoints = input->GetNumberOfPoints();
  double vOld[3], vNew[3];

  for (i = 0; i < noOfPoints; ++i){
    input->GetPoint(i, vOld);

    for (j = 0; j < 3; ++j){
      vNew[j] = factor * (vOld[j] + shift[j]);
    }

    input->GetPoints()->SetPoint(i, vNew);
  }
}

// -----------------------------------------------------------------------------
double meanRadius(vtkPolyData* input, double*cog)
{
  int i, noOfPoints;
  double point[3];
  double r, rSum = 0.0;

  noOfPoints = input->GetNumberOfPoints();

  for (i = 0; i < noOfPoints; i++){
    input->GetPoint (i, point);
    r = sqrt((point[0]-cog[0])*(point[0]-cog[0]) +
             (point[1]-cog[1])*(point[1]-cog[1]) +
             (point[2]-cog[2])*(point[2]-cog[2]));
    rSum += r;
  }

  return rSum / noOfPoints;
}

// -----------------------------------------------------------------------------
void AreaWeightedLaplacianSmoothing(vtkPolyData *input,
                                    int    noOfIterations,
                                    double relaxationFactor,
                                    double smoothnessThreshold = -1.0,
                                    bool   trackingOn = false)
{
  int i, j, k;
  double E_H2, area;

  double currPos[3];
  unsigned short noOfCells = 0;
  vtkIdType* cells = NULL;
  vtkTriangle* triangle = NULL;
  double totalArea = 0;
  double update[3];
  vtkIdList* ptIds = NULL;
  double v1[3], v2[3], v3[3], centre[3];
  double triangleArea = 0;
  double dx, dy, dz;
  double dist, val;
  double *normal;
  double h2norm;

  double cogOld[3];
  double cogNew[3];
  double radiusOld, radiusNew;
  double shift[3];
  double scaleFactor;

  int noOfPoints = input->GetNumberOfPoints();
  vtkDataArray *normals = input->GetPointData()->GetNormals();

  // Output points
  double *pts = new double[3 * noOfPoints];

  // Sum of signed distances traveled by a point
  vtkSmartPointer<vtkFloatArray> dists;
  if (trackingOn) {
    dists = vtkSmartPointer<vtkFloatArray>::New();
    dists->SetName("smoothingDists");
    dists->SetNumberOfComponents(1);
    dists->SetNumberOfTuples(noOfPoints);
    dists->FillComponent(0, .0);
    input->GetPointData()->AddArray(dists);
  }

  // Smoothing iterations
  for (i = 0; i <= noOfIterations; ++i) {
    if (verbose) cout << "iteration  " << i << " ";

    // Estimate \int H^2 dA by multiplying E(H^2) with Area.
    E_H2 = hSquareRobustMean(input);
    area = surfaceArea(input);

    // The L_2 norm using the Tosun formulation (MedIA 2004)
    h2norm = sqrt(E_H2 * area / 4.0 / M_PI);
    if (h2norm < smoothnessThreshold) break;
    if (verbose > 1) cout << h2norm << endl;

    getCoG(input, cogOld);
    radiusOld = meanRadius(input, cogOld);

    // Loop over surface.
    for (j = 0; j < noOfPoints; ++j){

      // Initialisation for current point.
      totalArea = 0;
      cells = NULL;

      update[0] = 0;
      update[1] = 0;
      update[2] = 0;

      // Store the current position of the node.
      input->GetPoint(j, currPos);

      // What cells does this node adjoin?
      input->GetPointCells(j, noOfCells, cells);
      if (cells == NULL) continue;

      for (k = 0; k < noOfCells; ++k){
        triangle = vtkTriangle::SafeDownCast(input->GetCell(cells[k]));
        if (triangle == NULL) continue;
        ptIds = triangle->GetPointIds();

        input->GetPoint(ptIds->GetId(0), v1);
        input->GetPoint(ptIds->GetId(1), v2);
        input->GetPoint(ptIds->GetId(2), v3);

        triangleArea = vtkTriangle::TriangleArea(v1, v2, v3);
        vtkTriangle::TriangleCenter(v1, v2, v3, centre);

        totalArea += triangleArea;

        update[0] += triangleArea * centre[0];
        update[1] += triangleArea * centre[1];
        update[2] += triangleArea * centre[2];
      }

      if (totalArea <= 0.0) {
        update[0] = currPos[0];
        update[1] = currPos[1];
        update[2] = currPos[2];
      } else {
      	update[0] /= totalArea;
      	update[1] /= totalArea;
      	update[2] /= totalArea;
      }

      dx = relaxationFactor * (update[0] - currPos[0]);
      dy = relaxationFactor * (update[1] - currPos[1]);
      dz = relaxationFactor * (update[2] - currPos[2]);

      pts[j*3  ] = currPos[0] + dx;
      pts[j*3+1] = currPos[1] + dy;
      pts[j*3+2] = currPos[2] + dz;

      if (dists) {
        dist = sqrt(dx*dx + dy*dy + dz*dz);
        normal = normals->GetTuple3(j);
        val = normal[0]*dx + normal[1]*dy + normal[2]*dz;
        if (val < 0) dist = -dist;
        dists->SetTuple1(j, dists->GetTuple1(j) + dist);
      }
    }

    for (j = 0; j < noOfPoints; ++j){
      input->GetPoints()->SetPoint(j, pts + j*3);
    }

    // update radius and centre of gravity
    getCoG(input, cogNew);
    radiusNew = meanRadius(input, cogNew);

    shift[0] = cogOld[0] - cogNew[0];
    shift[1] = cogOld[1] - cogNew[1];
    shift[2] = cogOld[2] - cogNew[2];

    scaleFactor = radiusOld / radiusNew;

    shiftAndScalePolyData(input, shift, scaleFactor);
    if (verbose) cout << endl;
  }

  if (verbose) {
    cout << "Final iterations : " << i << endl;
    cout << "Final L_2 norm of H^2 (threshold) : " << h2norm << " (" << smoothnessThreshold << ")" << endl;
  }

  delete[] pts;
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char **argv)
{
  REQUIRES_POSARGS(2);

  const char *input_name  = POSARG(1);
  const char *output_name = POSARG(2);

  int noOfIterations = 1;
  if (NUM_POSARGS > 2) {
    if (!FromString(POSARG(3), noOfIterations)) {
      cerr << "Invalid no. of iterations argument: " << POSARG(3) << endl;
      exit(1);
    }
  }

  double relaxationFactor = 1.0;
  if (NUM_POSARGS > 3) {
    if (!FromString(POSARG(4), relaxationFactor)) {
      cerr << "Invalid relaxation factor argument: " << POSARG(4) << endl;
      exit(1);
    }
  }

  if (NUM_POSARGS > 4) {
    PrintHelp(EXECNAME);
    exit(1);
  }

  // Read the input mesh
  vtkSmartPointer<vtkPolyData> polydata = ReadPolyData(input_name);

  // Parse remaining arguments
  irtkPolyDataSmoothing::WeightFunction weighting = irtkPolyDataSmoothing::Gaussian;
  irtkPolyDataSmoothing::ArrayNames     scalar_names;
  const char *tensor_name    = NULL;
  const char *e1_name        = NULL;
  const char *e2_name        = NULL;
  bool   smooth_points       = false;
  bool   area_weighted       = false;
  double smoothnessThreshold = -1.0;
  double sigma               = .0;
  double sigma2              = .0;
  bool   trackingOn          = false;
  bool   adjacent_only       = false;

  for (ALL_OPTIONS) {
    if      (OPTION("-points")) smooth_points = true;
    else if (OPTION("-scalars")) {
      if (HAS_ARGUMENT) {
        do {
          scalar_names.push_back(ARGUMENT);
        } while (HAS_ARGUMENT);
      } else {
        for (int i = 0; i < polydata->GetPointData()->GetNumberOfArrays(); ++i) {
          int type = polydata->GetPointData()->IsArrayAnAttribute(i);
          if (polydata->GetPointData()->GetArrayName(i) &&
              type != vtkDataSetAttributes::NORMALS &&
              type != vtkDataSetAttributes::TCOORDS &&
              type != vtkDataSetAttributes::PEDIGREEIDS &&
              type != vtkDataSetAttributes::GLOBALIDS &&
              type != vtkDataSetAttributes::EDGEFLAG) {
            scalar_names.push_back(polydata->GetPointData()->GetArrayName(i));
          }
        }
      }
    }
    else if (OPTION("-iterations")) noOfIterations = atoi(ARGUMENT);
    else if (OPTION("-relaxation")) relaxationFactor = atof(ARGUMENT);
    else if (OPTION("-combinatorial")) weighting = irtkPolyDataSmoothing::Combinatorial;
    else if (OPTION("-distance") || OPTION("inversedistance")) {
      weighting = irtkPolyDataSmoothing::InverseDistance;
      if (HAS_ARGUMENT) {
        const char *arg = ARGUMENT;
        if (!FromString(arg, sigma)) {
          cerr << "Invalid -inversedistance offset argument: " << arg << endl;
          exit(1);
        }
      }
    }
    else if (OPTION("-gaussian")) {
      weighting = irtkPolyDataSmoothing::Gaussian;
      if (HAS_ARGUMENT) {
        const char *arg = ARGUMENT;
        if (!FromString(arg, sigma)) {
          cerr << "Invalid -gaussian standard deviation argument: " << arg << endl;
          exit(1);
        }
      }
    }
    // -anisotropic
    // -anisotropic <sigma>
    // -anisotropic <sigma> <tensor_name>
    // -anisotropic <sigma> <sigma2>
    // -anisotropic <sigma> <sigma2> <e2_name>
    // -anisotropic <sigma> <sigma2> <e1_name> <e2_name>
    else if (OPTION("-anisotropic")) {
      weighting = irtkPolyDataSmoothing::AnisotropicGaussian;
      if (HAS_ARGUMENT) {
        const char *arg = ARGUMENT;
        if (!FromString(arg, sigma)) {
          cerr << "Invalid -anisotropic standard deviation argument: " << arg << endl;
          exit(1);
        }
      }
      if (HAS_ARGUMENT) {
        const char *arg = ARGUMENT;
        if (FromString(arg, sigma2)) {
          if (HAS_ARGUMENT) {
            e2_name = ARGUMENT;
            if (HAS_ARGUMENT) {
              e1_name = e2_name;
              e2_name = ARGUMENT;
            } else {
              if (polydata->GetPointData()->HasArray(irtkPolyDataCurvature::MINIMUM_DIRECTION)) {
                e1_name = irtkPolyDataCurvature::MINIMUM_DIRECTION;
              }
            }
          } else {
            if (polydata->GetPointData()->HasArray(irtkPolyDataCurvature::MINIMUM_DIRECTION)) {
              e1_name = irtkPolyDataCurvature::MINIMUM_DIRECTION;
            }
            if (polydata->GetPointData()->HasArray(irtkPolyDataCurvature::MAXIMUM_DIRECTION)) {
              e2_name = irtkPolyDataCurvature::MAXIMUM_DIRECTION;
            }
          }
        } else {
          sigma2      = sigma;
          tensor_name = ARGUMENT;
        }
      } else {
        if (polydata->GetPointData()->HasArray(irtkPolyDataCurvature::TENSOR)) {
          tensor_name = irtkPolyDataCurvature::TENSOR;
        }
      }
    }
    else if (OPTION("-exclnode") || OPTION("-adjacent")) adjacent_only = true;
    else if (OPTION("-areaweighted")) smooth_points = area_weighted = true;
    else if (OPTION("-track")) trackingOn = true;
    else if (OPTION("-threshold")) smoothnessThreshold = atof(ARGUMENT);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (!smooth_points && !area_weighted && scalar_names.empty()) {
    smooth_points = true;
    area_weighted = true;
  }

  if (verbose) {
    cout << "Input:          " << input_name << endl;
    cout << "Output:         " << output_name << endl;
    cout << "Iterations:     " << noOfIterations << endl;
    cout << "Relax factor:   " << relaxationFactor << endl;
    cout << "Smooth points:  " << ToString(smooth_points) << endl;
    cout << "Smooth scalars:";
    if (scalar_names.empty()) cout << " No";
    for (size_t i = 0; i < scalar_names.size(); ++i) cout << " " << scalar_names[i];
    cout << endl;
  }

  // Normals are required by Tosun's method
  if (area_weighted) {
    vtkSmartPointer<vtkPolyDataNormals> filter = vtkSmartPointer<vtkPolyDataNormals>::New();
    filter->SplittingOff();
    SetVTKInput(filter, polydata);
    filter->Update();
    polydata = filter->GetOutput();
  }

  // Build links
  polydata->BuildLinks();

  // Smooth node positions and/or scalar data
  if (area_weighted) {
    if (!scalar_names.empty()) {
      cerr << "-areaweighted smoothing of scalar data not implemented" << endl;
      exit(1);
    }
    AreaWeightedLaplacianSmoothing(polydata, noOfIterations, relaxationFactor, smoothnessThreshold, trackingOn);
  } else {
    if (smooth_points && trackingOn) {
      cerr << "-track only implemented for -areaweighted Laplacian smoothing" << endl;
      exit(1);
    }
    irtkPolyDataSmoothing smoother;
    smoother.Input(polydata);
    smoother.NumberOfIterations(noOfIterations);
    smoother.Lambda(relaxationFactor);
    smoother.Sigma(-sigma); // negative: multiple of avg. edge length
    smoother.MaximumDirectionSigma(-sigma2);
    smoother.Weighting(weighting);
    if (tensor_name) smoother.GeometryTensorName(tensor_name);
    if (e1_name    ) smoother.MinimumDirectionName(e1_name);
    if (e1_name    ) smoother.MaximumDirectionName(e2_name);
    smoother.SmoothPoints(smooth_points);
    smoother.SmoothArrays(scalar_names);
    smoother.AdjacentValuesOnly(adjacent_only);
    smoother.Run();
    polydata = smoother.Output();
  }

  // Recompute normals if node positions changed
  if (smooth_points && polydata->GetPointData()->GetNormals() != NULL) {
    if (verbose) cout << "\nRecalculating normals...", cout.flush();
    vtkSmartPointer<vtkPolyDataNormals> filter = vtkSmartPointer<vtkPolyDataNormals>::New();
    filter->SplittingOff();
    SetVTKInput(filter, polydata);
    filter->Update();
    polydata = filter->GetOutput();
    if (verbose) cout << " done" << endl;
  }

  // Save output mesh
  return WritePolyData(output_name, polydata) ? 0 : 1;
}


#else // HAS_VTK

#include <iostream>
using namespace std;
int main(int argc, char *argv[])
{
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
  exit(1);
}

#endif // HAS_VTK
