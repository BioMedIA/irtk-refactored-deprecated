#if (defined HAS_VTK)

#ifdef WIN32
#include <map>
#else
#include <map>
#endif

#include <irtkImage.h>
#include <irtkImageFunction.h>
#include <irtkEuclideanDistanceTransform.h>

#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkIntArray.h>

char *mask_name = NULL, *output_name = NULL;
char *input_surface_name = NULL;

void usage()
{
  cerr << "\t" << endl;
  cerr << "\tUsage:  polydatamask [mask] [surfaceIn] [surfaceOut] <options>" << endl;
  cerr << "\t" << endl;
  cerr << "\tRemove points from mesh [surfaceIn] that are near to the " << endl;
  cerr << "\tboundary of a structure that is labelled in the binary image [mask].  " << endl;
  cerr << "\tWrite the result to the mesh surfaceOut]." << endl;
  cerr << "\t" << endl;
  cerr << "\tOptions:" << endl;
  cerr << "\t-threshold t    The distance from the mask boundary for which" << endl;
  cerr << "\t                the surface points should be removed. The default" << endl;
  cerr << "\t                is 0 which means that points on the boundary or " << endl;
  cerr << "\t                inside the mask are removed.  Setting the " << endl;
  cerr << "\t                threshold to higher positive values removes points " << endl;
  cerr << "\t                further away from the boundary i.e. more points are " << endl;
  cerr << "\t                removed.  Negative values means that the removed are" << endl;
  cerr << "\t                inside the boundary." << endl;
  cerr << "\t" << endl;
  cerr << "\t" << endl;
  exit(0);
}

int main(int argc, char **argv)
{

  if (argc < 4){
    usage();
  }

  int i;
  bool ok;
  int noOfPoints;
  int newPtCount;
  double p[3];
  int* deletedPt;
  int *newPtIds;
  int newId;
  int j;
  int noOfFaces, newFaceCount;
  vtkIdType npts = 0;
  vtkIdType *oldPtId;
  int deleteFace;
  double surfaceBounds[6];
  double xmin, xmax, ymin, ymax, zmin, zmax;
  double val;
  double x, y, z;
  double threshold;
  threshold = 0.0;

  ////////////////////////////////////////////////
  // Main arguments.
  mask_name = argv[1];
  argc--;
  argv++;
  input_surface_name = argv[1];
  argc--;
  argv++;
  output_name = argv[1];
  argc--;
  argv++;

  /////////////////////////////////////////////////
  // Remaining arguments.
  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-threshold") == 0)){
      argc--;
      argv++;
      threshold = atof(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false){
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
    }
  }

  /////////////////////////////////////////
  // Distance transform.

  cerr << "Generating distance transform." << endl;

  irtkRealImage* dmap = new irtkRealImage(mask_name);

  irtkRealImage* dmapOut = new irtkRealImage(mask_name);
  irtkRealImage* dmapIn = new irtkRealImage(mask_name);

  irtkEuclideanDistanceTransform<irtkRealPixel> *edt = NULL;
  edt = new irtkEuclideanDistanceTransform<irtkRealPixel>
    (irtkEuclideanDistanceTransform<irtkRealPixel>::irtkDistanceTransform3D);

  irtkRealPixel* ptr2dmap, *ptr2dmapIn, *ptr2dmapOut;

  // dmap image actually contains mask (for now).
  // Threshold mask.
  ptr2dmap = dmap->GetPointerToVoxels();
  for (i = 0; i < dmap->GetNumberOfVoxels(); ++i){
    if (*ptr2dmap > 0){
      *ptr2dmap = 1;
    } else {
      *ptr2dmap = 0;
    }
    ++ptr2dmap;
  }

  edt->SetInput(dmap);
  edt->SetOutput(dmapOut);
  edt->Run();

  // Invert mask.
  ptr2dmap = dmap->GetPointerToVoxels();
  for (i = 0; i < dmap->GetNumberOfVoxels(); ++i){
    if (*ptr2dmap > 0){
      *ptr2dmap = 0;
    } else {
      *ptr2dmap = 1;
    }
    ++ptr2dmap;
  }

  edt->SetInput(dmap);
  edt->SetOutput(dmapIn);
  edt->Run();

  // Now we can calculate the actual dmap.
  ptr2dmap = dmap->GetPointerToVoxels();
  ptr2dmapIn = dmapIn->GetPointerToVoxels();
  ptr2dmapOut = dmapOut->GetPointerToVoxels();
  for (i = 0; i < dmap->GetNumberOfVoxels(); ++i){
    *ptr2dmap = sqrt(*ptr2dmapOut) - sqrt(*ptr2dmapIn);
    ++ptr2dmap;
    ++ptr2dmapIn;
    ++ptr2dmapOut;
  }

  ////////////////////////////////////////
  // Create interpolator.
  irtkImageFunction *interp = NULL;
  interp = new irtkLinearInterpolateImageFunction;
  interp->SetInput(dmap);
  interp->Initialize();

  /////////////////////////////////////////
  vtkPolyData *surface = vtkPolyData::New();

  cerr << "Reading surface ... " << endl;
  // Read surface
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(input_surface_name);
  reader->Update();
  surface = reader->GetOutput();

  //////////////////////////////////
  //Check that surface does not go outside FOV of image.

  surface->ComputeBounds();
  surface->GetBounds(surfaceBounds);

  xmin = surfaceBounds[0];
  xmax = surfaceBounds[1];
  ymin = surfaceBounds[2];
  ymax = surfaceBounds[3];
  zmin = surfaceBounds[4];
  zmax = surfaceBounds[5];

  dmap->WorldToImage(xmin, ymin, zmin);
  dmap->WorldToImage(xmax, ymax, zmax);
  if (xmin < -0.5 || xmax > dmap->GetX()-0.5 ||
      ymin < -0.5 || ymax > dmap->GetY()-0.5 ||
      zmin < -0.5 || zmax > dmap->GetZ()-0.5){
        cerr << "Surface outside bounds of image." << endl;
        exit(1);
  }

  //////////////////////////////////

  vtkPoints *oldPoints = vtkPoints::New();
  vtkPoints *newPoints = vtkPoints::New();

  noOfPoints= surface->GetNumberOfPoints();

  oldPoints = surface->GetPoints();

  newPtIds = new int[noOfPoints];
  deletedPt = new int[noOfPoints];

  for (i = 0; i < noOfPoints; ++i){
    deletedPt[i] = 0;
  }

  newPtCount = 0;
  for (i = 0; i < noOfPoints; ++i){

    oldPoints->GetPoint(i, p);
    x = p[0];
    y = p[1];
    z = p[2];

    dmap->WorldToImage(x, y, z);
    val = interp->Evaluate(x, y, z);

    if (val < threshold){
      continue;
    }

    ++newPtCount;

  }

  newPoints->SetNumberOfPoints(newPtCount);

  newPtCount = 0;
  for (i = 0; i < noOfPoints; ++i){

    oldPoints->GetPoint(i, p);
    x = p[0];
    y = p[1];
    z = p[2];

    dmap->WorldToImage(x, y, z);
    val = interp->Evaluate(x, y, z);

    if (val < threshold){
      deletedPt[i] = 1;
    } else {
      // Can keep the point.
      newPoints->InsertPoint(newPtCount, p);
      newPtIds[i] = newPtCount;
      ++newPtCount;
    }
  }

  cerr << "Points in original surface: " << noOfPoints << endl;
  cerr << "Points copied over        : " << newPtCount << endl;

  /////////////////////////////////////////////
  vtkCellArray* oldFaces = vtkCellArray::New();
  vtkCellArray* newFaces = vtkCellArray::New();


  oldFaces = surface->GetPolys();

  noOfFaces = oldFaces->GetNumberOfCells();

  newFaces->Initialize();
  newFaces->Allocate(noOfFaces, 0);
  newFaces->Reset();

  oldFaces->InitTraversal();
  newFaces->InitTraversal();
  newFaceCount = 0;

  for (i = 0; i < noOfFaces; ++i){
    oldFaces->GetNextCell(npts, oldPtId);

    deleteFace = 0;
    for (j = 0; j < npts; ++j){
      if(deletedPt[oldPtId[j]] == 1){
        deleteFace = 1;
      }
    }

    if (deleteFace == 1)
      continue;

    newFaces->InsertNextCell(npts);
    for (j = 0; j < npts; ++j){
      newId = newPtIds[oldPtId[j]];
      newFaces->InsertCellPoint(newId);
    }
    newFaces->UpdateCellCount(npts);
    ++newFaceCount;
  }
  newFaces->Squeeze();

  cerr << "Faces in original surface : " << noOfFaces << endl;
  cerr << "Faces copied over         : " << newFaceCount << endl;


  /////////////////////////////////////////////////////

  vtkPolyData* output = vtkPolyData::New();
  output->SetPoints(newPoints);
  output->SetPolys(newFaces);

  cerr << "Writing surface ... " << endl;
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
#if VTK_MAJOR_VERSION >= 6
  writer->SetInputData(output);
#else
  writer->SetInput(output);
#endif
  writer->SetFileName(output_name);
  writer->SetFileTypeToBinary();
  writer->Write();

  delete [] deletedPt;
  delete [] newPtIds;

  return 0;

}

/*

// Returns a floating point 0/1 image to show a label.  The image needs to
// be floating point so that it can later be used in a distance map filter.
irtkRealImage *getBinaryLabelImage(irtkGreyImage *labelImage, short label){

  int i, noOfVoxels;
  irtkRealPixel *ptr2mask;
  irtkGreyPixel *ptr2labels;
  irtkRealImage *mask = new irtkRealImage(*labelImage);
  ptr2mask = mask->GetPointerToVoxels();
  ptr2labels = labelImage->GetPointerToVoxels();

  noOfVoxels = labelImage->GetNumberOfVoxels();

  for (i = 0; i < noOfVoxels; ++i, ++ptr2mask, ++ptr2labels){
    if (*ptr2labels == label){
      *ptr2mask = 1;
    } else {
      *ptr2mask = 0;
    }
  }
  return mask;
}

void usage()
{
  cerr << "Usage:  polydataassignscalars [labelImage] [surfaceIn] [surfaceOut] <options>" << endl;
  cerr << "" << endl;
  cerr << "Assign scalars to the vertices of surfaceIn.  The scalar assigned" << endl;
  cerr << "is the label of the nearest voxel in labelImage to the vertex." << endl;
  cerr << "" << endl;
  cerr << "Write the result to surfaceOut." << endl;
  cerr << "Options: " << endl;
  cerr << "" << endl;

  exit(1);
}

int main(int argc, char **argv)
{
  int i, ok, noOfPoints;
  double surfaceBounds[6];
  double xmin, xmax, ymin, ymax, zmin, zmax;
  double pt[3];
  int val, zeroCount;
  int noOfVoxels, noOfLabels, currLabel;

  if (argc < 4){
    usage();
  }

  // Parse image
  input_image_name  = argv[1];
  argc--;
  argv++;
  input_surface_name  = argv[1];
  argc--;
  argv++;
  output_name  = argv[1];
  argc--;
  argv++;

  // Parse remaining arguments
  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-XXX") == 0)){
      argc--;
      argv++;
      // do stuff and maybe argv++ etc.
      ok = true;
    }
    if (!ok){
      cerr << "Cannot parse argument " << argv[1] << endl;
      usage();
    }
  }

  cerr << "Reading image ..." << endl;
  irtkGreyImage *labelImage = new irtkGreyImage(input_image_name);

  vtkPolyData *surface = vtkPolyData::New();

  cerr << "Reading surface ... " << endl;
  // Read surface
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(input_surface_name);
  reader->Update();
  surface = reader->GetOutput();

  noOfPoints= surface->GetNumberOfPoints();

  //cerr << "Making space for scalars ..., no. of points =  " << noOfPoints << endl;
  vtkIntArray *scalars = vtkIntArray::New();
  scalars->SetNumberOfComponents(1);
  scalars->SetNumberOfTuples(noOfPoints);

  //cerr << "Creating interpolator " << endl;
  irtkImageFunction *interp = NULL;
  interp = new irtkNearestNeighborInterpolateImageFunction;
  interp->SetInput(labelImage);
  interp->Initialize();

  //Check that surface does not go outside fov of label image.
  surface->ComputeBounds();
  surface->GetBounds(surfaceBounds);
  xmin = surfaceBounds[0];
  xmax = surfaceBounds[1];
  ymin = surfaceBounds[2];
  ymax = surfaceBounds[3];
  zmin = surfaceBounds[4];
  zmax = surfaceBounds[5];
  cerr << "Bounds of surface : ";
  cerr << "(" << xmin << ", " << ymin << ", " << zmin << ") and ";
  cerr << "(" << xmax << ", " << ymax << ", " << zmax << ")" << endl;

  labelImage->WorldToImage(xmin, ymin, zmin);
  labelImage->WorldToImage(xmax, ymax, zmax);
  cerr << "In image coords : ";
  cerr << "(" << xmin << ", " << ymin << ", " << zmin << ") and ";
  cerr << "(" << xmax << ", " << ymax << ", " << zmax << ")" << endl;
  if (xmin < -0.5 || xmax > labelImage->GetX()-0.5 ||
      ymin < -0.5 || ymax > labelImage->GetY()-0.5 ||
      zmin < -0.5 || zmax > labelImage->GetZ()-0.5){
        cerr << "Surface outside bounds of image." << endl;
        exit(1);
  }

  irtkRealImage *currLabelMask;
  irtkGreyImage *dilatedLabels;
  irtkRealImage *minDmap;
  irtkRealImage *currDmap;

  // Identify the number of distinct labels in the image.
  irtkRealPixel *ptr2dmap;
  irtkRealPixel *ptr2minDmap;
  irtkGreyPixel *ptr2label;
  irtkGreyPixel *ptr2dilatedLabel;

  countMap labelCount;

  // Count up different labels so we can identify the number of distinct labels.
  noOfVoxels = labelImage->GetNumberOfVoxels();
  ptr2label  = labelImage->GetPointerToVoxels();
  for (i = 0; i < noOfVoxels; ++i, ++ptr2label){
    if (*ptr2label > 0){
      ++labelCount[*ptr2label];
    }
  }

  cerr << "No. of Voxels          : " << noOfVoxels << endl;
  cerr << "Label Counts " << endl;
  for (iter = labelCount.begin(); iter != labelCount.end(); ++iter){
    cerr << iter->first << "\t" << iter->second << endl;
  }
  noOfLabels = labelCount.size();
  cerr << "No. of distinct labels : " << noOfLabels << endl;

  // Using the distance maps.
  minDmap  = new irtkRealImage(*labelImage);
  currDmap = new irtkRealImage(*labelImage);

  // Note that the dilated labels are initialised to the given label image.
  // I.e. the original labels are left alone and we seek to assign labels to
  // the zero voxels based on closest labeled voxels.
  dilatedLabels = new irtkGreyImage(*labelImage);

  // Initialise the minimum distance map.
  ptr2minDmap = minDmap->GetPointerToVoxels();
  for (i = 0; i < noOfVoxels; ++i, ++ptr2minDmap){
    *ptr2minDmap = FLT_MAX;
  }

  // Single distance transform filter for all labels.
  irtkEuclideanDistanceTransform<irtkRealPixel> *edt = NULL;
  edt = new irtkEuclideanDistanceTransform<irtkRealPixel>
    (irtkEuclideanDistanceTransform<irtkRealPixel>::irtkDistanceTransform3D);

  cerr << "Finding distance maps ..." << endl;
  cerr << "Current label : " << endl;
  for (iter = labelCount.begin(); iter != labelCount.end(); ++iter){
    currLabel = iter->first;
    cerr << "  " << currLabel << endl;

    // There is a new operator used for currLabelMask in the following function call.
    currLabelMask = getBinaryLabelImage(labelImage, currLabel);

    edt->SetInput(currLabelMask);
    edt->SetOutput(currDmap);
    edt->Run();

    ptr2minDmap      = minDmap->GetPointerToVoxels();
    ptr2dmap         = currDmap->GetPointerToVoxels();
    ptr2label        = labelImage->GetPointerToVoxels();
    ptr2dilatedLabel = dilatedLabels->GetPointerToVoxels();

    for (i = 0; i < noOfVoxels; ++i, ++ptr2minDmap, ++ptr2dmap, ++ptr2label, ++ptr2dilatedLabel){
      if (*ptr2label == 0 && *ptr2dmap < *ptr2minDmap){
        *ptr2minDmap = *ptr2dmap;
        *ptr2dilatedLabel = currLabel;
      }
    }

    // Tidy up.
    delete currLabelMask;
  }

  cerr << "Assigning scalars using dilated labels ... " << endl;
  interp->SetInput(dilatedLabels);
  interp->Initialize();
  zeroCount = 0;
  for (i = 0; i < noOfPoints; ++i){
    surface->GetPoint(i, pt);
    dilatedLabels->WorldToImage(pt[0], pt[1], pt[2]);

    val = (int) round(interp->Evaluate(pt[0], pt[1], pt[2]));
    scalars->SetTuple1(i,val);

    //if (val == 0){
    //  ++zeroCount;
    //}
  }

  //cerr << "Zero count : " << zeroCount << " = " << 100.0 * zeroCount / ((double) noOfPoints) << "%" << endl;
  cerr << "Updating surface ... " << endl;
  scalars->Modified();
  scalars->SetName("Labels");
  surface->GetPointData()->AddArray(scalars);
  surface->Update();

  cerr << "Writing surface ... " << endl;
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput(surface);
  writer->SetFileName(output_name);
  writer->SetFileTypeToBinary();
  writer->Write();

  return 0;
}


*/


#else
int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif

