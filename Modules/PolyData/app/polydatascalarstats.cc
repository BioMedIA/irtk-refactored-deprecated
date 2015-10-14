/* ========================================================
 * Get the statistics on scalars associated with a polydata
 * object.
 *
 * Paul Aljabar, July 2010.
======================================================== */

#if (defined HAS_VTK)

#include <irtkImage.h>

#include <vtkPolyData.h>
#include <vtkTriangle.h>
#include <vtkIdList.h>
#include <vtkPolyDataReader.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>

char *input_name = NULL;
char *scalar_name = NULL;
char *mask_name = NULL;

void usage()
{
  cerr << " " << endl;
  cerr << " Usage: polydatascalarstats [in] <options>" << endl;
  cerr << " " << endl;
  cerr << " Options: " << endl;
  cerr << " -q            Give the same output as normal but just the " << endl;
  cerr << "               numbers on a space separated line." << endl;
  cerr << " -name [name]  Name of scalars for which stats are required." << endl;
  cerr << " -mask [name]  Name of scalars to use as a mask, stats only" << endl;
  cerr << "               calculated where mask value is > 0." << endl;
  cerr << " " << endl;
  cerr << " " << endl;

  exit(1);
}

int main(int argc, char **argv)
{
  bool ok;
  double sum, sumSq, sumAbs;
  int i;
  int noOfPoints, count;
  double minVal, maxVal;
  double mean, var, sd, val;
  double meanSq, meanAbs, varAbs, sdAbs;
  int quiet = false;

  // For weighting by cell area.
  vtkTriangle *triangle;
  vtkIdType *cells;
  unsigned short noOfCells;
  int k;
  vtkIdList *ptIds;
  double v1[3], v2[3], v3[3], triangleArea;
  double totalArea;
  totalArea = 0.0;
  double integral = 0.0;
  int	nonTriangleFaces = 0;

  if (argc < 2){
    usage();
  }

  // Parse image
  input_name  = argv[1];
  argc--;
  argv++;

  // Parse remaining arguments
  while (argc > 1){
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-q") == 0)){
      argc--;
      argv++;
      quiet = true;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-name") == 0)){
       argc--;
       argv++;
       scalar_name = argv[1];
       argc--;
       argv++;
       ok = true;
     }
    if ((ok == false) && (strcmp(argv[1], "-mask") == 0)){
       argc--;
       argv++;
       mask_name = argv[1];
       argc--;
       argv++;
       ok = true;
     }
     if (!ok){
      cerr << "Cannot parse argument " << argv[1] << endl;
      usage();
    }
  }

  vtkPolyData *input = vtkPolyData::New();

  // Read surface
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName(input_name);
  reader->Update();
  input = reader->GetOutput();
  input->BuildCells();
  input->BuildLinks();

  noOfPoints= input->GetNumberOfPoints();

  vtkFloatArray *scalars;// = vtkFloatArray::New();
  int ind;

  if (scalar_name == NULL){
    scalars = (vtkFloatArray*) input->GetPointData()->GetArray(0);
    if (scalars == NULL){
      cerr << "No scalars available." << endl;
      exit(1);
    }
  } else {
    scalars = (vtkFloatArray*) input->GetPointData()->GetArray(scalar_name, ind);

    if (ind == -1 || scalars == NULL){
      cerr << "Scalars unavailable with name " << scalar_name << endl;
      exit(1);
    }
  }

  int *mask = new int[noOfPoints];
  for (i = 0; i < noOfPoints; ++i){
  	mask[i] = 1;
  }

  if (mask_name != NULL){
    vtkFloatArray *mask_scalars = (vtkFloatArray*) input->GetPointData()->GetArray(mask_name, ind);
    if (ind == -1 || mask_scalars == NULL){
      cerr << "Masking scalars unavailable with name " << mask_name << endl;
      exit(1);
    }
    if (mask_scalars->GetNumberOfComponents() > 1){
    	cerr << "Masking scalars " << mask_scalars->GetName() << " has more than one component." << endl;
    	exit(1);
    }

    for (i = 0; i < noOfPoints; ++i){
      if (mask_scalars->GetTuple1(i) > 0){
      	continue;
      }
      mask[i] = 0;
    }

  }

  if (scalars->GetNumberOfComponents() > 1){
  	cerr << "Scalars " << scalars->GetName() << " has more than one component." << endl;
  	exit(1);
  }

  sum = 0.0;
  sumSq = 0.0;
  sumAbs = 0.0;
  minVal = FLT_MAX;
  maxVal = -1 * minVal;

  count = 0;

  for (i = 0; i < noOfPoints; ++i){

  	if (mask[i] <= 0){
  		continue;
  	}
  	++count;

    val = scalars->GetTuple1(i);

    sum += val;
    sumSq += val*val;
    sumAbs += fabs(val);
    if (minVal > val)
      minVal = val;
    if (maxVal < val)
      maxVal = val;

    input->GetPointCells(i, noOfCells, cells);

    if ( cells == NULL )
      continue;

    for (k = 0; k < noOfCells; ++k){
      triangle = vtkTriangle::SafeDownCast(input->GetCell(cells[k]));

      if ( triangle != NULL ){
        ptIds = triangle->GetPointIds();

        input->GetPoint(ptIds->GetId(0), v1);
        input->GetPoint(ptIds->GetId(1), v2);
        input->GetPoint(ptIds->GetId(2), v3);

        triangleArea = vtkTriangle::TriangleArea(v1, v2, v3) / 3.0;
        integral += triangleArea * val;
        totalArea += triangleArea;

      } else {
      	++nonTriangleFaces;
      }

    }
  }

  if (count < 1){
  	cerr << "Zero points remain after masking, exiting " << endl;
  	exit(0);
  }

  mean    = sum / ((double) count);
  meanSq  = sumSq / ((double) count);
  meanAbs = sumAbs / ((double) count);

  var = meanSq - (mean*mean);
  sd  = sqrt(var);

  varAbs = meanSq - (meanAbs*meanAbs);
  sdAbs  = sqrt(varAbs);

  // Normalisation based on the area of a sphere. Becomes an L2 norm if input stat is a squared measurement.
  double normedIntegral = sqrt(integral / 4.0 / M_PI);

  if (quiet){
    cout << noOfPoints;
    cout << " " << count;
    cout << " " << mean;
    cout << " " << meanSq;
    cout << " " << sd;
    cout << " " << meanAbs;
    cout << " " << sdAbs;
    cout << " " << minVal << " " << maxVal;
    cout << " " << totalArea << endl;
    cout << " " << integral;
    cout << " " << normedIntegral << endl;
  } else {
    cout << "Scalar name   " << scalars->GetName() << endl;
    cout << "No of pts     " << noOfPoints << endl;
    cout << "After masking " << count << endl;
    cout << "Mean          " << mean << endl;
    cout << "Mean Sq       " << meanSq << endl;
    cout << "S.D.          " << sd << endl;
    cout << "Mean(abs)     " << meanAbs << endl;
    cout << "S.D(abs)      " << sdAbs << endl;
    cout << "Min/Max       " << minVal << " " << maxVal << endl;
    cout << "Area          " << totalArea << endl;
    cout << "Area int      " << integral << endl;
    cout << "sqrt(Area int / 4 pi) " << normedIntegral << endl;
    if (nonTriangleFaces > 0)
    	cout << "Non Triangle Faces " << nonTriangleFaces << endl;
    cout << "" << "" << endl;
  }

  delete [] mask;
  return 0;
}


#else
int main( int argc, char *argv[] ){
  cerr << argv[0] << " needs to be compiled with the VTK library " << endl;
}
#endif

