#include <irtkCommon.h>

#include <vtkMath.h>
#include <vtkPolyData.h>
#include <vtkCellLocator.h>
#include <vtkKdTreePointLocator.h>

using namespace std;


int main(int argc, char *argv[])
{
  if (argc != 3) {
    cerr << "usage: " << argv[0] << " <target> <source>" << endl;
    return 1;
  }

  vtkSmartPointer<vtkAbstractPointLocator> closest_point;
  vtkSmartPointer<vtkAbstractCellLocator>  closest_cell;
  vtkSmartPointer<vtkPolyData>             target, source;

  // read target polydata
  target = ReadPolyData(argv[1]);
  if (target->GetNumberOfPoints() == 0) {
    cerr << "Target has no points!" << endl;
    return 1;
  }

  // read source polydata
  source = ReadPolyData(argv[2]);
  if (source->GetNumberOfPoints() == 0) {
    cerr << "Source has no points!" << endl;
    return 1;
  }
 
  // find closest source points
  int subId;
  vtkIdType j, cellId;
  double a[3], b[3];

  double sum_point_dist = .0;
  double min_point_dist = +numeric_limits<double>::infinity();
  double max_point_dist = -numeric_limits<double>::infinity();
  vector<double> point_dist(target->GetNumberOfPoints(), .0);
 
  double sum_cell_dist = .0;
  double min_cell_dist = +numeric_limits<double>::infinity();
  double max_cell_dist = -numeric_limits<double>::infinity();
  vector<double> cell_dist(target->GetNumberOfPoints(), .0);

  closest_point = vtkSmartPointer<vtkKdTreePointLocator>::New();
  closest_point->SetDataSet(source);
  closest_point->BuildLocator();
  if (source->GetNumberOfCells() > 0) {
    closest_cell = vtkSmartPointer<vtkCellLocator>::New();
    closest_cell->SetDataSet(source);
    closest_cell->BuildLocator();
  }
 
  for (vtkIdType i = 0; i < target->GetNumberOfPoints(); ++i) {
    target->GetPoint(i, a);
    j = closest_point->FindClosestPoint(a);
    source->GetPoint(j, b);

    point_dist[i] = sqrt(vtkMath::Distance2BetweenPoints(a, b));
    sum_point_dist += point_dist[i];
    if (point_dist[i] < min_point_dist) min_point_dist = point_dist[i];
    if (point_dist[i] > max_point_dist) max_point_dist = point_dist[i];
 
    if (closest_cell == NULL) continue;
    closest_cell->FindClosestPoint(a, b, cellId, subId, cell_dist[i]);
    cell_dist[i] = sqrt(cell_dist[i]);
    sum_cell_dist += cell_dist[i];
    if (cell_dist[i] < min_cell_dist) min_cell_dist = cell_dist[i];
    if (cell_dist[i] > max_cell_dist) max_cell_dist = cell_dist[i];
  }
 
  closest_point = NULL;
  closest_cell  = NULL;
 
  // calculate Hausdorff distance
  double hausdorff_dist = max(max_cell_dist, max_point_dist);
  if (target->GetNumberOfCells() > 0) {
    double dist, max_target_dist = -numeric_limits<double>::infinity();
    closest_cell = vtkSmartPointer<vtkCellLocator>::New();
    closest_cell->SetDataSet(target);
    closest_cell->BuildLocator();
    for (vtkIdType i = 0; i < source->GetNumberOfPoints(); ++i) {
      source->GetPoint(i, a);
      closest_cell->FindClosestPoint(a, b, cellId, subId, dist);
      dist = sqrt(dist);
      if (dist > max_target_dist) max_target_dist = dist;
    }
    hausdorff_dist = max(hausdorff_dist, max_target_dist);
    closest_cell = NULL;
  }
 
  // print distances and summary
  if (verbose) {
    cout << "ID,PointDistance,CellDistance" << endl;
    for (vtkIdType i = 0; i < target->GetNumberOfPoints(); ++i) {
      cout << (i+1) << "," << point_dist[i] << "," << cell_dist[i] << endl;
    }
    cout << endl;
  }
  cout << "Point distance: min=" << min_point_dist << ", max=" << max_point_dist << ", avg=" << sum_point_dist / target->GetNumberOfPoints() << endl;
  cout << "Cell  distance: min=" << min_cell_dist  << ", max=" << max_cell_dist  << ", avg=" << sum_cell_dist  / target->GetNumberOfPoints() << endl;
  cout << endl;
  cout << "Hausdorff distance: " << hausdorff_dist << endl;
  return 0;
}
