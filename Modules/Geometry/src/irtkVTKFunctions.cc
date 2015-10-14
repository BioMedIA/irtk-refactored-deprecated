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

#include <irtkGeometry.h>

#ifdef HAS_VTK

void GetConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, int seed, vtkSmartPointer<vtkIdList> connectedVertices)
{

	//get all cells that vertex 'seed' is a part of
	vtkSmartPointer<vtkIdList> cellIdList =
		vtkSmartPointer<vtkIdList>::New();
	mesh->GetPointCells(seed, cellIdList);

	//cout << "There are " << cellIdList->GetNumberOfIds() << " cells that use point " << seed << endl;
		//loop through all the cells that use the seed point
		for(vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++)
		{

			vtkCell* cell = mesh->GetCell(cellIdList->GetId(i));
			//cout << "The cell has " << cell->GetNumberOfEdges() << " edges." << endl;

			//if the cell doesn't have any edges, it is a line
			if(cell->GetNumberOfEdges() <= 0)
			{
				continue;
			}

			for(vtkIdType e = 0; e < cell->GetNumberOfEdges(); e++)
			{
				vtkCell* edge = cell->GetEdge(e);

				vtkIdList* pointIdList = edge->GetPointIds();
				//cout << "This cell uses " << pointIdList->GetNumberOfIds() <<" points" << endl;
				/*
				for(vtkIdType p = 0; p < pointIdList->GetNumberOfIds(); p++)
				{
				cout << "Edge " << i << " uses point " << pointIdList->GetId(p) << endl;
				}
				*/
				if(pointIdList->GetId(0) == seed || pointIdList->GetId(1) == seed)
				{
					if(pointIdList->GetId(0) == seed)
					{
						connectedVertices->InsertUniqueId(pointIdList->GetId(1));
					}
					else
					{
						connectedVertices->InsertUniqueId(pointIdList->GetId(0));
					}
				}
			}
		}
		//cout << "There are " << connectedVertices->GetNumberOfIds() << " points connected to point " << seed << endl;
}

#endif
