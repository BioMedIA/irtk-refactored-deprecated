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

#include <algorithm>

// =============================================================================
// Miscellaneous functions
// =============================================================================

// -----------------------------------------------------------------------------
double median(float q, float p0, float p1)
{
    float m;
    m = p0>p1?p0:p1;
    m = q > m? m:q;
    return m;
}

// -----------------------------------------------------------------------------
void sort2(int index, float *array1, float *array2)
{
    int* indices = new int[index];
    float* array1_store = new float[index];
    float* array2_store = new float[index];

    for (int i = 0; i < index; i++) {
       // create array with indices to be sorted
       indices[i] = i;
       // keep a back up of the original arrays
       array1_store[i] = array1[i];
       array2_store[i] = array2[i];
    }

    // sort the first array and keep the index order
    class sort_indices
    {
    private:
        float* _array;

    public:
        sort_indices(float* array) : _array(array) {}
        bool operator()(int i, int j) { return _array[i] < _array[j]; }
    };

    std::sort(indices, indices + index, sort_indices(array1));

    // update the arrays
    for (int i = 0; i < index; i++) {
       array1[i+1] = array1_store[indices[i]];
       array2[i+1] = array2_store[indices[i]];
    }

    delete[] indices;
    delete[] array1_store;
    delete[] array2_store;
}

// -----------------------------------------------------------------------------
/// weight must be normalized to sum = 1 before entering here.
/// implemented fast weighted median selector algorithm from Y. Li and S. Osher.
double weightedmedian(int index, double lambda, double f, float *neighbors, float *weight)
{
    int i,pos;
    double p0,p1,csum;

    // sort the value according to neighbors
    sort2(index-1,neighbors,weight);
    
    csum = 1;

    if(lambda > 0){

        p0 = f+(2.0*csum/lambda);
        // find new median
        for(i = 1; i < index; i++){
             csum -= 2.0*weight[i];
             p1 = f+(2.0*csum/lambda);
             p0 = median(p0,p1,neighbors[i]);
        }

        return p0;
    }else{
        // 1/lambda = infinite, use a simple version.
        pos = 0;

        for(i = 1; i <= index; i++){
            if(csum > 0) pos++;
            if(csum < 0) pos--;
            if(i != index){
                csum = csum - 2.0*weight[i];
            }
        }
        return neighbors[(index+pos)/2];
    }
}
