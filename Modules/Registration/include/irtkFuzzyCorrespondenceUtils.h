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

#include <irtkSparseMatrix.h>


namespace irtkFuzzyCorrespondenceUtils {


/// Type of fuzzy correspondences weight matrix
typedef irtkGenericSparseMatrix<float> WeightMatrix;

/// Sinkhorn-Knopp algorithm for column and row normalization
///
/// \param[in,out] weight  Fuzzy correspondence weight matrix.
/// \param[in]     m       Number of columns to normalize.
/// \param[in]     n       Number of rows to normalize.
/// \param[in]     maxit   Maximum number of iterations.
/// \param[in]     maxerr  Maximum squared error. If non-positive, the maximum
///                        number of iterations is performed. This may be
///                        considerably faster than always checking the error.
void NormalizeWeights(WeightMatrix &weight, int m = -1, int n = -1,
                      int maxit = 10, double maxerr = .0);


} // namespace irtkFuzzyCorrespondenceUtils
