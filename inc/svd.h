/*
 * This file is part of Susa.

 * Susa is free software: you can redistribute it and/or modify
 * it under the terms of the Lesser GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * at your option) any later version.

 * Susa is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * Lesser GNU General Public License for more details.

 * You should have received a copy of the Lesser GNU General Public License
 * along with Susa.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file svd.h
 * @brief Singular Value Decomposition (SVD)
 * @author Behrooz, Kamary Aliabadi
 * @version 1.0.0
 *
 * @defgroup LALG Linear Algebra
 */

#ifndef SVD_H
#define SVD_H

namespace susa {


/**
 * @brief Singular Value Decomposition (SVD) 
 *
 * Takes a m by n <i>real</i> matrix a and decomposes it into udv, where u,v are
 * left and right orthogonal transformation matrices, and d is a 
 * diagonal matrix of singular values.
 *
 *
 * A = U * W * V
 *
 * @param mat_arg_a The input matrix A that will be replaced by U
 * @param mat_arg_w The diagonal values of W (eigenvalues of A) 
 * @param mat_arg_v The right orthogonal transformation matrix
 *
 * @ingroup LALG
*/
int svd(matrix <float> &mat_arg_a, matrix <float> &mat_arg_w, matrix <float> &mat_arg_v);

} // NAMESPACE SUSA
#endif // SVD_H
