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
 * @author Behrooz Kamary
 *
 * @defgroup LALG Linear Algebra
 */

#ifndef SUSA_SVD_H
#define SUSA_SVD_H

namespace susa {


/**
 * @brief Singular Value Decomposition (SVD)
 *
 * Takes a m by n <i>real</i> matrix and decomposes it into U * S * V, where U and V are
 * left and right orthogonal transformation matrices, and S is diagonal matrix of singular values.
 *
 *
 * A = U * S * V
 *
 * @param mat_arg_a The input matrix A
 * @param mat_arg_u The left orthogonal transformation matrix U
 * @param mat_arg_s The diagonal values of the eignevalue matrix S (eigenvalues of A)
 * @param mat_arg_v The right orthogonal transformation matrix V
 *
 * @ingroup LALG
 */
int svd(const matrix <float> &mat_arg_a,
              matrix <float> &mat_arg_u,
              matrix <float> &mat_arg_s,
              matrix <float> &mat_arg_v);

} // NAMESPACE SUSA
#endif // SUSA_SVD_H
