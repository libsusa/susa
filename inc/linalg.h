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
 * @file linalg.h
 * @brief Basic matrix and linear algebra operations on Susa types.
 * @author Behrooz, Aliabadi
 * @version 1.0.0
 *
 * @defgroup LALG Linear Algebra
 */

#ifndef LINALG_H
#define LINALG_H

namespace susa {


/**
 * @brief Multiplication
 *
 * @param mat_argl
 * @param mat_argr
 * @return Multiplication of the two input matrices
 * @ingroup LALG
 */
template <class T> matrix <T> mult( const matrix <T> &mat_argl,const matrix <T> &mat_argr ) {
    // To make it faster this method has been declared as 'friend' of
    // 'matrix' template class so 'mult' apears two times in 'matrix.h'
    // and two times in 'linalg.h'
    unsigned int uint_rows = mat_argl.no_rows();
    unsigned int uint_cols = mat_argr.no_cols();

    matrix <T> mat_ret;

    if (mat_argl._matrix == NULL || mat_argr._matrix == NULL) return mat_ret;

    if (mat_argl.no_cols() == mat_argr.no_rows()) {
        mat_ret = matrix <T> (uint_rows, uint_cols);
        for (unsigned int uint_row = 0; uint_row < uint_rows; uint_row++) {
            for (unsigned int uint_col = 0; uint_col < uint_cols; uint_col++) {
                for (unsigned int uint_index = 0; uint_index < mat_argl.no_cols(); uint_index++)
                    //mat_ret(uint_row, uint_col) += mat_argl(uint_row,uint_index) * mat_argr(uint_index,uint_col);
                    mat_ret._matrix[uint_row + uint_col * mat_ret.uint_rows] += mat_argl._matrix[uint_row + uint_index * mat_argl.uint_rows]
                            * mat_argr._matrix[uint_index + uint_col * mat_argr.uint_rows];
            }
        }
    } else {
        std::cout << std::endl << "[mult()] Matrices' dimensions mismatch.";
#ifdef _SUSA_TERMINATE_ON_ERROR
        exit(1);
#endif
    }
    return mat_ret;
}

/**
 * @brief Transpose
 *
 * @param mat_arg
 * @return Returns transpose of the input matrices
 * @ingroup LALG
 */
template <class T> matrix <T> transpose(const matrix <T> &mat_arg) {
    unsigned int uint_rows = mat_arg.no_rows();
    unsigned int uint_cols = mat_arg.no_cols();
    matrix <T> mat_ret;

    if (mat_arg._matrix == NULL) {
        return mat_ret;
    }
    mat_ret = matrix <T> (uint_cols, uint_rows);
    for (unsigned int uint_col = 0; uint_col < uint_cols; uint_col++) {
        for (unsigned int uint_row = 0; uint_row < uint_rows; uint_row++) {
            //mat_ret(uint_col,uint_row) = mat_arg(uint_row,uint_col);
            mat_ret._matrix[uint_col + uint_row * mat_ret.uint_rows] = mat_arg._matrix[uint_row + uint_col * mat_arg.uint_rows];
        }
    }

    return mat_ret;
}


/**
 * @brief Determinant
 *
 * @param mat_arg
 * @return Returns determinant of the input matrices
 * @ingroup LALG
 */
template <class T> double det(const matrix <T> &mat_arg) {
    double dbl_ret = 0.0;

    unsigned int uint_rows = mat_arg.no_rows();
    unsigned int uint_cols = mat_arg.no_cols();

    if (uint_rows == uint_cols) {
        if (uint_rows == 1) { // 1 x 1
            dbl_ret = mat_arg(0);
        } else if (uint_rows == 2) { // 2 x 2
            dbl_ret = mat_arg(0, 0) * mat_arg(1, 1) - mat_arg(1, 0) * mat_arg(0, 1);
        } else { // 3 x 3 and more
            for (unsigned int  uint_col = 0; uint_col < uint_cols; uint_col++) {
                matrix <T> mat_minor = mat_arg.shrink(0, uint_col);
                if (mat_minor.size() == mat_arg.size()) {
                    // If for any reason something goes wrong in the minor() it returns the matrix without any action.
                    // Therefore size of the input argument and the returned argument would the same. This case cause infinite
                    // iterations. To stop such a unwanted loop this if() returns from method.
                    // THIS METHOD IS NOT STABLE BUT WORKS.
                    std::cout << "[det()] Infinite loop avoided.";
                    return dbl_ret;
                }
                dbl_ret +=  (-2 * (int)(uint_col%2) + 1) * mat_arg(0, uint_col) * det(mat_minor);
            }
        }
    }

    return dbl_ret;
}



/**
 * @brief Concatenation
 *
 * @param mat_argl
 * @param mat_argr
 * @return Returns concatenation of the two input matrices
 * @ingroup LALG
 */
template <class T> matrix <T> concat(const matrix <T> &mat_argl, const matrix <T> &mat_argr) {
    // should be changed 20101224
    matrix <T> mat_ret;

    if ((mat_argl.no_cols() == mat_argr.no_cols() && mat_argl.no_rows() != mat_argr.no_rows()) ||
            (mat_argl.no_cols() == mat_argr.no_cols() && mat_argl.no_rows() == mat_argr.no_rows())) { // Case I

        mat_ret = matrix <T> (mat_argl.no_rows() + mat_argr.no_rows(), mat_argl.no_cols());

        for (unsigned int uint_row = 0; uint_row < mat_argl.no_rows(); uint_row++) {
            for (unsigned int uint_col = 0; uint_col < mat_argl.no_cols(); uint_col++) {
                mat_ret(uint_row,uint_col) = mat_argl(uint_row,uint_col);
            }
        }

        for (unsigned int uint_row = 0; uint_row < mat_argr.no_rows(); uint_row++) {
            for (unsigned int uint_col = 0; uint_col < mat_argr.no_cols(); uint_col++) {
                mat_ret(mat_argl.no_rows() + uint_row,uint_col) = mat_argr(uint_row,uint_col);
            }
        }
    } else if (mat_argl.no_cols() != mat_argr.no_cols() && mat_argl.no_rows() == mat_argr.no_rows()) { // Case II
        mat_ret = matrix <T> (mat_argl.no_rows(), mat_argl.no_cols() + mat_argr.no_cols());

        for (unsigned int uint_row = 0; uint_row < mat_argl.no_rows(); uint_row++) {
            for (unsigned int uint_col = 0; uint_col < mat_argl.no_cols(); uint_col++) {
                mat_ret(uint_row,uint_col) = mat_argl(uint_row,uint_col);
            }
        }

        for (unsigned int uint_row = 0; uint_row < mat_argr.no_rows(); uint_row++) {
            for (unsigned int uint_col = 0; uint_col < mat_argr.no_cols(); uint_col++) {
                mat_ret(uint_row,mat_argr.no_cols() + uint_col) = mat_argr(uint_row,uint_col);
            }
        }
    } else { // Case III
        std::cout << std::endl << "[concat()] To concatenate two matrices either columns, rows or both should have equal size.";
    }

    return mat_ret;
}

/**
 * @brief Flip left/right
 *
 * @param mat_arg
 * @return Returns flip (left/right) a matrix
 * @ingroup LALG
 */
template <class T> matrix <T> fliplr(const matrix <T> &mat_arg) {
    matrix <T> mat_ret(mat_arg.no_rows(), mat_arg.no_cols());

    if ( mat_arg.no_cols() == 1 || mat_arg.no_rows() == 1) {
        if (mat_arg.no_cols() == 1) std::cout << std::endl << "[fliplr()] Warning : flipud() was done.";
        for (unsigned int uint_i = 0; uint_i < mat_arg.size(); uint_i++) mat_ret(uint_i) = mat_arg(mat_arg.size() - uint_i - 1);
    } else if (mat_arg.no_cols() > 1 || mat_arg.no_rows() > 1) {
        for (unsigned int uint_row = 0; uint_row < mat_arg.no_rows(); uint_row++) {
            for (unsigned int uint_col = 0; uint_col < mat_arg.no_cols(); uint_col++) mat_ret(uint_row,uint_col) = mat_arg(uint_row,mat_arg.no_cols() - uint_col - 1);
        }
    }

    return mat_ret;
}

/**
 * @brief Flip up/down
 *
 * @param mat_arg
 * @return Returns flip (up/down) of the input matrix
 * @ingroup LALG
 */
template <class T> matrix <T> flipud(const matrix <T> &mat_arg) {
    matrix <T> mat_ret(mat_arg.no_rows(), mat_arg.no_cols());

    if ( mat_arg.no_cols() == 1 || mat_arg.no_rows() == 1) {
        if (mat_arg.no_rows() == 1) std::cout << std::endl << "[flipud()] Warning : fliplr() was done.";
        for (unsigned int uint_i = 0; uint_i < mat_arg.size(); uint_i++) mat_ret(uint_i) = mat_arg(mat_arg.size() - uint_i - 1);
    } else if (mat_arg.no_cols() > 1 || mat_arg.no_rows() > 1) {
        for (unsigned int uint_col = 0; uint_col < mat_arg.no_cols(); uint_col++) {
            for (unsigned int uint_row = 0; uint_row < mat_arg.no_rows(); uint_row++) mat_ret(uint_row,uint_col) = mat_arg(mat_arg.no_rows()  - uint_row - 1,uint_col);
        }
    }

    return mat_ret;
}


} // NAMESPACE SUSA

#endif //LINALG_H
