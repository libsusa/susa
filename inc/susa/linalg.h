/*
 * This file is part of Susa.
 *
 * Susa is free software: you can redistribute it and/or modify
 * it under the terms of the Lesser GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * at your option) any later version.
 *
 * Susa is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * Lesser GNU General Public License for more details.
 *
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

#ifndef SUSA_LINALG_H
#define SUSA_LINALG_H

namespace susa
{


/**
 * @brief Matrix multiplication
 *
 * @param mat_argl the left-hand-side matrix
 * @param mat_argr the right-hand-side matrix
 * @return multiplication of the two input matrices
 * @ingroup LALG
 */
template <class T> matrix <T> matmul(const matrix <T> &mat_argl, const matrix <T> &mat_argr)
{
    // To make it faster this method has been declared as 'friend' of
    // 'matrix' template class so 'matmul' apears two times in 'matrix.h'
    // and two times in 'linalg.h'
    size_t sizet_rows = mat_argl.no_rows();
    size_t sizet_cols = mat_argr.no_cols();

    matrix <T> mat_ret;

    if (mat_argl._matrix == NULL || mat_argr._matrix == NULL)
    {
        SUSA_ABORT("one or both matrices are empty.");
        return mat_ret;
    }

    if (mat_argl.no_cols() == mat_argr.no_rows())
    {
        mat_ret = matrix <T> (sizet_rows, sizet_cols, 0);
        for (size_t sizet_row = 0; sizet_row < sizet_rows; sizet_row++)
        {
            for (size_t sizet_col = 0; sizet_col < sizet_cols; sizet_col++)
            {
                for (size_t sizet_index = 0; sizet_index < mat_argl.no_cols(); sizet_index++)
                {
                    mat_ret._matrix[sizet_row + sizet_col * mat_ret.sizet_rows] +=
                            mat_argl._matrix[sizet_row + sizet_index * mat_argl.sizet_rows]
                            * mat_argr._matrix[sizet_index + sizet_col * mat_argr.sizet_rows];
                }
            }
        }
    }
    else
    {
        SUSA_ABORT("the matrices' dimensions mismatch.");
    }

    return mat_ret;
}

/**
 * @brief Transpose operator
 *
 * @param mat_arg the input matrix
 * @return returns transpose of the input matrices
 * @ingroup LALG
 */
template <class T> matrix <T> transpose(const matrix <T> &mat_arg)
{
    size_t sizet_rows = mat_arg.no_rows();
    size_t sizet_cols = mat_arg.no_cols();
    matrix <T> mat_ret;

    SUSA_ASSERT(mat_arg._matrix != NULL);

    if (mat_arg._matrix == NULL)
    {
        return mat_ret;
    }

    mat_ret = matrix <T> (sizet_cols, sizet_rows);

    for (size_t sizet_col = 0; sizet_col < sizet_cols; sizet_col++)
    {
        for (size_t sizet_row = 0; sizet_row < sizet_rows; sizet_row++)
        {
            mat_ret._matrix[sizet_col + sizet_row * mat_ret.sizet_rows] =
            mat_arg._matrix[sizet_row + sizet_col * mat_arg.sizet_rows];
        }
    }

    return mat_ret;
}


/**
 * @brief Determinant
 *
 * @param mat_arg
 * @return returns determinant of the input matrices
 * @ingroup LALG
 */
template <class T> double det(const matrix <T> &mat_arg)
{
    double dbl_ret = 0.0;

    size_t sizet_rows = mat_arg.no_rows();
    size_t sizet_cols = mat_arg.no_cols();

    if (sizet_rows == sizet_cols)
    {
        if (sizet_rows == 1)
        { // 1 x 1
            dbl_ret = mat_arg(0);
        }
        else if (sizet_rows == 2)
        { // 2 x 2
            dbl_ret = mat_arg(0, 0) * mat_arg(1, 1) - mat_arg(1, 0) * mat_arg(0, 1);
        }
        else
        { // 3 x 3 and more
            for (size_t  sizet_col = 0; sizet_col < sizet_cols; sizet_col++)
            {
                matrix <T> mat_minor = mat_arg.shrink(0, sizet_col);
                if (mat_minor.size() == mat_arg.size())
                {
                    // If for any reason something goes wrong in the minor() it returns the matrix without any action.
                    // Therefore size of the input argument and the returned argument would the same. This case cause infinite
                    // iterations. To stop such an unwanted loop this if() returns from the method.
                    // THIS METHOD IS NOT STABLE BUT WORKS.
                    SUSA_LOG_ERR("falling into an infinite loop avoided.");
                    return dbl_ret;
                }
                dbl_ret +=  (-2 * (int)(sizet_col % 2) + 1) * mat_arg(0, sizet_col) * det(mat_minor);
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
 * @return returns concatenation of the two input matrices
 * @ingroup LALG
 */
template <class T> matrix <T> concat(const matrix <T> &mat_argl, const matrix <T> &mat_argr)
{
    // should be changed 20101224
    matrix <T> mat_ret;

    if ((mat_argl.no_cols() == mat_argr.no_cols() && mat_argl.no_rows() != mat_argr.no_rows()) ||
            (mat_argl.no_cols() == mat_argr.no_cols() && mat_argl.no_rows() == mat_argr.no_rows()))
    { // Case I

        mat_ret = matrix <T> (mat_argl.no_rows() + mat_argr.no_rows(), mat_argl.no_cols());

        for (size_t sizet_row = 0; sizet_row < mat_argl.no_rows(); sizet_row++)
        {
            for (size_t sizet_col = 0; sizet_col < mat_argl.no_cols(); sizet_col++)
            {
                mat_ret(sizet_row,sizet_col) = mat_argl(sizet_row,sizet_col);
            }
        }

        for (size_t sizet_row = 0; sizet_row < mat_argr.no_rows(); sizet_row++)
        {
            for (size_t sizet_col = 0; sizet_col < mat_argr.no_cols(); sizet_col++)
            {
                mat_ret(mat_argl.no_rows() + sizet_row,sizet_col) = mat_argr(sizet_row,sizet_col);
            }
        }
    }
    else if (mat_argl.no_cols() != mat_argr.no_cols() && mat_argl.no_rows() == mat_argr.no_rows())
    { // Case II
        mat_ret = matrix <T> (mat_argl.no_rows(), mat_argl.no_cols() + mat_argr.no_cols());

        for (size_t sizet_row = 0; sizet_row < mat_argl.no_rows(); sizet_row++)
        {
            for (size_t sizet_col = 0; sizet_col < mat_argl.no_cols(); sizet_col++)
            {
                mat_ret(sizet_row,sizet_col) = mat_argl(sizet_row,sizet_col);
            }
        }

        for (size_t sizet_row = 0; sizet_row < mat_argr.no_rows(); sizet_row++)
        {
            for (size_t sizet_col = 0; sizet_col < mat_argr.no_cols(); sizet_col++)
            {
                mat_ret(sizet_row,mat_argr.no_cols() + sizet_col) = mat_argr(sizet_row,sizet_col);
            }
        }
    }
    else
    { // Case III
        SUSA_LOG_ERR("to concatenate two matrices either columns or rows, or both should have equal sizes.");
    }

    return mat_ret;
}

/**
 * @brief Flip left/right
 *
 * @param mat_arg the input matrix
 * @return returns flip (left/right) a matrix
 * @ingroup LALG
 */
template <class T> matrix <T> fliplr(const matrix <T> &mat_arg)
{
    matrix <T> mat_ret(mat_arg.no_rows(), mat_arg.no_cols());

    if ( mat_arg.no_cols() == 1 || mat_arg.no_rows() == 1)
    {
        if (mat_arg.no_cols() == 1) SUSA_LOG_INF("fliplr() has been replaced by flipud().");

        for (size_t sizet_i = 0; sizet_i < mat_arg.size(); sizet_i++)
        {
          mat_ret(sizet_i) = mat_arg(mat_arg.size() - sizet_i - 1);
        }
    }
    else if (mat_arg.no_cols() > 1 || mat_arg.no_rows() > 1)
    {
        for (size_t sizet_row = 0; sizet_row < mat_arg.no_rows(); sizet_row++)
        {
            for (size_t sizet_col = 0; sizet_col < mat_arg.no_cols(); sizet_col++)
            {
              mat_ret(sizet_row,sizet_col) = mat_arg(sizet_row,mat_arg.no_cols() - sizet_col - 1);
            }
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
template <class T> matrix <T> flipud(const matrix <T> &mat_arg)
{
    matrix <T> mat_ret(mat_arg.no_rows(), mat_arg.no_cols());

    if ( mat_arg.no_cols() == 1 || mat_arg.no_rows() == 1)
    {

        if (mat_arg.no_rows() == 1) SUSA_LOG_INF("flipud() has been replaced by fliplr().");

        for (size_t sizet_i = 0; sizet_i < mat_arg.size(); sizet_i++)
        {
          mat_ret(sizet_i) = mat_arg(mat_arg.size() - sizet_i - 1);
        }
    }
    else if (mat_arg.no_cols() > 1 || mat_arg.no_rows() > 1)
    {
        for (size_t sizet_col = 0; sizet_col < mat_arg.no_cols(); sizet_col++)
        {
            for (size_t sizet_row = 0; sizet_row < mat_arg.no_rows(); sizet_row++)
            {
              mat_ret(sizet_row,sizet_col) = mat_arg(mat_arg.no_rows()  - sizet_row - 1, sizet_col);
            }
        }
    }

    return mat_ret;
}

/**
 * @brief generates a square identity matrix
 *
 * @param sizet_n matrix dimention
 * @return square identity matrix
 * @ingroup LALG
 */
template <class T> matrix<T> eye(size_t sizet_n)
{
    matrix<T> mat_ret(sizet_n, sizet_n);

    for (size_t sizet_col = 0; sizet_col < sizet_n; sizet_col++)
    {
        for (size_t sizet_row = 0; sizet_row < sizet_n; sizet_row++)
        {
            if (sizet_col == sizet_row)
            {
                mat_ret(sizet_row, sizet_col) = 1;
            }
            else
            {
                mat_ret(sizet_row, sizet_col) = 0;
            }
        }
    }

    return mat_ret;
}

}      // NAMESPACE SUSA
#endif //SUSA_LINALG_H
