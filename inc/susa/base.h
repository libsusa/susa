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
 * along with Susa. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file base.h
 * @brief Basic operations on STL and Susa types (declaration and definition).
 *
 * @author Behrooz Kamary
 *
 * @defgroup Math Basic Mathematics
 */

#ifndef SUSA_BASE_H
#define SUSA_BASE_H

#include <type_traits>

namespace susa {

/**
 * @brief Minimum value indices
 *
 * @param mat_arg the input matrix
 * @return returns minimum value indices in the input matrix
 * @ingroup Math
 */
template <typename T, typename Allocator> matrix <size_t> min(const matrix <T, Allocator> &mat_arg);

/**
 * @brief Maximum
 *
 * @param mat_arg
 * @return Returns maximum value in the input matrix
 * @ingroup Math
 */
template <typename T, typename Allocator> matrix <size_t> max(const matrix <T, Allocator> &mat_arg);

/**
 * @brief Differential
 *
 * @param vec_arg STL vector
 * @return returns discrete differential of the input vector
 * @ingroup Math
 */
template <typename T, typename Allocator> std::vector <T> diff(const std::vector <T> &vec_arg);

/**
 * @brief Sum
 *
 * @param mat_arg input matrix
 * @return Returns sum of values in the input matrix
 * @ingroup Math
 */
template <typename T, typename Allocator> matrix <T, Allocator> sum(const matrix <T, Allocator> &mat_arg);

/**
 * @brief Sum
 *
 * @param vec_arg input STL vector
 * @return Returns sum of values in the input vector
 * @ingroup Math
 */
template <typename T> T sum(const std::vector <T> &vec_arg);

/**
 * @brief Mean
 *
 * @param mat_arg input matrix
 * @return Returns mean of values in the input vector
 * @ingroup Math
 */
template <typename T, typename Allocator> matrix <T, Allocator> mean(const matrix <T, Allocator> &mat_arg);

/**
 * @brief Norm
 *
 * @param mat_arg input matrix
 * @return compute 2-norm or square norm of the input matrix
 * @ingroup Math
 */
template <typename T, typename Allocator>
matrix <T, Allocator> norm(const matrix <T, Allocator> &mat_arg);

/**
 * @brief Real
 *
 * @param mat_arg input complex matrix
 * @return Returns real part of the complex input vector
 * @ingroup Math
 */
template <typename T, template <typename> typename Allocator>
matrix <T, Allocator<T>> real(const matrix <std::complex<T>, Allocator<std::complex<T>>> &mat_arg);

/**
 * @brief Imaginary
 *
 * @param mat_arg input complex matrix
 * @return Returns imaginary part of the complex input vector
 * @ingroup Math
 */
template <typename T, template <typename> typename Allocator>
matrix <T, Allocator<T>> imag(const matrix <std::complex<T>, Allocator<std::complex<T>>> &mat_arg);

/**
 * @brief Magnitude
 *
 * @param mat_arg input complex matrix
 * @return compute the absolute square of the complex input matrix elements
 * @ingroup Math
 */
template <typename T, template <typename> typename Allocator>
matrix <T, Allocator<T>> mag(const matrix <std::complex<T>, Allocator<std::complex<T>>> &mat_arg);

/**
 * @brief Absolute
 *
 * @param mat_arg input matrix
 * @return compute the absolute value i.e. the magnitude of the input matrix elements
 * @ingroup Math
 */
template <typename T, typename Allocator>
susa::matrix <T, Allocator> abs(const susa::matrix <T, Allocator> &mat_arg);

/**
 * @brief Absolute
 *
 * @param vec_arg input STL vector
 * @return compute the absolute value of the input vector
 * @ingroup Math
 */
template <typename T, typename Allocator>
std::vector <T> abs(const std::vector <T> &vec_arg);

/**
 * @brief Conjugate
 *
 * @param mat_arg input complex matrix
 * @return Returns conjugate of the complex input vector
 * @ingroup Math
 */
template <typename T, typename Allocator>
matrix <std::complex <T>, Allocator> conj(const matrix <std::complex <T>, Allocator> &mat_arg);

/**
 * @brief Sign
 *
 * @param mat_arg input matrix
 * @return Returns sign of the input
 * @ingroup Math
 */
template <typename T, typename Allocator> matrix <char> sign(const matrix <T, Allocator> &mat_arg);

/**
 * @brief Sign
 *
 * @param T_arg input matrix
 * @return Returns sign of the input
 * @ingroup Math
 */
template <typename T, typename Allocator> char sign(T T_arg);

/**
 * @brief Log
 *
 * @param mat_arg input matrix
 * @return Returns natural logarithm of input
 * @ingroup Math
 */
template <typename T, typename Allocator> matrix <T, Allocator> log(const matrix <T, Allocator> &mat_arg);

/**
 * @brief Pow
 *
 * @param uint_b base argument
 * @param uint_p power argument
 * @ingroup Math
 */
unsigned int pow(unsigned int uint_b, unsigned int uint_p);

/**
 * @brief Pow
 *
 * @param uint_b base argument
 * @param uint_p power argument
 * @ingroup Math
 */
int pow(int int_b, unsigned int uint_p);

/**
 * @brief log2
 *
 * @param uint_arg unsigned input argument
 * @ingroup Math
 */
template <typename T, typename std::enable_if_t<std::is_unsigned<T>::value,T>* = nullptr>
T log2(T int_arg);

/**
 * @brief Normal Cumulative Distribution Function (CDF)
 *
 * @param x input argument
 * @ingroup Math
 */
double normcdf(const double x);

/**
 * @brief qfunc
 *
 * @param x input argument
 * @ingroup Math
 */
double qfunc(const double x);

/**
 * @brief Modular operation
 *
 * @return lint_a "mod" lint_mod
 * @ingroup Math
 */
long int mod(long int lint_a, long int lint_mod);

/**
 * @brief Round
 *
 * @param  dbl_arg The input floating point number
 * @param  int_decimal Number of decimals
 * @return rounds the floating number dbl_arg up to int_decimal numbers
 *
 * @ingroup Math
 */
double round(const double dbl_arg, int int_decimal = -1);


/**
 * @brief Round
 *
 * @param  mat_arg the input susa::matrix<double>
 * @param  int_decimal number of decimals
 * @return rounds the matrix mat_arg elements up to int_decimal numbers
 *
 * @ingroup Math
 */
matrix <double> round(const matrix <double> &mat_arg, int int_decimal = -1);

// Implementations


template <typename T, typename std::enable_if_t<std::is_unsigned<T>::value,T>*>
T log2(T int_arg)
{
    T ret(0);

    while ((int_arg >> ret) != 0)
    {
        ret++;
    }

    return (ret - 1);
}

template <typename T> std::vector <T> diff(std::vector <T> &vec_arg)
{
    std::vector <T> vec_diff(vec_arg.size() - 1, T(0));

    for (size_t i = 0; i < vec_arg.size() - 1; i++)
    {
        vec_diff[i] = vec_arg[i + 1] - vec_arg[i];
    }

    return (vec_diff);
}

template <typename T> T sum(const std::vector<T> &vec_arg)
{
    T       ret(0);
    size_t  size = vec_arg.size();
    for (size_t i = 0; i < size; i++)
    {
        ret += vec_arg[i];
    }

    return ret;
}

template <typename T, typename Allocator>
matrix <T, Allocator> sum(const matrix <T, Allocator> &mat_arg)
{
    matrix <T, Allocator> mat_ret;

    if (mat_arg.no_rows() == 1 || mat_arg.no_cols() == 1)
    {
        mat_ret = matrix <T, Allocator> (1,1, T(0));
        for (size_t uint_i = 0; uint_i < mat_arg.size(); uint_i++)
            mat_ret(0) += mat_arg(uint_i);
    }
    else if (mat_arg.no_rows() > 1 && mat_arg.no_cols() > 1)
    {
        mat_ret = matrix <T, Allocator> (1, mat_arg.no_cols(), T(0));
        for (size_t uint_col = 0; uint_col < mat_arg.no_cols(); uint_col++)
        {
            for (size_t uint_row = 0; uint_row < mat_arg.no_rows(); uint_row++)
                mat_ret(uint_col) += mat_arg(uint_row, uint_col);
        }
    }

    return mat_ret;
}

template <typename T, typename Allocator>
matrix <T, Allocator> mean(const matrix <T, Allocator> &mat_arg)
{
    matrix <T, Allocator>   mat_ret;
    T                       T_ret(0);

    if (mat_arg.no_rows() == 1 || mat_arg.no_cols() == 1)
    {
        for (size_t uint_i = 0; uint_i < mat_arg.size(); uint_i++)
        {
            T_ret += mat_arg(uint_i);
        }
        T_ret /= mat_arg.size();
        mat_ret = matrix <T, Allocator> (1, 1, T_ret);
    }
    else if (mat_arg.no_rows() > 1 && mat_arg.no_cols() > 1)
    {
        size_t no_rows = mat_arg.no_rows();
        mat_ret = matrix <T, Allocator> (1, mat_arg.no_cols());
        for (size_t uint_col = 0; uint_col < mat_arg.no_cols(); uint_col++)
        {
            for (size_t uint_row = 0; uint_row < no_rows; uint_row++)
            {
                T_ret += mat_arg(uint_row, uint_col);
            }

            T_ret               /= no_rows;
            mat_ret(uint_col)   = T_ret;
            T_ret               = T(0);
        }
    }

    return mat_ret;
}

template <typename T>
std::vector <T> abs(const std::vector <T> &vec_arg)
{
    std::vector <T> vec_abs(vec_arg.size(),0);
    for (size_t i = 0; i < vec_arg.size(); i++) vec_abs[i] = abs(vec_arg[i]);
    return vec_abs;
}

template <typename T, typename Allocator>
susa::matrix <T, Allocator> abs(const susa::matrix <T, Allocator> &mat_arg)
{
    susa::matrix <T, Allocator>     mat_ret(mat_arg.shape());
    size_t size = mat_arg.size();
    for (size_t uint_i = 1; uint_i < size; uint_i++)
    {
        mat_ret(uint_i) = abs(mat_arg(uint_i));
    }

    return mat_ret;
}

template <typename T, typename Allocator>
matrix <size_t> min(const matrix <T, Allocator> &mat_arg)
{
    matrix <size_t>     mat_ret;
    T                   T_min;

    if (mat_arg.is_vector())
    {
        mat_ret = matrix <size_t> (1, 1, 0);
        T_min = mat_arg(0);
        size_t size = mat_arg.size();
        for (size_t uint_i = 1; uint_i < size; uint_i++)
        {
            if (mat_arg(uint_i) < T_min)
            {
                T_min = mat_arg(uint_i);
                mat_ret(0) = uint_i;
            }
        }
    }
    else
    {
        mat_ret = matrix <size_t> (1, mat_arg.no_cols(), 0);
        for (size_t uint_col = 0; uint_col < mat_arg.no_cols(); uint_col++)
        {
            T_min = mat_arg(0, uint_col);
            for (size_t uint_row = 1; uint_row < mat_arg.no_rows(); uint_row++)
            {
                if (mat_arg(uint_row, uint_col) < T_min)
                {
                    T_min = mat_arg(uint_row, uint_col);
                    mat_ret(uint_col) = uint_row;
                }
            }
        }
    }

    return mat_ret;
}

template <typename T, typename Allocator> matrix <size_t> max(const matrix <T, Allocator> &mat_arg)
{
    matrix <size_t> mat_ret;
    T               T_max;

    if (mat_arg.is_vector())
    {
        mat_ret = matrix <size_t> (1, 1, 0);
        T_max = mat_arg(0);
        for (size_t uint_i = 1; uint_i < mat_arg.size(); uint_i++)
        {
            if (mat_arg(uint_i) > T_max)
            {
                T_max = mat_arg(uint_i);
                mat_ret(0) = uint_i;
            }
        }
    }
    else if (mat_arg.no_rows() > 1 && mat_arg.no_cols() > 1)
    {
        mat_ret = matrix <size_t> (1, mat_arg.no_cols(), 0);
        for (size_t uint_col = 0; uint_col < mat_arg.no_cols(); uint_col++)
        {
            T_max = mat_arg(0, uint_col);
            for (size_t uint_row = 1; uint_row < mat_arg.no_rows(); uint_row++)
            {
                if (mat_arg(uint_row, uint_col) > T_max)
                {
                    T_max = mat_arg(uint_row,uint_col);
                    mat_ret(uint_col) = uint_row;
                }
            }
        }
    }

    return mat_ret;
}

template <typename T, typename Allocator>
matrix <std::complex <T>, Allocator> conj(const matrix <std::complex <T>, Allocator> &mat_arg)
{
    matrix <std::complex <T>, Allocator> mat_ret(mat_arg.shape());

    for (size_t uint_i = 0; uint_i < mat_arg.size(); uint_i++)
    {
        mat_ret(uint_i) = std::conj(mat_arg(uint_i));
    }
    return mat_ret;
}

template <typename T, typename Allocator>
matrix <T, Allocator> norm(const matrix <T, Allocator> &mat_arg)
{
    matrix <T, Allocator> mat_ret;

    if (mat_arg.is_vector())
    {
        mat_ret = matrix <T, Allocator> (1,1, T(0));
        for (size_t uint_i = 0; uint_i < mat_arg.size(); uint_i++)
        {
            mat_ret(0) += mat_arg(uint_i) * mat_arg(uint_i);
        }
        mat_ret(0) = std::sqrt(mat_ret(0));
    }
    else
    {
        mat_ret = matrix <T, Allocator> (1, mat_arg.no_cols(), T(0));
        for (size_t uint_col = 0; uint_col < mat_arg.no_cols(); uint_col++)
        {
            for (size_t uint_row = 0; uint_row < mat_arg.no_rows(); uint_row++)
            {
                mat_ret(uint_col) += mat_arg(uint_row,uint_col) * mat_arg(uint_row,uint_col);
            }
            mat_ret(uint_col) = std::sqrt(mat_ret(uint_col));
        }
    }
    return mat_ret;
}

template <typename T, template <typename> typename Allocator>
matrix <T, Allocator<T>> real(const matrix <std::complex<T>, Allocator<std::complex<T>>> &mat_arg)
{
    matrix <T, Allocator<T>> mat_ret(mat_arg.shape());

    for (size_t uint_i = 0; uint_i < mat_arg.size(); uint_i++)
    {
        mat_ret(uint_i) = mat_arg(uint_i).real();
    }

    return mat_ret;
}

template <typename T, template <typename> typename Allocator>
matrix <T, Allocator<T>> imag(const matrix <std::complex<T>, Allocator<std::complex<T>>> &mat_arg)
{
    matrix <T, Allocator<T>> mat_ret(mat_arg.shape());

    for (size_t uint_i = 0; uint_i < mat_arg.size(); uint_i++)
    {
        mat_ret(uint_i) = mat_arg(uint_i).imag();
    }

    return mat_ret;
}

template <typename T, template <typename> typename Allocator>
matrix <T, Allocator<T>> mag(const matrix <std::complex<T>, Allocator<std::complex<T>>> &mat_arg)
{
    matrix <T, Allocator<T>> mat_ret(mat_arg.shape());

    for (size_t uint_i = 0; uint_i < mat_arg.size(); uint_i++)
    {
        //mat_ret(uint_i) = (mat_arg(uint_i) * std::conj(mat_arg(uint_i))).real();
        mat_ret(uint_i) = std::norm(mat_arg(uint_i));
    }
    return mat_ret;
}

template <typename T, typename Allocator> matrix <char> sign(const matrix <T, Allocator> &mat_arg)
{
    matrix <T, Allocator> mat_ret(mat_arg.shape());

    for (size_t uint_i = 0; uint_i < mat_arg.size(); uint_i++)
    {
        if (mat_arg(uint_i) != 0) mat_ret(uint_i) = mat_arg(uint_i) > 0 ? 1 : -1;
        else mat_ret(uint_i) = 0;
    }

    return mat_ret;
}

template <typename T> char sign(T T_arg)
{
    if (T_arg > 0) return 1;
    else if (T_arg < 0) return -1;
    return 0;
}

template <typename T, typename Allocator> matrix <T, Allocator> log(const matrix <T, Allocator> &mat_arg)
{
    matrix <T, Allocator> mat_ret(mat_arg.shape());

    for (size_t uint_i = 0; uint_i < mat_arg.size(); uint_i++)
    {
        mat_ret(uint_i) = std::log(mat_arg(uint_i));
    }

    return mat_ret;
}


}       // NAMESPACE SUSA
#endif  // SUSA_BASE_H
