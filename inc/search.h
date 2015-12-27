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
 * @file search.h
 * @brief Search routines
 * This file containes special search routines.
 * @author Behrooz Kamary Aliabadi
 * @version 1.0.0
 *
 * @defgroup Search Search
 */

#ifndef SEARCH_H
#define SEARCH_H

#include "debug.h"

namespace susa {
/**
 * @brief Finds the <i>uint_num</i> least elements of the input vector.
 *
 * @param mat_arg Input vector
 * @param uint_num Number of elements to be find before the routine terminates
 * @return Returns a vector that contains indeces of elements in the input vector that have least values.
 *
 * @ingroup Search
*/
template <class T> matrix <unsigned int> select_least(const matrix <T> &mat_arg, unsigned int uint_num);

/**
 * @brief Finds the <i>uint_num</i> least elements of the input vector among a selected number of
 elements that their indeces are specified by <i>mat_limited_index</i>.
 *
 * @param mat_arg Input vector
 * @param uint_num Number of elements to be find before the routine terminates
 * @return Returns a vector that contains indeces of elements in the input vector that have least values.
 *
 * @ingroup Search
*/
template <class T> matrix <unsigned int> select_limited_least(const matrix <T> &mat_arg, const matrix <unsigned int> &mat_limited_index, unsigned int uint_num);

template <class T> matrix <unsigned int> select_most(const matrix <T> &mat_arg, unsigned int uint_num);


// Implementation
template <class T> matrix <unsigned int> select_least(const matrix <T> &mat_arg, unsigned int uint_num) {

    unsigned int uint_size = mat_arg.size();
    matrix <T> mat_ret = mat_arg;
    matrix <unsigned int> mat_index(uint_num,1);

    if (uint_size < uint_num) {
        std::cout << std::endl << "[select_least()] number of elements to be sorted is bigger than the matrix size.";
        return mat_index;
    }

    unsigned int uint_min_index;
    T T_min;

    for( unsigned int uint_i = 0; uint_i < uint_num; uint_i++) {
        uint_min_index = uint_i;
        T_min = mat_ret(uint_i);
        for (unsigned int uint_j = 0; uint_j < uint_size; uint_j++) {
            if (mat_ret(uint_j) < T_min) {
                uint_min_index = uint_j;
                T_min = mat_ret(uint_j);
            }
        }

        mat_index(uint_i) = uint_min_index;

        mat_ret(uint_min_index) = std::numeric_limits<T>::max();
    }

    return mat_index;
}

template <class T> matrix <unsigned int> select_limited_least(const matrix <T> &mat_arg, const matrix <unsigned int> &mat_limited_index, unsigned int uint_num) {

    unsigned int uint_size = mat_arg.size();
    matrix <T> mat_ret = mat_arg;
    matrix <unsigned int> mat_index(uint_num,1);

    if (uint_size < uint_num) {
        std::cout << std::endl << "[select_limited_least()] number of elements to be sorted is bigger than the matrix size.";
        return mat_index;
    }

    unsigned int uint_min_index;
    T T_min;

    for( unsigned int uint_i = 0; uint_i < uint_num; uint_i++) {
        uint_min_index = uint_i;
        T_min = mat_ret(uint_i);
        for (unsigned int uint_j = 0; uint_j < mat_limited_index.size(); uint_j++) {
            if (mat_ret(mat_limited_index(uint_j)) < T_min) {
                uint_min_index = mat_limited_index(uint_j);
                T_min = mat_ret(mat_limited_index(uint_j));
            }
        }

        mat_index(uint_i) = uint_min_index;

        mat_ret(uint_min_index) = std::numeric_limits<T>::max();
    }

    return mat_index;
}

template <class T> matrix <unsigned int> select_most(const matrix <T> &mat_arg, unsigned int uint_num) {

    unsigned int uint_size = mat_arg.size();
    matrix <T> mat_ret = mat_arg;
    matrix <unsigned int> mat_index(uint_num,1);

    if (uint_size < uint_num) {
        std::cout << std::endl << "[select_most()] number of elements to be sorted is bigger than the matrix size.";
        return mat_index;
    }

    unsigned int uint_max_index;
    T T_max;

    for( unsigned int uint_i = 0; uint_i < uint_num; uint_i++) {
        uint_max_index = uint_i;
        T_max = mat_ret(uint_i);
        for (unsigned int uint_j = 0; uint_j < uint_size; uint_j++) {
            if (mat_ret(uint_j) > T_max) {
                uint_max_index = uint_j;
                T_max = mat_ret(uint_j);
            }
        }

        mat_index(uint_i) = uint_max_index;

        mat_ret(uint_max_index) = std::numeric_limits<T>::min();
    }

    return mat_index;
}
} // NAMESPACE SUSA
#endif // SEARCH_H
