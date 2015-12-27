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
 * @file signal.h
 * @brief Signal processing operations on Susa types
 * @author Behrooz, Kamary Aliabadi
 * @version 1.0.0
 *
 * @defgroup Signal Signal Processing
 *
 * @todo <i>filter</i> and <i>conv</i> do not support all input matrices' sizes. it looks complicated to solve.
 */

#ifndef SIGNAL_H
#define SIGNAL_H

namespace susa {

/**
 * @brief Upsample
 *
 * @param mat_arg Input matrix
 * @param uint_u Upsampling rate
 * @return Upsampled matrix
 *
 * @ingroup Signal
 */
template <class T> matrix <T> upsample(const matrix <T> &mat_arg, unsigned int uint_u);

/**
 * @brief Downsample
 *
 * @param mat_arg Input matrix
 * @param uint_d Downsampling rate
 * @return Downsampled matrix
 *
 * @ingroup Signal
 */
template <class T> matrix <T> downsample(const matrix <T> &mat_arg, unsigned int uint_d);

/**
 * @brief Filter
 *
 * y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb) - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
 *
 * @param mat_arg_b Moving Average (MA) coefficients
 * @param mat_arg_a Autoregressive (AR) coefficients.
 * @param mat_arg_x Input data to be filtered.
 * @param uint_length Additional tail of zeros that can be appended to the input.
 * @return Returns filtered matrix
 *
 * @ingroup Signal
 */
template <class T, class TT> matrix <T> filter(const matrix <TT> &mat_arg_b, const matrix <TT> &mat_arg_a, const matrix <T> &mat_arg_x, unsigned int uint_length = 0);

/**
 * @brief Convolution
 *
 * @param mat_arg_a
 * @param mat_arg_b
 * @return Returns convolution of two input vectors
 *
 * @ingroup Signal
 */
template <class T> matrix <T> conv(const matrix <T> &mat_arg_a, const matrix <T> &mat_arg_b);

/**
 * @brief Constructs a convolution matrix from the impulse response <i>mat_arg</i>
 *
 * @param mat_arg input vector (impulse response)
 * @param uint_size
 * @return The convolution matrix
 *
 * @ingroup Signal
 */
template <class T> matrix <T> convmtx(const matrix <T> &mat_arg, unsigned int uint_size);

/**
 * @brief Constructs a toeplitz matrix from the input vector <i>mat_arg</i>
 *
 * @param mat_arg input vector
 * @return The toeplitz matrix
 *
 * @ingroup Signal
 */
template <class T> matrix <T> toeplitz(const matrix <T> &mat_arg);

// Implementation

// UPSAMPLE

template <class T> matrix <T> upsample(const matrix <T> &mat_arg, unsigned int uint_u) {
    matrix <T> mat_tmp(uint_u * mat_arg.no_rows(),mat_arg.no_cols());

    for (unsigned int uint_col = 0; uint_col < mat_arg.no_cols(); uint_col++) {
        for (unsigned int uint_row = 0; uint_row < mat_arg.no_rows(); uint_row++) {
            mat_tmp(uint_row * uint_u,uint_col) = mat_arg(uint_row, uint_col);
        }
    }
    return mat_tmp;
}

// DOWNSAMPLE

template <class T> matrix <T> downsample(const matrix <T> &mat_arg, unsigned int uint_d) {
    matrix <T> mat_tmp(mat_arg.no_rows()/uint_d,mat_arg.no_cols());

    for (unsigned int uint_col = 0; uint_col < mat_arg.no_cols(); uint_col++) {
        for (unsigned int uint_row = 0; uint_row < mat_arg.no_rows()/uint_d; uint_row++) {
            mat_tmp(uint_row,uint_col) = mat_arg(uint_row * uint_d,uint_col);
        }
    }
    return mat_tmp;
}

// FILTER

template <class T, class TT> matrix <T> filter( const matrix <TT> &mat_arg_b, const matrix <TT> &mat_arg_a, const matrix <T> &mat_arg_x, unsigned int uint_length) {

    unsigned int uint_response_length = 0;
    matrix <T> mat_y;

    int int_b_size = (int)mat_arg_b.size();
    int int_a_size = (int)mat_arg_a.size();

    T x;
    T y;

    if (mat_arg_x.no_cols() > 1 && mat_arg_x.no_rows() > 1) {
        uint_response_length = mat_arg_x.no_rows() + uint_length;
        mat_y = matrix <T> (uint_response_length,mat_arg_x.no_cols());
        for (unsigned int uint_col = 0; uint_col < mat_arg_x.no_cols(); uint_col++) {
            for (int i=0 ; i< uint_response_length; i++) {
                for (int k=0; k < int_b_size; k++) {
                    x = ((i - k) < (int)mat_arg_x.no_rows() && ((i - k) > -1)) ? mat_arg_x(i - k,uint_col) : 0;
                    mat_y(i,uint_col) += x * mat_arg_b(k);
                }
                for (int k=1; k < int_a_size; k++) {
                    y = ((i - k) < (int)mat_y.no_rows() && ((i - k) > -1)) ? mat_y(i - k,uint_col) : 0;
                    mat_y(i,uint_col) -= y * mat_arg_a(k);
                }
            }
        }
    } else if (mat_arg_x.no_cols() == 1 && mat_arg_x.no_rows() > 1) {
        uint_response_length = mat_arg_x.size() + uint_length;
        mat_y = matrix <T> (uint_response_length,1); // For column vector
        for (int i = 0; i< uint_response_length; i++) {
            for (int k = 0; k < int_b_size; k++) {
                x = ((i - k) < (int)mat_arg_x.size() && ((i - k) > -1)) ? mat_arg_x(i - k) : 0;
                mat_y(i) += x * mat_arg_b(k);
            }
            for (int k = 1; k < int_a_size; k++) {
                y = ((i - k) < (int)mat_y.size() && ((i - k) > -1)) ? mat_y(i - k) : 0;
                mat_y(i) -= y * mat_arg_a(k);
            }
        }
    } else if (mat_arg_x.no_cols() > 1 && mat_arg_x.no_rows() == 1) {
        uint_response_length = mat_arg_x.size() + uint_length;
        mat_y = matrix <T> (1,uint_response_length); // For row vector
        for (int i=0 ; i< uint_response_length; i++) {
            for (int k=0; k < int_b_size; k++) {
                x = ((i - k) < mat_arg_x.size() && ((i - k) > -1)) ? mat_arg_x(i - k) : 0;
                mat_y(i) += x * mat_arg_b(k);
            }
            for (int k=1; k < int_a_size; k++) {
                y = ((i - k) < mat_y.size() && ((i - k) > -1)) ? mat_y(i - k) : 0;
                mat_y(i) -= y * mat_arg_a(k);
            }
        }
    }
    return mat_y;
}


template <class T> matrix <T> conv(const matrix <T> &mat_arg_a, const matrix <T> &mat_arg_b) {

    matrix <T> mat_ret;
    unsigned int uint_length;

    if ((mat_arg_b.no_rows() > 1 && mat_arg_b.no_cols() == 1) || (mat_arg_b.no_cols() > 1 && mat_arg_b.no_rows() == 1)) { // mat_arg_b is a vector

        uint_length = mat_arg_b.size() - 1;
        mat_ret = filter(mat_arg_b, matrix <T> (1,1,1),mat_arg_a, uint_length);

    } else if ((mat_arg_a.no_rows() > 1 && mat_arg_a.no_cols() == 1) || (mat_arg_a.no_cols() > 1 && mat_arg_a.no_rows() == 1)) { // mat_arg_a is a vector

        uint_length = mat_arg_a.size() - 1;
        mat_ret = filter(mat_arg_a, matrix <T> (1,1,1),mat_arg_b, uint_length);

    } else if (mat_arg_b.no_cols() == 1 && mat_arg_b.no_rows() == 1) {
        mat_ret = mat_arg_b(0) * mat_arg_a;
    } else if (mat_arg_a.no_cols() == 1 && mat_arg_a.no_rows() == 1) {
        mat_ret = mat_arg_a(0) * mat_arg_b;
    } else if ((mat_arg_b.no_cols() > 1 && mat_arg_b.no_rows() > 1) && (mat_arg_a.no_cols() > 1 && mat_arg_a.no_rows() > 1) && (mat_arg_a.no_cols() == mat_arg_b.no_cols())) {
        for (unsigned int uint_col = 0; uint_col < mat_arg_a.no_cols(); uint_col++) {
            uint_length = mat_arg_a.no_rows() - 1;
            // I need a [set_row/set_col] to complete this part 2009-02-22
        }
    } else {
        std::cout << std::endl << "[conv()] not supported codition.";
    }

    return mat_ret;
}

template <class T> matrix <T> convmtx(const matrix <T> &mat_arg, unsigned int uint_n) {
    unsigned int uint_m = mat_arg.size();
    matrix <T> mat_ret;
    if (mat_arg.no_cols() == 1 && mat_arg.no_rows() >= 1) { // It is a column vector
        mat_ret = matrix <T> (uint_m + uint_n - 1, uint_n);
        for (unsigned int uint_col = 0; uint_col < uint_n; uint_col++) {
            for (unsigned int uint_row = 0; uint_row < uint_col; uint_row++) {
                mat_ret(uint_row,uint_col) = 0;
            }

            for (unsigned int uint_row = uint_col; uint_row < uint_m + uint_col; uint_row++) {
                mat_ret(uint_row,uint_col) = mat_arg(uint_row - uint_col);
            }

            for (unsigned int uint_row = uint_m + uint_col ; uint_row < uint_m + uint_n - 1; uint_row++) {
                mat_ret(uint_row,uint_col) = 0;
            }
        }
    } else if (mat_arg.no_rows() == 1 && mat_arg.no_cols() >= 1) { // It is a row vector
        mat_ret = matrix <T> (uint_n, uint_m + uint_n - 1);
        for (unsigned int uint_row = 0; uint_row < uint_n; uint_row++) {
            for (unsigned int uint_col = 0; uint_col < uint_row; uint_col++) {
                mat_ret(uint_row,uint_col) = 0;
            }

            for (unsigned int uint_col = uint_row; uint_col < uint_m + uint_row; uint_col++) {
                mat_ret(uint_row,uint_col) = mat_arg(uint_col - uint_row);
            }

            for (unsigned int uint_col = uint_m + uint_row ; uint_col < uint_m + uint_n - 1; uint_col++) {
                mat_ret(uint_row,uint_col) = 0;
            }
        }
    } else { // It is not a vector !
        std::cout << std::endl << "[convmtx()] The input argument is not a vector";
    }
    return mat_ret;
}

template <class T> matrix <T> toeplitz(const matrix <T> &mat_arg) {
    unsigned int uint_size = mat_arg.size();
    matrix <T> mat_ret(uint_size,uint_size);
    for ( unsigned int uint_row = 0; uint_row < uint_size; uint_row++) {
        for ( unsigned int uint_col = 0; uint_col < uint_size; uint_col++) {
            mat_ret(uint_row,uint_col) = mat_arg(std::abs((long int)(uint_col - uint_row)));
        }
    }
    return mat_ret;
}

} // NAMESPACE SUSA
#endif // SIGNAL_H
