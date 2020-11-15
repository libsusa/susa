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
 * @file fft.h
 * @brief Fast Fourier Transform (FFT) (declaration and definition).
 * @author Behrooz Kamary
 * @version 1.0.0
 */

#ifndef SUSA_FFT_H
#define SUSA_FFT_H

namespace susa {

/**
 * @brief The Fast Fourier Transform (FFT) class.
 *
 * @ingroup Signal
 */
template <typename T, typename Enable = void>
class fft;

template <typename T>
class fft <T, typename std::enable_if_t<std::is_floating_point<T>::value>>
{
  public:
    /**
    * @brief Fast Fourier Transform (FFT) using the radix-2 algorithm
    *
    * @todo it only supports 2^N input sizes. it shall cope with any vector size.
    * @param mat_arg input matrix
    * @return returns a matrix
    */
    matrix <std::complex<T>> radix2(const matrix <T>& mat_arg)
    {
        matrix <std::complex<T>> mat_ret;

        if (mat_arg.no_rows() == 1 || mat_arg.no_cols() == 1)
        {
            mat_ret = vector_radix2(mat_arg);
        }
        else if (mat_arg.no_rows() > 1 && mat_arg.no_cols() > 1)
        {
            mat_ret = matrix <std::complex<T>> (mat_arg.shape());
            for (size_t sz_col = 0; sz_col < mat_arg.no_cols(); sz_col++)
            {
                mat_ret.set_col(sz_col, vector_radix2(mat_arg.col(sz_col)));
            }
        }

        return mat_ret;
    }


    /**
    * @brief Discrete Fourier Transform (DFT)
    *
    * @param mat_arg input matrix
    * @return returns a matrix
    */
    matrix <std::complex<T>> dft(const matrix <T>& mat_arg)
    {
        matrix <std::complex<T>> mat_ret;

        if (mat_arg.no_rows() == 1 || mat_arg.no_cols() == 1)
        {
            mat_ret = vector_dft(mat_arg);
        }
        else if (mat_arg.no_rows() > 1 && mat_arg.no_cols() > 1)
        {
            mat_ret = matrix <std::complex<T>> (mat_arg.shape());
            for (size_t sz_col = 0; sz_col < mat_arg.no_cols(); sz_col++)
            {
                mat_ret.set_col(sz_col, vector_dft(mat_arg.col(sz_col)));
            }
        }

        return mat_ret;
    }


  private:

    matrix <std::complex<T>> vector_radix2(const matrix <T>& mat_arg)
    {
        size_t size = mat_arg.size();
        SUSA_ASSERT(is_power_of_two(size));

        matrix <std::complex<T>> mat_bit_rev        = bit_reverse(convert_to_complex(mat_arg));
        size_t sz_stage                             = susa::log2(size);

        for (size_t stage = 1; stage <= sz_stage; stage++)
        {
            size_t m = std::pow(2,stage);
            std::complex<T> w_m = std::exp(std::complex<T>(0,(T)-2 * PI / m));
            for (size_t k = 0; k < size; k += m)
            {
                std::complex<T> w(1,0);
                for (size_t j = 0; j < m/2; j++)
                {
                    std::complex<T> t           = w * mat_bit_rev(k + j + m/2);
                    std::complex<T> u           = mat_bit_rev(k + j);
                    mat_bit_rev(k + j)          = u + t;
                    mat_bit_rev(k + j + m/2)    = u - t;
                    w                           = w * w_m;
                }
            }
        }

        return mat_bit_rev;
    }

    matrix <std::complex<T>> vector_dft(const matrix <T>& mat_arg)
    {
        size_t size = mat_arg.size();
        matrix <std::complex<T>> mat_ret(size,1);

        for (size_t k = 0; k < size; k++)
        {
            for (size_t n = 0; n < size; n++)
            {
                mat_ret(k) += mat_arg(n) * std::exp(std::complex<T>(0, (T)-2 * PI * n * k / size));
            }
        }

        return mat_ret;
    }

    template <typename U>
    matrix <U> bit_reverse(matrix <U> mat_arg)
    {
        size_t  n       = mat_arg.size();
        size_t  n_bits  = susa::log2(n);
        U       U_swap  = 0;
        size_t  sz_brev;


        for (size_t indx=1; indx < n; indx++)
        {
            sz_brev = 0;
            for (size_t j=0; j < n_bits; j++)
            {
                sz_brev |= (( indx >> j) & 1U) << (n_bits - j - 1);
            }
            if (sz_brev > indx)
            {
                U_swap = mat_arg(indx);
                mat_arg(indx) = mat_arg(sz_brev);
                mat_arg(sz_brev) = U_swap;
            }
        }
        return mat_arg;
    }

    matrix <std::complex<T>> convert_to_complex(const matrix <T>& mat_arg)
    {
        size_t size = mat_arg.size();
        matrix <std::complex<T>> ret(size,1);

        for (size_t indx = 0; indx < size; indx++)
        {
            ret(indx) = std::complex<T>(mat_arg(indx),0);
        }

        return ret;
    }

    bool is_power_of_two(size_t arg)
    {
        return arg && (!(arg & (arg - 1)));
    }
};


} // NAMESPACE SUSA
#endif // SUSA_FFT_H
