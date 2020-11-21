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
 */

#ifndef SUSA_FFT_H
#define SUSA_FFT_H

namespace susa {

/**
 * @brief The Fast Fourier Transform (FFT) class.
 *
 * @ingroup Signal
 */
template <typename T, typename Allocator = std::allocator<T>, typename Enable = void>
class fft;

/**
 * @brief The Fast Fourier Transform (FFT) class.
 *
 * @ingroup Signal
 */
template <typename T, typename Allocator>
class fft <T, Allocator, typename std::enable_if_t<std::is_floating_point<T>::value>>
{
    private:
        using value_allocator = typename get_allocator<T, Allocator>::type;
        using complex_allocator = typename get_allocator<std::complex<T>, Allocator>::type;
    public:
  
    /**
    * @brief Fast Fourier Transform (FFT) using the radix-2 algorithm
    *
    * @todo it only supports 2^N input sizes. it shall cope with any vector size.
    * @param mat_arg input matrix of real numbers
    * @return returns a matrix of complex numbers
    */
    matrix <std::complex<T>, complex_allocator> radix2(const matrix <T, Allocator>& mat_arg)
    {
        matrix <std::complex<T>, complex_allocator> mat_ret;

        if (mat_arg.is_vector())
        {
            mat_ret = vector_radix2(mat_arg);
        }
        else
        {
            mat_ret = matrix <std::complex<T>, complex_allocator> (mat_arg.shape());
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
    * @param mat_arg input matrix of real numbers
    * @return returns a matrix of complex numbers
    */
    matrix <std::complex<T>, complex_allocator> dft(const matrix <T, Allocator>& mat_arg)
    {
        matrix <std::complex<T>, complex_allocator> mat_ret;

        if (mat_arg.is_vector())
        {
            mat_ret = vector_dft(mat_arg);
        }
        else
        {
            mat_ret = matrix <std::complex<T>, complex_allocator> (mat_arg.shape());
            for (size_t sz_col = 0; sz_col < mat_arg.no_cols(); sz_col++)
            {
                mat_ret.set_col(sz_col, vector_dft(mat_arg.col(sz_col)));
            }
        }

        return mat_ret;
    }


  private:

    matrix <std::complex<T>, complex_allocator> vector_radix2(const matrix <T, Allocator>& mat_arg)
    {
        size_t          size        = mat_arg.size();
        SUSA_ASSERT(is_power_of_two(size));
        size_t          sz_stage    = susa::log2(size);
        matrix <std::complex<T>, complex_allocator> mat_brev = convert_to_complex(mat_arg);

        bit_reverse(mat_brev);

        for (size_t stage = 1; stage <= sz_stage; stage++)
        {
            size_t m = std::pow(2,stage);
            std::complex<T> w_m = std::exp(std::complex<T>(0,(T)-2 * PI / m));
            for (size_t k = 0; k < size; k += m)
            {
                std::complex<T> w(1,0);
                for (size_t j = 0; j < m/2; j++)
                {
                    std::complex<T> t           = w * mat_brev(k + j + m/2);
                    std::complex<T> u           = mat_brev(k + j);
                    mat_brev(k + j)             = u + t;
                    mat_brev(k + j + m/2)       = u - t;
                    w                           = w * w_m;
                }
            }
        }

        return mat_brev;
    }

    matrix <std::complex<T>, complex_allocator> vector_dft(const matrix <T, Allocator>& mat_arg)
    {
        size_t size = mat_arg.size();
        matrix <std::complex<T>, complex_allocator> mat_ret(mat_arg.shape(), std::complex <T> (0,0));

        for (size_t k = 0; k < size; k++)
        {
            for (size_t n = 0; n < size; n++)
            {
                mat_ret(k) += mat_arg(n) * std::exp(std::complex<T>(0, (T)-2 * PI * n * k / size));
            }
        }

        return mat_ret;
    }

    void bit_reverse(matrix <std::complex<T>, complex_allocator>& mat_arg)
    {
        size_t              sz_size     = mat_arg.size();
        size_t              n_bits      = susa::log2(sz_size);
        std::complex<T>     U_swap      = 0;
        size_t              sz_brev;


        for (size_t indx=1; indx < sz_size; indx++)
        {
            sz_brev = 0;
            for (size_t j=0; j < n_bits; j++)
            {
                sz_brev |= (( indx >> j) & 1U) << (n_bits - j - 1);
            }
            if (sz_brev > indx)
            {
                U_swap              = mat_arg(indx);
                mat_arg(indx)       = mat_arg(sz_brev);
                mat_arg(sz_brev)    = U_swap;
            }
        }
    }


    matrix <std::complex<T>, complex_allocator> convert_to_complex(const matrix <T, Allocator>& mat_arg)
    {
        size_t size = mat_arg.size();
        matrix <std::complex<T>, complex_allocator> ret(mat_arg.shape());

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
