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
 * @file fft.cpp
 * @brief Unit Test Suit
 * @author Behrooz Kamary
 * @version 1.0.0
 */

#include "test.h"

int main(void)
{

{
    susa::fft <double, susa::allocator_log<double>>      fourier;
    susa::matrix <double, susa::allocator_log<double>>   mat_input("[4.,3.,2.,1.,0.,1.,2.,3.]");
    auto mat_output = fourier.radix2(mat_input);
    auto sum = std::abs(susa::sum(mat_output)(0));
    SUSA_TEST_EQ_DOUBLE(sum, 32.0, "fast fourier transform (FFT) with Radix-2 (symetric reals)");
}

{
    // test with the default allocator
    susa::fft<double> fourier;
    susa::matrix <double> mat_input("[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]");

    auto mat_output_dft = fourier.dft(mat_input);
    auto mat_output_fft = fourier.radix2(mat_input);
    bool bool_dft_fft_eq = true;

    for (size_t indx = 0; indx < mat_input.size(); indx++)
    {
        auto diff = std::abs(mat_output_dft(indx) - mat_output_fft(indx));
        if (diff > 1.0e-10) bool_dft_fft_eq = false;
    }

    SUSA_TEST_EQ(bool_dft_fft_eq, true, "fast fourier transform (FFT) with Radix-2 (asymetric reals)");
}

{
    // test with a custom allocator
    susa::fft<double, susa::allocator_log<double>> fourier;
    susa::matrix <double, susa::allocator_log<double>> mat_input("[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]");

    auto mat_output_dft = fourier.dft(mat_input);
    auto mat_output_fft = fourier.radix2(mat_input);
    bool bool_dft_fft_eq = true;

    for (size_t indx = 0; indx < mat_input.size(); indx++)
    {
        auto diff = std::abs(mat_output_dft(indx) - mat_output_fft(indx));
        if (diff > 1.0e-10) bool_dft_fft_eq = false;
    }

    SUSA_TEST_EQ(bool_dft_fft_eq, true, "fast fourier transform (FFT) with Radix-2 (asymetric reals)");
}

    SUSA_TEST_EQ(susa::memory_tacker::instance().read(), 0, "susa::array memory leak with susa allocator");

    SUSA_TEST_PRINT_STATS();

    return (uint_failed);
}