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
 * @file modulation.cpp
 * @brief These methods implement digital base-band communication schemes.
 *
 *
 * @author Behrooz Kamary Aliabadi
 * @version 1.0.0
 */


#include "susa.h"

namespace susa {

unsigned int qam::log2(unsigned int uint_x)
{
    unsigned int r = 0;

    while( (uint_x >> r) != 0)
    {
        r++;
    }
    return (r - 1);
}

qam::qam(unsigned int uint_m)
{

    uint_bps = log2(uint_m);
    double dbl_p = sqrt(uint_m);
    unsigned int  uint_l = (unsigned int)(uint_m/dbl_p);

    this->uint_m = uint_m;

    // Constellation generation
    mat_s    = matrix < std::complex <double> > (uint_l, uint_l, std::complex <double>(0,0));
    mat_axis = matrix < std::complex <double> > ((unsigned int)dbl_p, 1, std::complex <double>(0,0));

    for (unsigned int i=0; i<dbl_p; i++)
        mat_axis(i) = std::complex <double> (i - (dbl_p - 1)/2,i - (dbl_p - 1)/2);

    for (unsigned int i=0; i<uint_l; i++)
    {
        for (unsigned int j=0; j<uint_l; j++)
        {
            mat_s(i,j) = std::complex <double> (mat_axis(i).real(),mat_axis(j).imag());
        }
    }

    // Energy computations
    matrix <double> mat_mag = mag(mat_s);
    matrix <double> mat_sum = sum(mat_mag);
    dbl_es = sum(mat_sum)(0)/uint_m;
    dbl_eb = dbl_es / uint_bps;
}

std::complex <double> qam::demodulate_symbol(std::complex <double> complex_arg)
{
    matrix <double> mat_dist(uint_m,1);
    for (unsigned int uint_i = 0; uint_i < uint_m; uint_i++)
        mat_dist(uint_i) = std::abs(complex_arg - mat_s(uint_i));
    return mat_s(min(mat_dist)(0));
}

double qam::get_noise_deviation(double dbl_arg)
{
    return std::sqrt( ( 0.5 * dbl_es/uint_bps ) * std::pow(10, - ( dbl_arg/10 )) );
}

qam::~qam() {}

} // NAMESPACE SUSA
