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
 * @file base.cpp
 * @brief Basic mathematical operations (definition).
 *
 * @author Behrooz Kamary
 */


#include <susa.h>

namespace susa {

double normcdf(const double x)
{
    // maximum absolute error : 7.5E-8

    double dbl_ret = 0;

    const double b1 =  0.319381530;
    const double b2 = -0.356563782;
    const double b3 =  1.781477937;
    const double b4 = -1.821255978;
    const double b5 =  1.330274429;
    const double p  =  0.2316419;
    const double c  =  0.39894228;

    if(x >= 0.0)
    {
        double t = 1.0 / ( 1.0 + p * x );
        dbl_ret = (1.0 - c * std::exp( -x * x / 2.0 ) * t *
                   ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
    }
    else
    {
        double t = 1.0 / ( 1.0 - p * x );
        dbl_ret = ( c * std::exp( -x * x / 2.0 ) * t *
                    ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
    }

    return dbl_ret;
}

double qfunc(const double x)
{
    return (1.0 - normcdf(x));
}

unsigned int pow(unsigned int uint_b, unsigned int uint_p)
{
    unsigned int uint_ret = 1;
    if (uint_p == 0) return 1;
    for (unsigned int uint_i = 0; uint_i < uint_p; uint_i++) uint_ret *= uint_b;
    return uint_ret;
}

int pow(int int_b, unsigned int uint_p)
{
    int int_ret = 1;
    if (uint_p == 0) return 1;
    for (unsigned int uint_i = 0; uint_i < uint_p; uint_i++) int_ret *= int_b;
    return int_ret;
}

long int mod(long int lint_a, long int lint_mod)
{
    return(lint_a - std::floor((float)lint_a / lint_mod) * lint_mod);
}

double round(const double dbl_arg, int int_decimal)
{
    //TODO: to be fixed for large doubles to prevent the overflow.
    double dbl_p = std::pow(10.0, int_decimal);

    if (int_decimal > 0)
    {
        return (std::round(dbl_arg * dbl_p) / dbl_p);
    }

    return std::round(dbl_arg);
}


matrix <double> round(const matrix <double> &mat_arg, int int_decimal)
{
    unsigned int uint_size = mat_arg.size();
    matrix <double> mat_ret(mat_arg.no_rows(), mat_arg.no_cols());
    for (unsigned int uint_indx = 0; uint_indx < uint_size; uint_indx++)
    {
        mat_ret(uint_indx) = susa::round(mat_arg(uint_indx), int_decimal);
    }

    return mat_ret;
}

uint64_t log2ull(uint64_t uint_arg)
{
    static const int tab64[64] = {
        63,  0, 58,  1, 59, 47, 53,  2,
        60, 39, 48, 27, 54, 33, 42,  3,
        61, 51, 37, 40, 49, 18, 28, 20,
        55, 30, 34, 11, 43, 14, 22,  4,
        62, 57, 46, 52, 38, 26, 32, 41,
        50, 36, 17, 19, 29, 10, 13, 21,
        56, 45, 25, 31, 35, 16,  9, 12,
        44, 24, 15,  8, 23,  7,  6,  5};

    uint_arg |= uint_arg >> 1;
    uint_arg |= uint_arg >> 2;
    uint_arg |= uint_arg >> 4;
    uint_arg |= uint_arg >> 8;
    uint_arg |= uint_arg >> 16;
    uint_arg |= uint_arg >> 32;
    return tab64[((uint64_t)((uint_arg - (uint_arg >> 1))*0x07EDD5E59A4E28C2)) >> 58];
}

uint32_t log2u(uint32_t uint_arg)
{
    static const int tab32[32] = {
        0,  9,  1, 10, 13, 21,  2, 29,
        11, 14, 16, 18, 22, 25,  3, 30,
        8, 12, 20, 28, 15, 17, 24,  7,
        19, 27, 23,  6, 26,  5,  4, 31};

    uint_arg |= uint_arg >> 1;
    uint_arg |= uint_arg >> 2;
    uint_arg |= uint_arg >> 4;
    uint_arg |= uint_arg >> 8;
    uint_arg |= uint_arg >> 16;
    return tab32[(uint32_t)(uint_arg*0x07C4ACDD) >> 27];
}

} // NAMESPACE SUSA
