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
 * @file base.cpp
 * @brief Basic mathematical operations on STL and Susa types.
 * @author Behrooz, Kamary Aliabadi
 * @version 1.0.0
 */


#include "../inc/susa.h"

namespace susa {

double ncdf(const double x) {
    // maximum absolute error : 7.5E-8

    double dbl_ret = 0;

    const double b1 =  0.319381530;
    const double b2 = -0.356563782;
    const double b3 =  1.781477937;
    const double b4 = -1.821255978;
    const double b5 =  1.330274429;
    const double p  =  0.2316419;
    const double c  =  0.39894228;

    if(x >= 0.0) {
        double t = 1.0 / ( 1.0 + p * x );
        dbl_ret = (1.0 - c * std::exp( -x * x / 2.0 ) * t *
                   ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
    } else {
        double t = 1.0 / ( 1.0 - p * x );
        dbl_ret = ( c * std::exp( -x * x / 2.0 ) * t *
                    ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
    }

    return dbl_ret;
}

double qfunc(const double x) {
    // maximum absolute error : 7.5E-8

    double dbl_ret = 0;

    const double b1 =  0.319381530;
    const double b2 = -0.356563782;
    const double b3 =  1.781477937;
    const double b4 = -1.821255978;
    const double b5 =  1.330274429;
    const double p  =  0.2316419;
    const double c  =  0.39894228;

    if(x >= 0.0) {
        double t = 1.0 / ( 1.0 + p * x );
        dbl_ret = (1.0 - c * std::exp( -x * x / 2.0 ) * t *
                   ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
    } else {
        double t = 1.0 / ( 1.0 - p * x );
        dbl_ret = ( c * std::exp( -x * x / 2.0 ) * t *
                    ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
    }

    return (1 - dbl_ret);
}

unsigned int pow(unsigned int uint_b, unsigned int uint_p) {
    unsigned int uint_ret = 1;
    if (uint_p == 0) return 1;
    for (unsigned int uint_i = 0; uint_i < uint_p; uint_i++) uint_ret *= uint_b;
    return uint_ret;
}

int pow(int int_b, unsigned int uint_p) {
    int int_ret = 1;
    if (uint_p == 0) return 1;
    for (unsigned int uint_i = 0; uint_i < uint_p; uint_i++) int_ret *= int_b;
    return int_ret;
}


long int mod(long int lint_a, long int lint_mod) {
    return(lint_a - std::floor((float)lint_a/lint_mod) * lint_mod);
}

} // NAMESPACE SUSA
