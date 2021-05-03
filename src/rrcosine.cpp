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
 * @file rrcosine.cpp
 * @brief Root-Raised Cosine (definition).
 * @author Behrooz Kamary
 */

#include <susa.h>

namespace susa {


rrcosine::rrcosine(unsigned dbl_l, double dbl_alpha, int int_ord)
{
    SUSA_ASSERT(int_ord > 0);
    SUSA_ASSERT_MESSAGE(int_ord % 2 == 1, "the filter order must be odd.");

    cmat_g_T = cmatrix <double> (int_ord, 1, std::complex<double>(0,0));
    double rc_norm = 0;

    for (int i = 0; i < int_ord; i++)
    {
        for (int m=-(int_ord - 1)/2; m <= (int_ord - 1)/2 ; m++)
        {
            cmat_g_T(i) += std::sqrt(xrc(dbl_l * m / int_ord, dbl_alpha)) *
                            exp(std::complex<double>(0, 2 * PI * m * (i - (int_ord - 1)/2)/int_ord));
        }
    }

    mat_g_T = real(cmat_g_T);

    for (unsigned int i=0; i<mat_g_T.size(); i++)
    {
        rc_norm += mat_g_T(i) * mat_g_T(cmat_g_T.size() - i - 1);
    }

    for (unsigned int i=0; i<mat_g_T.size(); i++)
    {
        mat_g_T(i) = mat_g_T(i) / std::sqrt(rc_norm);
    }
}

// Private methods
double rrcosine::xrc(double dbl_f, double dbl_alpha)
{
    if (std::abs(dbl_f) > ((1 + dbl_alpha) * 0.5))
        return 0;
    else if (std::abs(dbl_f) > ((1 - dbl_alpha)*0.5))
        return 0.5 * (1 + cos((PI / dbl_alpha) * (std::abs(dbl_f) - (1 - dbl_alpha) * 0.5)));
    else return 1;
}

// Public methods
matrix <double> rrcosine::get()
{
    return mat_g_T;
}

}
