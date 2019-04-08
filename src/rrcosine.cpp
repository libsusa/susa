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
 * @version 1.0.0
 */

#include <susa.h>

namespace susa {

// RRCosine

RRCosine::RRCosine(double dFs, double dFd, double dAlpha, int intN)
{
    cmat_g_T = matrix <std::complex <double>> (intN,1,std::complex<double>(0,0));
    double rc_norm = 0;

    if (intN % 2 == 1)
    {
        for (int i=0; i<intN; i++)
        {
            for ( int m=-(intN - 1)/2; m <= (intN - 1)/2 ; m++ )
            {
                cmat_g_T(i) += sqrt(xrc(dFs*m/intN,dAlpha,dFd))*exp(std::complex<double>(0,2*PI*m*(i - (intN - 1)/2)/intN));
            }
        }

        mat_g_T = real(cmat_g_T);

        for (unsigned int i=0; i<mat_g_T.size(); i++)
        {
            rc_norm += mat_g_T(i) * mat_g_T(cmat_g_T.size() - i - 1);
        }

        for (unsigned int i=0; i<mat_g_T.size(); i++)
        {
            mat_g_T(i) = mat_g_T(i) / sqrt(rc_norm);
        }
    }
    else
    {
        // TODO: fix this odd case
        SUSA_ABORT("the filter order must be odd.");
    }
}

matrix <double> RRCosine::get()
{
    return mat_g_T;
}

// Private methods

double RRCosine::xrc(double dF, double dAlpha, double dT)
{
    if (std::abs(dF) > (1 + dAlpha)/(2*dT))
        return 0;
    else if (std::abs(dF) > ((1 - dAlpha)/(2*dT)))
        return (dT/2.0)*(1+cos((PI*dT/dAlpha)*(std::abs(dF)-(1-dAlpha)/(2*dT))));
    else return dT;

}



// Methodes
template <class T> int antipodal(T T_arg)
{
    if (T_arg == 0) return -1;
    else return 1;
}

}
