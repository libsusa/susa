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
 * @file rrcosine.h
 * @brief Root-Raised Cosine
 * @author Behrooz Kamary Aliabadi
 * @version 1.0.0
 *
 * @defgroup Communications
 */

#ifndef SUSA_RRCOSINE_H
#define SUSA_RRCOSINE_H

namespace susa
{

/**
 * @brief Root-Raised Cosine wrapper class.
 *
 * @ingroup Communications
 **/
class RRCosine {
  private:

    double xrc(double, double, double);

    matrix < std::complex <double> > cmat_g_T;
    matrix <double>  mat_g_T;

  public:

    RRCosine();

    RRCosine(double, double, double, int);

    matrix <double> get();
};


template <class T> int antipodal(T T_arg);
}       // NAMESPACE SUSA
#endif  // SUSA_RRCOSINE_H
