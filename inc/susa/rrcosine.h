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
 * @brief Root-Raised Cosine (declaration).
 * @author Behrooz Kamary
 *
 * @defgroup Communications
 */

#ifndef SUSA_RRCOSINE_H
#define SUSA_RRCOSINE_H

namespace susa
{
    

/**
 * @brief Root-Raised Cosine Filter
 *
 * This class generates Root-Raised Cosine FIR filter taps.
 * The filter taps are computed through frequency sampling method.
 *
 * @ingroup Communications
 */
class rrcosine
{
  private:

    double xrc(double, double);

    cmatrix<double>     cmat_g_T;
    matrix <double>     mat_g_T;

  public:

    /**
     * @brief Constructor
     *
     * @param dbl_l the filter interpolation (upsampling) factor L = Fs/Fd.
     * @param dbl_alpha the roll-off factor that determines the width of the transition band of the filter.
     * @param int_ord the FIR filter order that must be an odd positive integer.
     *
     */
    rrcosine(unsigned dbl_l, double dbl_alpha, int uint_ord);

    susa::matrix <double> get();
};

}       // NAMESPACE SUSA
#endif  // SUSA_RRCOSINE_H
