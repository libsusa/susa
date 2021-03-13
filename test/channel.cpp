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
 * @file channel.cpp
 * @brief Unit Test Suit
 * @author Behrooz Kamary
 */

#include "test.h"

int main(void)
{
    // ignore the tail of the bits
    const unsigned len = 15;

    // the actual transmitted bits
    susa::matrix<double> mat_bits("-1 -1 -1 -1 1 1 -1 1 -1 1 1 -1 1 -1 -1 1 -1 1");
    mat_bits.resize(18,1);

    // pulse amplitude (pam) vector
    susa::matrix<double> mat_pam("-1 1");

    // channel taps
    susa::matrix<double> mat_taps("1 -0.9 0.5 0.4");

    susa::channel<double> ch(mat_taps, mat_pam);

    susa::matrix<double> mat_coded = filter(mat_taps, susa::matrix<double> (1, 1, 1), mat_bits);

    susa::matrix<double> mat_decoded = ch.decode_bcjr(mat_coded, 10);

    susa::matrix<double> mat_decoded_bits(mat_decoded.shape());

    for (unsigned int uint_i = 0; uint_i < mat_decoded.size(); uint_i++)
    {
        mat_decoded_bits(uint_i) = mat_decoded(uint_i) > 1 ? 1 : -1;
    }

    susa::matrix <double> mat_decoded_mlse = ch.decode_mlse(mat_coded);

    SUSA_TEST_EQ(mat_bits.left(len), mat_decoded_bits.left(len), "BCJR channel equalizer.");
    SUSA_TEST_EQ(mat_bits.left(len), mat_decoded_mlse.left(len), "MLSE (Viterbi) channel equalizer.");
    SUSA_TEST_PRINT_STATS();

    return (uint_failed);
}
