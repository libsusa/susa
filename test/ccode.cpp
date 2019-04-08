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
 * @file ccode.cpp
 * @brief Unit Test Suit
 * @author Behrooz Kamary
 * @version 1.0.0
 */
#include <susa.h>

#include "test.h"

using namespace std;
using namespace susa;


int main(void)
{

    {
        susa::ccode state(2, 1, 2);
        state.set_generator(7, 0);
        state.set_generator(5, 1);

        matrix<uint8_t> mat_bits("1 0 1 1 0");

        matrix<uint8_t> mat_expected("1 1 1 0 0 0 0 1 0 1");

        matrix<uint8_t> mat_coded = state.encode(mat_bits);

        SUSA_TEST_EQ(mat_expected, transpose(mat_coded), "Convolutional Code Encoder");
    }

    {
        susa::ccode state(2, 1, 2);
        state.set_generator(7, 0);
        state.set_generator(5, 1);

        matrix<uint8_t> mat_bits("1 0 1 0 0 0");

        matrix<uint8_t> mat_expected("1 1 1 0 0 0 1 0 1 1 0 0");

        matrix<uint8_t> mat_coded = state.encode(mat_bits);

        SUSA_TEST_EQ(mat_expected, transpose(mat_coded), "Convolutional Code Encoder");
    }

    {
        susa::ccode state(2, 1, 2);
        state.set_generator(7, 0);
        state.set_generator(5, 1);

        matrix<uint8_t> mat_bits("0 1 0 1 1 1 0 0 1 0 1 0 0 0 1 0 0");

        matrix<uint8_t> mat_expected("0 0 1 1 1 0 0 0 0 1 1 0 0 1 1 1 1 1 1 0 0 0 1 0 1 1 0 0 1 1 1 0 1 1");

        matrix<uint8_t> mat_coded = state.encode(mat_bits);

        SUSA_TEST_EQ(mat_expected, transpose(mat_coded), "Convolutional Code Encoder");
    }

    SUSA_TEST_PRINT_STATS();

    return (uint_failed);
}