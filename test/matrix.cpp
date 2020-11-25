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
 * @file matrix.cpp
 * @brief Unit Test Suit
 * @author Behrooz Kamary
 * @version 1.0.0
 */

#include "test.h"

int main(void)
{
    susa::matrix <float> mat_m("[0 1 1; 2 3 2; 1 3 2; 4 2 2]");
    auto mat_left = mat_m.left(3);
    SUSA_TEST_EQ(mat_m, mat_left, "full matrix left");

    auto mat_right = mat_m.right(3);
    SUSA_TEST_EQ(mat_m, mat_right, "full matrix right");

    susa::matrix <float> mat_l("[0 1 ; 2 3 ; 1 3 ; 4 2 ]");
    mat_left = mat_m.left(2);
    SUSA_TEST_EQ(mat_l, mat_left, "matrix left");

    susa::matrix <float> mat_r("[1 1;3 2;3 2;2 2]");
    mat_right = mat_m.right(2);
    SUSA_TEST_EQ(mat_r, mat_right, "matrix right");

    auto mat_mid = mat_m.mid(1,2);
    SUSA_TEST_EQ(mat_mid, mat_right, "matrix mid");

    SUSA_TEST_PRINT_STATS();

    return (uint_failed);
}
