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
 * @file arithmetics.cpp
 * @brief Unit Test Suit
 * @author Behrooz Kamary
 */

#include "test.h"
#include <susa/arithmetics.h>

int main(void)
{

    double dbla = -457.103;
    susa::fixed_point <int64_t, 20,30> fpnuma(dbla);
    SUSA_TEST_EQ_DOUBLE(dbla, (double)fpnuma, "coversion to/from a floating point");

    double dblb = 74.166;
    susa::fixed_point <int64_t, 20,30> fpnumb(dblb);
    SUSA_TEST_EQ_DOUBLE(dblb, (double)fpnumb, "coversion to/from a floating point");


    double dbl_add = dbla + dblb;
    auto fp_add = fpnuma + fpnumb;
    SUSA_TEST_EQ_DOUBLE(dbl_add, (double)fp_add, "addition");


    double dbl_sub = dbla - dblb;
    auto fp_sub = fpnuma - fpnumb;
    SUSA_TEST_EQ_DOUBLE(dbl_sub, (double)fp_sub, "subtraction");


    double dbl_mul = dbla * dblb;
    auto fp_mul = fpnuma * fpnumb;
    SUSA_TEST_EQ_DOUBLE(dbl_mul, (double)fp_mul, "multiplication");


    double dbl_div = dbla / dblb;
    auto fp_div = fpnuma / fpnumb;
    SUSA_TEST_EQ_DOUBLE(dbl_div, (double)fp_div, "division");

    SUSA_TEST_PRINT_STATS();

    return (uint_failed);
}
