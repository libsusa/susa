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
 * @file signal.cpp
 * @brief Unit Test Suit
 * @author Behrooz Kamary Aliabadi
 * @version 1.0.0
 */

#include "test.h"

int main(int argc, char const *argv[])
{
    susa::matrix<int> mat_a("[1 1 1]");
    susa::matrix<int> mat_exp("[1 2 3 2 1]");
    susa::matrix<int> mat_res = conv(mat_a, mat_a);
    SUSA_TEST_EQ(mat_res, mat_exp, "convolution.");
    SUSA_TEST_PRINT_STATS();
}