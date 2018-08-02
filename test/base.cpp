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
 * @file base.cpp
 * @brief Unit Test Suit
 * @author Behrooz Kamary Aliabadi
 * @version 1.0.0
 */

#include "test.h"

int main(int argc, char const *argv[])
{
    SUSA_TEST_EQ(susa::round(42.163574), 42.0, "Round a double with no decimals.");
    SUSA_TEST_EQ(susa::round(42.163574, 2), 42.16, "Round a double with two decimals.");
    SUSA_TEST_EQ_DOUBLE(susa::normcdf(0.5), 0.6915, "Normal Cumulative Distribution Function.");
    SUSA_TEST_EQ(susa::round(susa::normcdf(0.5), 4), 0.6915, "Normal Cumulative Distribution Function.");

    susa::rng mt_rng(35748);
    susa::matrix <float> mat_a = mt_rng.rand(20,10);
    std::stringstream ss;
    ss << round(mat_a, 2);
    susa::matrix <double> mat_b(ss.str());
    SUSA_TEST_EQ(mat_b, round(mat_a, 2), "std::string to susa::matrix conversion.");

    SUSA_TEST_PRINT_STATS();

    return (uint_failed);
}
