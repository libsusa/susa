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
 * @author Behrooz Aliabadi
 * @version 1.0.0
 */

#include "test.h"

int main(int argc, char const *argv[])
{
    SUSA_TEST_EQ_DOUBLE(susa::normcdf(0.5), 0.6915, "Normal Cumulative Distribution Function");

    std::cout << std::endl << " -----------------";
    std::cout << std::endl << " NUMBER OF FAILED TESTS(" << uint_failed <<")";
    std::cout << std::endl << " NUMBER OF PASSED TESTS(" << uint_passed <<")";
    std::cout << std::endl << " TOTAL NUMBER OF TESTS (" << uint_total <<")";
    std::cout << std::endl;

    return (uint_failed);
}
