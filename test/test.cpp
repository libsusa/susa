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
 * @file test.cpp
 * @brief Unit Test Suit
 * @author Behrooz Aliabadi
 * @version 1.0.0
 */

#include "test.h"

int main(int argc, char const *argv[]) {

  std::cout << std::endl << " --- SUSA UNIT TEST SUIT ---";

  susa::matrix <double> mat_a("[1 2.3 -3.4;8 4.5 1.2;9.1 3 -5]");
  SUSA_TEST_EQ(mat_a(1,1),    4.5, "matrix parser");
  SUSA_TEST_EQ(mat_a(1),      8,   "matrix parser");
  SUSA_TEST_EQ(mat_a(4),      4.5, "matrix parser");
  SUSA_TEST_EQ(mat_a(2,2),    -5,  "matrix parser");

  SUSA_TEST_EQ_DOUBLE(susa::normcdf(0.5),0.6915, "Normal Cumulative Distribution Function");

  std::cout << std::endl << " -----------------";
  std::cout << std::endl << "NUMBER OF FAILED TESTS(" << uint_failed <<")";
  std::cout << std::endl << "NUMBER OF PASSED TESTS(" << uint_passed <<")";
  std::cout << std::endl << "TOTAL NUMBER OF TESTS (" << uint_total <<")";

  return (uint_failed);
}
