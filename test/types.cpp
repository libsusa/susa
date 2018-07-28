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
 * @file types.cpp
 * @brief Unit Test Suit
 * @author Behrooz Kamary Aliabadi
 * @version 1.0.0
 */

#include "test.h"
#include <utility>

int main(int argc, char const *argv[])
{


  susa::matrix <double> mat_a("[1 2.3 -3.4;8 4.5 1.2;9.1 3 -5]");
  SUSA_TEST_EQ(mat_a(1,1),    4.5, "matrix parser.");
  SUSA_TEST_EQ(mat_a(1),        8, "matrix parser.");
  SUSA_TEST_EQ(mat_a(4),      4.5, "matrix parser.");
  SUSA_TEST_EQ(mat_a(2,2),     -5, "matrix parser.");

  mat_a(1,1) = 6.6;
  SUSA_TEST_EQ(mat_a(1,1),    6.6, "matrix element assignment.");

  susa::matrix <int> mat_b("[1; 2; 3]");
  SUSA_TEST_EQ(mat_b.no_rows(), 3, "matrix parser rows.");
  SUSA_TEST_EQ(mat_b.no_cols(), 1, "matrix parser cols.");
  SUSA_TEST_EQ(mat_b.size(), 3, "matrix parser size.");

  susa::matrix <int> mat_e = mat_b;
  susa::matrix <int> mat_d (std::move(mat_b));
  SUSA_TEST_EQ(mat_d, mat_e, "move semantic.");

  susa::matrix<int> mat_f("2;2;2");
  SUSA_TEST_EQ(mat_f * mat_e, mat_e * 2, "scalar and matrix multiplication.");

  susa::matrix <int> mat_g = mat_a;
  SUSA_TEST_EQ(mat_g.no_rows(), mat_a.no_rows(), "matrix copy assignment rows.");
  SUSA_TEST_EQ(mat_g.no_cols(), mat_a.no_cols(), "matrix copy assignment cols.");
  SUSA_TEST_EQ(mat_g.size(), mat_a.size(), "matrix copy assignment cols.");

  susa::matrix <int> mat_h (mat_a);
  SUSA_TEST_EQ(mat_h.no_rows(), mat_a.no_rows(), "matrix copy constructor rows.");
  SUSA_TEST_EQ(mat_h.no_cols(), mat_a.no_cols(), "matrix copy constructor cols.");
  SUSA_TEST_EQ(mat_h.size(), mat_a.size(), "matrix copy constructor cols.");


  susa::array <int> arr_a({21,6,5,15,43});
  arr_a(2,4,3,0,1) = 55;
  arr_a(12,4,3,5,1) = 32;
  susa::array <int> arr_b = arr_a;

  SUSA_TEST_EQ(55, arr_b(2,4,3,0,1), "array copy assignment.");
  SUSA_TEST_EQ(32, arr_b(12,4,3,5,1), "array copy assignment.");


  SUSA_TEST_PRINT_STATS();

  return (uint_failed);
}
