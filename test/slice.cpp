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
 * @file slice.cpp
 * @brief Unit Test Suit
 * @author Behrooz Kamary
 */

#include "test.h"

int main(void)
{
  susa::matrix <int, susa::allocator_log<int>> mat_m("[0 1 1;"
                                                      "2 3 2;"
                                                      "1 3 2;"
                                                      "4 2 2]");
  {
  susa::matrix <int, susa::allocator_log<int>> mat_dagger("[0 1;"
                                                           "2 2;"
                                                           "4 2]");
  susa::slice <int, susa::slice_en::DAGGER, susa::allocator_log<int>> slc_m(mat_m, 2, 1);
  SUSA_TEST_EQ(mat_dagger, slc_m, "matrix and slice equality");
  }

  {
  susa::matrix <int, susa::allocator_log<int>> mat_dagger("[0 1;"
                                                           "2 3;"
                                                           "1 3]");
  susa::slice <int, susa::slice_en::DAGGER, susa::allocator_log<int>> slc_m(mat_m, 3, 2);
  SUSA_TEST_EQ(mat_dagger, slc_m, "matrix and slice equality");
  }

  SUSA_TEST_PRINT_STATS();

  return (uint_failed);
}
