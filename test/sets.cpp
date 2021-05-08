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
 * @file sets.cpp
 * @brief Unit Test Suit
 * @author Behrooz Kamary
 */

#include "test.h"
#include <susa/sets.h>

int main(void)
{

  susa::bitset<> bset(10);
  bset.set(6);
  SUSA_TEST_EQ(bset.exists(6), true, "bitset exists");
  SUSA_TEST_EQ(bset.exists(5), false, "bitset exists");
  SUSA_TEST_EQ(bset.exists(7), false, "bitset exists");


  SUSA_TEST_PRINT_STATS();

  return (uint_failed);
}
