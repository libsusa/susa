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
 * @file modulation.cpp
 * @brief Unit Test Suit
 * @author Behrooz Kamary
 */

#include "test.h"

int main(void)
{
  const unsigned  M = 64;
  const unsigned  l = 2000 * susa::log2u(M);
  susa::qam       mapper(M);
  susa::bitset<>  bset(l);
  bset.reset();


  susa::rng rnd(5786);
  auto binary = rnd.bernoulli(l);


  for (size_t index = 0; index < l; index++)
  {
    if (binary(index) > 0) bset.set(index);
  }

  auto mapped   = mapper.mod_bits(bset);
  auto bset_hat = mapper.demod_to_bits(mapped);


  SUSA_TEST_EQ(bset, bset_hat, "Gray coded modulation and demodulation");


  SUSA_TEST_PRINT_STATS();

  return (uint_failed);
}
