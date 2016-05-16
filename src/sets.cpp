/*
 * This file is part of Susa.
 *
 * Susa is free software: you can redistribute it and/or modify
 * it under the terms of the Lesser GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * at your option) any later version.
 *
 * Susa is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * Lesser GNU General Public License for more details.
 *
 * You should have received a copy of the Lesser GNU General Public License
 * along with Susa.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file set.cpp
 * @brief Basic mathematical operations on STL and Susa types.
 * @author Behrooz, Aliabadi
 * @version 1.0.0
 */


#include <susa.h>

namespace susa
{

  index_set::index_set(unsigned int uint_num_nodes)
  {
    uint_set_size = uint_num_nodes;

    nbytes = uint_num_nodes / 8;
    nbytes += (uint_num_nodes % 8) ? 1 : 0;

    try
    {
      entity = new char[nbytes];
    }
    catch ( std::bad_alloc ex)
    {
      SUSA_ASSERT_MESSAGE(false, "memory allocation exception.");
      std::exit(EXIT_FAILURE); // assert may be disabled but we exit anyway !
    }
  }

  index_set::~index_set()
  {
    delete[] entity;
    entity = nullptr;
  }

  void index_set::add(unsigned int index)
  {

    if (index >= uint_set_size) return;

    unsigned int byte = index / 8;
    unsigned int bit  = index % 8;


    if (byte > nbytes) return;
    entity[byte] |= (0x01 << bit);
  }

  void index_set::remove(unsigned int index)
  {
    if (exists(index))
    {
      unsigned int byte = index / 8;
      unsigned int bit  = index % 8;
      if (byte > nbytes) return;
      entity[byte] ^= (0x01 << bit); 
    }
  }

  unsigned int index_set::pop()
  {
    for (unsigned int uint_index = 0; uint_index < uint_set_size; uint_index++)
    {
      if (exists(uint_index))
      {
        remove(uint_index);
        return uint_index;
      }
    }

    return uint_set_size;
  }

  void index_set::push(unsigned int index)
  {
    add(index);
  }

  bool index_set::exists(unsigned int index)
  {
    if (index >= uint_set_size) return false;
    unsigned int byte = index / 8;
    unsigned int bit  = index % 8;

    if (byte > nbytes) return false;
    return ((entity[byte] >> bit) & 0x01);
  }

  void index_set::add_all()
  {
    for (unsigned int index = 0; index < nbytes; index++)
    {
      unsigned int byte = index / 8;
      unsigned int bit  = index % 8;
      entity[byte] |= (0x01 << bit);
    }
  }

  void index_set::remove_all()
  {
    for (unsigned int indx = 0; indx < nbytes; indx++) entity[indx] = 0x00;
  }

  bool index_set::is_not_empty()
  {
    char empty = 0x00;
    for (unsigned int indx = 0; indx < nbytes; indx++) empty |= entity[indx];
    return empty;
  }

} // NAMESPACE SUSA
