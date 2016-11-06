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
 * @file sets.cpp
 * @brief Set types and operations.
 * @author Behrooz, Aliabadi
 * @version 1.0.0
 */


#include <susa.h>

namespace susa
{

  iset::iset(unsigned int uint_num_nodes)
  {
    uint_set_size = uint_num_nodes;

    nbytes  = uint_set_size / 8;
    nbytes += (uint_set_size % 8) ? 1 : 0;

    this->allocate(nbytes);
  }

  void iset::add(unsigned int index)
  {

    if (index >= uint_set_size) return;

    unsigned int byte = index / 8;
    unsigned int bit  = index % 8;


    if (byte > nbytes) return;
    this->_matrix[byte] |= (0x01 << bit);
  }

  void iset::remove(unsigned int index)
  {
    if (exists(index))
    {
      unsigned int byte = index / 8;
      unsigned int bit  = index % 8;
      if (byte > nbytes) return;
      this->_matrix[byte] ^= (0x01 << bit);
    }
  }

  unsigned int iset::pop()
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

  void iset::push(unsigned int index)
  {
    add(index);
  }

  bool iset::exists(unsigned int index)
  {
    if (index >= uint_set_size) return false;
    unsigned int byte = index / 8;
    unsigned int bit  = index % 8;

    if (byte > nbytes) return false;
    return ((this->_matrix[byte] >> bit) & 0x01);
  }

  void iset::add_all()
  {
    for (unsigned int index = 0; index < nbytes; index++)
    {
      unsigned int byte = index / 8;
      unsigned int bit  = index % 8;
      this->_matrix[byte] |= (0x01 << bit);
    }
  }

  void iset::remove_all()
  {
    for (unsigned int indx = 0; indx < nbytes; indx++) this->_matrix[indx] = 0x00;
  }

  bool iset::is_not_empty()
  {
    char empty = 0x00;
    for (unsigned int indx = 0; indx < nbytes; indx++) empty |= this->_matrix[indx];
    return empty;
  }

} // NAMESPACE SUSA
