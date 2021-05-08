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
 * @file sets.h
 * @brief The <i>set</i> type (declaration and definition).
 * @author Behrooz Kamary
 *
 */

#ifndef SUSA_SETS_H
#define SUSA_SETS_H

#include <susa/debug.h>
#include <susa/memory.h>

namespace susa
{

/**
* @brief The <i>bitset</i> class.
* <i>bitset</i> is a set type aimed to store
* a set of bits in a quick (using single memory allocation)
* and memory efficient manner. You may alternatively use <i>std::bitset</i>
* where the set size is in hand at compile time (no runtime initialization).
*
* @ingroup TYPES
*
*/
template <typename Allocator = std::allocator<uintmax_t>>
class bitset
{

  public :

    //! constructor
    bitset(size_t sizet_size);

    //! constructor
    bitset(const bitset& bset_arg);

    //! destructor
    ~bitset();

    /**
     * @brief set a single bit at the specified index
     *
     * @param index the bit index
     */
    void set(size_t index);

    void reset(size_t index);

    size_t pop();

    void push(size_t index);

    bool exists(size_t index);

    void set();

    void reset();

    bool any();

    bitset& operator=(const bitset& bset_arg);

  private:

    Allocator   alloc;
    size_t      uint_set_size;
    unsigned    nbits;
    size_t      nblocks;
    uintmax_t*  _matrix;

};

template <typename Allocator>
bitset<Allocator>::bitset(size_t sizet_size)
: alloc()
, uint_set_size(sizet_size)
{

  nbits    = std::numeric_limits<uintmax_t>::digits;
  nblocks  = uint_set_size / nbits;
  nblocks += (uint_set_size % nbits) ? 1 : 0;

  try
  {
    _matrix = alloc.allocate(nblocks);
  }
  catch(const std::bad_alloc& e)
  {
    _matrix = nullptr;
    nblocks = 0;
    SUSA_LOG_ERR("memory allocation failed");
  }

  reset();
}

template <typename Allocator>
bitset<Allocator>::~bitset()
{
  if (_matrix != nullptr)
  {
    alloc.deallocate(_matrix, nblocks);
  }
}

template <typename Allocator>
bitset<Allocator>::bitset(const bitset& bset_arg)
{
  nblocks = bset_arg.nblocks;

  if (alloc.allocate(nblocks))
  {
    std::memcpy(_matrix, bset_arg._matrix, nblocks * sizeof(uintmax_t));
  }
}

template <typename Allocator>
bitset<Allocator>& bitset<Allocator>::operator=(const bitset<Allocator>& bset_arg)
{
  if (_matrix != nullptr && nblocks != 0)
  {
    alloc.deallocate(_matrix, nblocks);
    _matrix = nullptr;
    nblocks = 0;
  }

  if (alloc.allocate(bset_arg.nblocks))
  {
    nblocks = bset_arg.nblocks;
    std::memcpy(_matrix, bset_arg._matrix, nblocks * sizeof(uintmax_t));
  }

  return *this;
}

template <typename Allocator>
void bitset<Allocator>::set(size_t index)
{
  if (index >= uint_set_size) return;

  unsigned int byte = index / nbits;
  unsigned int bit  = index % nbits;

  if (byte > nblocks) return;
  _matrix[byte] |= (0x01 << bit);
}

template <typename Allocator>
void bitset<Allocator>::reset(size_t index)
{
  if (exists(index))
  {
    unsigned int byte = index / nbits;
    unsigned int bit  = index % nbits;
    if (byte > nblocks) return;
    _matrix[byte] ^= (0x01 << bit);
  }
}

template <typename Allocator>
size_t bitset<Allocator>::pop()
{
  for (unsigned int uint_index = 0; uint_index < uint_set_size; uint_index++)
  {
    if (exists(uint_index))
    {
      reset(uint_index);
      return uint_index;
    }
  }

  return uint_set_size;
}

template <typename Allocator>
void bitset<Allocator>::push(size_t index)
{
  set(index);
}

template <typename Allocator>
bool bitset<Allocator>::exists(size_t index)
{
  if (index >= uint_set_size) return false;
  unsigned int byte = index / nbits;
  unsigned int bit  = index % nbits;

  if (byte > nblocks) return false;
  return ((_matrix[byte] >> bit) & 0x01);
}

template <typename Allocator>
void bitset<Allocator>::set()
{
  for (unsigned int index = 0; index < nblocks; index++)
  {
    unsigned int byte = index / nbits;
    unsigned int bit  = index % nbits;
    _matrix[byte] |= (0x01 << bit);
  }
}

template <typename Allocator>
void bitset<Allocator>::reset()
{
  for (size_t indx = 0; indx < nblocks; indx++) this->_matrix[indx] = 0x00;
}

template <typename Allocator>
bool bitset<Allocator>::any()
{
  unsigned char empty = 0x00;
  for (size_t indx = 0; indx < nblocks; indx++) empty |= this->_matrix[indx];
  return (empty != 0);
}

}       // NAMESPACE SUSA
#endif  // SUSA_SETS_H
