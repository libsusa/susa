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
template <typename T = uintmax_t, typename Allocator = std::allocator<T>>
class bitset
{

  public :

    //! constructor
    bitset(size_t sizet_size);

    //! constructor
    bitset(const bitset& bset_arg);

    bitset() = delete;

    //! destructor
    ~bitset();

    /**
     * @brief set a single bit at the specified index
     *
     * @param index the bit index to be set
     */
    void set(size_t index);

    void reset(size_t index);

    size_t pop();

    void push(size_t index);

    bool exists(size_t index) const;

    void set();

    void reset();

    bool any();

    size_t size() const {return sz_size;}

    bitset&       operator=(const bitset& bset_arg);
    bool&         operator ()(size_t index);
    bool          operator ()(size_t index) const;

    friend bool   operator==(const bitset& bset_left, const bitset& bset_right)
    {
      if (bset_left.size() != bset_right.size()) return false;
      for (size_t index = 0; index < bset_left.nblocks; index++)
      {
        if (bset_left._matrix[index] != bset_right._matrix[index]) return false;
      }
      return true;
    }

  friend std::ostream &operator<<(std::ostream& out_stream, const bitset& bset_arg)
  {
    out_stream << "[ ";
    for (size_t index = 0; index < bset_arg.sz_size; index++)
    {
      out_stream << (bset_arg.exists(index) ? 1 : 0) << " ";
    }
    out_stream << "]";

    return out_stream;
  }

  private:

    bool allocate(size_t sz_blocks);

    Allocator   alloc;
    size_t      sz_size;
    unsigned    nbits;
    size_t      nblocks;
    T*          _matrix;

};

template <typename T, typename Allocator>
bitset<T, Allocator>::bitset(size_t sizet_size)
: alloc()
, sz_size(sizet_size)
{

  static_assert(std::is_unsigned<T>::value);

  nbits    = std::numeric_limits<T>::digits;
  nblocks  = sz_size / nbits;
  nblocks += (sz_size % nbits) ? 1 : 0;

  SUSA_ASSERT_MESSAGE(allocate(nblocks), "memory allocation failed.");

  reset();
}

template <typename T, typename Allocator>
bool bitset<T, Allocator>::allocate(size_t sz_blocks)
{
  try
  {
    _matrix = alloc.allocate(sz_blocks);
  }
  catch(const std::bad_alloc& e)
  {
    _matrix = nullptr;
    nblocks = 0;
    SUSA_LOG_ERR("memory allocation failed");
    return false;
  }

  return true;
}

template <typename T, typename Allocator>
bitset<T, Allocator>::~bitset()
{
  if (_matrix != nullptr)
  {
    alloc.deallocate(_matrix, nblocks);
  }
}

template <typename T, typename Allocator>
bitset<T, Allocator>::bitset(const bitset& bset_arg)
{
  nblocks         = bset_arg.nblocks;
  nbits           = bset_arg.nbits;
  sz_size         = bset_arg.sz_size;
  if (allocate(nblocks))
  {
    std::memcpy(_matrix, bset_arg._matrix, nblocks * sizeof(T));
  }
}

template <typename T, typename Allocator>
bitset<T, Allocator>& bitset<T, Allocator>::operator=(const bitset<T, Allocator>& bset_arg)
{
  if (_matrix != nullptr && nblocks != 0)
  {
    alloc.deallocate(_matrix, nblocks);
    _matrix = nullptr;
    nblocks = 0;
  }

  if (allocate(bset_arg.nblocks))
  {
    nblocks = bset_arg.nblocks;
    std::memcpy(_matrix, bset_arg._matrix, nblocks * sizeof(uintmax_t));
  }

  return *this;
}

template <typename T, typename Allocator>
void bitset<T, Allocator>::set(size_t index)
{
  SUSA_ASSERT_MESSAGE(index < sz_size, "index out of range");
  if (index >= sz_size) return;

  unsigned int byte = index / nbits;
  unsigned int bit  = index % nbits;

  _matrix[byte] |= (0x01 << bit);
}

template <typename T, typename Allocator>
void bitset<T, Allocator>::reset(size_t index)
{
  SUSA_ASSERT_MESSAGE(index < sz_size, "index out of range");
  if (index >= sz_size) return;

  if (exists(index))
  {
    unsigned int byte = index / nbits;
    unsigned int bit  = index % nbits;
    _matrix[byte] ^= (0x01 << bit);
  }
}

template <typename T, typename Allocator>
size_t bitset<T, Allocator>::pop()
{
  for (unsigned int uint_index = 0; uint_index < sz_size; uint_index++)
  {
    if (exists(uint_index))
    {
      reset(uint_index);
      return uint_index;
    }
  }

  return sz_size;
}

template <typename T, typename Allocator>
void bitset<T, Allocator>::push(size_t index)
{
  set(index);
}

template <typename T, typename Allocator>
bool bitset<T, Allocator>::exists(size_t index) const
{
  SUSA_ASSERT_MESSAGE(index < sz_size, "index out of range");
  if (index >= sz_size) return false;

  unsigned int byte = index / nbits;
  unsigned int bit  = index % nbits;

  return ((_matrix[byte] & (0x01 << bit)) > 0);
}

template <typename T, typename Allocator>
void bitset<T, Allocator>::set()
{
  for (unsigned int index = 0; index < nblocks; index++)
  {
    unsigned int byte = index / nbits;
    unsigned int bit  = index % nbits;
    _matrix[byte] |= (0x01 << bit);
  }
}

template <typename T, typename Allocator>
void bitset<T, Allocator>::reset()
{
  for (size_t indx = 0; indx < nblocks; indx++) _matrix[indx] = 0x00;
}

template <typename T, typename Allocator>
bool bitset<T, Allocator>::any()
{
  unsigned char empty = 0x00;
  for (size_t indx = 0; indx < nblocks; indx++) empty |= _matrix[indx];
  return (empty != 0);
}

template <typename T, typename Allocator>
bool bitset<T, Allocator>::operator()(size_t index) const
{
  return exists(index);
}

}       // NAMESPACE SUSA
#endif  // SUSA_SETS_H
