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
 * @file array.h
 * @brief The multi-dimensional array type (declaration and definition).
 *
 * This file contains the <i>array template class</i>. This is a non-intrinsic
 * data type of Susa. Use this class if you need multi-dimensional arrays.
 * This class handles the cloning when the assignment operator is used.
 * It also has the necessary memory management mechanisms to release
 * the allocated memory.
 *
 * @defgroup TYPES Types
 *
 * @author Behrooz Kamary
 */

#ifndef SUSA_ARRAY_H
#define SUSA_ARRAY_H

#include <susa/debug.h>
#include <susa/memory.h>

namespace susa
{
  /**
   * @brief <i>array</i> type.
   * An array is a multidimensional container type that should be used
   * to represent three or more dimensions.
   *
   * @ingroup TYPES
   *
   */
  template <typename T, typename Allocator = std::allocator<T>>
  class array
  {

  public:
    //! Constructor
    array();

    /**
     * @brief Constructor
     * @param list initialization list
     */
    array(std::initializer_list<size_t> list);

    //! Constructor
    array(array &&arg);

    //! Copy constructor
    array(const array &arg);

    //! Copy assignment constructor
    array &operator=(const array &arg);

    //! Destructor
    ~array() noexcept;

    T get(std::initializer_list<size_t> list) const;

    T get(const std::vector<size_t> &list) const;

    /**
     * @brief Clone the data by pointer that may be read from the disk.
     * The method assumes that you have loaded a compatible binary data
     * into the memory that has the right format.
     *
     * @param data pointer to the data
     * @param length length of the data
     */
    void clone(T *data, size_t length);

    template <typename... Args>
    T get(size_t uint_elem, Args... uint_args);

    //! () operator to set or get elements
    template <typename... Args>
    T operator()(size_t uint_elem, Args... uint_args) const
    {
      return get(uint_elem, uint_args...);
    }

    template <typename... Args>
    T &operator()(size_t uint_elem, Args... uint_args)
    {
      size_t uint_index = index(uint_elem, uint_args...);

      SUSA_ASSERT_MESSAGE(uint_index < uint_total, "the element index is out of range.");

      if (uint_index < uint_total && this->_matrix != NULL)
      {
        return this->_matrix[uint_index];
      }

      return *T_default;
    }

    template <typename... Args>
    size_t get_raw_index(size_t uint_elem, Args... uint_args)
    {
      return index(uint_elem, uint_args...);
    }

  private:
    std::vector<size_t> vec_dims;
    size_t uint_total;
    size_t uint_num_dims;
    T *T_default;

    std::vector<size_t> vec_get;
    size_t uint_get_step;
    T T_get_value;

    Allocator alloc;
    T *_matrix;

    size_t index(size_t uint_elem);
    size_t index(const std::vector<size_t> &list);
    template <typename... Args>
    size_t index(size_t uint_elem, Args... uint_args);
  };

  // Implementations
  template <typename T, typename Allocator>
  array<T, Allocator>::array()
  : uint_total(0)
  , T_default(new T)
  , uint_get_step(0)
  , alloc()
  , _matrix(nullptr)
  {
  }

  template <typename T, typename Allocator>
  array<T, Allocator>::array(std::initializer_list<size_t> list)
  : vec_dims(list)
  , T_default(new T)
  , alloc()
  {
    uint_total = 1;
    for (auto dim_size : vec_dims)
    {
      uint_total *= dim_size;
    }

    uint_num_dims = vec_dims.size();
    uint_get_step = 0;
    vec_get = std::vector<size_t>(uint_num_dims);

    _matrix = alloc.allocate(uint_total);
  }

  template <typename T, typename Allocator>
  void array<T, Allocator>::clone(T *data, size_t length)
  {

    if (length != uint_total)
    {
      if (uint_total != 0 && _matrix != nullptr)
      {
        alloc.deallocate(_matrix, uint_total);
      }
      _matrix       = alloc.allocate(uint_total);
      uint_total    = length;
    }

    std::memset(_matrix, 0x00, sizeof(T) * length);
    std::memcpy(_matrix, data, sizeof(T) * length);
  }

  template <typename T, typename Allocator>
  array<T, Allocator>::array(array &&arg)
  : alloc(std::move(arg.alloc))
  {

    uint_total      = arg.uint_total;
    vec_dims        = arg.vec_dims;
    uint_num_dims   = vec_dims.size();
    vec_get         = std::vector<size_t>(uint_num_dims);
    uint_get_step   = 0;

    T_default          = arg.T_default;
    _matrix         = arg._matrix;

    arg.T_default      = nullptr;
    arg._matrix     = nullptr;
    arg.uint_total  = 0;
  }

  template <typename T, typename Allocator>
  array<T, Allocator>::array(const array &arg)
  : T_default(new T)
  , alloc(arg.alloc)
  {

    uint_total    = arg.uint_total;
    vec_dims      = arg.vec_dims;
    uint_num_dims = vec_dims.size();
    vec_get       = std::vector<size_t>(uint_num_dims);
    uint_get_step = 0;
    _matrix       = alloc.allocate(uint_total);
    std::memcpy(_matrix, arg._matrix, uint_total * sizeof(T));
  }

  template <typename T, typename Allocator>
  array<T, Allocator> &array<T, Allocator>::operator=(const array &arg)
  {
    if (uint_total != 0 && _matrix != nullptr)
    {
      alloc.deallocate(_matrix, uint_total);
    }
    uint_total    = arg.uint_total;
    vec_dims      = arg.vec_dims;
    uint_num_dims = vec_dims.size();
    vec_get       = std::vector<size_t>(uint_num_dims);
    uint_get_step = 0;

    _matrix       = alloc.allocate(uint_total);
    std::memcpy(_matrix, arg._matrix, uint_total * sizeof(T));

    return *this;
  }

  template <typename T, typename Allocator>
  array<T, Allocator>::~array() noexcept
  {
    if (T_default != nullptr) delete T_default;
    if (_matrix != nullptr) alloc.deallocate(_matrix, uint_total);
  }

  template <typename T, typename Allocator>
  T array<T, Allocator>::get(std::initializer_list<size_t> list) const
  {

    SUSA_ASSERT(_matrix != nullptr);

    size_t uint_elem      = 0;
    size_t uint_dim_count = 0;
    size_t uint_factor    = 1;

    for (auto dim_size : list)
    {

      uint_factor = 1;
      for (size_t indx = 0; indx < uint_dim_count; ++indx)
      {
        uint_factor *= vec_dims[indx];
      }

      uint_elem += dim_size * uint_factor;
      uint_dim_count++;
    }

    SUSA_ASSERT_MESSAGE(uint_elem < uint_total, "the element index is out of range.");

    if (uint_elem < uint_total && _matrix != nullptr)
    {
      return _matrix[uint_elem];
    }

    return *T_default;
  }

  template <typename T, typename Allocator>
  T array<T, Allocator>::get(const std::vector<size_t> &list) const
  {

    SUSA_ASSERT(_matrix != nullptr);

    size_t uint_index = index(list);

    SUSA_ASSERT_MESSAGE(uint_index < uint_total, "the element index is out of range.");

    if (uint_index < uint_total && _matrix != nullptr)
    {
      return _matrix[uint_index];
    }

    return *T_default;
  }

  template <typename T, typename Allocator>
  template <typename... Args>
  T array<T, Allocator>::get(size_t uint_elem, Args... uint_args)
  {

    SUSA_ASSERT(_matrix != nullptr);

    size_t uint_index = index(uint_elem, uint_args...);

    SUSA_ASSERT_MESSAGE(uint_index < uint_total, "the element index is out of range.");

    if (uint_index < uint_total && _matrix != nullptr)
    {
      return _matrix[uint_index];
    }

    return *T_default;
  }

  template <typename T, typename Allocator>
  size_t array<T, Allocator>::index(const std::vector<size_t> &list)
  {

    size_t uint_elem = 0;
    size_t uint_dim_count = 0;
    size_t uint_factor = 1;

    for (auto dim_size : list)
    {

      uint_factor = 1;
      for (size_t indx = 0; indx < uint_dim_count; ++indx)
      {
        uint_factor *= vec_dims[indx];
      }

      uint_elem += dim_size * uint_factor;
      uint_dim_count++;
    }

    SUSA_ASSERT_MESSAGE(uint_elem < uint_total, "the element index is out of range.");

    return uint_elem;
  }

  template <typename T, typename Allocator>
  size_t array<T, Allocator>::index(size_t uint_elem)
  {

    size_t ret;

    if (uint_get_step == 0)
    {
      vec_get       = std::vector<size_t>(uint_num_dims);
      uint_get_step = 0;

      SUSA_ASSERT_MESSAGE(uint_elem < uint_total, "the element index is out of range.");

      if (uint_elem < uint_total) return uint_elem;

      return 0;
    }
    else
    {

      SUSA_ASSERT_MESSAGE(uint_get_step < uint_num_dims, "the number of arguments exceeded the number of dimensions.");

      vec_get[uint_get_step] = uint_elem;
      uint_get_step++;

      ret = index(vec_get);
    }

    vec_get = std::vector<size_t>(uint_num_dims);
    uint_get_step = 0;

    return ret;
  }

  template <typename T, typename Allocator>
  template <typename... Args>
  size_t array<T, Allocator>::index(size_t uint_elem, Args... uint_args)
  {

    SUSA_ASSERT_MESSAGE(uint_get_step < uint_num_dims, "the number of arguments exceeded the number of dimensions.");

    vec_get[uint_get_step] = uint_elem;
    uint_get_step++;
    return index(uint_args...);
  }

} // namespace susa
#endif // SUSA_ARRAY_H
