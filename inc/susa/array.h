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
 * @brief The multi-dimensional array type definition and declaration.
 *
 * This file contains the <i>array template class</i>. This is a non-intrinsic
 * data type of Susa. Use this class if you need multi-dimensional arrays.
 * This class handles the cloning when the assignment operator is used.
 * It also has the necessary memeory management mechanisms to release
 * the allocated memory.
 *
 * @defgroup TYPES Types
 *
 * @author Behrooz Kamary Aliabadi
 * @version 1.0.0
 */

#ifndef SUSA_ARRAY_H
#define SUSA_ARRAY_H

#include "debug.h"

namespace susa
{
  /**
   * @brief The <i>array</i> class.
   *
   * @ingroup TYPES
   *
   */
  template <class T> class array
  {

    public:
      /** Constructor
       * @param initialization list (since C++11)
       */
      array( std::initializer_list<unsigned int> list );

#ifdef HAS_MOVE_SEMANTICS
      array( array&& arg );
#endif

      ~array();

      T get( std::initializer_list<unsigned int> list ) const;
      T get(const std::vector<unsigned int>& list ) const;


      /** Clone
       * @brief Clone the data by pointer that may be read from the disk
       *
       * @param pointer to the data
       */
      void clone(T* data);

      template <typename... Args> T get(unsigned int uint_elem, Args... uint_args);

      //! () operator to set or get elements
      template <typename... Args> T operator ()( unsigned int uint_elem, Args... uint_args ) const
      {
        return get(uint_elem, uint_args...);
      }

      template <typename... Args> T &operator ()( unsigned int uint_elem, Args... uint_args )
      {
        unsigned int uint_index = index(uint_elem, uint_args...);

        SUSA_ASSERT_MESSAGE(uint_index < uint_total, "the element index is out of range.");

        if (uint_index < uint_total && _matrix != NULL)
        {
          return _matrix[uint_index];
        }

        return *T_fake;
      }

      template <typename... Args> unsigned int get_raw_index(unsigned int uint_elem, Args... uint_args)
      {
        return index(uint_elem, uint_args...);
      }

    private:
      void initialize(unsigned int bulk_size);

      std::vector<unsigned int> vec_dims;
      unsigned int uint_total;
      unsigned int uint_num_dims;
      T* _matrix;
      T* T_fake;

      unsigned int index( unsigned int uint_elem );
      unsigned int index(const std::vector<unsigned int>& list );
      template <typename... Args> unsigned int index( unsigned int uint_elem, Args... uint_args );

      std::vector<unsigned int> vec_get;
      unsigned int uint_get_step;
      T T_get_value;
  };

  // Implementations

  template <class T> array<T>::array( std::initializer_list<unsigned int> list )
    : vec_dims(list)
      , T_fake(new T)
  {
    uint_total = 1;
    for (auto dim_size : vec_dims)
    {
      uint_total *= dim_size;
    }

    uint_num_dims = vec_dims.size();
    uint_get_step = 0;
    vec_get = std::vector<unsigned int> (uint_num_dims);

    initialize(uint_total);
  }

  template <class T> void array<T>::clone(T* data)
  {
    std::memset(_matrix, 0x00, sizeof(T) * uint_total);
    std::memcpy(_matrix, data, sizeof(T) * uint_total);
  }

#ifdef HAS_MOVE_SEMANTICS
  template <class T> array<T>::array( array&& arg )
  {

    uint_total = arg.uint_total;
    vec_dims = arg.vec_dims;
    uint_num_dims = vec_dims.size();
    vec_get = std::vector<unsigned int>(uint_num_dims);
    uint_get_step = 0;

    _matrix = arg._matrix;
    arg._matrix = nullptr;

  }
#endif

  template <class T> array<T>::~array()
  {
    delete[] _matrix;
    delete T_fake;
  }

  template <class T> void array<T>::initialize( unsigned int bulk_size )
  {
    try
    {
      _matrix = new T [bulk_size];
    }
    catch ( std::bad_alloc ex)
    {
      SUSA_ASSERT_MESSAGE(false, "memory allocation exception.");
      std::exit(EXIT_FAILURE);
    }
  }

  template <class T> T array <T>::get( std::initializer_list<unsigned int> list ) const
  {

    SUSA_ASSERT(_matrix != NULL);

    unsigned int uint_elem      = 0;
    unsigned int uint_dim_count = 0;
    unsigned int uint_factor    = 1;

    for (auto dim_size : list) {

      uint_factor = 1;
      for (unsigned int indx = 0; indx < uint_dim_count; ++indx)
      {
        uint_factor *= vec_dims[indx];
      }

      uint_elem += dim_size * uint_factor;
      uint_dim_count++;
    }

    SUSA_ASSERT_MESSAGE(uint_elem < uint_total, "the element index is out of range.");

    if (uint_elem < uint_total && _matrix != NULL)
    {
      return _matrix[uint_elem];
    }

    return *T_fake;
  }

  template <class T> T array <T>::get(const std::vector<unsigned int>& list ) const
  {

    SUSA_ASSERT(_matrix != NULL);

    unsigned int uint_index = index(list);

    SUSA_ASSERT_MESSAGE(uint_index < uint_total, "the element index is out of range.");

    if (uint_index < uint_total && _matrix != NULL)
    {
      return _matrix[uint_index];
    }

    return *T_fake;
  }

  template <class T>
    template <typename... Args>
    T array<T>::get(unsigned int uint_elem, Args... uint_args)
    {

      SUSA_ASSERT(_matrix != NULL);

      unsigned int uint_index = index(uint_elem, uint_args...);

      SUSA_ASSERT_MESSAGE(uint_index < uint_total, "the element index is out of range.");

      if (uint_index < uint_total && _matrix != NULL)
      {
        return _matrix[uint_index];
      }

      return *T_fake;
    }

  template <class T> unsigned int array <T>::index(const std::vector<unsigned int>& list )
  {

    unsigned int uint_elem = 0;
    unsigned int uint_dim_count = 0;
    unsigned int uint_factor = 1;

    for (auto dim_size : list)
    {

      uint_factor = 1;
      for (unsigned int indx = 0; indx < uint_dim_count; ++indx)
      {
        uint_factor *= vec_dims[indx];
      }

      uint_elem += dim_size * uint_factor;
      uint_dim_count++;
    }

    return uint_elem;
  }

  template <class T> unsigned int array<T>::index(unsigned int uint_elem)
  {

    unsigned int ret;

    if (uint_get_step == 0)
    {
      vec_get = std::vector<unsigned int>(uint_num_dims);
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

    vec_get = std::vector<unsigned int>(uint_num_dims);
    uint_get_step = 0;

    return ret;
  }

  template <class T>
    template <typename... Args>
    unsigned int array<T>::index(unsigned int uint_elem, Args... uint_args)
    {

      SUSA_ASSERT_MESSAGE(uint_get_step < uint_num_dims, "the number of arguments exceeded the number of dimensions.");

      vec_get[uint_get_step] = uint_elem;
      uint_get_step++;
      return index(uint_args...);
    }

}      // NAMESPACE SUSA
#endif // SUSA_ARRAY_H
