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
 * @author Behrooz Kamary
 * @version 1.0.0
 */

#ifndef SUSA_ARRAY_H
#define SUSA_ARRAY_H

#include <susa/debug.h>
#include <susa/memory.h>

namespace susa
{
  /**
   * @brief The <i>array</i> class.
   * An array is a multidemsional container that should be used
   * to represent three or more dimensions.
   * 
   * @ingroup TYPES
   *
   */
  template <class T> class array : public susa::memory <T>
  {

    public:

      //! Constructor
      array();

      /**
       * @brief Constructor
       * @param list initialization list
       */
      array (std::initializer_list<size_t> list);

      //! Constructor
      array(array&& arg);

      //! Copy constructor
      array(const array& arg);
      
      //! Copy assignment constructor
      array& operator=(const array <T>& arg);
      
      //! Destructor
      ~array() noexcept;

      T get(std::initializer_list<size_t> list) const;
    
      T get(const std::vector<size_t>& list) const;


      /**
       * @brief Clone the data by pointer that may be read from the disk.
       * The method assumes that you have loaded a compatible binary data
       * into the memory that has the right format.
       *
       * @param data pointer to the data
       */
      void clone(T* data);

      template <typename... Args> T get(size_t uint_elem, Args... uint_args);

      //! () operator to set or get elements
      template <typename... Args> T operator ()( size_t uint_elem, Args... uint_args ) const
      {
        return get(uint_elem, uint_args...);
      }

      template <typename... Args> T &operator ()( size_t uint_elem, Args... uint_args )
      {
        size_t uint_index = index(uint_elem, uint_args...);

        SUSA_ASSERT_MESSAGE(uint_index < uint_total, "the element index is out of range.");

        if (uint_index < uint_total && this->_matrix != NULL)
        {
          return this->_matrix[uint_index];
        }

        return *T_fake;
      }

      template <typename... Args> size_t get_raw_index(size_t uint_elem, Args... uint_args)
      {
        return index(uint_elem, uint_args...);
      }

    private:

      std::vector<size_t> vec_dims;
      size_t              uint_total;
      size_t              uint_num_dims;
      T*                  T_fake;

      std::vector<size_t> vec_get;
      size_t              uint_get_step;
      T                   T_get_value;

      size_t index(size_t uint_elem);
      size_t index(const std::vector<size_t>& list);
      template <typename... Args> size_t index(size_t uint_elem, Args... uint_args);
  };

  // Implementations
  template <class T> array<T>::array()
  : susa::memory<T>()
  , T_fake(new T)
  {
      uint_get_step = 0;
  }

  template <class T> array<T>::array (std::initializer_list<size_t> list)
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
    vec_get       = std::vector<size_t> (uint_num_dims);

    this->allocate(uint_total);
  }

  template <class T> void array<T>::clone(T* data)
  {
    std::memset(this->_matrix, 0x00, sizeof(T) * uint_total);
    std::memcpy(this->_matrix, data, sizeof(T) * uint_total);
  }

  template <class T> array<T>::array(array&& arg)
  : susa::memory <T> (std::move(arg))
  {

    uint_total        = arg.uint_total;
    vec_dims          = arg.vec_dims;
    uint_num_dims     = vec_dims.size();
    vec_get           = std::vector<size_t>(uint_num_dims);
    uint_get_step     = 0;

    T_fake            = arg.T_fake;
    arg.T_fake        = nullptr;

  }

  template <class T> array<T>::array(const array& arg)
  : susa::memory<T>(arg)
  , T_fake(new T)
  {
    uint_total        = arg.uint_total;
    vec_dims          = arg.vec_dims;
    uint_num_dims     = vec_dims.size();
    vec_get           = std::vector<size_t>(uint_num_dims);
    uint_get_step     = 0;
  }
  
   template <class T> array<T>& array<T>::operator=(const array <T>& arg)
   {
     
      susa::memory<T>::operator=(arg);

      uint_total        = arg.uint_total;
      vec_dims          = arg.vec_dims;
      uint_num_dims     = vec_dims.size();
      vec_get           = std::vector<size_t>(uint_num_dims);
      uint_get_step     = 0;

      return *this;
   }

  template <class T> array<T>::~array() noexcept
  {
    delete T_fake;
  }

  template <class T> T array <T>::get( std::initializer_list<size_t> list ) const
  {

    SUSA_ASSERT(this->_matrix != NULL);

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

    if (uint_elem < uint_total && this->_matrix != NULL)
    {
      return this->_matrix[uint_elem];
    }

    return *T_fake;
  }

  template <class T> T array <T>::get(const std::vector<size_t>& list ) const
  {

    SUSA_ASSERT(this->_matrix != NULL);

    size_t uint_index = index(list);

    SUSA_ASSERT_MESSAGE(uint_index < uint_total, "the element index is out of range.");

    if (uint_index < uint_total && this->_matrix != NULL)
    {
      return this->_matrix[uint_index];
    }

    return *T_fake;
  }

  template <class T>
  template <typename... Args>
  T array<T>::get(size_t uint_elem, Args... uint_args)
  {

    SUSA_ASSERT(this->_matrix != NULL);

    size_t uint_index = index(uint_elem, uint_args...);

    SUSA_ASSERT_MESSAGE(uint_index < uint_total, "the element index is out of range.");

    if (uint_index < uint_total && this->_matrix != NULL)
    {
      return this->_matrix[uint_index];
    }

    return *T_fake;
  }

  template <class T> size_t array <T>::index (const std::vector<size_t>& list)
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

    return uint_elem;
  }

  template <class T> size_t array<T>::index(size_t uint_elem)
  {

    size_t ret;

    if (uint_get_step == 0)
    {
      vec_get = std::vector<size_t>(uint_num_dims);
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

  template <class T>
  template <typename... Args>
  size_t array<T>::index(size_t uint_elem, Args... uint_args)
  {

    SUSA_ASSERT_MESSAGE(uint_get_step < uint_num_dims, "the number of arguments exceeded the number of dimensions.");

    vec_get[uint_get_step] = uint_elem;
    uint_get_step++;
    return index(uint_args...);
  }

}      // NAMESPACE SUSA
#endif // SUSA_ARRAY_H
