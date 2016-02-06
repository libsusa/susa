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
 * along with Susa.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file mmatrix.h
 * @brief The multi-dimensional matrix class.
 * This file contains the <i>mmatrix template class</i>. This is a non-intrinsic
 * data type of Susa. Use this class if you need multi dimensional matrices.
 * This class handles the cloning of the matrices when the assignment operator is used.
 * It has the necessary memeory management mechanisms to release the allocated memory.
 * @author Behrooz Kamary Aliabadi
 * @version 1.0.0
 */

#ifndef MMATRIX_H
#define MMATRIX_H

#include "debug.h"

namespace susa
{

template <class T> class mmatrix
{

  public:
    /** Constructor
     * @param initialization list (since C++11)
     */
    mmatrix( std::initializer_list<unsigned int> list );

#ifdef HAS_MOVE_SEMANTICS
    mmatrix( mmatrix&& arg );
#endif

  ~mmatrix();

  T get( std::initializer_list<unsigned int> list ) const;
  T get(const std::vector<unsigned int>& list ) const;

  template <typename... Args> T get(unsigned int uint_elem, Args... uint_args);

  //! () operator to set or get elements
  template <typename... Args> T operator ()( unsigned int uint_elem, Args... uint_args ) const
  {
    return get(uint_elem, uint_args...);
  }

  template <typename... Args> T &operator ()( unsigned int uint_elem, Args... uint_args )
  {
    unsigned int uint_index = index(uint_elem, uint_args...);

    SUSA_ASSERT_MESSAGE(uint_index < uint_total, "the element index is out of range");

    if (uint_index < uint_total && _matrix != NULL)
     {
      return _matrix[uint_index];
    }

    return *T_fake;
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

template <class T> mmatrix<T>::mmatrix( std::initializer_list<unsigned int> list )
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

#ifdef HAS_MOVE_SEMANTICS
template <class T> mmatrix<T>::mmatrix( mmatrix&& arg )
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

template <class T> mmatrix<T>::~mmatrix()
{
  delete[] _matrix;
  delete T_fake;
}

template <class T> void mmatrix<T>::initialize( unsigned int bulk_size )
{
  try {
      _matrix = new T [bulk_size];
  } catch ( std::bad_alloc ex) {
      std::exit(EXIT_FAILURE);
  }
}

template <class T> T mmatrix <T>::get( std::initializer_list<unsigned int> list ) const
{

  SUSA_ASSERT(_matrix != NULL);

  unsigned int uint_elem = 0;
  unsigned int uint_dim_count = 0;
  unsigned int uint_factor = 1;

  for (auto dim_size : list) {

    uint_factor = 1;
    for (unsigned int indx = 0; indx < uint_dim_count; ++indx) {
      uint_factor *= vec_dims[indx];
    }

    uint_elem += dim_size * uint_factor;
    uint_dim_count++;
  }

  SUSA_ASSERT_MESSAGE(uint_elem < uint_total, "the element index is out of range");

  if (uint_elem < uint_total && _matrix != NULL) {
    return _matrix[uint_elem];
  }

  return *T_fake;
}

template <class T> T mmatrix <T>::get(const std::vector<unsigned int>& list ) const
{

  SUSA_ASSERT(_matrix != NULL);

  unsigned int uint_index = index(list);

  SUSA_ASSERT_MESSAGE(uint_index < uint_total, "the element index is out of range");

  if (uint_index < uint_total && _matrix != NULL) {
    return _matrix[uint_index];
  }

  return *T_fake;
}

template <class T>
template <typename... Args>
T mmatrix<T>::get(unsigned int uint_elem, Args... uint_args)
{

  SUSA_ASSERT(_matrix != NULL);

  unsigned int uint_index = index(uint_elem, uint_args...);

  SUSA_ASSERT_MESSAGE(uint_index < uint_total, "the element index is out of range");

  if (uint_index < uint_total && _matrix != NULL) {
    return _matrix[uint_index];
  }

  return *T_fake;
}

template <class T> unsigned int mmatrix <T>::index(const std::vector<unsigned int>& list )
{

    unsigned int uint_elem = 0;
    unsigned int uint_dim_count = 0;
    unsigned int uint_factor = 1;

    for (auto dim_size : list) {

      uint_factor = 1;
      for (unsigned int indx = 0; indx < uint_dim_count; ++indx) {
        uint_factor *= vec_dims[indx];
      }

      uint_elem += dim_size * uint_factor;
      uint_dim_count++;
    }

    return uint_elem;
}

template <class T> unsigned int mmatrix<T>::index(unsigned int uint_elem)
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
unsigned int mmatrix<T>::index(unsigned int uint_elem, Args... uint_args)
{

  SUSA_ASSERT_MESSAGE(uint_get_step < uint_num_dims, "the number of arguments exceeded the number of dimensions.");

  vec_get[uint_get_step] = uint_elem;
  uint_get_step++;
  return index(uint_args...);
}

}      // NAMESPACE SUSA
#endif // MMATRIX_H
