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
 * @file matrix.h
 * @brief The matrix type definition and declaration.
 *
 * This file contains the <i>matrix template class</i>. This is a non-intrinsic data type that is
 * used in Susa. Use this class if you need two dimensional matrices or vectors. It has been
 * designed such that the matrix elements can be accessed as vector elements. This class handles
 * the cloning of the matrices when the assignment operator is used.
 * It has the necessary memeory management mechanisms to release the allocated memory.
 *
 * @author Behrooz Aliabadi
 * @version 1.0.0
 */

#ifndef SUSA_MATRIX_H
#define SUSA_MATRIX_H

#include <debug.h>
#include <memory.h>

namespace susa
{

void pre_parser(std::string& str_string);

template <class T> class matrix;

template <class T>  matrix <T> operator-( const matrix <T>&, T );
template <class T>  matrix <T> operator-( T, const matrix <T>& );
template <class T>  matrix <T> operator-( const matrix <T>&, const matrix <T>& );

template <class T>  matrix <T> operator+( const matrix <T>&, T );
template <class T>  matrix <T> operator+( T, const matrix <T>& );
template <class T>  matrix <T> operator+( const matrix <T>&,const matrix <T>& );

template <class T>  matrix <T> operator*( const matrix <T>&, T );
template <class T>  matrix <T> operator*( T, const matrix <T>& );
template <class T>  matrix <T> operator*( const matrix <T>&,const matrix <T>& );

template <class T>  matrix <T> operator/( const matrix <T>&, T );
template <class T>  matrix <T> operator/( T, const matrix <T>& );
template <class T>  matrix <T> operator/( const matrix <T>&,const matrix <T>& );

template <class T> std::ostream &operator<<(std::ostream &, const matrix <T> &);
template <class T> std::istream &operator>>(std::istream &, matrix <T> &);


template <class T> matrix <T> matmul(const matrix <T> &mat_argl,const matrix <T> &mat_argr);
template <class T> matrix <T> transpose(const matrix <T> &mat_arg);

/**
 * @brief The <i>matrix</i> class.
 *
 * @ingroup TYPES
 *
 */
template <class T> class matrix :
  public susa::memory <T>
{
  private:
    unsigned int uint_rows;
    unsigned int uint_cols;

    T* T_fake;

    // UTILITY METHODS

    // This parses the strings and initialize the matrix elements.
    // The string initialization can be done using a constructor method
    // and an overloaded assignment operator (=) method
    void parser(std::string str_string);

  public:
    //! The default constructor
    matrix();

    /**
     * Constructor
     *
     * @param uint_rows Number of rows
     * @param uint_cols Number of columns
     * @param Tinitial Initial value of all rows and columns.
     */
    matrix( unsigned int uint_rows, unsigned int uint_cols, T Tinitial );

    /**
     * Constructor
     *
     * @param uint_rows Number of rows
     * @param uint_cols Number of columns.
     */
    matrix( unsigned int uint_rows, unsigned int uint_cols );

    //! Copy constructor
    matrix(const matrix <T> &mat_arg);

#ifdef HAS_MOVE_SEMANTICS
    //! Copy Constructor for rvalues
    matrix(matrix <T> &&mat_arg);
#endif
    //! Destructor
    ~matrix();

    /**
     * Constructor
     *
     * @param str_string string representation of the matrix
     */
    matrix(std::string str_string);

    //! Returns the value of a specific (row, column)
    T get( unsigned int uint_row, unsigned int uint_col ) const;

    //! Returns the value of a specific (elem)
    T get( unsigned int uint_elem ) const;


    //! Returns the number of columns
    unsigned int no_cols() const;

    //! Returns the number of rows
    unsigned int no_rows() const;

    //! Returns true if the matrix is square
    bool is_square() const;

    /** 
     * @brief Returns true if it represents a vector
     *
     * A susa::matrix <T> instance can be considered a <i>vector</i>
     * if and only if it has a single row or a single column.
     */
    bool is_vector() const;

    /** 
     * @brief Returns true if it represents a scalar
     *
     * A susa::matrix <T> instance can be considered a <i>scalar</i>
     * if and only if it has a single row and a single column.
     */
    bool is_scalar() const;

    //! Returns the size of matrix
    void set_all(T T_arg);

    //! Returns the indicated row
    matrix <T> row(unsigned int uint_row) const;

    //! Returns the indicated column
    matrix <T> col(unsigned int uint_col) const;

    //! Returns a matrix after eliminating the specified column and the row
    matrix <T> shrink(unsigned int uint_row, unsigned int uint_col) const;

    //! Resize the matrix to the newly specified dimensions
    bool resize(unsigned int uint_row, unsigned int uint_col);

    //! Considers the matrix object as a vector and return left side of that vector
    matrix <T> left(unsigned int uint_left) const;

    //! Considers the matrix object as a vector and return right side of that vector
    matrix <T> right(unsigned int uint_right) const;

    //! Considers the matrix object as a vector and return mid part of that vector
    matrix <T> mid(unsigned int uint_begin, unsigned int uint_end) const;

    //! Element wise Assignment by Addition operator
    matrix <T> operator+=( const matrix <T> &mat_arg );

    //! Element wise Assignment by Subtraction
    matrix <T> operator-=( const matrix <T> &mat_arg );

    //! Element wise Assignment operator
    matrix <T>& operator=( const matrix <T> &mat_arg );

    //! Element wise Assignment operator
    matrix <T> operator=( const matrix <T> &mat_arg ) const;

    //! Element wise Assignment operator
    matrix <T>& operator=( std::string str_string );

    //! Element wise Subtraction operator
    friend matrix <T> operator-<>( const matrix <T> &mat_argl, T T_arg);

    //! Element wise Subtraction operator
    friend matrix <T> operator-<>( T T_arg, const matrix <T> &mat_argr );

    //! Element wise Subtraction operator
    friend matrix <T> operator-<>( const matrix <T> &mat_argl, const matrix <T> &mat_argr);

    //! Element wise Addition operator
    friend matrix <T> operator+<>( const matrix <T> &mat_argl, T T_arg );

    //! Element wise Addition operator
    friend matrix <T> operator+<>( T T_arg, const matrix <T> &mat_argr );

    //! Element wise Addition operator
    friend matrix <T> operator+<>( const matrix <T> &mat_argl, const matrix <T> &mat_argr);

    //! Element wise Multiplication  operator
    friend matrix <T> operator*<>( const matrix <T> &mat_argl, T T_arg);

    //! Element wise Multiplication  operator
    friend matrix <T> operator*<>( T T_arg, const matrix <T> &mat_argr );

    //! Element wise Multiplication  operator
    friend matrix <T> operator*<>( const matrix <T> &mat_argl, const matrix <T> &mat_argr );


    //! Element wise Division operator
    friend matrix <T> operator/<>( const matrix <T> &mat_argl, T T_arg);

    //! Element wise Division operator
    friend matrix <T> operator/<>( T T_arg, const matrix <T> &mat_argr );

    //! Element wise Division operator
    friend matrix <T> operator/<>( const matrix <T> &mat_argl, const matrix <T> &mat_argr);

    //! Output stream
    friend std::ostream &operator<< <>(std::ostream &outStream, const matrix <T> &mat_arg);

    //! Input stream
    friend std::istream &operator>> <>(std::istream &inStream, matrix <T> &mat_arg);

    //! () operator to set or get elements
    T &operator ()( unsigned int uint_row, unsigned int uint_col );

    T operator ()( unsigned int uint_row, unsigned int uint_col ) const;

    //! () operator to set or get elements
    T &operator ()( unsigned int uint_elem);

    T operator ()( unsigned int uint_elem) const;

    //! Typecasting
    operator matrix <double> ();

    operator matrix <float> ();

    operator matrix <int> ();

    operator matrix <char> ();

    operator matrix <std::complex <double> > ();

    operator matrix <std::complex <float> > ();

    operator matrix <std::complex <int> > ();

    operator matrix <std::complex <char> > ();

    // Friend methods are used to speed up matrix operations.
    // Because they can access '_matrix' directly.
    // They are defined in 'linalg.h'
    friend matrix <T> matmul <> ( const matrix <T> &mat_argl, const matrix <T> &mat_argr);
    friend matrix <T> transpose <> (const matrix <T> &mat_arg);

    friend bool operator!=( const susa::matrix <T> &mat_argl, const susa::matrix <T> &mat_argr)
    {
      for (unsigned int uint_indx = 0; uint_indx < mat_argl.uint_objects; uint_indx++)
      {
          if (mat_argl._matrix[uint_indx] != mat_argr._matrix[uint_indx]) return true;
      }

      return false;
    }

    friend bool operator==( const susa::matrix <T> &mat_argl, const susa::matrix <T> &mat_argr)
    {
      for (unsigned int uint_indx = 0; uint_indx < mat_argl.uint_objects; uint_indx++)
      {
          if (mat_argl._matrix[uint_indx] != mat_argr._matrix[uint_indx]) return false;
      }

      return true;
    }
};

// Constructor and Destructor

template <class T> matrix <T>::matrix()
: T_fake(new T)
{
  uint_rows = 0;
  uint_cols = 0;
}

template <class T> matrix <T>::matrix(unsigned int uint_rows, unsigned int uint_cols, T Tinitial)
: T_fake(new T)
{

  SUSA_ASSERT(uint_cols > 0 && uint_rows > 0);

  this->uint_rows = uint_rows < 2 ? 1 : uint_rows;
  this->uint_cols = uint_cols < 2 ? 1 : uint_cols;

  this->allocate(this->uint_rows * this->uint_cols);

  for (unsigned int uint_index = 0; uint_index < this->uint_objects; uint_index++)
  {
      this->_matrix[uint_index] = Tinitial;
  }

}

template <class T> matrix <T>::matrix( unsigned int uint_rows, unsigned int uint_cols )
: T_fake(new T)
{
    SUSA_ASSERT(uint_cols > 0 && uint_rows > 0);

    this->uint_rows = uint_rows < 2 ? 1 : uint_rows;
    this->uint_cols = uint_cols < 2 ? 1 : uint_cols;

    this->allocate(this->uint_rows * this->uint_cols);
}

template <class T> matrix <T>::matrix(const matrix <T> &mat_arg)
: T_fake(new T)
{
    unsigned int uint_size = mat_arg.uint_rows * mat_arg.uint_cols;

    if (uint_size != 0)
    {
        this->allocate(uint_size);

        uint_rows = mat_arg.uint_rows;
        uint_cols = mat_arg.uint_cols;

        for (unsigned int uint_index = 0; uint_index < uint_size; uint_index++)
        {
            this->_matrix[uint_index] = mat_arg._matrix[uint_index];
        }
    }
    else
    {
        uint_rows = 0;
        uint_cols = 0;
        this->deallocate();
    }
}

#ifdef HAS_MOVE_SEMANTICS

template <class T> matrix <T>::matrix(matrix&& mat_arg) {

  uint_rows          = mat_arg.uint_rows;
  uint_cols          = mat_arg.uint_cols;
  this->uint_objects = mat_arg.uint_objects;
  this->uint_bytes   = mat_arg.uint_bytes;

  this->_matrix      = mat_arg._matrix;
  mat_arg._matrix    = nullptr;
}

#endif

template <class T> matrix <T>::matrix(std::string str_string)
: T_fake(new T)
{
  uint_rows     = 0;
  uint_cols     = 0;
  this->_matrix = NULL;

  parser(str_string);
}

template <class T> matrix <T>::~matrix()
{
    delete T_fake;
}

// Public methods

template <class T> T matrix <T>::get( unsigned int uint_row, unsigned int uint_col ) const
{

  SUSA_ASSERT(this->_matrix != NULL);

  SUSA_ASSERT_MESSAGE(uint_row < uint_rows && uint_col < uint_cols, "one or more indices are out of range.");

  return this->_matrix[uint_row + uint_col * uint_rows];

}

template <class T> T matrix <T>::get( unsigned int uint_elem ) const
{

  SUSA_ASSERT(this->_matrix != NULL);

  SUSA_ASSERT_MESSAGE(uint_elem < this->uint_objects, "the element index is out of range.");

  return this->_matrix[uint_elem];

}

template <class T> unsigned int  matrix <T>::no_cols() const
{
    return uint_cols;
}

template <class T> unsigned int  matrix <T>::no_rows() const
{
    return uint_rows;
}

template <class T> bool matrix <T>::is_square() const
{
    return (uint_rows == uint_cols);
}

template <class T> bool matrix <T>::is_vector() const
{
    return ((uint_rows == 1 && uint_cols > 1 ) || ( uint_rows > 1 && uint_cols == 1));
}

template <class T> bool matrix <T>::is_scalar() const
{
    return (uint_rows == 1 && uint_cols == 1);
}

template <class T> void  matrix <T>::set_all(T T_arg)
{
    for (unsigned int uint_counter = 0; uint_counter < this->uint_objects; uint_counter++) {
      this->_matrix[uint_counter] = T_arg;
    }
}

template <class T> matrix <T> matrix <T>::row(unsigned int uint_row) const
{
    matrix <T> mat_ret(1, uint_cols);

    if (uint_row < uint_rows)
    {
        for (unsigned int uint_i = 0; uint_i < uint_cols; uint_i++) {
          mat_ret(uint_i) = this->_matrix[uint_i * uint_rows + uint_row];
        }
    }

    return mat_ret;
}

template <class T> matrix <T> matrix <T>::col(unsigned int uint_col) const
{
    matrix <T> mat_ret(uint_rows,1);

    if (uint_col < uint_cols)
    {
        for (unsigned int uint_i = 0; uint_i < uint_rows; uint_i++)
        {
          mat_ret(uint_i) = this->_matrix[uint_col * uint_rows + uint_i];
        }
    }

    return mat_ret;
}

template <class T> matrix <T> matrix <T>::shrink(unsigned int uint_elim_row, unsigned int uint_elim_col) const
{

    matrix <T> mat_ret;

    unsigned int uint_new_col;
    unsigned int uint_new_row;

    //TODO: it should be verified; this condition may be incorrect
    SUSA_ASSERT_MESSAGE(uint_cols > 1 && uint_rows > 1
      && uint_elim_row >= 0 && uint_elim_col >= 0,
      "the input arguments error.");

    if (uint_cols > 1 && uint_rows > 1 && uint_elim_row >= 0 && uint_elim_col >= 0)
    {
        mat_ret = matrix <T> (uint_rows - 1,uint_cols - 1);
    }
    else
    {
        return *this;
    }


    SUSA_ASSERT_MESSAGE(uint_elim_col < uint_cols && uint_elim_row < uint_rows,
      "the input arguments exceed matrix size.");

    if (uint_elim_col < uint_cols && uint_elim_row < uint_rows)
    {
        for ( unsigned int uint_row = 0; uint_row < uint_rows; uint_row++ )
        {
            for ( unsigned int uint_col = 0; uint_col < uint_cols; uint_col++ )
            {
                if (uint_col != uint_elim_col || uint_row != uint_elim_row)
                {
                    uint_new_row = uint_row > uint_elim_row ? (uint_row - 1) : uint_row;
                    uint_new_col = uint_col > uint_elim_col ? (uint_col - 1) : uint_col;
                    if (uint_new_row < (uint_rows - 1) && uint_new_col < (uint_cols - 1))
                    {
                        mat_ret(uint_new_row, uint_new_col) = this->_matrix[uint_col * uint_rows + uint_row];
                    }
                }
            }
        }
    }
    else
    {
        return *this;
    }

    return mat_ret;
}

template <class T> bool matrix <T>::resize( unsigned int uint_rows, unsigned int uint_cols )
{
    SUSA_ASSERT(uint_cols > 0 && uint_rows > 0);

    this->uint_rows = uint_rows < 2 ? 1 : uint_rows;
    this->uint_cols = uint_cols < 2 ? 1 : uint_cols;

    this->allocate(this->uint_rows * this->uint_cols);

    // The design is to exit if the memory allocation
    // in susa::memory<T> fails. Thus the returned value
    // is not taken into account even if it is being tested.
    return true;
}

template <class T> matrix <T> matrix <T>::left(unsigned int uint_left) const
{
    matrix <T> mat_ret;

    SUSA_ASSERT(this->uint_objects >= uint_left);

    if (this->uint_objects >= uint_left )
    {
        if (uint_rows == 1 && uint_cols != 1) mat_ret = matrix <T> (1,uint_left);
        else if (uint_rows != 1 && uint_cols == 1) mat_ret = matrix <T> (uint_left,1);
        else if (uint_rows != 1 && uint_cols != 1) mat_ret = matrix <T> (uint_left,1);

        for (unsigned int uint_i = 0; uint_i < uint_left; uint_i++)
        {
          mat_ret(uint_i) = this->_matrix[uint_i];
        }
    }

    return mat_ret;
}


template <class T> matrix <T> matrix <T>::right(unsigned int uint_right) const
{
    matrix <T> mat_ret;

    SUSA_ASSERT(this->uint_objects >= uint_right);

    if (this->uint_objects >= uint_right )
    {
        if (uint_rows == 1 && uint_cols != 1) mat_ret = matrix <T> (1,uint_right);
        else if (uint_rows != 1 && uint_cols == 1) mat_ret = matrix <T> (uint_right,1);
        else if (uint_rows != 1 && uint_cols != 1) mat_ret = matrix <T> (uint_right,1);

        for (unsigned int uint_i = this->uint_objects; uint_i > (this->uint_objects - uint_right); uint_i--)
        {
            mat_ret(uint_right - this->uint_objects + uint_i - 1) = this->_matrix[uint_i - 1];
        }
    }

    return mat_ret;
}


template <class T> matrix <T> matrix <T>::mid(unsigned int uint_begin, unsigned int uint_end) const
{
    matrix <T> mat_ret;
    unsigned int uint_mid = uint_end - uint_begin;

    SUSA_ASSERT(this->uint_objects > uint_begin && this->uint_objects > uint_end && uint_end > uint_begin);

    if (this->uint_objects > uint_begin && this->uint_objects > uint_end && uint_end > uint_begin)
    {
        if (uint_rows == 1 && uint_cols != 1) mat_ret = matrix <T> (1,uint_mid + 1);
        else if (uint_rows != 1 && uint_cols == 1) mat_ret = matrix <T> (uint_mid + 1,1);
        else if (uint_rows != 1 && uint_cols != 1) mat_ret = matrix <T> (uint_mid + 1,1);

        for (unsigned int uint_i = uint_begin; uint_i <= uint_end; uint_i++)
        {
          mat_ret(uint_i - uint_begin) = this->_matrix[uint_i];
        }
    }

    return mat_ret;
}

// Operators

//  ()
template <class T> T &matrix<T>::operator ()( unsigned int uint_row, unsigned int uint_col )
{

  SUSA_ASSERT(this->_matrix != NULL);

  SUSA_ASSERT_MESSAGE(uint_row < uint_rows && uint_col < uint_cols, "one or more indices are out of range.");

  if (uint_row < uint_rows && uint_col < uint_cols && this->_matrix != NULL)
  {
      return this->_matrix[uint_row + uint_col * uint_rows];
  }

  return *T_fake;
}

template <class T> T matrix<T>::operator ()( unsigned int uint_row, unsigned int uint_col ) const
{
  return get(uint_row, uint_col);
}

template <class T> T &matrix<T>::operator ()( unsigned int uint_elem )
{

    SUSA_ASSERT(this->_matrix != NULL);

    SUSA_ASSERT_MESSAGE(uint_elem < this->uint_objects, "the index is out of range.");

    if (uint_elem < uint_cols * uint_rows && this->_matrix != NULL) {
      return this->_matrix[uint_elem];
    }

    return *T_fake;
}

template <class T> T matrix<T>::operator ()( unsigned int uint_elem ) const
{
  return get(uint_elem);
}

// Typecasting
template <class T> matrix<T>::operator matrix <double> ()
{

    matrix <double> mat_ret(uint_rows, uint_cols);

    for (unsigned int uint_elem = 0; uint_elem < this->uint_objects; uint_elem++) {
        mat_ret(uint_elem) = (double)this->_matrix[uint_elem];
    }

    return mat_ret;
}

template <class T> matrix<T>::operator matrix <float> ()
{

    matrix <float> mat_ret(uint_rows, uint_cols);

    for (unsigned int uint_elem = 0; uint_elem < this->uint_objects; uint_elem++) {
        mat_ret(uint_elem) = (float)this->_matrix[uint_elem];
    }
    return mat_ret;
}

template <class T> matrix<T>::operator matrix <int> ()
{

    matrix <int> mat_ret(uint_rows, uint_cols);

    for (unsigned int uint_elem = 0; uint_elem < this->uint_objects; uint_elem++) {
        mat_ret(uint_elem) = (int)this->_matrix[uint_elem];
    }
    return mat_ret;
}


template <class T> matrix<T>::operator matrix <char> ()
{

    matrix <char> mat_ret(uint_rows, uint_cols);

    for (unsigned int uint_elem = 0; uint_elem < this->uint_objects; uint_elem++)
    {
        mat_ret(uint_elem) = (char)this->_matrix[uint_elem];
    }
    return mat_ret;
}


template <class T> matrix<T>::operator matrix <std::complex <double> > ()
{

    matrix <std::complex <double> > mat_ret(uint_rows, uint_cols);

    for (unsigned int uint_elem = 0; uint_elem < this->uint_objects; uint_elem++)
    {
        mat_ret(uint_elem) = std::complex <double> ((double)this->_matrix[uint_elem].real(), (double)this->_matrix[uint_elem].imag());
    }
    return mat_ret;
}

template <class T> matrix<T>::operator matrix <std::complex <float> > ()
{

    matrix <std::complex <float> > mat_ret(uint_rows, uint_cols);

    for (unsigned int uint_elem = 0; uint_elem < this->uint_objects; uint_elem++)
    {
        mat_ret(uint_elem) = std::complex <float> ((float)this->_matrix[uint_elem].real(), (float)this->_matrix[uint_elem].imag());
    }
    return mat_ret;
}

template <class T> matrix<T>::operator matrix <std::complex <int> > ()
{

    matrix <std::complex <int> > mat_ret(uint_rows, uint_cols);

    for (unsigned int uint_elem = 0; uint_elem < this->uint_objects; uint_elem++)
    {
        mat_ret(uint_elem) = std::complex <int> ((int)this->_matrix[uint_elem].real(), (int)this->_matrix[uint_elem].imag());
    }
    return mat_ret;
}



template <class T> matrix<T>::operator matrix <std::complex <char> > ()
{

    matrix <std::complex <char> > mat_ret(uint_rows, uint_cols);

    for (unsigned int uint_elem = 0; uint_elem < this->uint_objects; uint_elem++)
    {
        mat_ret(uint_elem) = std::complex <char> ((char)this->_matrix[uint_elem].real(), (char)this->_matrix[uint_elem].imag());
    }
    return mat_ret;
}


//  +
template <class T> matrix<T> operator+(const matrix <T> &mat_argl, const matrix <T> &mat_argr)
{

    matrix <T> mat_temp(mat_argl.uint_rows, mat_argl.uint_cols);
    unsigned int uint_size = mat_argl.uint_rows * mat_argl.uint_cols;

    SUSA_ASSERT_MESSAGE((mat_argl.uint_rows == mat_argr.uint_rows)
      && (mat_argl.uint_cols == mat_argr.uint_cols),
      "the matrices have different sizes.");

    if ((mat_argl.uint_rows == mat_argr.uint_rows) && (mat_argl.uint_cols == mat_argr.uint_cols))
    {
        for ( unsigned int uint_index = 0; uint_index < uint_size; uint_index++ )
        {
            mat_temp._matrix[uint_index] = mat_argl._matrix[uint_index] + mat_argr._matrix[uint_index];
        }
    }

    return mat_temp;
}

template <class T> matrix<T> operator+( const matrix <T> &mat_argl, T T_arg )
{

    matrix <T> mat_temp(mat_argl.uint_rows, mat_argl.uint_cols);
    unsigned int uint_size = mat_argl.uint_rows * mat_argl.uint_cols;

    for ( unsigned int uint_index = 0; uint_index < uint_size; uint_index++ )
    {
        mat_temp._matrix[uint_index] = mat_argl._matrix[uint_index] + T_arg;
    }

    return mat_temp;
}

template <class T> matrix<T> operator+( T T_arg, const matrix <T> &mat_argr )
{

    matrix <T> mat_temp(mat_argr.uint_rows, mat_argr.uint_cols);
    unsigned int uint_size = mat_argr.uint_rows * mat_argr.uint_cols;

    for ( unsigned int uint_index = 0; uint_index < uint_size; uint_index++ )
    {
        mat_temp._matrix[uint_index] = mat_argr._matrix[uint_index] + T_arg;
    }

    return mat_temp;
}



//  -
template <class T> matrix<T> operator-( const matrix <T> &mat_argl, const matrix <T> &mat_argr)
{

    matrix <T> mat_temp(mat_argl.uint_rows, mat_argl.uint_cols);
    unsigned int uint_size = mat_argl.uint_rows * mat_argl.uint_cols;

    SUSA_ASSERT_MESSAGE((mat_argl.uint_rows == mat_argr.uint_rows)
      && (mat_argl.uint_cols == mat_argr.uint_cols),
      "the matrices have different sizes.");

    if ((mat_argl.uint_rows == mat_argr.uint_rows) && (mat_argl.uint_cols == mat_argr.uint_cols))
    {
        for ( unsigned int uint_index = 0; uint_index < uint_size; uint_index++ )
        {
            mat_temp._matrix[uint_index] = mat_argl._matrix[uint_index] - mat_argr._matrix[uint_index];
        }
    }

    return mat_temp;
}

template <class T> matrix<T> operator-( const matrix <T> &mat_argl, T T_arg )
{

    matrix <T> mat_temp(mat_argl.uint_rows, mat_argl.uint_cols);
    unsigned int uint_size = mat_argl.uint_rows * mat_argl.uint_cols;

    for ( unsigned int uint_index = 0; uint_index < uint_size; uint_index++ )
    {
        mat_temp._matrix[uint_index] = mat_argl._matrix[uint_index] - T_arg;
    }

    return mat_temp;
}

template <class T> matrix<T> operator-( T T_arg, const matrix <T> &mat_argr )
{

    matrix <T> mat_temp(mat_argr.uint_rows, mat_argr.uint_cols);
    unsigned int uint_size = mat_argr.uint_rows * mat_argr.uint_cols;

    for ( unsigned int uint_index = 0; uint_index < uint_size; uint_index++ )
    {
        mat_temp._matrix[uint_index] = T_arg - mat_argr._matrix[uint_index];
    }

    return mat_temp;
}

//  *
template <class T> matrix<T> operator*( const matrix <T> &mat_argl, const matrix <T> &mat_argr)
{

    matrix <T> mat_temp(mat_argl.uint_rows, mat_argl.uint_cols);
    unsigned int uint_size = mat_argl.uint_rows * mat_argl.uint_cols;

    SUSA_ASSERT_MESSAGE((mat_argl.uint_rows == mat_argr.uint_rows)
      && (mat_argl.uint_cols == mat_argr.uint_cols),
      "the matrices have different sizes.");

    if ((mat_argl.uint_rows == mat_argr.uint_rows) && (mat_argl.uint_cols == mat_argr.uint_cols))
    {
        for ( unsigned int uint_index = 0; uint_index < uint_size; uint_index++ )
        {
            mat_temp._matrix[uint_index] = mat_argl._matrix[uint_index] * mat_argr._matrix[uint_index];
        }
    }

    return mat_temp;
}

template <class T> matrix<T> operator*( const matrix <T> &mat_argl, T T_arg )
{

    matrix <T> mat_temp(mat_argl.uint_rows, mat_argl.uint_cols);
    unsigned int uint_size = mat_argl.uint_rows * mat_argl.uint_cols;

    for ( unsigned int uint_index = 0; uint_index < uint_size; uint_index++ )
    {
        mat_temp._matrix[uint_index] = mat_argl._matrix[uint_index] * T_arg;
    }

    return mat_temp;
}

template <class T> matrix<T> operator*( T T_arg, const matrix <T> &mat_argr )
{

    matrix <T> mat_temp(mat_argr.uint_rows, mat_argr.uint_cols);
    unsigned int uint_size = mat_argr.uint_rows * mat_argr.uint_cols;

    for ( unsigned int uint_index = 0; uint_index < uint_size; uint_index++ )
    {
        mat_temp._matrix[uint_index] = T_arg * mat_argr._matrix[uint_index];
    }

    return mat_temp;
}

//  /
template <class T> matrix<T> operator/( const matrix <T> &mat_argl, const matrix <T> &mat_argr)
{

    matrix <T> mat_temp(mat_argl.uint_rows, mat_argl.uint_cols);
    unsigned int uint_size = mat_argl.uint_rows * mat_argl.uint_cols;

    SUSA_ASSERT_MESSAGE((mat_argl.uint_rows == mat_argr.uint_rows)
      && (mat_argl.uint_cols == mat_argr.uint_cols),
      "the matrices have different sizes.");

    if ((mat_argl.uint_rows == mat_argr.uint_rows) && (mat_argl.uint_cols == mat_argr.uint_cols))
    {
        for ( unsigned int uint_index = 0; uint_index < uint_size; uint_index++ )
        {
            mat_temp._matrix[uint_index] = mat_argl._matrix[uint_index] / mat_argr._matrix[uint_index];
        }
    }

    return mat_temp;
}

template <class T> matrix<T> operator/( const matrix <T> &mat_argl, T T_arg ) {

    matrix <T> mat_temp(mat_argl.uint_rows, mat_argl.uint_cols);
    unsigned int uint_size = mat_argl.uint_rows * mat_argl.uint_cols;

    for ( unsigned int uint_index = 0; uint_index < uint_size; uint_index++ ) {
        mat_temp._matrix[uint_index] = mat_argl._matrix[uint_index] / T_arg;
    }

    return mat_temp;
}

template <class T> matrix<T> operator/( T T_arg, const matrix <T> &mat_argr ) {

    matrix <T> mat_temp(mat_argr.uint_row, mat_argr.uint_col);
    //
    // NOT IMPLEMENTED YET
    //
    SUSA_ASSERT_MESSAGE(false, "NOT IMPLEMENTED YET.");
    return mat_temp;
}

//  -=
template <class T> matrix<T> matrix <T>::operator-=(const matrix <T> &mat_arg)
{

    SUSA_ASSERT(this->_matrix != NULL);

    SUSA_ASSERT_MESSAGE((mat_arg.uint_rows == uint_rows)
      && (mat_arg.uint_cols == uint_cols),
      "the matrices have different sizes.");

    if ((mat_arg.uint_rows == uint_rows) && (mat_arg.uint_cols == uint_cols) && this->_matrix != NULL)
    {
        for (unsigned int uint_index = 0; uint_index < this->uint_objects; uint_index++)
        {
            this->_matrix[uint_index] -= mat_arg._matrix[uint_index];
        }
    }

    return *this;
}

//  +=
template <class T> matrix<T> matrix <T>::operator+=(const matrix <T> &mat_arg)
{

    SUSA_ASSERT(this->_matrix != NULL);

    SUSA_ASSERT_MESSAGE((mat_arg.uint_rows == uint_rows)
      && (mat_arg.uint_cols == uint_cols),
      "the matrices have different sizes.");

    if ((mat_arg.uint_rows == uint_rows) && (mat_arg.uint_cols == uint_cols) && this->_matrix != NULL)
    {
        for (unsigned int uint_index = 0; uint_index < this->uint_objects; uint_index++)
        {
            this->_matrix[uint_index] += mat_arg._matrix[uint_index];
        }
    }

    return *this;
}

//  =
template <class T> matrix<T>& matrix <T>::operator=( const matrix <T> &mat_arg )
{

    unsigned int uint_size = mat_arg.uint_cols * mat_arg.uint_rows;

    if (this != &mat_arg)
    {
        if (uint_size != 0)
        {

            if (this->uint_objects != uint_size)
            {
                this->allocate(uint_size);
            }

            uint_rows = mat_arg.uint_rows;
            uint_cols = mat_arg.uint_cols;

            for (unsigned int uint_index = 0; uint_index < uint_size; uint_index++)
            {
                this->_matrix[uint_index] = mat_arg._matrix[uint_index];
            }
        }
        else
        {
            uint_rows = 0;
            uint_cols = 0;
            this->deallocate();
        }
    }


    return *this;
}

template <class T> matrix<T> matrix <T>::operator=( const matrix <T> &mat_arg ) const
{

    unsigned int uint_size = mat_arg.uint_cols * mat_arg.uint_rows;

    if (this != &mat_arg)
    {
        if (uint_size != 0)
        {

            if ((uint_rows * uint_cols) != uint_size)
            {
                this->allocate(uint_size);
            }

            uint_rows = mat_arg.uint_rows;
            uint_cols = mat_arg.uint_cols;

            for (unsigned int uint_index = 0; uint_index < uint_size; uint_index++)
            {
                this->_matrix[uint_index] = mat_arg._matrix[uint_index];
            }
        }
        else
        {
            uint_rows = 0;
            uint_cols = 0;

            this->deallocate();

        }
    }


    return *this;
}

template <class T> matrix<T>& matrix <T>::operator=( std::string str_string )
{

    parser(str_string);
    return *this;
}


template <class T> std::ostream &operator<<(std::ostream &outStream, const matrix <T> &mat_arg) {
    outStream << "[";
    for (unsigned int uint_row = 0; uint_row < mat_arg.uint_rows; uint_row++) {
        for (unsigned int uint_col = 0; uint_col < mat_arg.uint_cols; uint_col++) {
            outStream << mat_arg._matrix[uint_row + uint_col * mat_arg.uint_rows];
            if (uint_col < mat_arg.uint_cols - 1) outStream << " ";
        }
        if (uint_row < mat_arg.uint_rows - 1) outStream << "\n ";
    }
    outStream << "]";
    return outStream;
}


template <class T> void matrix <T>::parser(std::string str_string)
{

    pre_parser(str_string);

    unsigned int uint_size = 0;           // linear size of the matrix (#rows * #columns)
    int int_length = str_string.length(); // length of the input string
    unsigned int uint_cols_ = 0;          // number of columns
    unsigned int uint_rows_ = 0;          // number of rows

    // This loop gets total number of rows and columns
    for (int int_i = 0; int_i < int_length; int_i++)
    {
        switch (str_string[int_i])
        {
        case ';':
            uint_rows_++;
            // since the number of rows are counted by number of space and the
            // spaces around ';' are cut, we also count a column here !
            uint_cols_++;
            break;

        case 0xA:
            uint_rows_++;
            break;

        case 0x20:
            uint_cols_++;
            break;

        default:
            break;
        };
    }

    uint_cols_++;
    uint_rows_++;
    uint_size = uint_cols_ * uint_rows_;

    SUSA_ASSERT_MESSAGE(uint_cols_ % uint_rows_ == 0,
      "the number of columns are not equal in each row.");

    if (uint_cols_ % uint_rows_ != 0)
    {
        return;
    }


    if (uint_size != 0)
    {

        if ((this->uint_objects) != uint_size)
        {
            this->allocate(uint_size);
        }

        uint_rows = uint_rows_;
        uint_cols = uint_cols_ / uint_rows_;


        std::stringstream ss_all(str_string);
        char* char_buff = (char *)std::malloc(4096);
        SUSA_ASSERT_MESSAGE(char_buff != nullptr, "memory allocation failed.");
        if (char_buff == nullptr) std::exit(EXIT_FAILURE);
        T T_tmp;

        for (unsigned int uint_row = 0; uint_row < uint_rows; uint_row++)
        {
            ss_all.getline(char_buff, 4096, ';');
            std::stringstream ssrow(char_buff);
            for (unsigned int uint_col = 0; uint_col < uint_cols; uint_col++)
            {
                ssrow.getline(char_buff, 4096, ' ');
                if (!(std::istringstream(char_buff) >> T_tmp)) T_tmp = 0;
                this->_matrix[uint_row + uint_col * uint_rows] = T_tmp;
            }
        }

        std::free(char_buff);

    }
    else
    {
        uint_rows = 0;
        uint_cols = 0;

        this->deallocate();
    }
}
}      // NAMESPACE SUSA

#endif // SUSA_MATRIX_H
