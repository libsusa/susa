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

#include <susa/debug.h>
#include <susa/memory.h>

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
    size_t  sizet_rows;
    size_t  sizet_cols;

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
     * @param sizet_rows Number of rows
     * @param sizet_cols Number of columns
     * @param Tinitial Initial value of all rows and columns.
     */
    matrix( size_t sizet_rows, size_t sizet_cols, T Tinitial );

    /**
     * Constructor
     *
     * @param sizet_rows Number of rows
     * @param sizet_cols Number of columns.
     */
    matrix( size_t sizet_rows, size_t sizet_cols );

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
    T get( size_t sizet_row, size_t sizet_col ) const;

    //! Returns the value of a specific (elem)
    T get( size_t sizet_elem ) const;


    //! Returns the number of columns
    size_t no_cols() const;

    //! Returns the number of rows
    size_t no_rows() const;

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
    matrix <T> row(size_t sizet_row) const;

    //! Returns the indicated column
    matrix <T> col(size_t sizet_col) const;

    //! Returns a matrix after eliminating the specified column and the row
    matrix <T> shrink(size_t sizet_row, size_t sizet_col) const;

    //! Resize the matrix to the newly specified dimensions
    bool resize(size_t sizet_row, size_t sizet_col);

    //! Considers the matrix object as a vector and return left side of that vector
    matrix <T> left(size_t sizet_left) const;

    //! Considers the matrix object as a vector and return right side of that vector
    matrix <T> right(size_t sizet_right) const;

    //! Considers the matrix object as a vector and return mid part of that vector
    matrix <T> mid(size_t sizet_begin, size_t sizet_end) const;

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
    T &operator ()( size_t sizet_row, size_t sizet_col );

    T operator ()( size_t sizet_row, size_t sizet_col ) const;

    //! () operator to set or get elements
    T &operator ()( size_t sizet_elem);

    T operator ()( size_t sizet_elem) const;

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
      for (size_t sizet_indx = 0; sizet_indx < mat_argl.sizet_objects; sizet_indx++)
      {
          if (mat_argl._matrix[sizet_indx] != mat_argr._matrix[sizet_indx]) return true;
      }

      return false;
    }

    friend bool operator==( const susa::matrix <T> &mat_argl, const susa::matrix <T> &mat_argr)
    {
      for (size_t sizet_indx = 0; sizet_indx < mat_argl.sizet_objects; sizet_indx++)
      {
          if (mat_argl._matrix[sizet_indx] != mat_argr._matrix[sizet_indx]) return false;
      }

      return true;
    }
};

// Constructor and Destructor

template <class T> matrix <T>::matrix()
: T_fake(new T)
{
  sizet_rows = 0;
  sizet_cols = 0;
}

template <class T> matrix <T>::matrix(size_t sizet_rows, size_t sizet_cols, T Tinitial)
: T_fake(new T)
{

  SUSA_ASSERT(sizet_cols > 0 && sizet_rows > 0);

  this->sizet_rows = sizet_rows < 2 ? 1 : sizet_rows;
  this->sizet_cols = sizet_cols < 2 ? 1 : sizet_cols;

  this->allocate(this->sizet_rows * this->sizet_cols);

  for (size_t sizet_index = 0; sizet_index < this->sizet_objects; sizet_index++)
  {
      this->_matrix[sizet_index] = Tinitial;
  }

}

template <class T> matrix <T>::matrix( size_t sizet_rows, size_t sizet_cols )
: T_fake(new T)
{
    SUSA_ASSERT(sizet_cols > 0 && sizet_rows > 0);

    this->sizet_rows = sizet_rows < 2 ? 1 : sizet_rows;
    this->sizet_cols = sizet_cols < 2 ? 1 : sizet_cols;

    this->allocate(this->sizet_rows * this->sizet_cols);
}

template <class T> matrix <T>::matrix(const matrix <T> &mat_arg)
: T_fake(new T)
{
    size_t sizet_size = mat_arg.sizet_objects;

    if (sizet_size != 0)
    {
        this->allocate(sizet_size);

        sizet_rows = mat_arg.sizet_rows;
        sizet_cols = mat_arg.sizet_cols;

        for (size_t sizet_index = 0; sizet_index < sizet_size; sizet_index++)
        {
            this->_matrix[sizet_index] = mat_arg._matrix[sizet_index];
        }
    }
    else
    {
        sizet_rows = 0;
        sizet_cols = 0;
        this->deallocate();
    }
}

#ifdef HAS_MOVE_SEMANTICS

template <class T> matrix <T>::matrix(matrix&& mat_arg) {

  sizet_rows          = mat_arg.sizet_rows;
  sizet_cols          = mat_arg.sizet_cols;
  this->sizet_objects = mat_arg.sizet_objects;
  this->sizet_bytes   = mat_arg.sizet_bytes;

  this->_matrix      = mat_arg._matrix;
  mat_arg._matrix    = nullptr;
}

#endif

template <class T> matrix <T>::matrix(std::string str_string)
: T_fake(new T)
{
  sizet_rows     = 0;
  sizet_cols     = 0;
  this->_matrix = NULL;

  parser(str_string);
}

template <class T> matrix <T>::~matrix()
{
    delete T_fake;
}

// Public methods

template <class T> T matrix <T>::get( size_t sizet_row, size_t sizet_col ) const
{

  SUSA_ASSERT(this->_matrix != NULL);

  SUSA_ASSERT_MESSAGE(sizet_row < sizet_rows && sizet_col < sizet_cols, "one or more indices is/are out of range.");

  return this->_matrix[sizet_row + sizet_col * sizet_rows];

}

template <class T> T matrix <T>::get( size_t sizet_elem ) const
{

  SUSA_ASSERT(this->_matrix != NULL);

  SUSA_ASSERT_MESSAGE(sizet_elem < this->sizet_objects, "the element index is out of range.");

  return this->_matrix[sizet_elem];

}

template <class T> size_t  matrix <T>::no_cols() const
{
    return sizet_cols;
}

template <class T> size_t  matrix <T>::no_rows() const
{
    return sizet_rows;
}

template <class T> bool  matrix <T>::is_square() const
{
    return (sizet_rows == sizet_cols);
}

template <class T> bool  matrix <T>::is_vector() const
{
    return ((sizet_rows == 1 && sizet_cols > 1 ) || ( sizet_rows > 1 && sizet_cols == 1));
}

template <class T> bool  matrix <T>::is_scalar() const
{
    return (sizet_rows == 1 && sizet_cols == 1);
}

template <class T> void  matrix <T>::set_all(T T_arg)
{
    for (size_t sizet_counter = 0; sizet_counter < this->sizet_objects; sizet_counter++)
    {
      this->_matrix[sizet_counter] = T_arg;
    }
}

template <class T> matrix <T> matrix <T>::row(size_t sizet_row) const
{
    matrix <T> mat_ret(1, sizet_cols);

    if (sizet_row < sizet_rows)
    {
        for (size_t sizet_i = 0; sizet_i < sizet_cols; sizet_i++)
        {
          mat_ret(sizet_i) = this->_matrix[sizet_i * sizet_rows + sizet_row];
        }
    }

    return mat_ret;
}

template <class T> matrix <T> matrix <T>::col(size_t sizet_col) const
{
    matrix <T> mat_ret(sizet_rows, 1);

    if (sizet_col < sizet_cols)
    {
        size_t sizet_fix = sizet_col * sizet_rows;

        for (size_t sizet_i = 0; sizet_i < sizet_rows; sizet_i++)
        {
          mat_ret(sizet_i) = this->_matrix[sizet_fix + sizet_i];
        }
    }

    return mat_ret;
}

template <class T> matrix <T> matrix <T>::shrink(size_t sizet_elim_row, size_t sizet_elim_col) const
{

    matrix <T> mat_ret;

    size_t sizet_new_col;
    size_t sizet_new_row;

    //TODO: it should be verified; this condition may be incorrect
    SUSA_ASSERT_MESSAGE(sizet_cols > 1 && sizet_rows > 1
      && sizet_elim_row >= 0 && sizet_elim_col >= 0,
      "the input arguments error.");

    if (sizet_cols > 1 && sizet_rows > 1 && sizet_elim_row >= 0 && sizet_elim_col >= 0)
    {
        mat_ret = matrix <T> (sizet_rows - 1,sizet_cols - 1);
    }
    else
    {
        return *this;
    }


    SUSA_ASSERT_MESSAGE(sizet_elim_col < sizet_cols && sizet_elim_row < sizet_rows,
      "the input arguments exceed matrix size.");

    if (sizet_elim_col < sizet_cols && sizet_elim_row < sizet_rows)
    {
        for ( size_t sizet_row = 0; sizet_row < sizet_rows; sizet_row++ )
        {
            for ( size_t sizet_col = 0; sizet_col < sizet_cols; sizet_col++ )
            {
                if (sizet_col != sizet_elim_col || sizet_row != sizet_elim_row)
                {
                    sizet_new_row = sizet_row > sizet_elim_row ? (sizet_row - 1) : sizet_row;
                    sizet_new_col = sizet_col > sizet_elim_col ? (sizet_col - 1) : sizet_col;
                    if (sizet_new_row < (sizet_rows - 1) && sizet_new_col < (sizet_cols - 1))
                    {
                        mat_ret(sizet_new_row, sizet_new_col) =
                        this->_matrix[sizet_col * sizet_rows + sizet_row];
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

template <class T> bool matrix <T>::resize( size_t sizet_rows, size_t sizet_cols )
{
    SUSA_ASSERT(sizet_cols > 0 && sizet_rows > 0);

    this->sizet_rows = sizet_rows < 2 ? 1 : sizet_rows;
    this->sizet_cols = sizet_cols < 2 ? 1 : sizet_cols;

    this->allocate(this->sizet_rows * this->sizet_cols);

    // The design is to exit if the memory allocation
    // in susa::memory<T> fails. Thus the returned value
    // is not taken into account even if it is being tested.
    return true;
}

template <class T> matrix <T> matrix <T>::left(size_t sizet_left) const
{
    matrix <T> mat_ret;

    SUSA_ASSERT(this->sizet_objects >= sizet_left);

    if (this->sizet_objects >= sizet_left )
    {
        if (sizet_rows == 1 && sizet_cols != 1) mat_ret = matrix <T> (1,sizet_left);
        else if (sizet_rows != 1 && sizet_cols == 1) mat_ret = matrix <T> (sizet_left,1);
        else if (sizet_rows != 1 && sizet_cols != 1) mat_ret = matrix <T> (sizet_left,1);

        for (size_t sizet_i = 0; sizet_i < sizet_left; sizet_i++)
        {
          mat_ret(sizet_i) = this->_matrix[sizet_i];
        }
    }

    return mat_ret;
}


template <class T> matrix <T> matrix <T>::right(size_t sizet_right) const
{
    matrix <T> mat_ret;

    SUSA_ASSERT(this->sizet_objects >= sizet_right);

    if (this->sizet_objects >= sizet_right )
    {
        if (sizet_rows == 1 && sizet_cols != 1) mat_ret = matrix <T> (1,sizet_right);
        else if (sizet_rows != 1 && sizet_cols == 1) mat_ret = matrix <T> (sizet_right,1);
        else if (sizet_rows != 1 && sizet_cols != 1) mat_ret = matrix <T> (sizet_right,1);

        for (size_t sizet_i = this->sizet_objects; sizet_i > (this->sizet_objects - sizet_right); sizet_i--)
        {
            mat_ret(sizet_right - this->sizet_objects + sizet_i - 1) = this->_matrix[sizet_i - 1];
        }
    }

    return mat_ret;
}


template <class T> matrix <T> matrix <T>::mid(size_t sizet_begin, size_t sizet_end) const
{
    matrix <T> mat_ret;
    size_t sizet_mid = sizet_end - sizet_begin;

    SUSA_ASSERT(this->sizet_objects > sizet_begin && this->sizet_objects > sizet_end && sizet_end > sizet_begin);

    if (this->sizet_objects > sizet_begin && this->sizet_objects > sizet_end && sizet_end > sizet_begin)
    {
        if (sizet_rows == 1 && sizet_cols != 1) mat_ret = matrix <T> (1,sizet_mid + 1);
        else if (sizet_rows != 1 && sizet_cols == 1) mat_ret = matrix <T> (sizet_mid + 1,1);
        else if (sizet_rows != 1 && sizet_cols != 1) mat_ret = matrix <T> (sizet_mid + 1,1);

        for (size_t sizet_i = sizet_begin; sizet_i <= sizet_end; sizet_i++)
        {
          mat_ret(sizet_i - sizet_begin) = this->_matrix[sizet_i];
        }
    }

    return mat_ret;
}

// Operators

//  ()
template <class T> T &matrix<T>::operator ()( size_t sizet_row, size_t sizet_col )
{

  SUSA_ASSERT(this->_matrix != NULL);

  SUSA_ASSERT_MESSAGE(sizet_row < sizet_rows && sizet_col < sizet_cols,
                      "one or more indices is/are out of range.");

  if (sizet_row < sizet_rows && sizet_col < sizet_cols && this->_matrix != NULL)
  {
      return this->_matrix[sizet_row + sizet_col * sizet_rows];
  }

  return *T_fake;
}

template <class T> T matrix<T>::operator ()( size_t sizet_row, size_t sizet_col ) const
{
  return get(sizet_row, sizet_col);
}

template <class T> T &matrix<T>::operator ()( size_t sizet_elem )
{

    SUSA_ASSERT(this->_matrix != NULL);

    SUSA_ASSERT_MESSAGE(sizet_elem < this->sizet_objects, "the index is out of range.");

    if (sizet_elem < sizet_cols * sizet_rows && this->_matrix != NULL) {
      return this->_matrix[sizet_elem];
    }

    return *T_fake;
}

template <class T> T matrix<T>::operator ()( size_t sizet_elem ) const
{
  return get(sizet_elem);
}

// Typecasting
template <class T> matrix<T>::operator matrix <double> ()
{

    matrix <double> mat_ret(sizet_rows, sizet_cols);

    for (size_t sizet_elem = 0; sizet_elem < this->sizet_objects; sizet_elem++)
    {
        mat_ret(sizet_elem) = (double)this->_matrix[sizet_elem];
    }

    return mat_ret;
}

template <class T> matrix<T>::operator matrix <float> ()
{

    matrix <float> mat_ret(sizet_rows, sizet_cols);

    for (size_t sizet_elem = 0; sizet_elem < this->sizet_objects; sizet_elem++)
    {
        mat_ret(sizet_elem) = (float)this->_matrix[sizet_elem];
    }
    return mat_ret;
}

template <class T> matrix<T>::operator matrix <int> ()
{

    matrix <int> mat_ret(sizet_rows, sizet_cols);

    for (size_t sizet_elem = 0; sizet_elem < this->sizet_objects; sizet_elem++) {
        mat_ret(sizet_elem) = (int)this->_matrix[sizet_elem];
    }
    return mat_ret;
}


template <class T> matrix<T>::operator matrix <char> ()
{

    matrix <char> mat_ret(sizet_rows, sizet_cols);

    for (size_t sizet_elem = 0; sizet_elem < this->sizet_objects; sizet_elem++)
    {
        mat_ret(sizet_elem) = (char)this->_matrix[sizet_elem];
    }
    return mat_ret;
}


template <class T> matrix<T>::operator matrix <std::complex <double> > ()
{

    matrix <std::complex <double> > mat_ret(sizet_rows, sizet_cols);

    for (size_t sizet_elem = 0; sizet_elem < this->sizet_objects; sizet_elem++)
    {
        mat_ret(sizet_elem) = std::complex <double> ((double)this->_matrix[sizet_elem].real(), (double)this->_matrix[sizet_elem].imag());
    }
    return mat_ret;
}

template <class T> matrix<T>::operator matrix <std::complex <float> > ()
{

    matrix <std::complex <float> > mat_ret(sizet_rows, sizet_cols);

    for (size_t sizet_elem = 0; sizet_elem < this->sizet_objects; sizet_elem++)
    {
        mat_ret(sizet_elem) = std::complex <float> ((float)this->_matrix[sizet_elem].real(), (float)this->_matrix[sizet_elem].imag());
    }
    return mat_ret;
}

template <class T> matrix<T>::operator matrix <std::complex <int> > ()
{

    matrix <std::complex <int> > mat_ret(sizet_rows, sizet_cols);

    for (size_t sizet_elem = 0; sizet_elem < this->sizet_objects; sizet_elem++)
    {
        mat_ret(sizet_elem) = std::complex <int> ((int)this->_matrix[sizet_elem].real(), (int)this->_matrix[sizet_elem].imag());
    }
    return mat_ret;
}



template <class T> matrix<T>::operator matrix <std::complex <char> > ()
{

    matrix <std::complex <char> > mat_ret(sizet_rows, sizet_cols);

    for (size_t sizet_elem = 0; sizet_elem < this->sizet_objects; sizet_elem++)
    {
        mat_ret(sizet_elem) = std::complex <char> ((char)this->_matrix[sizet_elem].real(), (char)this->_matrix[sizet_elem].imag());
    }
    return mat_ret;
}


//  +
template <class T> matrix<T> operator+(const matrix <T> &mat_argl, const matrix <T> &mat_argr)
{

    matrix <T> mat_temp(mat_argl.sizet_rows, mat_argl.sizet_cols);
    size_t sizet_size = mat_argl.sizet_rows * mat_argl.sizet_cols;

    SUSA_ASSERT_MESSAGE((mat_argl.sizet_rows == mat_argr.sizet_rows)
      && (mat_argl.sizet_cols == mat_argr.sizet_cols),
      "the matrices have different sizes.");

    if ((mat_argl.sizet_rows == mat_argr.sizet_rows) && (mat_argl.sizet_cols == mat_argr.sizet_cols))
    {
        for ( size_t sizet_index = 0; sizet_index < sizet_size; sizet_index++ )
        {
            mat_temp._matrix[sizet_index] = mat_argl._matrix[sizet_index] + mat_argr._matrix[sizet_index];
        }
    }

    return mat_temp;
}

template <class T> matrix<T> operator+( const matrix <T> &mat_argl, T T_arg )
{

    matrix <T> mat_temp(mat_argl.sizet_rows, mat_argl.sizet_cols);
    size_t sizet_size = mat_argl.sizet_rows * mat_argl.sizet_cols;

    for ( size_t sizet_index = 0; sizet_index < sizet_size; sizet_index++ )
    {
        mat_temp._matrix[sizet_index] = mat_argl._matrix[sizet_index] + T_arg;
    }

    return mat_temp;
}

template <class T> matrix<T> operator+( T T_arg, const matrix <T> &mat_argr )
{

    matrix <T> mat_temp(mat_argr.sizet_rows, mat_argr.sizet_cols);
    size_t sizet_size = mat_argr.sizet_rows * mat_argr.sizet_cols;

    for ( size_t sizet_index = 0; sizet_index < sizet_size; sizet_index++ )
    {
        mat_temp._matrix[sizet_index] = mat_argr._matrix[sizet_index] + T_arg;
    }

    return mat_temp;
}



//  -
template <class T> matrix<T> operator-( const matrix <T> &mat_argl, const matrix <T> &mat_argr)
{

    matrix <T> mat_temp(mat_argl.sizet_rows, mat_argl.sizet_cols);
    size_t sizet_size = mat_argl.sizet_rows * mat_argl.sizet_cols;

    SUSA_ASSERT_MESSAGE((mat_argl.sizet_rows == mat_argr.sizet_rows)
      && (mat_argl.sizet_cols == mat_argr.sizet_cols),
      "the matrices have different sizes.");

    if ((mat_argl.sizet_rows == mat_argr.sizet_rows) && (mat_argl.sizet_cols == mat_argr.sizet_cols))
    {
        for ( size_t sizet_index = 0; sizet_index < sizet_size; sizet_index++ )
        {
            mat_temp._matrix[sizet_index] = mat_argl._matrix[sizet_index] - mat_argr._matrix[sizet_index];
        }
    }

    return mat_temp;
}

template <class T> matrix<T> operator-( const matrix <T> &mat_argl, T T_arg )
{

    matrix <T> mat_temp(mat_argl.sizet_rows, mat_argl.sizet_cols);
    size_t sizet_size = mat_argl.sizet_rows * mat_argl.sizet_cols;

    for ( size_t sizet_index = 0; sizet_index < sizet_size; sizet_index++ )
    {
        mat_temp._matrix[sizet_index] = mat_argl._matrix[sizet_index] - T_arg;
    }

    return mat_temp;
}

template <class T> matrix<T> operator-( T T_arg, const matrix <T> &mat_argr )
{

    matrix <T> mat_temp(mat_argr.sizet_rows, mat_argr.sizet_cols);
    size_t sizet_size = mat_argr.sizet_rows * mat_argr.sizet_cols;

    for ( size_t sizet_index = 0; sizet_index < sizet_size; sizet_index++ )
    {
        mat_temp._matrix[sizet_index] = T_arg - mat_argr._matrix[sizet_index];
    }

    return mat_temp;
}

//  *
template <class T> matrix<T> operator*( const matrix <T> &mat_argl, const matrix <T> &mat_argr)
{

    matrix <T> mat_temp(mat_argl.sizet_rows, mat_argl.sizet_cols);
    size_t sizet_size = mat_argl.sizet_rows * mat_argl.sizet_cols;

    SUSA_ASSERT_MESSAGE((mat_argl.sizet_rows == mat_argr.sizet_rows)
      && (mat_argl.sizet_cols == mat_argr.sizet_cols),
      "the matrices have different sizes.");

    if ((mat_argl.sizet_rows == mat_argr.sizet_rows) && (mat_argl.sizet_cols == mat_argr.sizet_cols))
    {
        for ( size_t sizet_index = 0; sizet_index < sizet_size; sizet_index++ )
        {
            mat_temp._matrix[sizet_index] = mat_argl._matrix[sizet_index] * mat_argr._matrix[sizet_index];
        }
    }

    return mat_temp;
}

template <class T> matrix<T> operator*( const matrix <T> &mat_argl, T T_arg )
{

    matrix <T> mat_temp(mat_argl.sizet_rows, mat_argl.sizet_cols);
    size_t sizet_size = mat_argl.sizet_rows * mat_argl.sizet_cols;

    for ( size_t sizet_index = 0; sizet_index < sizet_size; sizet_index++ )
    {
        mat_temp._matrix[sizet_index] = mat_argl._matrix[sizet_index] * T_arg;
    }

    return mat_temp;
}

template <class T> matrix<T> operator*( T T_arg, const matrix <T> &mat_argr )
{

    matrix <T> mat_temp(mat_argr.sizet_rows, mat_argr.sizet_cols);
    size_t sizet_size = mat_argr.sizet_rows * mat_argr.sizet_cols;

    for ( size_t sizet_index = 0; sizet_index < sizet_size; sizet_index++ )
    {
        mat_temp._matrix[sizet_index] = T_arg * mat_argr._matrix[sizet_index];
    }

    return mat_temp;
}

//  /
template <class T> matrix<T> operator/( const matrix <T> &mat_argl, const matrix <T> &mat_argr)
{

    matrix <T> mat_temp(mat_argl.sizet_rows, mat_argl.sizet_cols);
    size_t sizet_size = mat_argl.sizet_rows * mat_argl.sizet_cols;

    SUSA_ASSERT_MESSAGE((mat_argl.sizet_rows == mat_argr.sizet_rows)
      && (mat_argl.sizet_cols == mat_argr.sizet_cols),
      "the matrices have different sizes.");

    if ((mat_argl.sizet_rows == mat_argr.sizet_rows) && (mat_argl.sizet_cols == mat_argr.sizet_cols))
    {
        for ( size_t sizet_index = 0; sizet_index < sizet_size; sizet_index++ )
        {
            mat_temp._matrix[sizet_index] = mat_argl._matrix[sizet_index] / mat_argr._matrix[sizet_index];
        }
    }

    return mat_temp;
}

template <class T> matrix<T> operator/( const matrix <T> &mat_argl, T T_arg ) {

    matrix <T> mat_temp(mat_argl.sizet_rows, mat_argl.sizet_cols);
    size_t sizet_size = mat_argl.sizet_rows * mat_argl.sizet_cols;

    for ( size_t sizet_index = 0; sizet_index < sizet_size; sizet_index++ ) {
        mat_temp._matrix[sizet_index] = mat_argl._matrix[sizet_index] / T_arg;
    }

    return mat_temp;
}

template <class T> matrix<T> operator/( T T_arg, const matrix <T> &mat_argr ) {

    matrix <T> mat_temp(mat_argr.sizet_row, mat_argr.sizet_col);
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

    SUSA_ASSERT_MESSAGE((mat_arg.sizet_rows == sizet_rows)
      && (mat_arg.sizet_cols == sizet_cols),
      "the matrices have different sizes.");

    if ((mat_arg.sizet_rows == sizet_rows) && (mat_arg.sizet_cols == sizet_cols) && this->_matrix != NULL)
    {
        for (size_t sizet_index = 0; sizet_index < this->sizet_objects; sizet_index++)
        {
            this->_matrix[sizet_index] -= mat_arg._matrix[sizet_index];
        }
    }

    return *this;
}

//  +=
template <class T> matrix<T> matrix <T>::operator+=(const matrix <T> &mat_arg)
{

    SUSA_ASSERT(this->_matrix != NULL);

    SUSA_ASSERT_MESSAGE((mat_arg.sizet_rows == sizet_rows)
      && (mat_arg.sizet_cols == sizet_cols),
      "the matrices have different sizes.");

    if ((mat_arg.sizet_rows == sizet_rows) && (mat_arg.sizet_cols == sizet_cols) && this->_matrix != NULL)
    {
        for (size_t sizet_index = 0; sizet_index < this->sizet_objects; sizet_index++)
        {
            this->_matrix[sizet_index] += mat_arg._matrix[sizet_index];
        }
    }

    return *this;
}

//  =
template <class T> matrix<T>& matrix <T>::operator=( const matrix <T> &mat_arg )
{

    size_t sizet_size = mat_arg.sizet_cols * mat_arg.sizet_rows;

    if (this != &mat_arg)
    {
        if (sizet_size != 0)
        {

            if (this->sizet_objects != sizet_size)
            {
                this->allocate(sizet_size);
            }

            sizet_rows = mat_arg.sizet_rows;
            sizet_cols = mat_arg.sizet_cols;

            for (size_t sizet_index = 0; sizet_index < sizet_size; sizet_index++)
            {
                this->_matrix[sizet_index] = mat_arg._matrix[sizet_index];
            }
        }
        else
        {
            sizet_rows = 0;
            sizet_cols = 0;
            this->deallocate();
        }
    }


    return *this;
}

template <class T> matrix<T> matrix <T>::operator=( const matrix <T> &mat_arg ) const
{

    size_t sizet_size = mat_arg.sizet_cols * mat_arg.sizet_rows;

    if (this != &mat_arg)
    {
        if (sizet_size != 0)
        {

            if ((sizet_rows * sizet_cols) != sizet_size)
            {
                this->allocate(sizet_size);
            }

            sizet_rows = mat_arg.sizet_rows;
            sizet_cols = mat_arg.sizet_cols;

            for (size_t sizet_index = 0; sizet_index < sizet_size; sizet_index++)
            {
                this->_matrix[sizet_index] = mat_arg._matrix[sizet_index];
            }
        }
        else
        {
            sizet_rows = 0;
            sizet_cols = 0;

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
    for (size_t sizet_row = 0; sizet_row < mat_arg.sizet_rows; sizet_row++) {
        for (size_t sizet_col = 0; sizet_col < mat_arg.sizet_cols; sizet_col++) {
            outStream << mat_arg._matrix[sizet_row + sizet_col * mat_arg.sizet_rows];
            if (sizet_col < mat_arg.sizet_cols - 1) outStream << " ";
        }
        if (sizet_row < mat_arg.sizet_rows - 1) outStream << "\n ";
    }
    outStream << "]";
    return outStream;
}


template <class T> void matrix <T>::parser(std::string str_string)
{

    pre_parser(str_string);

    size_t sizet_size = 0;           // linear size of the matrix (#rows * #columns)
    int int_length = str_string.length(); // length of the input string
    size_t sizet_cols_ = 0;          // number of columns
    size_t sizet_rows_ = 0;          // number of rows

    // This loop gets total number of rows and columns
    for (int int_i = 0; int_i < int_length; int_i++)
    {
        switch (str_string[int_i])
        {
        case ';':
            sizet_rows_++;
            // since the number of rows are counted by number of space and the
            // spaces around ';' are cut, we also count a column here !
            sizet_cols_++;
            break;

        case 0xA:
            sizet_rows_++;
            break;

        case 0x20:
            sizet_cols_++;
            break;

        default:
            break;
        };
    }

    sizet_cols_++;
    sizet_rows_++;
    sizet_size = sizet_cols_ * sizet_rows_;

    SUSA_ASSERT_MESSAGE(sizet_cols_ % sizet_rows_ == 0,
      "the number of columns are not equal in each row.");

    if (sizet_cols_ % sizet_rows_ != 0)
    {
        return;
    }


    if (sizet_size != 0)
    {

        if ((this->sizet_objects) != sizet_size)
        {
            this->allocate(sizet_size);
        }

        sizet_rows = sizet_rows_;
        sizet_cols = sizet_cols_ / sizet_rows_;


        std::stringstream ss_all(str_string);
        char* char_buff = (char *)std::malloc(4096);
        SUSA_ASSERT_MESSAGE(char_buff != nullptr, "memory allocation failed.");
        if (char_buff == nullptr) std::exit(EXIT_FAILURE);
        T T_tmp;

        for (size_t sizet_row = 0; sizet_row < sizet_rows; sizet_row++)
        {
            ss_all.getline(char_buff, 4096, ';');
            std::stringstream ssrow(char_buff);
            for (size_t sizet_col = 0; sizet_col < sizet_cols; sizet_col++)
            {
                ssrow.getline(char_buff, 4096, ' ');
                if (!(std::istringstream(char_buff) >> T_tmp)) T_tmp = 0;
                this->_matrix[sizet_row + sizet_col * sizet_rows] = T_tmp;
            }
        }

        std::free(char_buff);

    }
    else
    {
        sizet_rows = 0;
        sizet_cols = 0;

        this->deallocate();
    }
}
}      // NAMESPACE SUSA

#endif // SUSA_MATRIX_H
