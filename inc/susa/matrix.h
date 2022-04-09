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
 * @author Behrooz Kamary
 */

#ifndef SUSA_MATRIX_H
#define SUSA_MATRIX_H


// refer to "Efficient Processing of Two-Dimensional Arrays with C or C++"
// by David I. Donato
// https://pubs.usgs.gov/tm/07/e01/tm7e1.pdf
// for comparison of different memory addressing techniques
// #define get_lindex(R, C)      R + C * sizet_rows

#define MAX_STR_LEN           4096

#include <type_traits>

#include <susa/debug.h>
#include <susa/memory.h>
#include <utility>

namespace susa
{

void pre_parser(std::string& str_string);

template <typename T, typename Allocator> class matrix;

template<class T, template <typename> typename Allocator = std::allocator>
using cmatrix = matrix<std::complex<T>, Allocator<std::complex<T>>>;

template <typename T, typename Allocator>  matrix <T, Allocator> operator-( const matrix <T, Allocator>&, T );
template <typename T, typename Allocator>  matrix <T, Allocator> operator-( T, const matrix <T, Allocator>& );
template <typename T, typename Allocator>  matrix <T, Allocator> operator-( const matrix <T, Allocator>&, const matrix <T, Allocator>& );

template <typename T, typename Allocator>  matrix <T, Allocator> operator+( const matrix <T, Allocator>&, T );
template <typename T, typename Allocator>  matrix <T, Allocator> operator+( T, const matrix <T, Allocator>& );
template <typename T, typename Allocator>  matrix <T, Allocator> operator+( const matrix <T, Allocator>&,const matrix <T, Allocator>& );

template <typename T, typename Allocator>  matrix <T, Allocator> operator*( const matrix <T, Allocator>&, T );
template <typename T, typename Allocator>  matrix <T, Allocator> operator*( T, const matrix <T, Allocator>& );
template <typename T, typename Allocator>  matrix <T, Allocator> operator*( const matrix <T, Allocator>&,const matrix <T, Allocator>& );

template <typename T, typename Allocator>  matrix <T, Allocator> operator/( const matrix <T, Allocator>&, T );
template <typename T, typename Allocator>  matrix <T, Allocator> operator/( T, const matrix <T, Allocator>& );
template <typename T, typename Allocator>  matrix <T, Allocator> operator/( const matrix <T, Allocator>&,const matrix <T, Allocator>& );

template <typename T, typename Allocator> std::ostream &operator<<(std::ostream &, const matrix <T, Allocator> &);
template <typename T, typename Allocator> std::istream &operator>>(std::istream &, matrix <T, Allocator> &);


template <typename T, typename Allocator> matrix <T, Allocator> matmul(const matrix <T, Allocator> &mat_argl,const matrix <T, Allocator> &mat_argr);
template <typename T, typename Allocator> matrix <T, Allocator> transpose(const matrix <T, Allocator> &mat_arg);

template <typename E>
class matrix_expression
{

public:

    double operator()(size_t sz_index) const
    {
      return static_cast<E const&>(*this)(sz_index);
    }

    size_t size() const
    {
       return static_cast<E const&>(*this).size();
    }

    const E& clone(void) const
    {
        return *static_cast<const E*>(this);
    }
};

template <typename OPR, typename LHS, typename RHS>
struct expr_binary : public matrix_expression<expr_binary<OPR, LHS, RHS>>
{
  const LHS&  lhs;
  const RHS&  rhs;
  OPR         opr;


  expr_binary(OPR&& functor, const LHS& lhs, const RHS& rhs)
  : lhs(lhs)
  , rhs(rhs)
  , opr(std::forward<OPR>(functor))
  {}

  auto operator()() const
  {
    return opr(lhs, rhs);
  }

  using result_t = decltype(std::declval<OPR>()(std::declval<const LHS&>(), std::declval<const RHS&>()));
  operator result_t() const
  {
    return this->operator()();
  }
};

struct expr_add
{
  template <class T, class U>
  auto operator()(const T& left, const U& right) const
  {
    return left + right;
  }
};

struct expr_sub
{
  template <class T, class U>
  auto operator()(const T& left, const U& right) const
  {
    return left + right;
  }
};

/**
 * @class matrix
 * @brief <i>matrix</i> type.
 * A matrix is a two dimensional array aimed to
 * be used with algebraic operators and functions.
 * A single row or a single column matrix is a <i>vector</i>.
 * The methods, operators and functions that take vectors as their
 * input parameters perform column-wise computation when they are passed
 * matrices (and not vectors based on the above definition). This is to follow
 * the algebraic convention where a matrix column is a vector in a space.
 *
 * @ingroup TYPES
 *
 */
template <typename T, typename Allocator = std::allocator <T>>
class matrix
: public matrix_expression<matrix<T, Allocator>>
{
  private:
    size_t          sizet_rows;
    size_t          sizet_cols;
    size_t          sizet_objects;
    T*              T_default;
    Allocator       alloc;
    T*              _matrix;

    // UTILITY METHODS

    // This parses the strings and initialize the matrix elements.
    // The string initialization can be done using a constructor method
    // and an overloaded assignment operator (=) method
    bool    parser(std::string str_string);

    bool    allocate(size_t sz_objects);

  public:
    //! Constructor
    matrix();

    /**
     * Constructor
     *
     * @param sizet_rows Number of rows
     * @param sizet_cols Number of columns
     * @param Tinitial Initial value of all elements.
     */
    matrix(size_t sizet_rows, size_t sizet_cols, T Tinitial);

    /**
     * Constructor
     * the elements are not initialized to speed up the instantiation.
     *
     * @param sizet_rows Number of rows
     * @param sizet_cols Number of columns.
     */
    matrix(size_t sizet_rows, size_t sizet_cols);

    /**
     * Constructor
     *
     * @param tshape the shape of the matrix as std::tuple
     */
    matrix(const std::tuple <size_t,size_t>& tshape);

    /**
     * Constructor
     *
     * @param tshape the shape of the matrix as std::tuple
     * @param Tinitial Initial value of all elements.
     */
    matrix(const std::tuple <size_t,size_t>& tshape, T Tinitial);

    //! Copy constructor
    matrix(const matrix <T, Allocator> &mat_arg);

    //! Copy Constructor for rvalues
    matrix(matrix <T, Allocator> &&mat_arg) noexcept;

    //! Destructor
    ~matrix() noexcept;

    /**
     * @brief Constructor
     *
     * @param str_string string representation of the matrix
     */
    matrix(std::string str_string);

    /**
     * @brief Constructor
     *
     * @param c_string literal string representation of the matrix
     */
    matrix(const char char_string[]);

    //! Returns the value of a specific (row, column)
    T get(size_t sizet_row, size_t sizet_col) const;

    //! Returns the value of a specific (elem)
    T get(size_t sizet_elem) const;

    //! Returns the linear index of a specific (row, column) pair
    constexpr size_t get_lindex(size_t sizet_row, size_t sizet_col) const
    {
        return (sizet_row + sizet_col * sizet_rows);
    }

    //! Returns the number of columns
    size_t no_cols() const;

    //! Returns the number of rows
    size_t no_rows() const;

    //! Returns the shape (number of rows and columns) as a tuple
    std::tuple <size_t,size_t> shape() const;

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

    /**
     * @brief set a row
     *
     * copy the input matrix elemnts into a row
     *
     * @param row the row index
     * @param mat_arg the input matrix
     */
    bool set_row(size_t row, const matrix <T, Allocator>& mat_arg);

    /**
     * @brief set a column
     * copy the input matrix elemnts into a column
     *
     * @param col the column index
     * @param mat_arg the input matrix
     */
    bool set_col(size_t col, const matrix <T, Allocator>& mat_arg);

    /**
     * @brief get a column as std::vector
     * copy the specified matrix column into a std::vector
     *
     * @param col the column index
     */
    std::vector <T> get_col_as_vec(size_t col);

    /**
     * @brief get a row as std::vector
     * copy the specified matrix row into a std::vector
     *
     * @param row the row index
     */
    std::vector <T> get_row_as_vec(size_t row);

    /**
     * @brief swap two columns of the matrix
     */
    void swap_cols(size_t col_a, size_t col_b);

    /**
     * @brief swap two rows of the matrix
     */
    void swap_rows(size_t row_a, size_t row_b);

    //! returns the size of matrix
    void set_all(T T_arg);

    //! Returns the indicated row
    matrix <T, Allocator> row(size_t sizet_row) const;

    //! Returns the indicated column
    matrix <T, Allocator> col(size_t sizet_col) const;

    //! Returns the number of allocated objects.
    size_t size() const;

    //! Returns a matrix after eliminating the specified column and the row
    matrix <T, Allocator> shrink(size_t sizet_row, size_t sizet_col) const;

    //! Resize the matrix to the newly specified dimensions
    bool resize(size_t sizet_row, size_t sizet_col);

    //! Reshape the matrix to the newly specified dimensions
    matrix <T, Allocator> reshape(size_t sizet_row, size_t sizet_col);

    //! Considers the matrix object as a vector and return left side of that vector
    matrix <T, Allocator> left(size_t sizet_left) const;

    //! Considers the matrix object as a vector and return right side of that vector
    matrix <T, Allocator> right(size_t sizet_right) const;

    //! Considers the matrix object as a vector and return mid part of that vector
    matrix <T, Allocator> mid(size_t sizet_begin, size_t sizet_end) const;

    //! Element wise Assignment by Addition operator
    matrix <T, Allocator> operator+=( const matrix <T, Allocator> &mat_arg );

    //! Element wise Assignment by Subtraction
    matrix <T, Allocator> operator-=( const matrix <T, Allocator> &mat_arg );

    //! Element wise Assignment operator
    matrix <T, Allocator>& operator=( const matrix <T, Allocator> &mat_arg );

    //! Element wise Assignment operator
    matrix <T, Allocator>& operator=( const std::string& str_string );

    //! Element wise Assignment operator
    matrix <T, Allocator>& operator=( const char char_string[] );

    //! Element wise Subtraction operator
    friend matrix <T, Allocator> operator-<>( const matrix <T, Allocator> &mat_argl, T T_arg);

    //! Element wise Subtraction operator
    friend matrix <T, Allocator> operator-<>( T T_arg, const matrix <T, Allocator> &mat_argr );

    //! Element wise Subtraction operator
    friend matrix <T, Allocator> operator-<>( const matrix <T, Allocator> &mat_argl, const matrix <T, Allocator> &mat_argr);

    //! Element wise Addition operator
    friend matrix <T, Allocator> operator+<>( const matrix <T, Allocator> &mat_argl, T T_arg );

    //! Element wise Addition operator
    friend matrix <T, Allocator> operator+<>( T T_arg, const matrix <T, Allocator> &mat_argr );

    //! Element wise Addition operator
    friend matrix <T, Allocator> operator+<>( const matrix <T, Allocator> &mat_argl, const matrix <T, Allocator> &mat_argr);

    //! Element wise Multiplication  operator
    friend matrix <T, Allocator> operator*<>( const matrix <T, Allocator> &mat_argl, T T_arg);

    //! Element wise Multiplication  operator
    friend matrix <T, Allocator> operator*<>( T T_arg, const matrix <T, Allocator> &mat_argr );

    //! Element wise Multiplication  operator
    friend matrix <T, Allocator> operator*<>( const matrix <T, Allocator> &mat_argl, const matrix <T, Allocator> &mat_argr );

    //! Element wise Division operator
    friend matrix <T, Allocator> operator/<>( const matrix <T, Allocator> &mat_argl, T T_arg);

    //! Element wise Division operator
    friend matrix <T, Allocator> operator/<>( const matrix <T, Allocator> &mat_argl, const matrix <T, Allocator> &mat_argr);

    //! Output stream
    friend std::ostream &operator<< <>(std::ostream &outStream, const matrix <T, Allocator> &mat_arg);

    //! Input stream
    friend std::istream &operator>> <>(std::istream &inStream, matrix <T, Allocator> &mat_arg);

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

    operator matrix <int8_t> ();

    operator matrix <uint8_t> ();

    operator matrix <std::complex <double> > ();

    operator matrix <std::complex <float> > ();

    operator matrix <std::complex <int> > ();

    operator matrix <std::complex <int8_t> > ();

    operator matrix <std::complex <uint8_t> > ();


    // Friend methods are used to speed up matrix operations.
    // Because they can access '_matrix' directly.
    // They are defined in 'linalg.h'
    friend matrix <T, Allocator> matmul <> (const matrix <T, Allocator> &mat_argl, const matrix <T, Allocator> &mat_argr);
    friend matrix <T, Allocator> transpose <> (const matrix <T, Allocator> &mat_arg);

    friend bool operator!=( const susa::matrix <T, Allocator> &mat_argl, const susa::matrix <T, Allocator> &mat_argr)
    {
        if (mat_argl.shape() != mat_argr.shape()) return true;

        for (size_t sizet_indx = 0; sizet_indx < mat_argl.sizet_objects; sizet_indx++)
        {
            if (mat_argl._matrix[sizet_indx] != mat_argr._matrix[sizet_indx]) return true;
        }

        return false;
    }

    friend bool operator==( const susa::matrix <T, Allocator> &mat_argl, const susa::matrix <T, Allocator> &mat_argr)
    {
        if (mat_argl.shape() != mat_argr.shape()) return false;

        for (size_t sizet_indx = 0; sizet_indx < mat_argl.sizet_objects; sizet_indx++)
        {
            if (mat_argl._matrix[sizet_indx] != mat_argr._matrix[sizet_indx]) return false;
        }

        return true;
    }
};

// Constructors and Destructor

template <typename T, typename Allocator> matrix <T, Allocator>::matrix()
: sizet_rows(0)
, sizet_cols(0)
, sizet_objects(0)
, T_default(new T)
, alloc()
, _matrix(nullptr)
{

}

template <typename T, typename Allocator>
matrix <T, Allocator>::matrix(size_t sizet_rows, size_t sizet_cols, T Tinitial)
: T_default(new T)
, alloc()
, _matrix(nullptr)
{

  SUSA_ASSERT(sizet_cols > 0 && sizet_rows > 0);

  this->sizet_rows = sizet_rows < 2 ? 1 : sizet_rows;
  this->sizet_cols = sizet_cols < 2 ? 1 : sizet_cols;
  sizet_objects = this->sizet_rows * this->sizet_cols;

  if (!allocate(sizet_objects)) return;

  for (size_t sizet_index = 0; sizet_index < sizet_objects; sizet_index++)
  {
    _matrix[sizet_index] = Tinitial;
  }
}

template <typename T, typename Allocator>
matrix <T, Allocator>::matrix(size_t sizet_rows, size_t sizet_cols)
: T_default(new T)
, alloc()
, _matrix(nullptr)
{
  SUSA_ASSERT(sizet_cols > 0 && sizet_rows > 0);
  this->sizet_rows = sizet_rows < 2 ? 1 : sizet_rows;
  this->sizet_cols = sizet_cols < 2 ? 1 : sizet_cols;
  sizet_objects = this->sizet_rows * this->sizet_cols;
  allocate(sizet_objects);
}

template <typename T, typename Allocator>
matrix <T, Allocator>::matrix(const matrix <T, Allocator> &mat_arg)
: T_default(new T)
, alloc(mat_arg.alloc)
, _matrix(nullptr)
{
  this->sizet_rows      = mat_arg.sizet_rows;
  this->sizet_cols      = mat_arg.sizet_cols;
  this->sizet_objects   = mat_arg.sizet_objects;

  if (allocate(sizet_objects))
  {
    std::memcpy(_matrix, mat_arg._matrix, sizet_objects * sizeof(T));
  }
}

template <typename T, typename Allocator>
matrix <T, Allocator>::matrix(const std::tuple<size_t, size_t>& tshape)
: T_default(new T)
, _matrix(nullptr)
{
  size_t sizet_rows;
  size_t sizet_cols;
  std::tie(sizet_rows, sizet_cols) = tshape;

  SUSA_ASSERT(sizet_cols > 0 && sizet_rows > 0);
  this->sizet_rows = sizet_rows < 2 ? 1 : sizet_rows;
  this->sizet_cols = sizet_cols < 2 ? 1 : sizet_cols;
  sizet_objects = this->sizet_rows * this->sizet_cols;
  allocate(sizet_objects);
}

template <typename T, typename Allocator>
matrix <T, Allocator>::matrix(const std::tuple<size_t, size_t>& tshape, T Tinitial)
: T_default(new T)
, _matrix(nullptr)
{
  size_t sizet_rows;
  size_t sizet_cols;
  std::tie(sizet_rows, sizet_cols) = tshape;

  SUSA_ASSERT(sizet_cols > 0 && sizet_rows > 0);
  this->sizet_rows    = sizet_rows < 2 ? 1 : sizet_rows;
  this->sizet_cols    = sizet_cols < 2 ? 1 : sizet_cols;
  sizet_objects       = this->sizet_rows * this->sizet_cols;

  if(!allocate(sizet_objects)) return;

  for (size_t sizet_index = 0; sizet_index < sizet_objects; sizet_index++)
  {
    _matrix[sizet_index] = Tinitial;
  }
}

template <typename T, typename Allocator>
matrix <T, Allocator>::matrix(matrix<T, Allocator>&& mat_arg) noexcept
{
  sizet_rows          = mat_arg.sizet_rows;
  sizet_cols          = mat_arg.sizet_cols;
  sizet_objects       = mat_arg.sizet_objects;
  T_default           = mat_arg.T_default;
  T_default           = mat_arg.T_default;
  _matrix             = mat_arg._matrix;
  mat_arg.T_default   = nullptr;
  mat_arg._matrix     = nullptr;
}


template <typename T, typename Allocator>
matrix <T, Allocator>::matrix(std::string str_string)
: sizet_rows(0)
, sizet_cols(0)
, sizet_objects(0)
, T_default(new T)
, alloc()
, _matrix(nullptr)
{
  parser(str_string);
}

template <typename T, typename Allocator>
matrix <T, Allocator>::matrix(const char char_string[])
: sizet_rows(0)
, sizet_cols(0)
, sizet_objects(0)
, T_default(new T)
, alloc()
, _matrix(nullptr)
{
  parser(std::string(char_string));
}

template <typename T, typename Allocator> matrix <T, Allocator>::~matrix() noexcept
{
    if (T_default != nullptr) delete T_default;
    if (_matrix != nullptr) alloc.deallocate(_matrix, sizet_objects);
}

// Public methods

template <typename T, typename Allocator>
T matrix <T, Allocator>::get(size_t sizet_row, size_t sizet_col) const
{
  SUSA_ASSERT(_matrix != nullptr);

  SUSA_ASSERT_MESSAGE(sizet_row < sizet_rows && sizet_col < sizet_cols, "one or more indices is/are out of range.");

  return _matrix[get_lindex(sizet_row, sizet_col)];
}

template <typename T, typename Allocator>
T matrix <T, Allocator>::get(size_t sizet_elem) const
{
  SUSA_ASSERT(_matrix != nullptr);

  SUSA_ASSERT_MESSAGE(sizet_elem < this->sizet_objects, "the element index is out of range.");

  return _matrix[sizet_elem];
}

template <typename T, typename Allocator> size_t inline matrix <T, Allocator>::no_cols() const
{
  return sizet_cols;
}

template <typename T, typename Allocator> size_t inline matrix <T, Allocator>::no_rows() const
{
  return sizet_rows;
}

template <typename T, typename Allocator>  std::tuple<size_t, size_t> matrix<T, Allocator>::shape() const
{
  return std::make_tuple(sizet_rows, sizet_cols);
}

template <typename T, typename Allocator> bool  matrix <T, Allocator>::is_square() const
{
  return (sizet_rows == sizet_cols);
}

template <typename T, typename Allocator> bool  matrix <T, Allocator>::is_vector() const
{
  return ((sizet_rows == 1 && sizet_cols > 1 ) || ( sizet_rows > 1 && sizet_cols == 1));
}

template <typename T, typename Allocator> bool  matrix <T, Allocator>::is_scalar() const
{
  return (sizet_rows == 1 && sizet_cols == 1);
}

template <typename T, typename Allocator> void  matrix <T, Allocator>::set_all(T T_arg)
{
  for (size_t sizet_counter = 0; sizet_counter < this->sizet_objects; sizet_counter++)
  {
    _matrix[sizet_counter] = T_arg;
  }
}

template <typename T, typename Allocator> matrix <T, Allocator> matrix <T, Allocator>::row(size_t sizet_row) const
{
    matrix <T, Allocator> mat_ret(1, sizet_cols);

    if (sizet_row < sizet_rows)
    {
        for (size_t sizet_i = 0; sizet_i < sizet_cols; sizet_i++)
        {
          mat_ret(sizet_i) = _matrix[get_lindex(sizet_row, sizet_i)];
        }
    }

    return mat_ret;
}

template <typename T, typename Allocator> matrix <T, Allocator> matrix <T, Allocator>::col(size_t sizet_col) const
{
    matrix <T, Allocator> mat_ret(sizet_rows, 1);

    if (sizet_col < sizet_cols)
    {
        for (size_t sizet_i = 0; sizet_i < sizet_rows; sizet_i++)
        {
          mat_ret(sizet_i) = _matrix[get_lindex(sizet_i, sizet_col)];
        }
    }

    return mat_ret;
}

template <typename T, typename Allocator>
size_t matrix <T, Allocator>::size() const
{
    return sizet_objects;
}

template <typename T, typename Allocator>
bool matrix <T, Allocator>::set_row(size_t row, const matrix <T, Allocator>& mat_arg)
{
  SUSA_ASSERT_MESSAGE((row < sizet_rows || mat_arg.sizet_objects >= sizet_cols), "dimention mismatch");
  if (row > sizet_rows || mat_arg.sizet_objects < sizet_cols) return false;

  for (size_t indx = 0; indx < sizet_cols; indx++)
  {
    _matrix[get_lindex(row,indx)] = mat_arg._matrix[indx];
  }

  return true;
}

template <typename T, typename Allocator>
bool matrix <T, Allocator>::set_col(size_t col, const matrix <T, Allocator>& mat_arg)
{
  SUSA_ASSERT_MESSAGE((col < sizet_cols || mat_arg.sizet_objects >= sizet_rows), "dimension mismatch");
  if (col > sizet_cols || mat_arg.sizet_objects < sizet_rows) return false;

  for (size_t indx = 0; indx < sizet_rows; indx++)
  {
    _matrix[get_lindex(indx,col)] = mat_arg._matrix[indx];
  }

  return true;
}

template <typename T, typename Allocator> std::vector <T> matrix <T, Allocator>::get_col_as_vec(size_t col)
{
  std::vector<T> ret(sizet_rows);
  for (size_t indx = 0; indx < sizet_rows; indx++)
  {
    ret[indx] = get(indx, col);
  }

  return ret;
}

template <typename T, typename Allocator> std::vector <T> matrix <T, Allocator>::get_row_as_vec(size_t row)
{
    std::vector<T> ret(sizet_cols);
    for (size_t indx = 0; indx < sizet_cols; indx++)
    {
        ret[indx] = get(row, indx);
    }

    return ret;
}

template <typename T, typename Allocator>
void matrix <T, Allocator>::swap_cols(size_t col_a, size_t col_b)
{
  T T_tmp;

  for (size_t row = 0; row < sizet_rows; row++)
  {
      T_tmp = _matrix[get_lindex(row,col_a)];
      _matrix[get_lindex(row,col_a)] = _matrix[get_lindex(row,col_b)];
      _matrix[get_lindex(row,col_b)] = T_tmp;
  }
}

template <typename T, typename Allocator> void matrix <T, Allocator>::swap_rows(size_t row_a, size_t row_b)
{
  T T_tmp;

  for (size_t col = 0; col < sizet_cols; col++)
  {
      T_tmp = _matrix[get_lindex(row_a,col)];
      _matrix[get_lindex(row_a,col)] = _matrix[get_lindex(row_b,col)];
      _matrix[get_lindex(row_b,col)] = T_tmp;
  }
}

template <typename T, typename Allocator>
matrix <T, Allocator> matrix <T, Allocator>::shrink(size_t sizet_elim_row, size_t sizet_elim_col) const
{

  matrix <T, Allocator> mat_ret;

  size_t sizet_new_col;
  size_t sizet_new_row;

  SUSA_ASSERT_MESSAGE(sizet_cols > 1 && sizet_rows > 1, "the input arguments error.");

  if (sizet_cols > 1 && sizet_rows > 1)
  {
    mat_ret = matrix <T, Allocator> (sizet_rows - 1,sizet_cols - 1);
  }
  else
  {
    return *this;
  }


  SUSA_ASSERT_MESSAGE(sizet_elim_col < sizet_cols && sizet_elim_row < sizet_rows, "the input arguments exceed matrix size.");

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
            _matrix[get_lindex(sizet_row,sizet_col)];
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

template <typename T, typename Allocator>
matrix <T, Allocator> matrix <T, Allocator>::reshape(size_t sizet_rows, size_t sizet_cols)
{
    SUSA_ASSERT((sizet_rows * sizet_cols) == sizet_objects);
    matrix <T, Allocator> mat_ret;
    if ((sizet_rows * sizet_cols) != sizet_objects) return mat_ret;
    mat_ret = *this;
    mat_ret.resize(sizet_rows, sizet_cols);
    return mat_ret;
}

template <typename T, typename Allocator>
bool matrix<T, Allocator>::resize(size_t sizet_rows, size_t sizet_cols)
{
  SUSA_ASSERT(sizet_cols > 0 && sizet_rows > 0);

  this->sizet_rows        = sizet_rows < 2 ? 1 : sizet_rows;
  this->sizet_cols        = sizet_cols < 2 ? 1 : sizet_cols;
  size_t sizet_objects    = this->sizet_rows * this->sizet_cols;

  if (sizet_objects != this->sizet_objects)
  {
    if (_matrix != nullptr)
    {
      alloc.deallocate(_matrix, this->sizet_objects);
      _matrix = nullptr;
    }
    if (allocate(sizet_objects))
    {
      this->sizet_objects = sizet_objects;
    }
  }

  return (_matrix != nullptr);
}

template <typename T, typename Allocator>
matrix <T, Allocator> matrix <T, Allocator>::left(size_t sizet_left) const
{
  matrix <T, Allocator> mat_ret;

  SUSA_ASSERT(_matrix != nullptr);
  if (_matrix == nullptr) return mat_ret;

  if (is_vector())
  {
    SUSA_ASSERT(sizet_objects >= sizet_left);

    if (sizet_rows == 1) mat_ret = matrix <T, Allocator> (1,sizet_left);
    else if (sizet_cols == 1) mat_ret = matrix <T, Allocator> (sizet_left,1);

    for (size_t sizet_i = 0; sizet_i < sizet_left; sizet_i++)
    {
      mat_ret(sizet_i) = _matrix[sizet_i];
    }
  }
  else
  {
    mat_ret = matrix <T, Allocator> (sizet_rows,sizet_left);
    for (size_t sizet_row = 0; sizet_row < sizet_rows; sizet_row++)
      for (size_t sizet_col = 0; sizet_col < sizet_cols && sizet_col < sizet_left; sizet_col++)
          mat_ret(sizet_row, sizet_col) = get(sizet_row,sizet_col);
  }

  return mat_ret;
}


template <typename T, typename Allocator>
matrix <T, Allocator> matrix <T, Allocator>::right(size_t sizet_right) const
{
  matrix <T, Allocator> mat_ret;

  SUSA_ASSERT(_matrix != nullptr);
  if (_matrix == nullptr) return mat_ret;

  if (is_vector())
  {
    SUSA_ASSERT(sizet_objects >= sizet_right);

    if (sizet_rows == 1) mat_ret = matrix <T, Allocator> (1,sizet_right);
    else if (sizet_cols == 1) mat_ret = matrix <T, Allocator> (sizet_right,1);

    size_t sz_start  = sizet_objects - sizet_right - 1;
    T*     _r_matrix = &_matrix[sz_start];
    for (size_t sizet_i = 0; sizet_i < sizet_right; sizet_i++)
    {
      mat_ret(sizet_i) = _r_matrix[sizet_i];
    }
  }
  else
  {
    size_t sz_start = sizet_cols - sizet_right;
    mat_ret = matrix <T, Allocator> (sizet_rows, sizet_right);
    for (size_t sizet_row = 0; sizet_row < sizet_rows; sizet_row++)
      for (size_t sizet_col = 0; sizet_col < sizet_cols && sizet_col < sizet_right; sizet_col++)
          mat_ret(sizet_row, sizet_col) = get(sizet_row, sz_start + sizet_col);
  }

  return mat_ret;
}


template <typename T, typename Allocator>
matrix <T, Allocator> matrix <T, Allocator>::mid(size_t sizet_begin, size_t sizet_end) const
{
  matrix <T, Allocator> mat_ret;

  SUSA_ASSERT(_matrix != nullptr);
  if (_matrix == nullptr) return mat_ret;
  SUSA_ASSERT(sizet_end > sizet_begin);
  if (sizet_end < sizet_begin) return mat_ret;

  size_t sz_width = sizet_end - sizet_begin + 1;

  if (is_vector())
  {
    SUSA_ASSERT(sizet_objects > sizet_begin && sizet_objects > sizet_end);
    if (sizet_begin > sizet_objects || sizet_end > sizet_objects) return mat_ret;

    if (sizet_rows == 1) mat_ret = matrix <T, Allocator> (1,sz_width);
    else if (sizet_cols == 1) mat_ret = matrix <T, Allocator> (sz_width,1);

    T*     _r_matrix = &_matrix[sizet_begin];
    for (size_t sizet_i = 0; sizet_i < sz_width; sizet_i++)
    {
      mat_ret(sizet_i) = _r_matrix[sizet_i];
    }
  }
  else
  {
    SUSA_ASSERT(sizet_cols > sizet_begin && sizet_cols > sizet_end);
    if (sizet_begin > sizet_cols || sizet_end > sizet_cols) return mat_ret;

    mat_ret = matrix <T, Allocator> (sizet_rows, sz_width);
    for (size_t sizet_row = 0; sizet_row < sizet_rows; sizet_row++)
      for (size_t sizet_col = 0; sizet_col < sizet_cols && sizet_col < sz_width; sizet_col++)
          mat_ret(sizet_row, sizet_col) = get(sizet_row, sizet_begin + sizet_col);
  }

  return mat_ret;
}

// Operators

//  ()
template <typename T, typename Allocator>
T& matrix<T, Allocator>::operator ()(size_t sizet_row, size_t sizet_col)
{
  SUSA_ASSERT(_matrix != nullptr);
  SUSA_ASSERT_MESSAGE(sizet_row < sizet_rows && sizet_col < sizet_cols, "one or more indices is/are out of range.");

  if (sizet_row < sizet_rows && sizet_col < sizet_cols && _matrix != nullptr)
  {
    return _matrix[get_lindex(sizet_row,sizet_col)];
  }

  return *T_default;
}

template <typename T, typename Allocator>
T matrix<T, Allocator>::operator ()(size_t sizet_row, size_t sizet_col) const
{
  return get(sizet_row, sizet_col);
}

template <typename T, typename Allocator>
T& matrix<T, Allocator>::operator ()(size_t sizet_elem)
{
  SUSA_ASSERT(_matrix != nullptr);
  SUSA_ASSERT_MESSAGE(sizet_elem < sizet_objects, "the index is out of range.");

  if (sizet_elem < sizet_objects && _matrix != nullptr)
  {
    return _matrix[sizet_elem];
  }

  return *T_default;
}

template <typename T, typename Allocator>
T matrix<T, Allocator>::operator ()(size_t sizet_elem) const
{
  return get(sizet_elem);
}

// Typecasting
template <typename T, typename Allocator>
matrix<T, Allocator>::operator matrix <double> ()
{
  matrix <double> mat_ret(sizet_rows, sizet_cols);

  for (size_t sizet_elem = 0; sizet_elem < sizet_objects; sizet_elem++)
  {
    mat_ret(sizet_elem) = (double)_matrix[sizet_elem];
  }

  return mat_ret;
}

template <typename T, typename Allocator>
matrix<T, Allocator>::operator matrix <float> ()
{
  matrix <float> mat_ret(sizet_rows, sizet_cols);

  for (size_t sizet_elem = 0; sizet_elem < sizet_objects; sizet_elem++)
  {
    mat_ret(sizet_elem) = (float)_matrix[sizet_elem];
  }
  return mat_ret;
}

template <typename T, typename Allocator>
matrix<T, Allocator>::operator matrix <int> ()
{
  matrix <int> mat_ret(sizet_rows, sizet_cols);

  for (size_t sizet_elem = 0; sizet_elem < sizet_objects; sizet_elem++)
  {
    mat_ret(sizet_elem) = (int)_matrix[sizet_elem];
  }
  return mat_ret;
}


template <typename T, typename Allocator>
matrix<T, Allocator>::operator matrix <int8_t> ()
{
  matrix <int8_t> mat_ret(sizet_rows, sizet_cols);

  for (size_t sizet_elem = 0; sizet_elem < sizet_objects; sizet_elem++)
  {
    mat_ret(sizet_elem) = (int8_t)_matrix[sizet_elem];
  }
  return mat_ret;
}

template <typename T, typename Allocator>
matrix<T, Allocator>::operator matrix <uint8_t> ()
{
  matrix <uint8_t> mat_ret(sizet_rows, sizet_cols);

  for (size_t sizet_elem = 0; sizet_elem < sizet_objects; sizet_elem++)
  {
    mat_ret(sizet_elem) = (uint8_t)_matrix[sizet_elem];
  }
  return mat_ret;
}

template <typename T, typename Allocator>
matrix<T, Allocator>::operator matrix <std::complex <double> > ()
{
  matrix <std::complex <double> > mat_ret(sizet_rows, sizet_cols);

  for (size_t sizet_elem = 0; sizet_elem < sizet_objects; sizet_elem++)
  {
    mat_ret(sizet_elem) = std::complex <double> ((double)_matrix[sizet_elem].real(), (double)_matrix[sizet_elem].imag());
  }
  return mat_ret;
}

template <typename T, typename Allocator>
matrix<T, Allocator>::operator matrix <std::complex <float> > ()
{
  matrix <std::complex <float> > mat_ret(sizet_rows, sizet_cols);

  for (size_t sizet_elem = 0; sizet_elem < sizet_objects; sizet_elem++)
  {
    mat_ret(sizet_elem) = std::complex <float> ((float)_matrix[sizet_elem].real(), (float)_matrix[sizet_elem].imag());
  }
  return mat_ret;
}

template <typename T, typename Allocator>
matrix<T, Allocator>::operator matrix <std::complex <int> > ()
{
  matrix <std::complex <int> > mat_ret(sizet_rows, sizet_cols);

  for (size_t sizet_elem = 0; sizet_elem < this->sizet_objects; sizet_elem++)
  {
    mat_ret(sizet_elem) = std::complex <int> ((int)_matrix[sizet_elem].real(), (int)_matrix[sizet_elem].imag());
  }
  return mat_ret;
}

template <typename T, typename Allocator>
matrix<T, Allocator>::operator matrix <std::complex <int8_t> > ()
{
  matrix <std::complex <int8_t> > mat_ret(sizet_rows, sizet_cols);

  for (size_t sizet_elem = 0; sizet_elem < this->sizet_objects; sizet_elem++)
  {
    mat_ret(sizet_elem) = std::complex <int8_t> ((int8_t)_matrix[sizet_elem].real(), (int8_t)_matrix[sizet_elem].imag());
  }
  return mat_ret;
}

template <typename T, typename Allocator>
matrix<T, Allocator>::operator matrix <std::complex <uint8_t> > ()
{
  matrix <std::complex <uint8_t> > mat_ret(sizet_rows, sizet_cols);

  for (size_t sizet_elem = 0; sizet_elem < this->sizet_objects; sizet_elem++)
  {
    mat_ret(sizet_elem) = std::complex <uint8_t> ((uint8_t)_matrix[sizet_elem].real(), (uint8_t)_matrix[sizet_elem].imag());
  }
  return mat_ret;
}

//  +
template <typename T, typename Allocator>
matrix<T, Allocator> operator+(const matrix <T, Allocator> &mat_argl, const matrix <T, Allocator> &mat_argr)
{
  return expr_binary(expr_add{}, mat_argl, mat_argr);
}
/*
{
  matrix <T, Allocator> mat_ret(mat_argl.sizet_rows, mat_argl.sizet_cols);
  size_t sizet_size = mat_argl.sizet_rows * mat_argl.sizet_cols;

  SUSA_ASSERT_MESSAGE((mat_argl.sizet_rows == mat_argr.sizet_rows)
    && (mat_argl.sizet_cols == mat_argr.sizet_cols),
    "the matrices have different sizes.");

  if ((mat_argl.sizet_rows == mat_argr.sizet_rows) && (mat_argl.sizet_cols == mat_argr.sizet_cols))
  {
      for ( size_t sizet_index = 0; sizet_index < sizet_size; sizet_index++ )
      {
        mat_ret._matrix[sizet_index] = mat_argl._matrix[sizet_index] + mat_argr._matrix[sizet_index];
      }
  }

  return mat_ret;
}
*/

template <typename T, typename Allocator>
matrix<T, Allocator> operator+(const matrix <T, Allocator> &mat_argl, T T_arg)
{
  matrix <T, Allocator> mat_ret(mat_argl.shape());
  size_t sizet_size = mat_argl.sizet_objects;

  for ( size_t sizet_index = 0; sizet_index < sizet_size; sizet_index++ )
  {
    mat_ret._matrix[sizet_index] = mat_argl._matrix[sizet_index] + T_arg;
  }

  return mat_ret;
}

template <typename T, typename Allocator>
matrix<T, Allocator> operator+(T T_arg, const matrix <T, Allocator> &mat_argr)
{
  matrix <T, Allocator> mat_ret(mat_argr.shape());
  size_t sizet_size = mat_argr.sizet_objects;

  for (size_t sizet_index = 0; sizet_index < sizet_size; sizet_index++)
  {
    mat_ret._matrix[sizet_index] = mat_argr._matrix[sizet_index] + T_arg;
  }

  return mat_ret;
}



//  -
template <typename T, typename Allocator>
matrix<T, Allocator> operator-( const matrix <T, Allocator>& mat_argl, const matrix <T, Allocator>& mat_argr)
{
  matrix <T, Allocator> mat_ret(mat_argl.shape());
  size_t sizet_size = mat_argl.sizet_objects;

  SUSA_ASSERT_MESSAGE((mat_argl.sizet_rows == mat_argr.sizet_rows)
    && (mat_argl.sizet_cols == mat_argr.sizet_cols),
    "the matrices have different sizes.");

  if ((mat_argl.sizet_rows == mat_argr.sizet_rows) && (mat_argl.sizet_cols == mat_argr.sizet_cols))
  {
    for ( size_t sizet_index = 0; sizet_index < sizet_size; sizet_index++ )
    {
      mat_ret._matrix[sizet_index] = mat_argl._matrix[sizet_index] - mat_argr._matrix[sizet_index];
    }
  }

  return mat_ret;
}

template <typename T, typename Allocator>
matrix<T, Allocator> operator-( const matrix <T, Allocator> &mat_argl, T T_arg )
{
  matrix <T, Allocator> mat_ret(mat_argl.shape());
  size_t sizet_size = mat_argl.sizet_objects;

  for (size_t sizet_index = 0; sizet_index < sizet_size; sizet_index++)
  {
    mat_ret._matrix[sizet_index] = mat_argl._matrix[sizet_index] - T_arg;
  }

  return mat_ret;
}

template <typename T, typename Allocator> matrix<T, Allocator> operator-( T T_arg, const matrix <T, Allocator> &mat_argr )
{
  matrix <T, Allocator> mat_ret(mat_argr.shape());
  size_t sizet_size = mat_argr.sizet_objects;

  for (size_t sizet_index = 0; sizet_index < sizet_size; sizet_index++)
  {
    mat_ret._matrix[sizet_index] = T_arg - mat_argr._matrix[sizet_index];
  }

  return mat_ret;
}

//  *
template <typename T, typename Allocator>
matrix<T, Allocator> operator*( const matrix <T, Allocator> &mat_argl, const matrix <T, Allocator> &mat_argr)
{
  matrix <T, Allocator> mat_ret(mat_argl.shape());
  size_t sizet_size = mat_argl.sizet_objects;

  SUSA_ASSERT_MESSAGE((mat_argl.sizet_rows == mat_argr.sizet_rows)
    && (mat_argl.sizet_cols == mat_argr.sizet_cols),
    "the matrices have different sizes.");

  if ((mat_argl.sizet_rows == mat_argr.sizet_rows) && (mat_argl.sizet_cols == mat_argr.sizet_cols))
  {
      for (size_t sizet_index = 0; sizet_index < sizet_size; sizet_index++)
      {
        mat_ret._matrix[sizet_index] = mat_argl._matrix[sizet_index] * mat_argr._matrix[sizet_index];
      }
  }

  return mat_ret;
}

template <typename T, typename Allocator>
matrix<T, Allocator> operator*( const matrix <T, Allocator> &mat_argl, T T_arg )
{
  matrix <T, Allocator> mat_ret(mat_argl.shape());
  size_t sizet_size = mat_argl.sizet_objects;

  for (size_t sizet_index = 0; sizet_index < sizet_size; sizet_index++)
  {
    mat_ret._matrix[sizet_index] = mat_argl._matrix[sizet_index] * T_arg;
  }

  return mat_ret;
}

template <typename T, typename Allocator>
matrix<T, Allocator> operator*( T T_arg, const matrix <T, Allocator> &mat_argr )
{
  matrix <T, Allocator> mat_ret(mat_argr.shape());
  size_t sizet_size = mat_argr.sizet_objects;

  for ( size_t sizet_index = 0; sizet_index < sizet_size; sizet_index++ )
  {
    mat_ret._matrix[sizet_index] = T_arg * mat_argr._matrix[sizet_index];
  }

  return mat_ret;
}

//  /
template <typename T, typename Allocator>
matrix<T, Allocator> operator/( const matrix <T, Allocator> &mat_argl, const matrix <T, Allocator> &mat_argr)
{
  matrix <T, Allocator> mat_ret(mat_argl.shape());
  size_t sizet_size = mat_argl.sizet_objects;

  SUSA_ASSERT_MESSAGE((mat_argl.sizet_rows == mat_argr.sizet_rows)
    && (mat_argl.sizet_cols == mat_argr.sizet_cols),
    "the matrices have different sizes.");

  if ((mat_argl.sizet_rows == mat_argr.sizet_rows) && (mat_argl.sizet_cols == mat_argr.sizet_cols))
  {
    for ( size_t sizet_index = 0; sizet_index < sizet_size; sizet_index++ )
    {
      mat_ret._matrix[sizet_index] = mat_argl._matrix[sizet_index] / mat_argr._matrix[sizet_index];
    }
  }

  return mat_ret;
}

template <typename T, typename Allocator>
matrix<T, Allocator> operator/( const matrix <T, Allocator> &mat_argl, T T_arg )
{

  matrix <T, Allocator> mat_ret(mat_argl.shape());
  size_t sizet_size = mat_argl.sizet_objects;

  for ( size_t sizet_index = 0; sizet_index < sizet_size; sizet_index++ )
  {
    mat_ret._matrix[sizet_index] = mat_argl._matrix[sizet_index] / T_arg;
  }

  return mat_ret;
}


//  -=
template <typename T, typename Allocator> matrix<T, Allocator> matrix <T, Allocator>::operator-=(const matrix <T, Allocator> &mat_arg)
{

    SUSA_ASSERT(_matrix != nullptr);

    SUSA_ASSERT_MESSAGE((mat_arg.sizet_rows == sizet_rows)
      && (mat_arg.sizet_cols == sizet_cols),
      "the matrices have different sizes.");

    if ((mat_arg.sizet_rows == sizet_rows) && (mat_arg.sizet_cols == sizet_cols) && _matrix != nullptr)
    {
        for (size_t sizet_index = 0; sizet_index < sizet_objects; sizet_index++)
        {
            _matrix[sizet_index] -= mat_arg._matrix[sizet_index];
        }
    }

    return *this;
}

//  +=
template <typename T, typename Allocator> matrix<T, Allocator> matrix <T, Allocator>::operator+=(const matrix <T, Allocator> &mat_arg)
{

  SUSA_ASSERT(_matrix != nullptr);
  if (_matrix == nullptr) return *this;

  SUSA_ASSERT_MESSAGE((mat_arg.sizet_rows == sizet_rows)
    && (mat_arg.sizet_cols == sizet_cols),
    "the matrices have different sizes.");

  if ((mat_arg.sizet_rows == sizet_rows) && (mat_arg.sizet_cols == sizet_cols))
  {
    for (size_t sizet_index = 0; sizet_index < sizet_objects; sizet_index++)
    {
      _matrix[sizet_index] += mat_arg._matrix[sizet_index];
    }
  }

  return *this;
}

//  =
template <typename T, typename Allocator>
matrix<T, Allocator>& matrix <T, Allocator>::operator=(const matrix <T, Allocator>& mat_arg)
{
  if (sizet_objects != 0 && _matrix != nullptr)
  {
    alloc.deallocate(_matrix, sizet_objects);
    _matrix         = nullptr;
    sizet_objects   = 0;
  }

  if (allocate(mat_arg.sizet_objects))
  {
    sizet_rows      = mat_arg.sizet_rows;
    sizet_cols      = mat_arg.sizet_cols;
    sizet_objects   = sizet_rows * sizet_cols;
    std::memcpy(_matrix, mat_arg._matrix, sizet_objects * sizeof(T));
  }

  return *this;
}

template <typename T, typename Allocator>
matrix<T, Allocator>& matrix <T, Allocator>::operator=(const std::string& str_string)
{
  parser(str_string);
  return *this;
}

template <typename T, typename Allocator>
matrix<T, Allocator>& matrix <T, Allocator>::operator=(const char char_string[])
{
  parser(std::string(char_string));
  return *this;
}

inline std::ostream& operator<<(std::ostream& os, char c)
{
  return std::is_signed<char>::value ? os << static_cast<int>(c): os << static_cast<unsigned int>(c);
}

inline std::ostream& operator<<(std::ostream& os, signed char c)
{
  return os << static_cast<int>(c);
}

inline std::ostream& operator<<(std::ostream& os, unsigned char c)
{
  return os << static_cast<unsigned int>(c);
}

inline std::ostream& operator<<(std::ostream& os, std::complex<int8_t> c)
{
  return os << "(" << static_cast<int>(c.real()) << "," << static_cast<int>(c.imag()) << ")";
}

inline std::ostream& operator<<(std::ostream& os, std::complex<uint8_t> c)
{
  return os << "(" << static_cast<unsigned int>(c.real()) << "," << static_cast<unsigned int>(c.imag()) << ")";
}

template <typename T, typename Allocator>
std::ostream &operator<<(std::ostream& out_stream, const matrix <T, Allocator>& mat_arg)
{
  out_stream << "[";
  for (size_t sizet_row = 0; sizet_row < mat_arg.sizet_rows; sizet_row++)
  {
    for (size_t sizet_col = 0; sizet_col < mat_arg.sizet_cols; sizet_col++)
    {
      out_stream << mat_arg(sizet_row, sizet_col);
      if (sizet_col < mat_arg.sizet_cols - 1) out_stream << " ";
    }
    if (sizet_row < mat_arg.sizet_rows - 1) out_stream << "\n ";
  }
  out_stream << "]";
  return out_stream;
}


template <typename T, typename Allocator>
bool matrix <T, Allocator>::parser(std::string str_string)
{
    pre_parser(str_string);

    size_t sizet_size  = 0;                     // linear size of the matrix (#rows * #columns)
    size_t sizet_cols_ = 0;                     // number of columns
    size_t sizet_rows_ = 0;                     // number of rows
    int int_length     = str_string.length();   // length of the input string


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

    SUSA_ASSERT_MESSAGE(sizet_cols_ % sizet_rows_ == 0, "the number of columns are not equal in each row.");
    if (sizet_cols_ % sizet_rows_ != 0) return false;

    if (sizet_size == 0)
    {
        sizet_rows      = 0;
        sizet_cols      = 0;
        sizet_objects   = 0;
        if (_matrix != nullptr) alloc.deallocate(_matrix, sizet_size);
        _matrix = nullptr;
        return false;
    }


    sizet_rows = sizet_rows_;
    sizet_cols = sizet_cols_ / sizet_rows_;

    sizet_size = sizet_cols * sizet_rows;

    if ((sizet_objects) != sizet_size)
    {
        if (sizet_objects != 0 && _matrix != nullptr)
        {
            alloc.deallocate(_matrix, sizet_objects);
            _matrix         = nullptr;
            sizet_objects   = 0;
        }

        if (allocate(sizet_size))
        {
            sizet_objects = sizet_size;
        }
    }

    std::stringstream ss_all(str_string);
    char char_buff[MAX_STR_LEN];

    std::conditional_t<std::is_same_v<T, uint8_t> || std::is_same_v<T, int8_t>, int, T> T_tmp;

    for (size_t sizet_row = 0; sizet_row < sizet_rows; sizet_row++)
    {
        ss_all.getline(char_buff, MAX_STR_LEN, ';');
        std::stringstream ssrow(char_buff);
        for (size_t sizet_col = 0; sizet_col < sizet_cols; sizet_col++)
        {
            ssrow.getline(char_buff, MAX_STR_LEN, 0x20);
            if (!(std::istringstream(char_buff) >> T_tmp))
            {
                SUSA_LOG_ERR("std::istringstream failed to parse the string");
                T_tmp = 0;
            }
            _matrix[get_lindex(sizet_row,sizet_col)] = T_tmp;
        }
    }

    return true;
}

template <typename T, typename Allocator>
bool matrix <T, Allocator>::allocate(size_t sz_objects)
{
  SUSA_ASSERT_MESSAGE(_matrix == nullptr, "not deallocated and may cause memory leak");
  try
  {
    _matrix = alloc.allocate(sz_objects);
  }
  catch(const std::bad_alloc& e)
  {
    _matrix       = nullptr;
    sizet_rows    = 0;
    sizet_cols    = 0;
    sizet_objects = 0;
    SUSA_LOG_ERR("memory allocation failed");
  }
  SUSA_ASSERT(_matrix != nullptr);
  return (_matrix != nullptr);
}
}      // NAMESPACE SUSA

#endif // SUSA_MATRIX_H
