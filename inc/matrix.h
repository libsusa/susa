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
 * @file matrix.h
 * @brief The matrix operation wrapper class.
 * This file contains the <i>matrix template class</i>. This is the only non-intrinsic data type that all
 * Susa methods work with.
 * @author Behrooz, Kamary Aliabadi
 * @version 1.0.0
 */

#ifndef MATRIX_H
#define MATRIX_H


namespace susa {

void pre_parser(std::string &str_string);

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

template <class T> matrix <T> mult( const matrix <T> &mat_argl,const matrix <T> &mat_argr);
template <class T> matrix <T> transpose(const matrix <T> &mat_arg);

/**
 * @brief Matrix template class
 *
 * @ingroup LALG
 *
 * @todo There is no feature to initialize/set elements of the matrix through an input string.
 */
template <class T> class matrix {
  private:
    unsigned int uint_rows;
    unsigned int uint_cols;

    unsigned int uint_id;
    unsigned int *uint_gc_index;
    static unsigned int uint_object_counter;


    T* _matrix;
    T Tfake;

    // UTILITY METHODS

    // This parses the strings and initialize the matrix elements.
    // The string initialization can be done using a constructor method
    // and an overloaded assignement operator (=) method
    void parser(std::string str_string);

  public:
    //! The default constructor
    matrix();

    /** The constructor
     * @param uint_row Number of rows
     * @param uint_col Number of columns
     * @param Tinitial Initial value of all rows and columns. If not set, initial value will be zero.
     */
    matrix( unsigned int uint_rows, unsigned int uint_cols, T Tinitial);

    /** The constructor
     * @param uint_row Number of rows
     * @param uint_col Number of columns.
     */
    matrix( unsigned int uint_rows, unsigned int uint_cols );

    //! The destructor
    ~matrix();

    //! Copy constructor
    matrix(const matrix <T> &mat_arg);

    matrix(std::string str_string);

    //! Returns the value of specific (row, column)
    T get( unsigned int uint_row, unsigned int uint_col ) const;

    //! Returns the value of specific (elem)
    T get( unsigned int uint_elem ) const;

    //! Returns the value of specific (elem)
    unsigned int get_id() const;

    //! Returns number of columns
    unsigned int no_cols() const;

    //! Returns number of rows
    unsigned int no_rows() const;

    //! Returns size of matrix
    unsigned int size() const;

    //! Returns size of matrix
    void set_all(T T_arg);

    //! Returns specified row
    matrix <T> row(unsigned int uint_row) const;

    //! Returns specified column
    matrix <T> col(unsigned int uint_col) const;

    //! Returns a matrix after eliminating the specified column and the row
    matrix <T> shrink(unsigned int uint_row, unsigned int uint_col) const;

    //! Considers the matrix object as a vector and return left side of that vector
    matrix <T> left(unsigned int uint_left) const;

    //! Considers the matrix object as a vector and return right side of that vector
    matrix <T> right(unsigned int uint_right) const;

    //! Considers the matrix object as a vector and return mid part of that vector
    matrix <T> mid(unsigned int uint_begin, unsigned int uint_end) const;

    //! Elementwise Assignment by Addition operator
    matrix <T> operator+=( const matrix <T> &mat_arg );

    //! Elementwise Assignment by Substraction
    matrix <T> operator-=( const matrix <T> &mat_arg );

    //! Elementwise Assignment operator
    matrix <T>& operator=( const matrix <T> &mat_arg );

    //! Elementwise Assignment operator
    matrix <T> operator=( const matrix <T> &mat_arg ) const;

    //! Elementwise Assignment operator
    matrix <T>& operator=( std::string str_string );

    //! Elementwise Substraction operator
    friend matrix <T> operator-<>( const matrix <T> &mat_argl, T T_arg);

    //! Elementwise Substraction operator
    friend matrix <T> operator-<>( T T_arg, const matrix <T> &mat_argr );

    //! Elementwise Substraction operator
    friend matrix <T> operator-<>( const matrix <T> &mat_argl, const matrix <T> &mat_argr);


    //! Elementwise Addition operator
    friend matrix <T> operator+<>( const matrix <T> &mat_argl, T T_arg );

    //! Elementwise Addition operator
    friend matrix <T> operator+<>( T T_arg, const matrix <T> &mat_argr );

    //! Elementwise Addition operator
    friend matrix <T> operator+<>( const matrix <T> &mat_argl, const matrix <T> &mat_argr);


    //! Elementwise Multiplication  operator
    friend matrix <T> operator*<>( const matrix <T> &mat_argl, T T_arg);

    //! Elementwise Multiplication  operator
    friend matrix <T> operator*<>( T T_arg, const matrix <T> &mat_argr );

    //! Elementwise Multiplication  operator
    friend matrix <T> operator*<>( const matrix <T> &mat_argl, const matrix <T> &mat_argr );


    //! Elementwise Division operator
    friend matrix <T> operator/<>( const matrix <T> &mat_argl, T T_arg);

    //! Elementwise Division operator
    friend matrix <T> operator/<>( T T_arg, const matrix <T> &mat_argr );

    //! Elementwise Division operator
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
    // This because they can access '_matrix' directly.
    // They are defined in 'linalg.h'
    friend matrix <T> mult<>( const matrix <T> &mat_argl,const matrix <T> &mat_argr);
    friend matrix <T> transpose<>(const matrix <T> &mat_arg);
};

// Static constants
template <class T> unsigned int matrix <T>::uint_object_counter = 100;

// Constructor and Destructor

template <class T> matrix <T>::matrix() {
    this->uint_rows = 0;
    this->uint_cols = 0;

    _matrix = NULL;

    Tfake = 0;

    uint_gc_index = new unsigned int;
    *uint_gc_index = 0;

    uint_id = uint_object_counter;
    uint_object_counter++;

}

template <class T> matrix <T>::matrix(unsigned int uint_rows, unsigned int uint_cols, T Tinitial) {
    this->uint_rows = uint_rows < 2 ? 1 : uint_rows;
    this->uint_cols = uint_cols < 2 ? 1 : uint_cols;

    try {
        _matrix = new T [this->uint_rows * this->uint_cols];
        for (unsigned int uint_counter = 0; uint_counter < this->uint_rows * this->uint_cols; uint_counter++) _matrix[uint_counter] = Tinitial;
    } catch ( std::bad_alloc ex) {
        std::cout << std::endl << "[matrix::matrix()] memory allocation exception.";
        exit(0);
    }

    Tfake = 0;

    uint_gc_index = new unsigned int;
    *uint_gc_index = 0;

    uint_id = uint_object_counter;
    uint_object_counter++;

}

template <class T> matrix <T>::matrix( unsigned int uint_rows, unsigned int uint_cols ) {
    this->uint_rows = uint_rows < 2 ? 1 : uint_rows;
    this->uint_cols = uint_cols < 2 ? 1 : uint_cols;

    try {
        _matrix = new T [this->uint_rows * this->uint_cols];
        for (unsigned int uint_counter = 0; uint_counter < this->uint_rows * this->uint_cols; uint_counter++) _matrix[uint_counter] = 0;
    } catch ( std::bad_alloc ex) {
        std::cout << std::endl << "[matrix::matrix()] memory allocation exception.";
        exit(0);
    }
    Tfake = 0;

    uint_gc_index = new unsigned int;
    *uint_gc_index = 0;

    uint_id = uint_object_counter;
    uint_object_counter++;

}

template <class T> matrix <T>::matrix(const matrix <T> &mat_arg) {
    unsigned int uint_size = mat_arg.uint_rows * mat_arg.uint_cols;

    if (uint_size != 0) {

        try {
            _matrix = new T [uint_size];
        } catch ( std::bad_alloc ex) {
            std::cout << std::endl << "[matrix::matrix()] memory allocation exception.";
            exit(0);
        }

        uint_rows = mat_arg.uint_rows;
        uint_cols = mat_arg.uint_cols;

        for (unsigned int uint_index = 0; uint_index < uint_size; uint_index++) {
            _matrix[uint_index] = mat_arg._matrix[uint_index];
        }
    } else {
        uint_rows = 0;
        uint_cols = 0;

        if (_matrix != NULL) {
            delete [] _matrix;
            _matrix = NULL;
        }
    }

    Tfake = 0;

    uint_gc_index = new unsigned int;
    *uint_gc_index = 0;

    uint_id = uint_object_counter;
    uint_object_counter++;

}

template <class T> matrix <T>::matrix(std::string str_string) {

    this->uint_rows = 0;
    this->uint_cols = 0;

    _matrix = NULL;

    parser(str_string);
    Tfake = 0;

    uint_gc_index = new unsigned int;
    *uint_gc_index = 0;

    uint_id = uint_object_counter;
    uint_object_counter++;
}


template <class T> matrix <T>::~matrix() {

    if (_matrix != NULL) {
        delete [] _matrix;
        _matrix = NULL;
    }

    if (uint_gc_index != NULL) {
        delete uint_gc_index;
        uint_gc_index = NULL;
    }

    uint_object_counter--;
}

// Public methods

template <class T> T matrix <T>::get( unsigned int uint_row, unsigned int uint_col ) const {

    if (uint_row < uint_rows && uint_col < uint_cols && _matrix != NULL)
        return _matrix[uint_row + uint_col * uint_rows];
    else {
        if (_matrix == NULL) std::cout << std::endl << "[matrix::get()] NULL memory allocation.";
        else {
            std::cout << std::endl <<"(" << uint_row << "," << uint_col << ")  "<<"[matrix::get()] matrix index(es) are out of range. - " << uint_id;
        }
        return Tfake;
    }

}

template <class T> T matrix <T>::get( unsigned int uint_elem ) const {

    if (uint_elem < uint_rows * uint_cols && _matrix != NULL)
        return _matrix[uint_elem];
    else {
        if (_matrix == NULL) std::cout << std::endl << "[matrix::get()] NULL memory allocation.";
        else {
            std::cout << std::endl <<"(" << uint_elem << ")  "<<"[matrix::get()] matrix index(es) are out of range. - " << uint_id;
        }
        return Tfake;
    }

}

template <class T> unsigned int matrix <T>::get_id() const {
    return uint_id;
}

template <class T> unsigned int  matrix <T>::no_cols() const {
    return uint_cols;
}

template <class T> unsigned int  matrix <T>::no_rows() const {
    return uint_rows;
}

template <class T> unsigned int  matrix <T>::size() const {
    return (uint_rows * uint_cols);
}

template <class T> void  matrix <T>::set_all(T T_arg) {
    for (unsigned int uint_counter = 0; uint_counter < this->uint_rows * this->uint_cols; uint_counter++) _matrix[uint_counter] = T_arg;
}

template <class T> matrix <T> matrix <T>::row(unsigned int uint_row) const {
    matrix <T> mat_ret(1,uint_cols);

    if (uint_row < uint_rows) {
        for (unsigned int uint_i = 0; uint_i < uint_cols; uint_i++) mat_ret(uint_i) = _matrix[uint_i*uint_rows + uint_row];
    }

    return mat_ret;
}

template <class T> matrix <T> matrix <T>::col(unsigned int uint_col) const {
    matrix <T> mat_ret(uint_rows,1);

    if (uint_col < uint_cols) {
        for (unsigned int uint_i = 0; uint_i < uint_rows; uint_i++) mat_ret(uint_i) = _matrix[uint_col*uint_rows + uint_i];
    }

    return mat_ret;
}

template <class T> matrix <T> matrix <T>::shrink(unsigned int uint_elim_row, unsigned int uint_elim_col) const {

    matrix <T> mat_ret;

    unsigned int uint_new_col;
    unsigned int uint_new_row;

    if (uint_cols > 1 && uint_rows > 1 && uint_elim_row >= 0 && uint_elim_col >= 0) {
        mat_ret = matrix <T> (uint_rows - 1,uint_cols - 1);
    } else {
        std::cout << std::endl << "[matrix::minor()] Input arguments error.";
        return *this;
    }

    if (uint_elim_col < uint_cols && uint_elim_row < uint_rows) {
        for ( unsigned int uint_row = 0; uint_row < uint_rows; uint_row++ ) {
            for ( unsigned int uint_col = 0; uint_col < uint_cols; uint_col++ ) {
                if (uint_col != uint_elim_col || uint_row != uint_elim_row) {
                    uint_new_row = uint_row > uint_elim_row ? (uint_row - 1) : uint_row;
                    uint_new_col = uint_col > uint_elim_col ? (uint_col - 1) : uint_col;
                    if (uint_new_row < (uint_rows - 1) && uint_new_col < (uint_cols - 1)) {
                        mat_ret(uint_new_row, uint_new_col) = _matrix[uint_col*uint_rows + uint_row];
                    }
                }
            }
        }
    } else {
        std::cout << std::endl << "[matrix::shrink()] Input arguments exceed matrix size.";
        return *this;
    }

    return mat_ret;
}

template <class T> matrix <T> matrix <T>::left(unsigned int uint_left) const {
    matrix <T> mat_ret;
    unsigned int uint_size = uint_rows * uint_cols;
    if ( uint_size >= uint_left ) {
        if (uint_rows == 1 && uint_cols != 1) mat_ret = matrix <T> (1,uint_left);
        else if (uint_rows != 1 && uint_cols == 1) mat_ret = matrix <T> (uint_left,1);
        else if (uint_rows != 1 && uint_cols != 1) mat_ret = matrix <T> (uint_left,1);

        for (unsigned int uint_i = 0; uint_i < uint_left; uint_i++) mat_ret(uint_i) = _matrix[uint_i];
    } else {
        std::cout << std::endl << "[matrix::left()] Error.";
    }

    return mat_ret;

}


template <class T> matrix <T> matrix <T>::right(unsigned int uint_right) const {
    matrix <T> mat_ret;
    unsigned int uint_size = uint_rows * uint_cols;
    if ( uint_size >= uint_right ) {
        if (uint_rows == 1 && uint_cols != 1) mat_ret = matrix <T> (1,uint_right);
        else if (uint_rows != 1 && uint_cols == 1) mat_ret = matrix <T> (uint_right,1);
        else if (uint_rows != 1 && uint_cols != 1) mat_ret = matrix <T> (uint_right,1);

        for (unsigned int uint_i = uint_size; uint_i > (uint_size - uint_right); uint_i--)
            mat_ret(uint_right - uint_size + uint_i - 1) = _matrix[uint_i - 1];
    } else {
        std::cout << std::endl << "[matrix::right()] Error.";
    }

    return mat_ret;
}


template <class T> matrix <T> matrix <T>::mid(unsigned int uint_begin, unsigned int uint_end) const {
    matrix <T> mat_ret;
    unsigned int uint_size = uint_rows * uint_cols;
    unsigned int uint_mid = uint_end - uint_begin;

    if ( uint_size > uint_begin && uint_size > uint_end && uint_end > uint_begin) {
        if (uint_rows == 1 && uint_cols != 1) mat_ret = matrix <T> (1,uint_mid + 1);
        else if (uint_rows != 1 && uint_cols == 1) mat_ret = matrix <T> (uint_mid + 1,1);
        else if (uint_rows != 1 && uint_cols != 1) mat_ret = matrix <T> (uint_mid + 1,1);

        for (unsigned int uint_i = uint_begin; uint_i <= uint_end; uint_i++) mat_ret(uint_i - uint_begin) = _matrix[uint_i];
    } else {
        std::cout << std::endl << "[matrix::mid()] Error.";
    }

    return mat_ret;
}

// Operators

//  ()
template <class T> T &matrix<T>::operator ()( unsigned int uint_row, unsigned int uint_col ) {

    if (uint_row < uint_rows && uint_col < uint_cols && _matrix != NULL)
        return _matrix[uint_row + uint_col * uint_rows];
    else {
        if (_matrix == NULL) std::cout << std::endl << "[matrix::operator ()] NULL memory allocation.";
        else {
            std::cout << std::endl <<"(" << uint_row << "," << uint_col << ")  "<<"[matrix::operator ()] matrix index(es) are out of range. - " << uint_id;
        }
        return Tfake;
    }
}

template <class T> T matrix<T>::operator ()( unsigned int uint_row, unsigned int uint_col ) const {

    if (uint_row < uint_rows && uint_col < uint_cols && _matrix != NULL)
        return _matrix[uint_row + uint_col * uint_rows];
    else {
        if (_matrix == NULL) std::cout << std::endl << "[matrix::operator () const] NULL memory allocation.";
        else {
            std::cout << std::endl <<"(" << uint_row << "," << uint_col << ")  "<<"[matrix::operator () const] matrix index(es) are out of range.";
        }
        return Tfake;
    }
}

template <class T> T &matrix<T>::operator ()( unsigned int uint_elem ) {

    if (uint_elem < uint_cols * uint_rows && _matrix != NULL)
        return _matrix[uint_elem];
    else {
        if (_matrix == NULL) std::cout << std::endl << "[matrix::operator ()] NULL memory allocation.";
        else {
            std::cout << std::endl <<"(" << uint_elem << ")  "<<"[matrix::operator ()] matrix index(es) are out of range. - " << uint_id;
        }
        return Tfake;
    }
}

template <class T> T matrix<T>::operator ()( unsigned int uint_elem ) const {

    if (uint_elem < uint_cols * uint_rows && _matrix != NULL)
        return _matrix[uint_elem];
    else {
        if (_matrix == NULL) std::cout << std::endl << "[matrix::operator () const] NULL memory allocation.";
        else {
            std::cout << std::endl <<"(" << uint_elem << ")  " << "[matrix::operator () const] matrix index(es) are out of range. - " << uint_id;
        }
        return Tfake;
    }
}

// Typecasting
template <class T> matrix<T>::operator matrix <double> () {

    matrix <double> mat_ret(uint_rows, uint_cols);

    unsigned int uint_size = uint_rows * uint_cols;

    for (unsigned int uint_elem = 0; uint_elem < uint_size; uint_elem++) {
        mat_ret(uint_elem) = (double)_matrix[uint_elem];
    }

    return mat_ret;
}

template <class T> matrix<T>::operator matrix <float> () {

    matrix <float> mat_ret(uint_rows, uint_cols);

    unsigned int uint_size = uint_rows * uint_cols;

    for (unsigned int uint_elem = 0; uint_elem < uint_size; uint_elem++) {
        mat_ret(uint_elem) = (float)_matrix[uint_elem];
    }
    return mat_ret;
}

template <class T> matrix<T>::operator matrix <int> () {

    matrix <int> mat_ret(uint_rows, uint_cols);

    unsigned int uint_size = uint_rows * uint_cols;

    for (unsigned int uint_elem = 0; uint_elem < uint_size; uint_elem++) {
        mat_ret(uint_elem) = (int)_matrix[uint_elem];
    }
    return mat_ret;
}


template <class T> matrix<T>::operator matrix <char> () {

    matrix <char> mat_ret(uint_rows, uint_cols);

    unsigned int uint_size = uint_rows * uint_cols;

    for (unsigned int uint_elem = 0; uint_elem < uint_size; uint_elem++) {
        mat_ret(uint_elem) = (char)_matrix[uint_elem];
    }
    return mat_ret;
}


template <class T> matrix<T>::operator matrix <std::complex <double> > () {

    matrix <std::complex <double> > mat_ret(uint_rows, uint_cols);

    unsigned int uint_size = uint_rows * uint_cols;

    for (unsigned int uint_elem = 0; uint_elem < uint_size; uint_elem++) {
        mat_ret(uint_elem) = std::complex <double> ((double)_matrix[uint_elem].real(), (double)_matrix[uint_elem].imag());
    }
    return mat_ret;
}

template <class T> matrix<T>::operator matrix <std::complex <float> > () {

    matrix <std::complex <float> > mat_ret(uint_rows, uint_cols);

    unsigned int uint_size = uint_rows * uint_cols;

    for (unsigned int uint_elem = 0; uint_elem < uint_size; uint_elem++) {
        mat_ret(uint_elem) = std::complex <float> ((float)_matrix[uint_elem].real(), (float)_matrix[uint_elem].imag());
    }
    return mat_ret;
}

template <class T> matrix<T>::operator matrix <std::complex <int> > () {

    matrix <std::complex <int> > mat_ret(uint_rows, uint_cols);

    unsigned int uint_size = uint_rows * uint_cols;

    for (unsigned int uint_elem = 0; uint_elem < uint_size; uint_elem++) {
        mat_ret(uint_elem) = std::complex <int> ((int)_matrix[uint_elem].real(), (int)_matrix[uint_elem].imag());
    }
    return mat_ret;
}



template <class T> matrix<T>::operator matrix <std::complex <char> > () {

    matrix <std::complex <char> > mat_ret(uint_rows, uint_cols);

    unsigned int uint_size = uint_rows * uint_cols;

    for (unsigned int uint_elem = 0; uint_elem < uint_size; uint_elem++) {
        mat_ret(uint_elem) = std::complex <char> ((char)_matrix[uint_elem].real(), (char)_matrix[uint_elem].imag());
    }
    return mat_ret;
}


//  +
template <class T> matrix<T> operator+(const matrix <T> &mat_argl, const matrix <T> &mat_argr) {

    matrix <T> mat_temp(mat_argl.uint_rows, mat_argl.uint_cols);
    unsigned int uint_size = mat_argl.uint_rows * mat_argl.uint_cols;

    if ((mat_argl.uint_rows == mat_argr.uint_rows) && (mat_argl.uint_cols == mat_argr.uint_cols)) {
        for ( unsigned int uint_index = 0; uint_index < uint_size; uint_index++ ) {
            mat_temp._matrix[uint_index] = mat_argl._matrix[uint_index] + mat_argr._matrix[uint_index];
        }
    } else {
        std::cout << "\n [matrix::operator +()] Matrices have different sizes";
    }
    return mat_temp;
}

template <class T> matrix<T> operator+( const matrix <T> &mat_argl, T T_arg ) {

    matrix <T> mat_temp(mat_argl.uint_rows, mat_argl.uint_cols);
    unsigned int uint_size = mat_argl.uint_rows * mat_argl.uint_cols;

    for ( unsigned int uint_index = 0; uint_index < uint_size; uint_index++ ) {
        mat_temp._matrix[uint_index] = mat_argl._matrix[uint_index] + T_arg;
    }

    return mat_temp;
}

template <class T> matrix<T> operator+( T T_arg, const matrix <T> &mat_argr ) {

    matrix <T> mat_temp(mat_argr.uint_rows, mat_argr.uint_cols);
    unsigned int uint_size = mat_argr.uint_rows * mat_argr.uint_cols;

    for ( unsigned int uint_index = 0; uint_index < uint_size; uint_index++ ) {
        mat_temp._matrix[uint_index] = mat_argr._matrix[uint_index] + T_arg;
    }

    return mat_temp;
}



//  -
template <class T> matrix<T> operator-( const matrix <T> &mat_argl, const matrix <T> &mat_argr) {

    matrix <T> mat_temp(mat_argl.uint_rows, mat_argl.uint_cols);
    unsigned int uint_size = mat_argl.uint_rows * mat_argl.uint_cols;

    if ((mat_argl.uint_rows == mat_argr.uint_rows) && (mat_argl.uint_cols == mat_argr.uint_cols)) {
        for ( unsigned int uint_index = 0; uint_index < uint_size; uint_index++ ) {
            mat_temp._matrix[uint_index] = mat_argl._matrix[uint_index] - mat_argr._matrix[uint_index];
        }
    } else {
        std::cout << "\n [matrix::operator -()] Matrices have different sizes";
    }
    return mat_temp;
}

template <class T> matrix<T> operator-( const matrix <T> &mat_argl, T T_arg ) {

    matrix <T> mat_temp(mat_argl.uint_rows, mat_argl.uint_cols);
    unsigned int uint_size = mat_argl.uint_rows * mat_argl.uint_cols;

    for ( unsigned int uint_index = 0; uint_index < uint_size; uint_index++ ) {
        mat_temp._matrix[uint_index] = mat_argl._matrix[uint_index] - T_arg;
    }

    return mat_temp;
}

template <class T> matrix<T> operator-( T T_arg, const matrix <T> &mat_argr ) {

    matrix <T> mat_temp(mat_argr.uint_rows, mat_argr.uint_cols);
    unsigned int uint_size = mat_argr.uint_rows * mat_argr.uint_cols;

    for ( unsigned int uint_index = 0; uint_index < uint_size; uint_index++ ) {
        mat_temp._matrix[uint_index] = T_arg - mat_argr._matrix[uint_index];
    }

    return mat_temp;
}

//  *
template <class T> matrix<T> operator*( const matrix <T> &mat_argl, const matrix <T> &mat_argr) {

    matrix <T> mat_temp(mat_argl.uint_rows, mat_argl.uint_cols);
    unsigned int uint_size = mat_argl.uint_rows * mat_argl.uint_cols;

    if ((mat_argl.uint_rows == mat_argr.uint_rows) && (mat_argl.uint_cols == mat_argr.uint_cols)) {
        for ( unsigned int uint_index = 0; uint_index < uint_size; uint_index++ ) {
            mat_temp._matrix[uint_index] = mat_argl._matrix[uint_index] * mat_argr._matrix[uint_index];
        }
    } else {
        std::cout << "\n [matrix::operator *()] Matrices have different sizes";
    }
    return mat_temp;
}

template <class T> matrix<T> operator*( const matrix <T> &mat_argl, T T_arg ) {

    matrix <T> mat_temp(mat_argl.uint_rows, mat_argl.uint_cols);
    unsigned int uint_size = mat_argl.uint_rows * mat_argl.uint_cols;

    for ( unsigned int uint_index = 0; uint_index < uint_size; uint_index++ ) {
        mat_temp._matrix[uint_index] = mat_argl._matrix[uint_index] * T_arg;
    }

    return mat_temp;
}

template <class T> matrix<T> operator*( T T_arg, const matrix <T> &mat_argr ) {

    matrix <T> mat_temp(mat_argr.uint_rows, mat_argr.uint_cols);
    unsigned int uint_size = mat_argr.uint_rows * mat_argr.uint_cols;

    for ( unsigned int uint_index = 0; uint_index < uint_size; uint_index++ ) {
        mat_temp._matrix[uint_index] = T_arg * mat_argr._matrix[uint_index];
    }

    return mat_temp;
}

//  /
template <class T> matrix<T> operator/( const matrix <T> &mat_argl, const matrix <T> &mat_argr) {

    matrix <T> mat_temp(mat_argl.uint_rows, mat_argl.uint_cols);
    unsigned int uint_size = mat_argl.uint_rows * mat_argl.uint_cols;

    if ((mat_argl.uint_rows == mat_argr.uint_rows) && (mat_argl.uint_cols == mat_argr.uint_cols)) {
        for ( unsigned int uint_index = 0; uint_index < uint_size; uint_index++ ) {
            mat_temp._matrix[uint_index] = mat_argl._matrix[uint_index] / mat_argr._matrix[uint_index];
        }
    } else {
        std::cout << "\n [matrix::operator /()] Matrices have different sizes";
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
    return mat_temp;
}
//  -=
template <class T> matrix<T> matrix <T>::operator-=(const matrix <T> &mat_arg) {

    unsigned int uint_size = uint_cols * uint_rows;;

    if ((mat_arg.uint_rows == uint_rows) && (mat_arg.uint_cols == uint_cols) && _matrix != NULL) {
        for (unsigned int uint_index = 0; uint_index < uint_size; uint_index++) {
            _matrix[uint_index] -= mat_arg._matrix[uint_index] ;
        }
    } else {
        std::cout << std::endl << "[matrix::operator -=] matrix dimension mismatch.";
    }
    return *this;
}

//  +=
template <class T> matrix<T> matrix <T>::operator+=(const matrix <T> &mat_arg) {

    unsigned int uint_size = uint_cols * uint_rows;;

    if ((mat_arg.uint_rows == uint_rows) && (mat_arg.uint_cols == uint_cols) && _matrix != NULL) {
        for (unsigned int uint_index = 0; uint_index < uint_size; uint_index++) {
            _matrix[uint_index] += mat_arg._matrix[uint_index] ;
        }
    } else {
        std::cout << std::endl << "[matrix::operator +=] matrix dimension mismatch.";
    }
    return *this;
}

//  =
template <class T> matrix<T>& matrix <T>::operator=( const matrix <T> &mat_arg ) {

    unsigned int uint_size = mat_arg.uint_cols * mat_arg.uint_rows;

    if (this != &mat_arg) {
        if (uint_size != 0) {

            if ((uint_rows * uint_cols) != uint_size) {

                if (_matrix != NULL) {
                    delete [] _matrix;
                    _matrix = NULL;
                }

                try {
                    _matrix = new T [uint_size];
                } catch ( std::bad_alloc ex) {
                    std::cout << std::endl << "[matrix::operator =] memory allocation exception.";
                    exit(0);
                }

            }

            uint_rows = mat_arg.uint_rows;
            uint_cols = mat_arg.uint_cols;

            for (unsigned int uint_index = 0; uint_index < uint_size; uint_index++) {
                _matrix[uint_index] = mat_arg._matrix[uint_index];
            }
        } else {
            uint_rows = 0;
            uint_cols = 0;

            if (_matrix != NULL) {
                delete [] _matrix;
                _matrix = NULL;
            }
        }
    }


    return *this;
}

template <class T> matrix<T> matrix <T>::operator=( const matrix <T> &mat_arg ) const {

    unsigned int uint_size = mat_arg.uint_cols * mat_arg.uint_rows;

    if (this != &mat_arg) {
        if (uint_size != 0) {

            if ((uint_rows * uint_cols) != uint_size) {

                if (_matrix != NULL) {
                    delete [] _matrix;
                    _matrix = NULL;
                }

                try {
                    _matrix = new T [uint_size];
                } catch ( std::bad_alloc ex) {
                    std::cout << std::endl << "[matrix::operator = const] memory allocation exception.";
                    exit(0);
                }

            }

            uint_rows = mat_arg.uint_rows;
            uint_cols = mat_arg.uint_cols;

            for (unsigned int uint_index = 0; uint_index < uint_size; uint_index++) {
                _matrix[uint_index] = mat_arg._matrix[uint_index];
            }
        } else {
            uint_rows = 0;
            uint_cols = 0;

            if (_matrix != NULL) {
                delete [] _matrix;
                _matrix = NULL;
            }
        }
    }


    return *this;
}

template <class T> matrix<T>& matrix <T>::operator=( std::string str_string ) {

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

template <class T> void matrix <T>::parser(std::string str_string) {

    pre_parser(str_string);

    unsigned int uint_size = 0; // linear size of the matrix (#rows * #columns)
    int int_length = str_string.length(); // length of the input string
    unsigned int uint_cols_ = 0; // number of columns
    unsigned int uint_rows_ = 0; // number of rows

    // This loop gets total number of rows and columns
    for (int int_i = 0; int_i < int_length; int_i++) {
        switch (str_string[int_i]) {
        case ';':
            uint_rows_++;
            // since the number of rows are counted by number of space and the
            // spaces around ';' are cut, we count a column here too !
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
    uint_rows_++;
    uint_cols_ = uint_cols_ + 1;
    uint_size = uint_cols_ * uint_rows_;

    if (uint_cols_ % uint_rows_ != 0) {
        std::cout << std::endl << "[matrix::parser] number of columns are not equal in each row." <<std::endl;
        return;
    }


    if (uint_size != 0) {

        if ((uint_rows * uint_cols) != uint_size) {

            if (_matrix != NULL) {
                delete [] _matrix;
                _matrix = NULL;
            }

            try {
                _matrix = new T [uint_size];
            } catch ( std::bad_alloc ex) {
                std::cout << std::endl << "[matrix::operator = const] memory allocation exception.";
                exit(0);
            }

        }

        uint_rows = uint_rows_;
        uint_cols = uint_cols_ / uint_rows_;


        std::stringstream ss_all(str_string);
        char* char_buff;
        char_buff = new char [4096];

        for (unsigned int uint_row = 0; uint_row < uint_rows; uint_row++) {
            ss_all.getline(char_buff,4096,';');
            std::stringstream ssrow(char_buff);
            for (unsigned int uint_col = 0; uint_col < uint_cols; uint_col++) {
                ssrow.getline(char_buff,4096,' ');
                std::stringstream(char_buff) >> _matrix[uint_row + uint_col * uint_rows];
            }
        }

    } else {
        uint_rows = 0;
        uint_cols = 0;

        if (_matrix != NULL) {
            delete [] _matrix;
            _matrix = NULL;
        }
    }
}
}      // NAMESPACE SUSA
#endif // MATRIX_H

