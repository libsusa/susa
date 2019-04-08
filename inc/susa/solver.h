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
 * @file solver.h
 * @brief Linear solvers and algorithms
 * @author Behrooz Kamary
 * @version 1.0.0
 *
 * @defgroup LALG
 */

#ifndef SUSA_SOLVER_H
#define SUSA_SOLVER_H

#include <susa.h>

namespace susa
{

/**
 * @brief LU decomposition algebraic methods and operations
 * @ingroup LALG
 */
template <class T> class lu 
{
  public:

    /**
     * @brief constructor
     *
     * @param mat_arg the input square matrix
     */
    lu(const matrix<T>& mat_arg);

    //! destructor
    ~lu();

    /**
     * @brief decomposition
     *
     * Doolittle algorithm implementation
     * @param mat_arg the input square matrix
     * @return true if the decomposition succeeds
     */
    bool                decompose();

    //! get the upper decomposed matrix
    const matrix<T>&    get_upper() const { return mat_upper; }

    //! get the lower decomposed matrix
    const matrix<T>&    get_lower() const { return mat_lower; }

  private:
    const matrix <T>&   mat_a;
    matrix <T>          mat_upper;
    matrix <T>          mat_lower;

    size_t              sizet_n;

};


/**
 * @brief LUP decomposition algebraic methods and operations
 * @ingroup LALG
 */
template <class T> class lup
{
  public:

    /**
     * @brief constructor
     *
     * @param mat_arg the input square matrix
     * @param dbl_tolerance the tolerance of pivoting
     */
    lup(const matrix<T>& mat_arg, double dbl_tolerance = 1e-4);

    //! destructor
    ~lup();

    /**
     * @brief decomposition
     *
     * Doolittle algorithm implementation
     * @param mat_arg the input square matrix
     * @return true if the decomposition succeeds
     */
    bool                decompose();

    //! compute the pivot and LU matrix
    bool                compute_pivot();

    /**
     * @brief solve a linear equation of the form Ax = b
     *
     * @param mat_b is the input vector
     * @return the solution to the linear equation
     */
    matrix <T>          solve(const matrix<T>& mat_b);

    //! get the linear pivot vector
    const matrix<T>&    get_pivot() const { return mat_pp; }

    //! get the LU matrix after calling lup::decompose or lup::compute_pivot
    const matrix<T>&    get_lu() const {return mat_lu;}

  private:
    matrix <T>          mat_lu;
    const matrix <T>&   mat_a;
    matrix <size_t>     mat_p;
    matrix <T>          mat_pp;

    double              dbl_tol;
    size_t              sizet_n;

};

template <class T> lu<T>::lu(const matrix<T>& mat_arg)
: mat_a(mat_arg)
, mat_upper(mat_arg.shape(), 0)
, mat_lower(mat_arg.shape(), 0)
, sizet_n(0)
{
    SUSA_ASSERT_MESSAGE(mat_arg.is_square(), "the matrix must be square");
}

template <class T> lu<T>::~lu()
{

}

template <class T> lup<T>::lup(const matrix<T>& mat_arg, double dbl_tolerance)
: mat_a(mat_arg)
, dbl_tol(dbl_tolerance)
, sizet_n(0)
{
    SUSA_ASSERT_MESSAGE(mat_arg.is_square(), "the matrix must be square");
}

template <class T> lup<T>::~lup()
{

}

template <class T> bool lu<T>::decompose()
{
    if (!mat_a.is_square()) return false;

    sizet_n = mat_a.no_cols();

    for (size_t i = 0; i < sizet_n; i++)
    {

        for (size_t k = i; k < sizet_n; k++)
        {

            T sum = 0;
            for (size_t j = 0; j < i; j++)
            {
                sum += mat_lower(i,j) * mat_upper(j,k);
            }

            mat_upper(i,k) = mat_a(i,k) - sum;
        }

        for (size_t k = i; k < sizet_n; k++)
        {
            if (i == k)
            {
                mat_lower(i,i) = 1;
            }
            else
            {
                T sum = 0;
                for (size_t j = 0; j < i; j++)
                {
                    sum += mat_lower(k,j) * mat_upper(j,i);
                }

                mat_lower(k,i) = (mat_a(k,i) - sum) / mat_upper(i,i);
            }
        }
    }

    return true;
}

template <class T> bool lup<T>::decompose()
{
    if (!mat_a.is_square()) return false;


    sizet_n     = mat_a.no_cols();
    mat_p       = matrix<size_t>(sizet_n + 1, 1);

    mat_lu = mat_a;

    for (size_t indx = 0; indx <= sizet_n; indx++)
    {
        mat_p(indx) = indx;
    }

    for (size_t sizet_col = 0; sizet_col < sizet_n; sizet_col++)
    {
        T T_max = 0;
        T T_abs = 0;
        size_t sizet_max_row = sizet_col;

        for (size_t sizet_row = sizet_col; sizet_row < sizet_n; sizet_row++)
        {
            T_abs = mat_a(sizet_row, sizet_col);
            if (T_abs < 0) T_abs *= -1.0;
            if ( T_abs > T_max)
            {
                T_max = T_abs;
                sizet_max_row = sizet_row;
            }
        }

        if (T_max < dbl_tol) return false;

        if (sizet_max_row != sizet_col)
        {
            std::swap(mat_p(sizet_max_row), mat_p(sizet_col));
            mat_lu.swap_rows(sizet_max_row, sizet_col);
            mat_p(sizet_n)++;
        }

        for (size_t sizet_row = sizet_col + 1; sizet_row < sizet_n; sizet_row++)
        {
            mat_lu(sizet_row, sizet_col) /= mat_lu(sizet_col, sizet_col);

            for (size_t sizet_col_in = sizet_col + 1; sizet_col_in < sizet_n; sizet_col_in++)
            {
                mat_lu(sizet_row, sizet_col_in) -= mat_lu(sizet_row, sizet_col) * mat_lu(sizet_col, sizet_col_in);
            }
        }
    }

    return true;
}

template <class T> matrix<T> lup<T>::solve(const matrix<T> &mat_b)
{
    matrix<T> mat_ret(sizet_n, 1, 0);

    if (sizet_n == 0 || !mat_b.is_vector() || mat_b.size() != sizet_n) return mat_ret;

    for (size_t sizet_row = 0; sizet_row < sizet_n; sizet_row++)
    {
        mat_ret(sizet_row) = mat_b(mat_p(sizet_row));

        for (size_t sizet_col = 0; sizet_col < sizet_row; sizet_col++)
        {
            mat_ret(sizet_row) -= mat_lu(sizet_row, sizet_col) * mat_ret(sizet_col);
        }
    }


    for (size_t sizet_row = sizet_n; sizet_row-- > 0;)
    {

        for (size_t sizet_col = sizet_row + 1; sizet_col < sizet_n; sizet_col++)
        {
            mat_ret(sizet_row) -= mat_lu(sizet_row, sizet_col) * mat_ret(sizet_col);
        }

        mat_ret(sizet_row) = mat_ret(sizet_row) / mat_lu(sizet_row, sizet_row);
    }

    return mat_ret;
}


template <class T> bool lup<T>::compute_pivot()
{
    sizet_n     = mat_a.no_cols();
    mat_lu      = mat_a;
    mat_pp      = eye <T> (sizet_n);

    size_t exchanges = 0;

    for (size_t i = 0; i < sizet_n; i++)
    {

        //pivot section
        T       Umax = 0;
        size_t  row  = 0;
        size_t  q    = 0;
        for (size_t r = i; r < sizet_n; r++)
        {
            T Uii   = mat_lu(r,i);
            q       = 0;
            while (q < i)
            {
                Uii -= mat_lu(r,q) * mat_lu(q,r);
                q++;
            }
            if (std::abs(Uii) > Umax)
            {
                Umax = std::abs(Uii);
                row = r;
            }
        }

        if (Umax < dbl_tol) return false;

        if (i != row)
        { //swap rows
            exchanges++;
            mat_lu.swap_rows(i,row);
            mat_pp.swap_rows(i,row);
        }

        size_t j = i;
        while (j < sizet_n)
        { //determine U across row i
            q = 0;
            while (q < i)
            {
                mat_lu(i,j) -= mat_lu(i,q) * mat_lu(q,j);
                q++;
            }
            j++;
        }
        j = i + 1;
        while (j < sizet_n)
        { //determine L down column i
            q = 0;
            while (q < i)
            {
                mat_lu(j,i) -= mat_lu(j,q) * mat_lu(q,i);
                q++;
            }
            mat_lu(j,i) = mat_lu(j,i) / mat_lu(i,i);
            j++;
        }
    }

    return true;
}

/**
 * @brief solve a set of linear equation
 *
 * solves a set of linear equation presented
 * in matrix form Ax = b where A is the coeficient
 * matrix and x is the solution vector.
 *
 * @param mat_a square coefficient matrix
 * @return solution vector
 * @ingroup LALG
 */
template <class T> matrix<T> linsolve(matrix<T> mat_a, matrix<T> mat_b)
{
    susa::lup <T> solver(mat_a);
    solver.decompose();
    return solver.solve(mat_b);
}

} // namespace susa

#endif // SUSA_LU_H