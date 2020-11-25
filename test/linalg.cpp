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
 * along with Susa. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file linalg.cpp
 * @brief Unit Test Suit
 * @author Behrooz Kamary
 */

#include "test.h"
#include <iomanip>

int main(void)
{
    {
    susa::matrix <float> mat_a("2 0 2; 0 1 0; 0 0 0");
    susa::matrix <float> mat_u;
    susa::matrix <float> mat_s;
    susa::matrix <float> mat_v;

    susa::svd(mat_a, mat_u, mat_s, mat_v);
    susa::matrix<int> expected("28284;0;10000");
    auto experiment = (susa::matrix<int>)(mat_s * 10000.0f);
    SUSA_TEST_EQ (experiment, expected, "Singular Value Decomposition (SVD) for a sample matrix.");

    mat_a = "1 2; 3 4; 5 6; 7 8";
    susa::svd(mat_a, mat_u, mat_s, mat_v);
    expected   = susa::matrix<int> ("6268;142690");
    experiment = (susa::matrix<int>)(mat_s * 10000.0f);
    SUSA_TEST_EQ (experiment, expected, "Singular Value Decomposition (SVD) for a sample matrix.");

    mat_a = "16,  2,  3, 13;5, 11, 10,  8;9,  7,  6, 12;4, 14, 15,  1";
    susa::svd(mat_a, mat_u, mat_s, mat_v);

    expected   = susa::matrix<int> ("-5000 6708 -4999 2236;\
                                     -5000 -2236 5000 6708;\
                                     -4999 2236 4999 -6708;\
                                     -5000 -6708 -4999 -2236");
    experiment = (susa::matrix<int>)(mat_u * 10000.0f);
    SUSA_TEST_EQ (experiment, expected, "Singular Value Decomposition (SVD) for a sample matrix.");

    expected = (susa::matrix<int>) (susa::matrix<float> ("[34;17.88854408;4.472135067;0]") * 10000.0f);
    experiment = (susa::matrix<int>)(mat_s * 10000.0f);
    SUSA_TEST_EQ (experiment, expected, "Singular Value Decomposition (SVD) for a sample matrix.");
    }

    {
    susa::matrix <int> mat_left ("5 6;3 2; 7 4;-4 8");
    susa::matrix <int> mat_right ("5 6 9 -3;-1 -2 0 1");
    susa::matrix <int> experiment = matmul(mat_left,mat_right);
    susa::matrix <int> expected   = susa::matrix <int> ("19 18 45 -9;13 14 27 -7;31 34 63 -17;-28 -40 -36 20");
    SUSA_TEST_EQ(experiment, expected, "matmul() matrix multiplication.");
    }

    {
    susa::matrix <int>      mat_a("2, -1, -2;-4, 6, 3;-4, -2, 8");
    susa::lu <int>          solver(mat_a);
    solver.decompose();

    susa::matrix <int> mat_res = susa::matmul(solver.get_lower(), solver.get_upper());
    SUSA_TEST_EQ(mat_res, mat_a, "matrix LU decomposition.");

    }

    {
    susa::matrix <float>      mat_a("2, 1, 1;-1, 1, -1;1, 2, 3");
    susa::matrix <float>      mat_b("2;3;-10");
    susa::lup <float>         solver(mat_a);
    SUSA_TEST_EQ(solver.decompose(), true, "LUP decomposition ret value");

    susa::matrix <float>      mat_res = solver.solve(mat_b);
    SUSA_TEST_EQ((susa::matmul(mat_a, mat_res) == mat_b), true, "LUP linear equation solution");
    }

    {
    susa::matrix <float>      mat_a("2,7,6,2;9,5,1,3;4,3,8,4;5,6,7,8");
    susa::matrix <float>      mat_p("0 1 0 0;1 0 0 0;0 0 1 0;0 0 0 1");
    susa::lup <float>         solver(mat_a);
    solver.decompose_alt();

    SUSA_TEST_EQ(solver.get_pivot(), mat_p, "compute pivot of LUP decomposition");

    susa::matrix <float>      mat_pa("9 5 1 3;2 7 6 2;4 3 8 4;5 6 7 8");
    susa::matrix <float>      mat_res = susa::matmul(solver.get_pivot(), mat_a);

    SUSA_TEST_EQ(mat_pa, mat_res, "compute pivot of LUP decomposition");

    susa::matrix <float>      mat_lu("9 5 1 3;"
                                     "0.2222222222222222 5.888888888888889 5.777777777777778 1.3333333333333335;"
                                     "0.4444444444444444 0.13207547169811318 6.7924528301886795 2.4905660377358494;"
                                     "0.5555555555555556 0.5471698113207547 0.48333333333333334 4.3999999999999995");

    susa::matrix <float>      mat_err_sum = susa::sum(susa::sum(mat_lu - solver.get_lu()));
    SUSA_TEST_EQ_DOUBLE(mat_err_sum(0), 0, "compute LUP decomposition");
    }

    {
    susa::matrix <float>      mat_a("1 0 2;-1 5 0; 0 3 -9");
    susa::matrix <float>      mat_res = susa::inv(mat_a);
    susa::matrix <float>      mat_exp("0.8824 -0.1176 0.1961;0.1765 0.1765 0.0392;0.0588 0.0588 -0.0980");

    susa::matrix <float>      mat_err_sum = susa::sum(susa::sum(mat_exp - mat_res));
    SUSA_TEST_EQ_DOUBLE(mat_err_sum(0), 0, "inverse of a square matrix");
    std::cout << susa::sum(susa::sum(mat_exp - mat_res)) << std::endl;
    }


    SUSA_TEST_PRINT_STATS();

    return (uint_failed);
}
