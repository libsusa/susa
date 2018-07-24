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
 * @author Behrooz Aliabadi
 * @version 1.0.0
 */

#include "test.h"
#include <iomanip>

int main(int argc, char const *argv[])
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

    susa::matrix <int> mat_left ("5 6;3 2; 7 4;-4 8");
    susa::matrix <int> mat_right ("5 6 9 -3;-1 -2 0 1");
    susa::matrix <int> experiment = matmul(mat_left,mat_right);
    susa::matrix <int> expected   = susa::matrix <int> ("19 18 45 -9;13 14 27 -7;31 34 63 -17;-28 -40 -36 20");
    SUSA_TEST_EQ(experiment, expected, "matmul() matrix multiplication.");

    SUSA_TEST_PRINT_STATS();

    return (uint_failed);
}
