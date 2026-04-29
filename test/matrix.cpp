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
 * @file matrix.cpp
 * @brief Unit Test Suit
 * @author Kamary
 */

#include "test.h"

int main(void)
{
    susa::matrix <int> mat_m("[0 1 1; 2 3 2; 1 3 2; 4 2 2]");
    auto mat_left = mat_m.left(3);
    SUSA_TEST_EQ(mat_m, mat_left, "full matrix left");

    auto mat_right = mat_m.right(3);
    SUSA_TEST_EQ(mat_m, mat_right, "full matrix right");

    susa::matrix <int> mat_l("[0 1 ; 2 3 ; 1 3 ; 4 2 ]");
    mat_left = mat_m.left(2);
    SUSA_TEST_EQ(mat_l, mat_left, "matrix left");

    susa::matrix <int> mat_r("[1 1;3 2;3 2;2 2]");
    mat_right = mat_m.right(2);
    SUSA_TEST_EQ(mat_r, mat_right, "matrix right");

    auto mat_mid = mat_m.mid(1,2);
    SUSA_TEST_EQ(mat_mid, mat_right, "matrix mid");

    susa::matrix <int> mat_a = "[0 1 1; 2 3 2; 1 3 2; 4 2 2]";
    SUSA_TEST_EQ(mat_a, mat_m, "literal string matrix parser");

    // ========== expr_add Tests ==========
    // Test expr_add with matrix + matrix
    susa::matrix<int> mat_a1("[1 2; 3 4]");
    susa::matrix<int> mat_a2("[5 6; 7 8]");
    susa::matrix<int> mat_add_result = mat_a1 + mat_a2;
    susa::matrix<int> mat_add_expected("[6 8; 10 12]");
    SUSA_TEST_EQ(mat_add_result, mat_add_expected, "expr_add: matrix + matrix");

    // Test expr_add with matrix + scalar
    susa::matrix<int> mat_s("[1 2; 3 4]");
    susa::matrix<int> mat_add_scalar_result = mat_s + 5;
    susa::matrix<int> mat_add_scalar_expected("[6 7; 8 9]");
    SUSA_TEST_EQ(mat_add_scalar_result, mat_add_scalar_expected, "expr_add: matrix + scalar");

    // Test expr_add with scalar + matrix
    susa::matrix<int> mat_add_scalar_result2 = 5 + mat_s;
    SUSA_TEST_EQ(mat_add_scalar_result2, mat_add_scalar_expected, "expr_add: scalar + matrix");

    // Test expr_add with zero
    susa::matrix<int> mat_zero("[0 0; 0 0]");
    susa::matrix<int> mat_add_zero = mat_s + mat_zero;
    SUSA_TEST_EQ(mat_add_zero, mat_s, "expr_add: matrix + zero matrix");

    // ========== expr_sub Tests ==========
    // Test expr_sub with matrix - matrix
    susa::matrix<int> mat_s1("[10 12; 14 16]");
    susa::matrix<int> mat_s2("[1 2; 3 4]");
    susa::matrix<int> mat_sub_result = mat_s1 - mat_s2;
    susa::matrix<int> mat_sub_expected("[9 10; 11 12]");
    SUSA_TEST_EQ(mat_sub_result, mat_sub_expected, "expr_sub: matrix - matrix");

    // Test expr_sub with matrix - scalar
    susa::matrix<int> mat_sub_scalar_result = mat_s1 - 5;
    susa::matrix<int> mat_sub_scalar_expected("[5 7; 9 11]");
    SUSA_TEST_EQ(mat_sub_scalar_result, mat_sub_scalar_expected, "expr_sub: matrix - scalar");

    // Test expr_sub with scalar - matrix
    susa::matrix<int> mat_scalar_sub_result = 20 - mat_s2;
    susa::matrix<int> mat_scalar_sub_expected("[19 18; 17 16]");
    SUSA_TEST_EQ(mat_scalar_sub_result, mat_scalar_sub_expected, "expr_sub: scalar - matrix");

    // Test expr_sub with zero
    susa::matrix<int> mat_sub_zero = mat_s2 - mat_zero;
    SUSA_TEST_EQ(mat_sub_zero, mat_s2, "expr_sub: matrix - zero matrix");

    // Test expr_sub result equals zero when matrices are same
    susa::matrix<int> mat_sub_self = mat_s2 - mat_s2;
    susa::matrix<int> mat_zero_result("[0 0; 0 0]");
    SUSA_TEST_EQ(mat_sub_self, mat_zero_result, "expr_sub: matrix - itself equals zero");

    // ========== Expression Template Combination Tests ==========
    // Test (a + b) - c
    susa::matrix<int> mat_b1("[1 1; 1 1]");
    susa::matrix<int> mat_b2("[2 2; 2 2]");
    susa::matrix<int> mat_b3("[1 1; 1 1]");
    susa::matrix<int> mat_combined = (mat_b1 + mat_b2) - mat_b3;
    susa::matrix<int> mat_combined_expected("[2 2; 2 2]");
    SUSA_TEST_EQ(mat_combined, mat_combined_expected, "expr: (matrix + matrix) - matrix");

    // Test (a - b) + c
    susa::matrix<int> mat_c1("[5 5; 5 5]");
    susa::matrix<int> mat_c2("[2 2; 2 2]");
    susa::matrix<int> mat_c3("[1 1; 1 1]");
    susa::matrix<int> mat_combined2 = (mat_c1 - mat_c2) + mat_c3;
    susa::matrix<int> mat_combined2_expected("[4 4; 4 4]");
    SUSA_TEST_EQ(mat_combined2, mat_combined2_expected, "expr: (matrix - matrix) + matrix");

    SUSA_TEST_PRINT_STATS();

    return (uint_failed);
}
