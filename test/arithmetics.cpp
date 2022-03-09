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
 * @file arithmetics.cpp
 * @brief Unit Test Suit
 * @author Behrooz Kamary
 */

#include "test.h"
#include <susa/fixed_point.h>
#include <susa/arithmetics.h>

int main(void)
{

    {
    // High Precision Fixed-Point
    double dbla = -457.103;
    susa::fixed_point <int64_t, 40,41> fpnuma(dbla);
    SUSA_TEST_EQ_DOUBLE4(dbla, (double)fpnuma, "coversion to/from a floating point");

    double dblb = 74.166;
    susa::fixed_point <int64_t, 40,41> fpnumb(dblb);
    SUSA_TEST_EQ_DOUBLE4(dblb, (double)fpnumb, "coversion to/from a floating point");


    double dbl_add = dbla + dblb;
    auto fp_add = fpnuma + fpnumb;
    SUSA_TEST_EQ_DOUBLE4(dbl_add, (double)fp_add, "addition");


    double dbl_sub = dbla - dblb;
    auto fp_sub = fpnuma - fpnumb;
    SUSA_TEST_EQ_DOUBLE4(dbl_sub, (double)fp_sub, "subtraction");


    double dbl_mul = dbla * dblb;
    auto fp_mul = fpnuma * fpnumb;
    SUSA_TEST_EQ_DOUBLE4(dbl_mul, (double)fp_mul, "multiplication");


    double dbl_div = dbla / dblb;
    auto fp_div = fpnuma / fpnumb;
    SUSA_TEST_EQ_DOUBLE4(dbl_div, (double)fp_div, "division");

    double dblc = 74.166;
    susa::fixed_point <int64_t, 40,41> fpnumc(dblc);

    SUSA_TEST_EQ((fpnuma > fpnumb), false, "less than operator for two floating points");
    SUSA_TEST_EQ((fpnuma >= fpnumb), false, "less than or equal operator for two floating points");
    SUSA_TEST_EQ((fpnuma < fpnumb), true, "less than operator for two floating points");
    SUSA_TEST_EQ((fpnuma <= fpnumb), true, "less than or equal operator for two floating points");
    SUSA_TEST_EQ((fpnumc == fpnumb), true, "equal operator for two floating points");
    SUSA_TEST_EQ((fpnumc != fpnumb), false, "not equal operator for two floating points");
    SUSA_TEST_EQ((fpnuma != fpnumb), true, "not equal operator for two floating points");

    }

    {
    // Low Precision Fixed-Point
    double dbla = -457.103;
    susa::fixed_point <int64_t, 20, 30> fpnuma(dbla);
    SUSA_TEST_EQ_DOUBLE4(dbla, (double)fpnuma, "coversion to/from a floating point");

    double dblb = 74.166;
    susa::fixed_point <int64_t, 20, 30> fpnumb(dblb);
    SUSA_TEST_EQ_DOUBLE4(dblb, (double)fpnumb, "coversion to/from a floating point");

    double dbl_add = dbla + dblb;
    auto fp_add = fpnuma + fpnumb;
    SUSA_TEST_EQ_DOUBLE4(dbl_add, (double)fp_add, "addition");


    double dbl_sub = dbla - dblb;
    auto fp_sub = fpnuma - fpnumb;
    SUSA_TEST_EQ_DOUBLE4(dbl_sub, (double)fp_sub, "subtraction");


    double dbl_mul = dbla * dblb;
    auto fp_mul = fpnuma * fpnumb;
    SUSA_TEST_EQ_DOUBLE2(dbl_mul, (double)fp_mul, "multiplication");


    double dbl_div = dbla / dblb;
    auto fp_div = fpnuma / fpnumb;
    SUSA_TEST_EQ_DOUBLE4(dbl_div, (double)fp_div, "division");

    double dblc = 74.166;
    susa::fixed_point <int64_t, 20, 30> fpnumc(dblc);

    SUSA_TEST_EQ((fpnuma > fpnumb), false, "less than operator for two floating points");
    SUSA_TEST_EQ((fpnuma >= fpnumb), false, "less than or equal operator for two floating points");
    SUSA_TEST_EQ((fpnuma < fpnumb), true, "less than operator for two floating points");
    SUSA_TEST_EQ((fpnuma <= fpnumb), true, "less than or equal operator for two floating points");
    SUSA_TEST_EQ((fpnumc == fpnumb), true, "equal operator for two floating points");
    SUSA_TEST_EQ((fpnumc != fpnumb), false, "not equal operator for two floating points");
    SUSA_TEST_EQ((fpnuma != fpnumb), true, "not equal operator for two floating points");
    }

    {
    // Low Precision Fixed-Point
    double dbla = -457.103;
    susa::fixed_point <int64_t, 20, 30> fpnuma(dbla);
    SUSA_TEST_EQ_DOUBLE4(dbla, (double)fpnuma, "coversion to/from a floating point");

    double dblb = -74.166;
    susa::fixed_point <int64_t, 20, 30> fpnumb(dblb);
    SUSA_TEST_EQ_DOUBLE4(dblb, (double)fpnumb, "coversion to/from a floating point");

    double dblc = 74.166;
    susa::fixed_point <int64_t, 20, 30> fpnumc(dblc);

    SUSA_TEST_EQ((fpnuma > fpnumb), false, "less than operator for two floating points");
    SUSA_TEST_EQ((fpnuma >= fpnumb), false, "less than or equal operator for two floating points");
    SUSA_TEST_EQ((fpnuma < fpnumb), true, "less than operator for two floating points");
    SUSA_TEST_EQ((fpnuma <= fpnumb), true, "less than or equal operator for two floating points");
    SUSA_TEST_EQ((fpnumc == fpnumb), false, "equal operator for two floating points");
    SUSA_TEST_EQ((fpnumc != fpnumb), true, "not equal operator for two floating points");
    SUSA_TEST_EQ((fpnumc == -fpnumb), true, "negation operator for two floating points");
    SUSA_TEST_EQ((fpnumc == fpnumb * -1), true, "multiply by minus one on the right-hand side");
    SUSA_TEST_EQ((fpnumc == -1 * fpnumb), true, "multiply by minus one on the left-hand side");
    }

    {
    unsigned uint_a = 432;
    unsigned uint_b = 832;
    auto     uint_r = susa::intmul(uint_a, uint_b);
    SUSA_TEST_EQ(std::get<0>(uint_r) + (std::get<1>(uint_r) << (std::numeric_limits<unsigned>::digits / 2)), 359424, "unsigned integer multiplication");
    }

    {
    unsigned uint_a = 65000;
    unsigned uint_b = 64000;
    auto     uint_r = susa::intmul(uint_a, uint_b);
    SUSA_TEST_EQ(std::get<0>(uint_r) + (std::get<1>(uint_r) << (std::numeric_limits<unsigned>::digits / 2)), 416e7, "unsigned integer multiplication");
    }

    {
      double dblnum = 0.5;
      susa::fixed_point<int64_t,30, 33> fpnum(dblnum);
      SUSA_TEST_EQ_DOUBLE4(dblnum, (double)fpnum, "ctor cast to double");
      SUSA_TEST_EQ_DOUBLE4(dblnum, (float)fpnum, "ctor cast to float");

      dblnum = -0.6;
      fpnum = dblnum;
      SUSA_TEST_EQ_DOUBLE4(dblnum, (double)fpnum, "cast to double");
      SUSA_TEST_EQ_DOUBLE4(dblnum, (float)fpnum, "cast to float");
    }

    {
      SUSA_TEST_EQ(susa::ffloor(4.5), 4, "susa fast floor implementation");
      SUSA_TEST_EQ(susa::ffloor(-4.5), -5, "susa fast floor implementation");
    }

    {
      unsigned uint_qube = 8563472;
      SUSA_TEST_EQ(susa::uisqrt(uint_qube), 2926, "susa square root for unsigned integers");
    }

    {
      susa::fixed_point<int32_t, 10, 21> fpnum(1.0);
      SUSA_TEST_EQ_DOUBLE3((double)susa::texp(fpnum, 20), 2.71828, "susa Taylor series e^x for fixed_point type");
    }

    {
      susa::cordic<double, 15> c;

      SUSA_TEST_EQ_DOUBLE4(c.sin(M_PI/ 6.0f), 0.5, "CORDIC sine calculation");
      SUSA_TEST_EQ_DOUBLE4(c.cos(M_PI/ 6.0f), 0.866025, "CORDIC cosine calculation");
    }

    {
      // Low Precision Fixed-Point
      susa::cordic<susa::fixed_point<int64_t, 30, 33>, 15> c;
      SUSA_TEST_EQ_DOUBLE4((double)c.sin(M_PI/ 6.0f), 0.5, "CORDIC sine calculation (fp)");
      SUSA_TEST_EQ_DOUBLE4((double)c.cos(M_PI/ 6.0f), 0.866025, "CORDIC cosine calculation (fp)");
    }

    SUSA_TEST_PRINT_STATS();

    return (uint_failed);
}
