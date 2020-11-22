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


//#define SUSA_NASSERT

#include <susa.h>

int main(void)
{

  susa::matrix <double> mat_a("[1 2.3 -3.4;8 4.5 1.2;9.1 3 -5]");

  susa::matrix <double> mat_b(3,3);
  mat_b(0,0) = 1;mat_b(0,1) = 2;mat_b(0,2) = 3;
  mat_b(1,0) = 4;mat_b(1,1) = 5;mat_b(1,2) = 6;
  mat_b(2,0) = 7;mat_b(2,1) = 8;mat_b(2,2) = 9;

  std::cout << std::endl << "mat_a            = " << std::endl << mat_a;
  std::cout << std::endl << "mat_b            = " << std::endl << mat_b;
  std::cout << std::endl << "transpose(mat_a) = " << std::endl << transpose(mat_a);
  std::cout << std::endl << "det(mat_a)       = " << susa::det(mat_a) << std::endl;
  std::cout << std::endl << "mat_a .* mat_b   = " << mat_a * mat_b;
  std::cout << std::endl << "mat_a * mat_b    = " << susa::matmul(mat_a,mat_b);

  susa::array <double> arr({21,6,5,15,43});
  arr(2,4,3,0,1) = 55.4;
  std::cout << std::endl << "arr.get({2,4,3,0,1})         = " << arr.get({2,4,3,0,1});
  std::cout << std::endl << "arr.get(2,4,3,0,1)           = " << arr.get(2,4,3,0,1);
  std::cout << std::endl << "arr(2,4,3,0,1)               = " << arr(2,4,3,0,1);
  std::cout << std::endl << "arr.get(21*6*5*15*4)         = " << arr.get(21*6*5*15*4);
  std::cout << std::endl << "arr.get_raw_index(2,4,3,0,1) = " << arr.get_raw_index(2,4,3,0,1);
  std::cout << std::endl;


  susa::matrix <int> mat_c("[1; 2; 3]");

  std::cout << "num rows  = " << mat_c.no_rows() << std::endl;
  std::cout << "num cols  = " << mat_c.no_cols() << std::endl;
  std::cout << "num elems = " << mat_c.size() << std::endl;
  std::cout << std::endl;

  SUSA_LOG_INF("logs can be disabled by defining SUSA_NDEBUG.");

  // this line demonstrates SUSA_ASSERT_MESSAGE macro usage.
  // to disable the message and the abortion of the software uncomment "#define SUSA_NASSERT"
  std::cout << std::endl << "mat_a = " << std::endl << mat_a(mat_a.size() + 1);

  return 0;
}
