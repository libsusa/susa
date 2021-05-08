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
 * @file arithmetics.h
 * @brief Integeral arithmetics operations (declaration and definition).
 * @author Behrooz Kamary
 *
 * @defgroup ARITHMETICS Arithmetics
 */

#ifndef SUSA_ARITHMETICS_H
#define SUSA_ARITHMETICS_H

#include <limits>
#include <cstdint>
#include <type_traits>

namespace susa {

/**
 * @brief detect the endianness of the machine at compile time.
 *
 * this enumeration class has been borrowed from C++20 standard.
 */
enum class endian
{
#ifdef _WIN32
    little = 0,
    big    = 1,
    native = little
#else
    little = __ORDER_LITTLE_ENDIAN__,
    big    = __ORDER_BIG_ENDIAN__,
    native = __BYTE_ORDER__
#endif
};

/**
 * @brief Greatest Common Divisor (GCD)
 *
 * @param arg_a
 * @param arg_b
 *
 * @ingroup ARITHMETICS
 *
 */
template <typename T>
constexpr T gcd(T arg_a, T arg_b)
{
    while (arg_a != arg_b)
    {
        if (arg_a > arg_b)
            arg_a = arg_a - arg_b;
        else
            arg_b = arg_b - arg_a;
    }

    return arg_a;
}


/**
 * @brief mapping an integeral matrix from one domain i.e. range to another.
 *
 * assume having an tegeral vector/matrix of intrinsic type X lying between
 * <i>arg_x_min</i> and <i>arg_x_max</i> and the target is to scale/map them
 * to another integeral vector/matrix of intrinsic type Y lying between
 * <i>arg_y_min</i> and <i>arg_y_max</i>. <i>scale</i> function does this mapping
 * and returns the result in a vector/matrix same as the input vector/matrix.
 *
 * @param arg_x_min
 * @param arg_x_max
 * @param arg_y_min
 * @param arg_y_max
 * @param mat_arg
 *
 * @ingroup ARITHMETICS
 *
 */
template <typename X, typename Y,  template <typename> typename Allocator>
matrix <Y, Allocator<Y>>
scale(X arg_x_min, X arg_x_max, Y arg_y_min, Y arg_y_max, matrix <X, Allocator<X>> mat_arg)
{
  static_assert(std::is_integral<X>::value);
  static_assert(std::is_integral<Y>::value);

  matrix <Y, Allocator<Y>> mat_ret(mat_arg.shape());

  size_t sz_length    = mat_arg.size();
  X x_range           = arg_x_max - arg_x_min;
  Y y_range           = arg_y_max - arg_y_min;
  X r                 = x_range / 2; // integer rounding correction
  r                  *= y_range < 0 ? -1 : 1;

  for (size_t indx = 0; indx < sz_length; indx++)
  {
    mat_ret(indx) = ((mat_arg(indx) - arg_x_min) * y_range + r) / x_range;
  }

  return mat_ret;
}

template <typename T>
std::tuple<T,T> intmul(const T& T_left, const T& T_right)
{
  static_assert(std::is_unsigned<T>::value);
  static_assert(susa::endian::native == susa::endian::little);
  const unsigned NB = std::numeric_limits<T>::digits / 2;
  const T mask_low  = (1 << NB) - 1;
  const T mask_high = ~mask_low;

  T left_low    = T_left  & mask_low;
  T left_high   = (T_left  & mask_high) >> NB;
  T right_low   = T_right & mask_low;
  T right_high  = (T_right & mask_high) >> NB;

  T low         = left_low * right_low + ((left_low * right_high) << NB);
  T high        = left_high * right_low + ((left_high * right_high) << NB);

  return std::make_tuple(low,high);
}

template <typename T>
constexpr bool mult_overflow(const T& arg_a, const T& arg_b)
{
    return ((arg_b >= 0) && (arg_a >= 0) && (arg_a > std::numeric_limits<T>::max() / arg_b))
        || ((arg_b < 0) && (arg_a < 0) && (arg_a < std::numeric_limits<T>::max() / arg_b));
}

template <typename T>
constexpr bool mult_underflow(const T& arg_a, const T& arg_b)
{
    return ((arg_b >= 0) && (arg_a < 0) && (arg_a < std::numeric_limits<T>::min() / arg_b))
        || ((arg_b < 0) && (arg_a >= 0) && (arg_a > std::numeric_limits<T>::min() / arg_b));
}

template <typename T>
constexpr bool div_overflow(const T& arg_a, const T& arg_b)
{
    return (arg_a == std::numeric_limits<T>::min()) && (arg_b == -1) && (arg_a != 0);
}

template <typename T>
constexpr bool add_overflow(const T& arg_a, const T& arg_b)
{
    return (arg_b >= 0) && (arg_a > std::numeric_limits<T>::max() - arg_b);
}

template <typename T>
constexpr bool add_underflow(const T& arg_a, const T& arg_b)
{
    return (arg_b < 0) && (arg_a < std::numeric_limits<T>::min() - arg_b);
}

template <typename T>
constexpr bool sub_overflow(const T& arg_a, const T& arg_b)
{
    return (arg_b < 0) && (arg_a > std::numeric_limits<T>::max() + arg_b);
}

template <typename T>
constexpr bool sub_underflow(const T& arg_a, const T& arg_b)
{
    return (arg_b >= 0) && (arg_a < std::numeric_limits<T>::min() + arg_b);
}

}      // NAMESPACE SUSA
#endif // SUSA_ARITHMETICS_H
