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
#include <cmath>
#include <type_traits>
#include <susa/fixed_point.h>

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
 * @brief square root of an unsigned integer
 *
 * @ingroup ARITHMETICS
 *
 */
template <typename T>
constexpr T isqrt(T T_arg)
{

  static_assert(std::is_unsigned<T>::value);

	T ret = T_arg >> 1;
	if (!ret) return T_arg;

  T interim = ( ret + T_arg / ret ) >> 1;
  while (interim < ret)
  {
    ret = interim;
    interim = (ret + T_arg / ret) >> 1;
  }

  return ret;
}

/**
 * @brief Euler's number exponent function
 *
 * @ingroup ARITHMETICS
 *
 */
template <typename B, unsigned char I, unsigned char F>
fixed_point<B,I,F> texp(fixed_point<B,I,F> fp_arg, unsigned int uint_tail = 10)
{

  fixed_point<B,I,F> ret(1.0f);

  for (unsigned i = uint_tail - 1; i > 0; i--)
  {
    ret = fp_arg * ret / i + 1u;
  }

  return ret;
}

/**
 * @brief fast floor function
 *
 * @ingroup ARITHMETICS
 *
 */
template<typename T>
constexpr T ffloor(T T_arg)
{
  static_assert(std::is_floating_point<T>::value);

  if (T_arg < 0) return ((int)(T_arg - 1));
  return (int)(T_arg);
}

/**
 * @brief fast cosine function
 *
 * @ingroup ARITHMETICS
 *
 */
template<typename T>
T fcos(T T_arg)
{
  static_assert(std::is_floating_point<T>::value);

  constexpr T tp = 1./(2.*M_PI);

  T_arg *= tp;
  T_arg -= T(.25) + ffloor(T_arg + T(.25));
  T_arg *= T(16.) * (std::abs(T_arg) - T(.5));

  // extra precision
  T_arg += T(.225) * T_arg * (std::abs(T_arg) - T(1.));

  return T_arg;
}

/**
 * @brief fast sine function
 *
 * @ingroup ARITHMETICS
 *
 */
template<typename T>
T fsin(T T_arg)
{
  return fcos(T_arg - M_PI_2);
}

/**
 * @brief CORDIC algorithm for calculating trigonometric functions
 *
 * @ingroup ARITHMETICS
 */
template <typename T, unsigned R = 128>
class cordic
{

private:
  std::vector<T>  vec_theta;
  std::vector<T>  vec_atan;
  std::vector<T>  vec_kval;
  double          dbl_k;

public:
  cordic()
  : vec_theta(R)
  , vec_atan(R)
  , vec_kval(R)
  , dbl_k(1.0)
  {
    // n.b. R <= F when T is susa::fixed_point<>
    for (unsigned uint_index = 0; uint_index < R; uint_index++)
    {
      dbl_k *= std::sqrt(1.0f / (1 + std::pow(2.0, -2.0 * uint_index)));
      vec_theta[uint_index] = (double)std::pow(2.0, -(double)uint_index);
      vec_atan[uint_index]  = (double)std::atan((double)vec_theta[uint_index]);
      vec_kval[uint_index]  = (double)dbl_k;
    }
  }

  /**
   * @brief sine
   */
  T sin(const T& T_arg)
  {
    auto vec_result = calculate(T_arg);
    return vec_result[0];
  }

  /**
   * @brief cosine
   */
  T cos(const T& T_arg)
  {
    auto vec_result = calculate(T_arg);
    return vec_result[1];
  }

  /**
   * @brief the natural exponential funtion
   *
   * @param T_arg imaginary number in the form jT_arg
   */
  std::complex<T> expj(const T& T_arg)
  {
    auto vec_result = calculate(T_arg);
    return std::complex<T>(vec_result[1], vec_result[0]);
  }

private:
  std::vector<T> calculate(const T& T_arg)
  {
    std::vector<T> vec_result(3);

    T T_beta(T_arg);
    T T_sin(0.0);
    T T_cos(1.0);
    T T_sin_prev, T_cos_prev;
    int sgn = 1;
    for (unsigned i = 0; i < R; i++)
    {
      if ( T_beta < 0) sgn = -1;
      else sgn = 1;

      T_sin_prev  = T_sin;
      T_cos_prev  = T_cos;
      T_sin       = sgn * T_cos_prev * vec_theta[i] + T_sin_prev;
      T_cos       = T_cos_prev - sgn * T_sin_prev * vec_theta[i];
      T_beta      = T_beta - sgn * vec_atan[i];
    }
    vec_result[0] = dbl_k * T_sin;
    vec_result[1] = dbl_k * T_cos;
    vec_result[2] = 0;

    return vec_result;
  }
};

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
