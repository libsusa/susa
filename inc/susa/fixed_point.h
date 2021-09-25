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
 * @file fixed_point.h
 * @brief Fixed-Point Number Type
 * @author Behrooz Kamary
 */

#ifndef SUSA_FIXEDPOINT_H
#define SUSA_FIXEDPOINT_H

#include <limits>
#include <cstdint>
#include <type_traits>

namespace susa {
template <typename B, unsigned char I, unsigned char F, bool Enable = ((I+F) > std::numeric_limits<B>::digits)> class fixed_point;

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F>
operator*(const fixed_point<B, I, F>& fp_lhs, const fixed_point<B, I, F>& fp_rhs)
{
  fixed_point <B, I, F> r(fp_lhs);
  return r *= fp_rhs;
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F>
operator/(const fixed_point<B, I, F>& fp_lhs, const fixed_point<B, I, F>& fp_rhs)
{
  fixed_point <B, I, F> r(fp_lhs);
  return r /= fp_rhs;
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F>
operator*(const fixed_point<B, I, F>& fp_lhs, double dbl_rhs)
{
  fixed_point <B, I, F> r(fp_lhs);
  return r *= fixed_point <B, I, F> (dbl_rhs);
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F>
operator*(double dbl_lhs, const fixed_point<B, I, F>& fp_rhs)
{
  fixed_point <B, I, F> r(dbl_lhs);
  return r *= fp_rhs;
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F>
operator*(const fixed_point<B, I, F>& fp_lhs, float flt_rhs)
{
  fixed_point <B, I, F> r(fp_lhs);
  return r *= fixed_point <B, I, F> (flt_rhs);
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F>
operator*(float flt_lhs, const fixed_point<B, I, F>& fp_rhs)
{
  fixed_point <B, I, F> r(flt_lhs);
  return r *= fp_rhs;
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F>
operator/(const fixed_point<B, I, F>& fp_lhs, unsigned uint_rhs)
{
  fixed_point <B, I, F> r(fp_lhs);
  return r /= uint_rhs;
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F>
operator/(unsigned uint_lhs, const fixed_point<B, I, F>& fp_rhs)
{
  fixed_point <B, I, F> r(uint_lhs, 0, false);
  return r /= fp_rhs;
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F>
operator+(const fixed_point<B, I, F>& fp_lhs, unsigned uint_rhs)
{
  fixed_point <B, I, F> r(fp_lhs);
  if (fp_lhs.is_negative())
    return r -= uint_rhs;
  else
    return r += uint_rhs;
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F>
operator+(unsigned uint_lhs, const fixed_point<B, I, F>& fp_rhs)
{
  fixed_point <B, I, F> r(uint_lhs, 0, false);
  if (fp_rhs.is_negative())
    return r -= fp_rhs;
  else
    return r += fp_rhs;
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F>
operator-(const fixed_point<B, I, F>& fp_lhs, unsigned uint_rhs)
{
  fixed_point <B, I, F> r(fp_lhs);
  if (fp_lhs.is_negative())
    return r += uint_rhs;
  else
    return r -= uint_rhs;
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F>
operator-(unsigned uint_lhs, const fixed_point<B, I, F>& fp_rhs)
{
  fixed_point <B, I, F> r(uint_lhs, 0, false);
  if (fp_rhs.is_negative())
    return r += fp_rhs;
  else
    return r -= fp_rhs;
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F>
operator*(const fixed_point<B, I, F>& fp_lhs, int int_rhs)
{
  fixed_point <B, I, F> r(fp_lhs);
  return r *= int_rhs;
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F>
operator*(int int_lhs, const fixed_point<B, I, F>& fp_rhs)
{

  fixed_point <B, I, F> r(std::abs(int_lhs), 0, (int_lhs > 0) ? false : true);
  return r *= fp_rhs;
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F>
operator+(const fixed_point<B, I, F>& fp_lhs, double dbl_rhs)
{
  fixed_point <B, I, F> r(fp_lhs);
  return r += fixed_point <B, I, F> (dbl_rhs);
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F>
operator+(double dbl_lhs, const fixed_point<B, I, F>& fp_rhs)
{
  fixed_point <B, I, F> r(dbl_lhs);
  return r += fp_rhs;
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F>
operator-(const fixed_point<B, I, F>& fp_lhs, double dbl_rhs)
{
  fixed_point <B, I, F> r(fp_lhs);
  return r -= fixed_point <B, I, F> (dbl_rhs);
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F>
operator-(double dbl_lhs, const fixed_point<B, I, F>& fp_rhs)
{
  fixed_point <B, I, F> r(dbl_lhs);
  return r -= fp_rhs;
}

template <typename B, unsigned char I, unsigned char F>
std::ostream &operator<<(std::ostream &out, const fixed_point<B, I, F>& fp_arg)
{
  return out << (double)fp_arg;
}

}

#endif
