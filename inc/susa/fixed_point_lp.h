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
 * @file fixed_point_lp.h
 * @brief Fixed-Point Number Type (Low Precision)
 * @author Behrooz Kamary
 */

#ifndef SUSA_FIXEDPOINT_LP_H
#define SUSA_FIXEDPOINT_LP_H

#include <limits>
#include <cstdint>
#include <type_traits>

#include "fixed_point.h"

namespace susa {
/**
 * @class fixed_point
 * @brief <i>fixed-point</i> type represents a fixed-point number that is defined by its base template B,
 * the number of integer bits I and the number of fraction bits F. I and F can take values up to the number of bits
 * in the base type. F can not be more than half of the number of bits in <b>uintmax_t</b>.
 *
 * this is the light <i>fixed_point</i> implementation where the number of bits in base type B is greater
 * than the number of required bits for integeral and fractional parts (I + F).
 *
 * @ingroup TYPES
 *
 */
template <typename B, unsigned char I, unsigned char F>
class fixed_point <B,I,F, false>
{

public:
  using UB = typename std::make_unsigned<B>::type;
  const unsigned NB   = std::numeric_limits<UB>::digits;
  const unsigned FH   = F / 2;
  const unsigned FR   = F % 2;

private:
  UB      integer;
  double  dbl_orig;

  UB      factor;
  UB      mask_int;
  UB      mask_frac;
  UB      mask_tot;
  bool    overflow;
  bool    negative;

public:
  fixed_point();
  fixed_point(double dbl_arg);
  fixed_point(UB arg_integer, UB arg_fraction, bool negative = false);
  fixed_point(const fixed_point& fp_arg);

  operator double() const;
  operator float() const;

  fixed_point  operator~() const
  {
    fixed_point ret(*this);
    ret.integer  = ~ret.integer;
    ret.negative = !ret.negative;
    return ret;
  }

  fixed_point  operator-() const
  {
    fixed_point ret(*this);
    ret.negative = !negative;
    ret.dbl_orig = -dbl_orig;
    return ret;
  }

  // prefix increment
  fixed_point& operator++()
  {
    integer++;
    return *this;
  }

  // postfix increment
  fixed_point operator++(int)
  {
    fixed_point ret(*this);
    ++(*this);
    return ret;
  }

  // prefix decrement
  fixed_point& operator--()
  {
      integer--;
      return *this;
  }

  // postfix decrement
  fixed_point operator--(int)
  {
      fixed_point old = *this;
      operator--();
      return old;
  }

  fixed_point& operator=(const fixed_point& fp_arg)
  {
    integer     = fp_arg.integer;
    negative    = fp_arg.negative;
    overflow    = fp_arg.overflow;
    dbl_orig    = fp_arg.dbl_orig;
    factor      = fp_arg.factor;
    mask_int    = fp_arg.mask_int;
    mask_frac   = fp_arg.mask_frac;
    mask_tot    = fp_arg.mask_tot;

    return *this;
  }

  fixed_point& operator=(double dbl_arg)
  {
    dbl_orig = dbl_arg;
    if (dbl_arg < 0)
    {
      negative = true;
      dbl_arg *= -1;
    }
    else
    {
      negative = false;
    }

    integer             = static_cast <UB> (dbl_arg);
    double  remainder   = dbl_arg - integer;
    integer           <<= F;
    integer            |= static_cast <UB> (std::round(remainder * factor)) & mask_frac;
    integer            &= mask_tot;

    return *this;
  }

  fixed_point& operator*=(int int_arg)
  {

    if (int_arg == 1) return *this;
    else if (int_arg == -1)
    {
      negative = !negative;
      return *this;
    }

    integer *= std::abs(int_arg);
    overflow = (integer & ~mask_tot) > 0 ? true : false;
    if (negative != (int_arg < 0)) negative = !negative;
    else negative = false;
    return *this;
  }

  fixed_point& operator>>=(unsigned amount)
  {
    integer <<= amount;
    integer &= mask_tot;
    return *this;
  }

  fixed_point& operator<<=(unsigned amount)
  {
    integer >>= amount;
    return *this;
  }

  fixed_point& operator*=(unsigned uint_arg)
  {
    integer *= uint_arg;
    overflow = (integer >> (I + F)) > 0 ? true : false;
    return *this;
  }

  fixed_point& operator+=(unsigned uint_arg)
  {
    UB operand    = uint_arg;
    operand     <<= F;
    integer += operand;
    overflow = (integer >> (I + F)) > 0 ? true : false;
    dbl_orig = (double)*this;
    return *this;
  }

  fixed_point& operator-=(unsigned uint_arg)
  {
    UB operand    = uint_arg;
    operand     <<= F;
    if (integer > operand)
      integer -= operand;
    else
    {
      integer = operand - integer;
      negative = true;
    }
    dbl_orig = (double)*this;
    return *this;
  }

  fixed_point& operator/=(unsigned uint_arg)
  {
    integer /= uint_arg;
    dbl_orig = (double)*this;
    return *this;
  }

  fixed_point& operator/=(const fixed_point& fp_den)
  {
    // keep the precision if the operand is unity
    if (fp_den.integer == (1ul << F))
    {
      negative  = (negative != fp_den.negative);
      dbl_orig  = (double)*this;
      return *this;
    }
    else if (integer == (1ul << F))
    {
      negative  = (negative != fp_den.negative);
      integer   = fp_den.integer;
      dbl_orig  = (double)*this;
      return *this;
    }

    UB rem = integer;
    integer = 0;
    while (rem > fp_den.integer)
    {
      rem -= fp_den.integer;
      integer++;
    }

    for (unsigned f = 0; f < F; f++)
    {
      rem <<= 1;
      integer <<= 1;

      while (rem > fp_den.integer)
      {
        integer++;
        rem -= fp_den.integer;
      }
    }

    if (negative != fp_den.is_negative()) negative = true;
    dbl_orig = (double)*this;
    return *this;
  }

  fixed_point& operator*=(const fixed_point& fp_rhs)
  {
    // keep the precision if the operand is unity
    if (fp_rhs.integer == (1ul << F))
    {
      negative  = (negative != fp_rhs.negative);
      dbl_orig  = (double)*this;
      return *this;
    }
    else if (integer == (1ul << F))
    {
      negative  = (negative != fp_rhs.negative);
      integer   = fp_rhs.integer;
      dbl_orig  = (double)*this;
      return *this;
    }

    integer   = (integer >> FH) * (fp_rhs.integer >> FH);
    integer >>= FR;
    overflow = (integer & ~mask_tot) > 0 ? true : false;
    integer &= mask_tot;
    negative = (negative != fp_rhs.negative);
    dbl_orig = (double)*this;
    return *this;
  }

  fixed_point& operator-=(const fixed_point& fp_rhs)
  {
    return this->operator+=(-fp_rhs);
  }

  fixed_point &operator+=(const fixed_point &fp_rhs)
  {
    if (negative == fp_rhs.negative)
    {
      integer += fp_rhs.integer;
    }
    else if (integer > fp_rhs.integer)
    {
      integer -= fp_rhs.integer;
    }
    else
    {
      integer = fp_rhs.integer - integer;
      negative = fp_rhs.negative;
    }

    overflow = ((integer & ~mask_tot) > 0);
    integer &= mask_tot;
    dbl_orig = (double)*this;
    return *this;
  }

  friend bool operator!=(const fixed_point& fp_lhs, const fixed_point& fp_rhs)
  {
    return ((fp_lhs.integer != fp_rhs.integer) || (fp_lhs.negative != fp_rhs.negative));
  }

  friend bool operator==(const fixed_point& fp_lhs, const fixed_point& fp_rhs)
  {
    return ((fp_lhs.integer == fp_rhs.integer) && (fp_lhs.negative == fp_rhs.negative));
  }

  friend bool operator>(const fixed_point& fp_lhs, const fixed_point& fp_rhs)
  {
    if (fp_lhs.integer > fp_rhs.integer) return !fp_lhs.negative;
    else if (fp_lhs.integer < fp_rhs.integer) return fp_lhs.negative;

    return !fp_lhs.negative;
  }

  friend bool operator<(const fixed_point& fp_lhs, const fixed_point& fp_rhs)
  {
    if (fp_lhs.negative != fp_rhs.negative)
    {
      if (fp_lhs.negative) return true;
      else return false;
    }

    if (fp_lhs.integer < fp_rhs.integer) return !fp_lhs.negative;
    else if (fp_lhs.integer > fp_rhs.integer) return fp_lhs.negative;

    return fp_lhs.negative;
  }

  friend bool operator<(const fixed_point& fp_lhs, int int_rhs)
  {
    unsigned uint_rhs = static_cast<unsigned>(int_rhs);

    if (fp_lhs.negative != (int_rhs < 0))
    {
      if (fp_lhs.negative) return true;
      else return false;
    }

    if (fp_lhs.integer < uint_rhs) return !fp_lhs.negative;
    else if (fp_lhs.integer > uint_rhs) return fp_lhs.negative;

    return false;
  }

  friend bool operator>=(const fixed_point& fp_lhs, const fixed_point& fp_rhs)
  {
    if (fp_lhs.integer > fp_rhs.integer) return !fp_lhs.negative;
    else if (fp_lhs.integer < fp_rhs.integer) return fp_lhs.negative;

    return !fp_lhs.negative;
  }

  friend bool operator<=(const fixed_point& fp_lhs, const fixed_point& fp_rhs)
  {
    if (fp_lhs.integer < fp_rhs.integer) return !fp_lhs.negative;
    else if (fp_lhs.integer > fp_rhs.integer) return fp_lhs.negative;

    return fp_lhs.negative;
  }

  //! returns the sign of fixed-point integer
  bool is_negative() const {return negative;}

  //! returns true if the latest arithmetic operation has overflowed
  bool is_overflow() const {return overflow;}
};

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F, false>::fixed_point()
: integer (0)
, dbl_orig(0)
, factor(1ull << F)
, mask_int (((1ull << I) - 1) << F)
, mask_frac ((1ull << F) - 1)
, mask_tot ((1ull << (I + F)) - 1)
, overflow (false)
, negative (false)
{
  static_assert(std::is_integral<B>::value);

  // required for all operations
  static_assert(std::numeric_limits<UB>::digits > (I + F));

}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F, false>::fixed_point(const fixed_point& fp_arg)
{
  static_assert(std::is_integral<B>::value);

  // required for all operations
  static_assert(std::numeric_limits<UB>::digits > (I + F));

  integer     = fp_arg.integer;
  dbl_orig    = fp_arg.dbl_orig;
  negative    = fp_arg.negative;
  overflow    = fp_arg.overflow;
  factor      = fp_arg.factor;
  mask_frac   = fp_arg.mask_frac;
  mask_int    = fp_arg.mask_int;
  mask_tot    = fp_arg.mask_tot;
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F, false>::fixed_point(UB arg_integer, UB arg_fraction, bool negative)
: integer(0)
, factor(1ull << F)
, mask_int (((1ull << I) - 1) << F)
, mask_frac ((1ull << F) - 1)
, mask_tot ((1ull << (I + F)) - 1)
, overflow (false)
, negative (negative)
{
  static_assert(std::is_integral<B>::value);

  // required for all operations
  static_assert(std::numeric_limits<UB>::digits > (I + F));

  integer |= ((arg_integer << F) & mask_int);
  integer |= (arg_fraction & mask_frac);
  integer &= mask_tot;
  dbl_orig = (double)*this;
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F, false>::fixed_point(double dbl_arg)
: dbl_orig(dbl_arg)
, factor(1ull << F)
, mask_int (((1ull << I) - 1) << F)
, mask_frac ((1ull << F) - 1)
, mask_tot ((1ull << (I + F)) - 1)
, overflow (false)
{
  static_assert(std::is_integral<B>::value);

  // required for all operations
  static_assert(std::numeric_limits<UB>::digits > (I + F));

  if (dbl_arg < 0)
  {
    negative = true;
    dbl_arg *= -1;
  }
  else
  {
    negative = false;
  }

  integer             = static_cast <UB> (dbl_arg);
  double  remainder   = dbl_arg - integer;
  integer           <<= F;
  integer            |= static_cast <UB> (std::round(remainder * factor)) & mask_frac;
  integer            &= mask_tot;
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F, false>::operator double() const
{
  return (negative ? -1.0 : 1.0) * (((integer) >> F) + (double)(integer & mask_frac) / factor);
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F, false>::operator float() const
{
  return (negative ? -1.0 : 1.0) * ((float)(integer >> F) + (float)(integer & mask_frac) / factor);
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F, false>
operator+(const fixed_point<B, I, F, false>& fp_lhs, const fixed_point<B, I, F, false>& fp_rhs)
{
  fixed_point <B, I, F> r(fp_lhs);
  return r += fp_rhs;
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F, false>
operator-(const fixed_point<B, I, F, false>& fp_lhs, const fixed_point<B, I, F, false>& fp_rhs)
{
  fixed_point <B, I, F> r(fp_lhs);
  return r -= fp_rhs;
}

}
#endif
