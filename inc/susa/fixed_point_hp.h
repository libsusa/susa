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
 * @file fixed_point_hp.h
 * @brief Fixed-Point Number Type (High Precision)
 * @author Behrooz Kamary
 */

#ifndef SUSA_FIXEDPOINT_HP_H
#define SUSA_FIXEDPOINT_HP_H

#include <limits>
#include <cstdint>
#include <type_traits>

#include "fixed_point.h"

namespace susa {

/**
 * @class fixed_point
 * @brief <i>fixed-point</i> type represents a fixed-point number that is defined by its base template B,
 * the number of integer bits I and the number of fraction bits F. I and F can take values up to the number of bits
 * in the base type.
 *
 * @ingroup TYPES
 *
 */
template <typename B, unsigned char I, unsigned char F>
class fixed_point <B,I,F, true>
{

public:

  using UB = typename std::make_unsigned<B>::type;
  const unsigned NB = std::numeric_limits<uintmax_t>::digits / 2;
  const unsigned FH   = F / 2;
  const unsigned FR   = F % 2;

private:

  UB      integer;
  UB      fraction;
  double  dbl_orig;

  UB      factor;
  UB      mask_int;
  UB      mask_frac;
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
    ret.fraction = ~ret.fraction;
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

  fixed_point& operator|=(const fixed_point& fp_rhs);
  fixed_point& operator&=(const fixed_point& fp_rhs);
  fixed_point& operator^=(const fixed_point& fp_rhs);
  fixed_point& operator%=(const fixed_point& fp_rhs);

  fixed_point& operator>>=(unsigned amount)
  {
    fraction >>= amount;

    if (F > amount)
    {
      fraction |= (integer << (F - amount));
    }
    else
    {
      fraction  = (integer >> (amount - F));
    }

    integer >>= amount;

    integer  &= mask_int;
    fraction &= mask_frac;
    return *this;
  }

  fixed_point& operator<<=(unsigned amount)
  {
    integer <<= amount;

    if (F > amount)
    {
      integer |= (fraction >> (F - amount));
    }
    else
    {
      integer = (fraction << (amount - F));
    }

    fraction <<= amount;

    integer  &= mask_int;
    fraction &= mask_frac;
    return *this;
  }

  // prefix increment
  fixed_point& operator++()
  {
    integer++;
    dbl_orig     = (double)*this;
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
      dbl_orig     = (double)*this;
      return *this;
  }

  // postfix decrement
  fixed_point operator--(int)
  {
      fixed_point old = *this;
      operator--();
      return old;
  }

  fixed_point& operator=(const fixed_point& fp_rhs)
  {
    integer     = fp_rhs.integer;
    fraction    = fp_rhs.fraction;
    negative    = fp_rhs.negative;
    dbl_orig    = fp_rhs.dbl_orig;
    factor      = fp_rhs.factor;
    mask_int    = fp_rhs.mask_int;
    mask_frac   = fp_rhs.mask_frac;
    overflow    = false;

    return *this;
  }

  fixed_point& operator+=(const fixed_point& fp_rhs)
  {
    fixed_point r(fp_rhs);
    fraction    += r.fraction;
    UB carry     = (fraction & ~mask_frac) >> F;
    fraction    &= mask_frac;
    integer     += r.integer;
    integer     += carry;
    overflow     = (integer >> I) > 0 ? true : false;
    integer     &= mask_int;
    dbl_orig     = (double)*this;
    return *this;
  }

  fixed_point& operator-=(const fixed_point& fp_rhs)
  {
    fixed_point r(fp_rhs);
    if (fraction < r.fraction)
    {
        fraction += (1ull << F);
        integer--;
    }
    fraction    -= r.fraction;
    UB carry     = (fraction & ~mask_frac) >> F;
    fraction    &= mask_frac;
    integer     -= r.integer;
    integer     -= carry;
    overflow     = (integer >> I) > 0 ? true : false;
    integer     &= mask_int;
    dbl_orig     = (double)*this;
    return *this;
  }

  fixed_point& operator*=(const fixed_point& fp_rhs)
  {
    // keep the precision if the operand is unity
    if (fp_rhs.integer == 1ul and fp_rhs.fraction == 0)
    {
      negative = (negative != fp_rhs.negative);
      dbl_orig = (double)*this;
      return *this;
    }
    else if (integer == 1ul and fraction == 0)
    {
      negative = (negative != fp_rhs.negative);
      integer  = fp_rhs.integer;
      fraction = fp_rhs.fraction;
      dbl_orig = (double)*this;
      return *this;
    }

    fixed_point r(fp_rhs);
    uintmax_t intint        = (uintmax_t)integer * r.integer;
    uintmax_t fracint       = (uintmax_t)fraction * r.integer;
    uintmax_t intfrac       = (uintmax_t)integer * r.fraction;
    uintmax_t fracfrac      = (uintmax_t)(fraction >> FH) * (r.fraction >> FH);
    fracfrac              >>= FR;

    integer  = intint + (fracint >> F) + (intfrac >> F) + ((fracfrac >> F));
    fraction = (fracint & mask_frac) + (intfrac & mask_frac) +  ((fracfrac >> 0) & mask_frac);
    overflow = (integer >> I) > 0 ? true : false;
    integer  &= mask_int;

    if (negative != r.is_negative()) negative = true;

    return *this;
  }

  fixed_point& operator/=(const fixed_point& fp_den)
  {
    // keep the precision if the operand is unity
    if (fp_den.integer == 1ul and fp_den.fraction == 0)
    {
      negative = (negative != fp_den.negative);
      dbl_orig = (double)*this;
      return *this;
    }
    else if (integer == 1ul and fraction == 0)
    {
      negative = (negative != fp_den.negative);
      integer  = fp_den.integer;
      fraction = fp_den.fraction;
      dbl_orig = (double)*this;
      return *this;
    }

    fixed_point rem(*this);
    rem.negative = false;
    integer      = 0;
    while (rem >= fp_den)
    {
      integer++;
      rem -= fp_den;
    }

    fraction = 0;

    for (unsigned f = 0; f < F; f++)
    {
      rem <<= 1;
      fraction <<= 1;
      while (rem >= fp_den)
      {
        fraction++;
        rem -= fp_den;
      }
    }

    if (negative != fp_den.is_negative()) negative = true;

    return *this;
  }

  friend bool operator!=(const fixed_point& fp_lhs, const fixed_point& fp_rhs)
  {
    return fp_lhs.fraction != fp_rhs.fraction || fp_lhs.integer != fp_rhs.integer;
  }

  friend bool operator==(const fixed_point& fp_lhs, const fixed_point& fp_rhs)
  {
    return fp_lhs.fraction == fp_rhs.fraction && fp_lhs.integer == fp_rhs.integer;
  }

  friend bool operator>(const fixed_point& fp_lhs, const fixed_point& fp_rhs)
  {
    if (fp_lhs.integer > fp_rhs.integer) return !fp_lhs.negative;
    else if (fp_lhs.integer < fp_rhs.integer) return fp_lhs.negative;
    else if (fp_lhs.fraction > fp_rhs.fraction) return !fp_lhs.negative;

    return fp_lhs.negative;
  }

  friend bool operator<(const fixed_point& fp_lhs, const fixed_point& fp_rhs)
  {
    if (fp_lhs.integer < fp_rhs.integer) return !fp_lhs.negative;
    else if (fp_lhs.integer > fp_rhs.integer) return fp_lhs.negative;
    else if (fp_lhs.fraction < fp_rhs.fraction) return !fp_lhs.negative;

    return fp_lhs.negative;
  }

  friend bool operator>=(const fixed_point& fp_lhs, const fixed_point& fp_rhs)
  {
    if (fp_lhs.integer > fp_rhs.integer) return !fp_lhs.negative;
    else if (fp_lhs.integer < fp_rhs.integer) return fp_lhs.negative;
    else if (fp_lhs.fraction >= fp_rhs.fraction) return !fp_lhs.negative;

    return fp_lhs.negative;
  }

  friend bool operator<=(const fixed_point& fp_lhs, const fixed_point& fp_rhs)
  {
    if (fp_lhs.integer < fp_rhs.integer) return !fp_lhs.negative;
    else if (fp_lhs.integer > fp_rhs.integer) return fp_lhs.negative;
    else if (fp_lhs.fraction <= fp_rhs.fraction) return !fp_lhs.negative;

    return fp_lhs.negative;
  }

  //! returns the sign of fixed-point integer
  bool is_negative() const {return negative;}

  //! returns true if the latest arithmetic operation has overflowed
  bool is_overflow() const {return overflow;}

};

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F, true>::fixed_point()
: integer (0)
, fraction (0)
, dbl_orig(0)
, factor(1ull << F)
, mask_int ((1ull << I) - 1)
, mask_frac ((1ull << F) - 1)
, overflow (false)
, negative (false)
{
  static_assert(std::is_integral<B>::value);

  // required for all operations
  static_assert(std::numeric_limits<UB>::digits > I);
  static_assert(std::numeric_limits<UB>::digits > F);

}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F, true>::fixed_point(const fixed_point& fp_arg)
{
  static_assert(std::is_integral<B>::value);

  // required for all operations
  static_assert(std::numeric_limits<UB>::digits > I);
  static_assert(std::numeric_limits<UB>::digits > F);

  integer     = fp_arg.integer;
  fraction    = fp_arg.fraction;
  dbl_orig    = fp_arg.dbl_orig;
  negative    = fp_arg.negative;
  factor      = fp_arg.factor;
  mask_frac   = fp_arg.mask_frac;
  mask_int    = fp_arg.mask_int;
  overflow    = false;
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F, true>::fixed_point(UB arg_integer, UB arg_fraction, bool negative)
: integer (arg_integer)
, fraction (arg_fraction)
, factor(1ull << F)
, mask_int ((1ull << I) - 1)
, mask_frac ((1ull << F) - 1)
, overflow (false)
, negative (negative)
{
  static_assert(std::is_integral<B>::value);

  // required for all operations
  static_assert(std::numeric_limits<UB>::digits > I);
  static_assert(std::numeric_limits<UB>::digits > F);

  dbl_orig = (double)*this;
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F, true>::fixed_point(double dbl_arg)
: dbl_orig(dbl_arg)
, factor(1ull << F)
, mask_int ((1ull << I) - 1)
, mask_frac ((1ull << F) - 1)
, overflow (false)
{
  static_assert(std::is_integral<B>::value);

  // required for all operations
  static_assert(std::numeric_limits<UB>::digits > I);
  static_assert(std::numeric_limits<UB>::digits > F);

  if (dbl_arg < 0)
  {
    negative = true;
    dbl_arg *= -1;
  }
  else
  {
    negative = false;
  }

  integer             = static_cast <UB> (dbl_arg + 0.5);
  double  remainder   = dbl_arg - integer;
  fraction            = (UB)(std::round(remainder * factor)) & mask_frac;
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F, true>::operator double() const
{
  return (negative ? -1.0 : 1.0) * (integer + (double)fraction / factor);
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F, true>::operator float() const
{
  return (negative ? -1.0 : 1.0) * (integer + (float)fraction / factor);
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F, true>
operator+(const fixed_point<B, I, F, true>& fp_lhs, const fixed_point<B, I, F, true>& fp_rhs)
{
  fixed_point <B, I, F> r(fp_lhs);
  if (fp_lhs.is_negative() != fp_rhs.is_negative())
    return r -= fp_rhs;
  else
    return r += fp_rhs;
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F, true>
operator-(const fixed_point<B, I, F, true>& fp_lhs, const fixed_point<B, I, F, true>& fp_rhs)
{
  fixed_point <B, I, F> r(fp_lhs);
  if (fp_lhs.is_negative() != fp_rhs.is_negative())
    return r += fp_rhs;
  else
    return r -= fp_rhs;
}

}      // NAMESPACE SUSA
#endif // SUSA_FIXEDPOINT_H
