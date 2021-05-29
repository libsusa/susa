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

    integer             = static_cast <UB> (dbl_arg + 0.5);
    double  remainder   = dbl_arg - integer;
    integer           <<= F;
    integer            |= static_cast <UB> (std::round(remainder * factor)) & mask_frac;

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
    if (negative != (int_arg < 0)) negative = true;
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
    // avoid loosing precision when multiplying by one
    if (integer == (1ul << F))
    {
      negative  = (negative != fp_rhs.negative);
      integer   = fp_rhs.integer;
      dbl_orig = (double)*this;
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
    if (integer > fp_rhs.integer)
      integer -= fp_rhs.integer;
    else
    {
      integer = fp_rhs.integer - integer;
      negative = true;
    }
    dbl_orig = (double)*this;
    return *this;
  }

  fixed_point& operator+=(const fixed_point& fp_rhs)
  {
    integer += fp_rhs.integer;
    overflow = (integer & ~mask_tot) > 0 ? true : false;
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

  integer             = static_cast <UB> (dbl_arg + 0.5);
  double  remainder   = dbl_arg - integer;
  integer           <<= F;
  integer            |= static_cast <UB> (std::round(remainder * factor)) & mask_frac;
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F, false>::operator double() const
{
  return (negative ? -1.0 : 1.0) * (((integer) >> F) + (double)(integer & mask_frac) / factor);
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F, false>::operator float() const
{
  return (negative ? -1.0 : 1.0) * (((integer) >> F) + (float)(integer & mask_frac) / factor);
}

/**
 * @class fixed_point
 * @brief <i>fixed-point</i> type represents a fixed-point number that is defined by its base template B,
 * the number of integer bits I and the number of fraction bits F. I and F can take values up to the number of bits
 * in the base type. F can not be more than half of the number of bits in <b>uintmax_t</b>.
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
    if (fp_rhs.integer == 1ul)
    {
      negative = (negative != fp_rhs.negative);
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
fixed_point<B, I, F>
operator+(const fixed_point<B, I, F>& fp_lhs, const fixed_point<B, I, F>& fp_rhs)
{
  fixed_point <B, I, F> r(fp_lhs);
  if (fp_lhs.is_negative() != fp_rhs.is_negative())
    return r -= fp_rhs;
  else
    return r += fp_rhs;
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F>
operator-(const fixed_point<B, I, F>& fp_lhs, const fixed_point<B, I, F>& fp_rhs)
{
  fixed_point <B, I, F> r(fp_lhs);
  if (fp_lhs.is_negative() != fp_rhs.is_negative())
    return r += fp_rhs;
  else
    return r -= fp_rhs;
}

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

}      // NAMESPACE SUSA
#endif // SUSA_FIXEDPOINT_H
