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
 * @file susa.h
 * @brief Fixed-Point Number Type and Arithmetics
 * @author Behrooz Kamary
 */

#ifndef SUSA_ARITHMETICS_H
#define SUSA_ARITHMETICS_H

#include <limits>
#include <cstdint>
#include <type_traits>

namespace susa {

template <typename B, unsigned char I, unsigned char F, typename Enable = void> class fixed_point;

template <typename B, unsigned char I, unsigned char F>
class fixed_point <B, I, F, typename std::enable_if<std::is_integral<B>::value>::type>
{

private:
    
    using UB = typename std::make_unsigned<B>::type;
    const unsigned NB = std::numeric_limits<uintmax_t>::digits / 2;

    double  dbl_orig;

    UB      integer;
    UB      fraction;

    UB      factor;
    UB      mask_int;
    UB      mask_frac;
    bool    overflow;
    bool    negative;

public:

    fixed_point();
    fixed_point(double dbl_arg);

    operator double() const;
    operator float() const;
    operator int() const;
    operator bool() const;

    bool operator!() const;
    fixed_point operator~() const;
    fixed_point operator-() const;
    fixed_point& operator|=(const fixed_point& arg_rhs);
    fixed_point& operator&=(const fixed_point& arg_rhs);
    fixed_point& operator^=(const fixed_point& arg_rhs);

    fixed_point& operator+=(const fixed_point& arg_rhs)
    {
        fixed_point r(arg_rhs);
        fraction += r.fraction;
        UB carry = (fraction & ~mask_frac) >> F;
        fraction &= mask_frac;
        integer += r.integer;
        integer += carry;
        overflow = (integer >> I) > 0;
        integer &= mask_int;
        return *this;
    }

    fixed_point& operator-=(const fixed_point& arg_rhs)
    {
        fixed_point r(arg_rhs);
        if (fraction < r.fraction)
        {
            fraction += (1ull << F);
            integer--;
        }
        fraction -= r.fraction;
        UB carry  = (fraction & ~mask_frac) >> F;
        fraction &= mask_frac;
        integer  -= r.integer;
        integer  -= carry;
        overflow = (integer >> I) > 0;
        integer  &= mask_int;
        return *this;
    }

    fixed_point& operator*=(const fixed_point& arg_rhs)
    {
        fixed_point r(arg_rhs);
        uintmax_t intint        = (uintmax_t)integer * r.integer;
        uintmax_t fracint       = (uintmax_t)fraction * r.integer;
        uintmax_t intfrac       = (uintmax_t)integer * r.fraction;
        uintmax_t fracfrac      = (uintmax_t)fraction * r.fraction;

        integer  = intint + (fracint >> F) + (intfrac >> F) + ((fracfrac >> F) >> F);
        fraction = (fracint & mask_frac) + (intfrac & mask_frac) +  ((fracfrac >> F) & mask_frac);
        overflow = (integer >> I) > 0;
        integer  &= mask_int;

        if (negative != r.is_negative()) negative = true;
        
        return *this;  
    }

    fixed_point& operator/=(const fixed_point& arg_rhs);
    fixed_point& operator%=(const fixed_point& arg_rhs);
    fixed_point& operator<<=(int amount);
    fixed_point& operator>>=(int amount);

    bool is_negative() const {return negative;}

};

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F, typename std::enable_if<std::is_integral<B>::value>::type>::fixed_point()
: dbl_orig(0)
, factor(1ull << F)
, mask_int ((1ull << I) - 1)
, mask_frac ((1ull << F) - 1)
, overflow (false)
, negative (false)
{
    integer     = 0;
    fraction    = 0;
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F,
typename std::enable_if<std::is_integral<B>::value>::type>::
fixed_point(double dbl_arg)
: dbl_orig(dbl_arg)
, factor(1ull << F)
, mask_int ((1ull << I) - 1)
, mask_frac ((1ull << F) - 1)
, overflow (false)
{
    // required for all operations
    static_assert(std::numeric_limits<UB>::digits > I);
    static_assert(std::numeric_limits<UB>::digits > F);

    // required for multiplication
    static_assert(std::numeric_limits<UB>::digits > 2 * F);

    if (dbl_arg < 0)
    {
        negative = true;
        dbl_arg *= -1;
    }
    else
    {
        negative = false;
    }
    integer             = static_cast <UB> (dbl_arg > 0 ? std::floor(dbl_arg) : std::ceil(dbl_arg));    
    double  remainder   = dbl_arg - integer;
    fraction            = (UB)(std::round(remainder * factor)) & mask_frac;
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F, typename std::enable_if<std::is_integral<B>::value>::type>::
operator double() const
{
    return (negative ? -1.0 : 1.0) * (integer + (double)fraction / factor);
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F, typename std::enable_if<std::is_integral<B>::value>::type>::
operator float() const
{
    return (negative ? -1.0 : 1.0) * (integer + (float)fraction / factor);
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F>
operator+(const fixed_point<B, I, F>& arg_lhs, const fixed_point<B, I, F>& arg_rhs)
{
    fixed_point <B, I, F> r(arg_lhs);
    if (arg_lhs.is_negative() != arg_rhs.is_negative())
        return r -= arg_rhs;
    else
        return r += arg_rhs;
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F>
operator-(const fixed_point<B, I, F>& arg_lhs, const fixed_point<B, I, F>& arg_rhs)
{
    fixed_point <B, I, F> r(arg_lhs);
    if (arg_lhs.is_negative() != arg_rhs.is_negative())
        return r += arg_rhs;
    else
        return r -= arg_rhs;
}

template <typename B, unsigned char I, unsigned char F>
fixed_point<B, I, F>
operator*(const fixed_point<B, I, F>& arg_lhs, const fixed_point<B, I, F>& arg_rhs)
{
    fixed_point <B, I, F> r(arg_lhs);
    return r *= arg_rhs;
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