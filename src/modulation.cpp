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
 * @file modulation.cpp
 * @brief digital base-band modulation and demodulation schemes (definiton).
 *
 * @author Behrooz Kamary
 */


#include "susa.h"

namespace susa {

matrix <int8_t> bpsk(const matrix <uint8_t>& mat_arg)
{
  matrix <int8_t> mat_ret(mat_arg.shape());
  size_t uint_length = mat_arg.size();

  for (size_t uint_index = 0; uint_index < uint_length; uint_index++)
  {
    mat_ret(uint_index) = (mat_arg(uint_index) != 0) ? 1 : -1;
  }

  return mat_ret;
}

matrix <int8_t> bpsk(const bitset<>& bset_arg)
{
  matrix <int8_t> mat_ret(bset_arg.size(), 1);
  size_t uint_length = bset_arg.size();

  for (size_t uint_index = 0; uint_index < uint_length; uint_index++)
  {
    mat_ret(uint_index) = (bset_arg(uint_index)) ? 1 : -1;
  }

  return mat_ret;
}

qam::qam(unsigned int uint_m)
: uint_m(uint_m)
{
    uint_bps                = log2u(uint_m);
    double          dbl_p   = std::sqrt(uint_m);
    unsigned int    uint_l  = (unsigned int)(uint_m/dbl_p);

    // constellation generation
    mat_s    = cmatrix <double> (uint_l, uint_l, std::complex <double>(0,0));
    mat_axis = matrix <double> ((unsigned int)dbl_p, 1, 0);

    for (int i = -(dbl_p - 1), j=0; i < dbl_p; j++, i += 2)
    {
      mat_axis(j) = i;
    }

    for (unsigned int i=0; i < uint_l; i++)
    {
      for (unsigned int j=0; j < uint_l; j++)
      {
        mat_s(i,j) = std::complex <double> (mat_axis(i), mat_axis(j));
      }
    }

    // energy computations
    matrix <double> mat_mag     = mag(mat_s);
    matrix <double> mat_sum     = sum(mat_mag);
    dbl_es                      = sum(mat_sum)(0)/uint_m;
    dbl_eb                      = dbl_es / uint_bps;
}

std::complex <double> qam::demod_sym(std::complex <double> complex_arg)
{
  matrix <double> mat_dist(uint_m,1);
  for (unsigned int uint_i = 0; uint_i < uint_m; uint_i++)
  {
    mat_dist(uint_i) = std::abs(complex_arg - mat_s(uint_i));
  }
  return mat_s(min(mat_dist)(0));
}

std::tuple<unsigned, unsigned> qam::demod_to_tuple(std::complex <double> complex_arg)
{
  std::tuple <unsigned,unsigned> ret;
  uint32_t uint_length = std::sqrt(uint_m);
  double dbl_dist_min = std::numeric_limits<double>::max();
  for (uint32_t uint_i = 0; uint_i < uint_length; uint_i++)
  {
    for (uint32_t uint_q = 0; uint_q < uint_length; uint_q++)
    {
      double dbl_dist = std::abs(complex_arg - mat_s(uint_i, uint_q));
      if (dbl_dist < dbl_dist_min)
      {
        ret = std::make_tuple(uint_i, uint_q);
        dbl_dist_min = dbl_dist;
      }
    }
  }
  return ret;
}

double qam::get_noise_deviation(double dbl_arg)
{
  return std::sqrt(0.5f * dbl_eb * std::pow(10, - (dbl_arg/10)));
}

cmatrix <double> qam::mod_bits(const bitset<>& bset_arg)
{
  SUSA_ASSERT_MESSAGE((bset_arg.size() % uint_bps) == 0, "number of bits must be a factor of log2(M)");
  size_t sz_num_symbols = bset_arg.size() / uint_bps;
  matrix <std::complex<double>> mat_ret (1, sz_num_symbols);

  for (size_t sz_sym = 0; sz_sym < sz_num_symbols; sz_sym++)
  {
    unsigned uint_i(0), uint_q(0);
    for (unsigned uint_bit = 0; uint_bit < (uint_bps / 2); uint_bit++)
    {
      if (bset_arg.exists(sz_sym * uint_bps + uint_bit)) uint_i |= (1 << uint_bit);
      if (bset_arg.exists(sz_sym * uint_bps + (uint_bps / 2) + uint_bit)) uint_q |= (1 << uint_bit);
    }

    uint_i = bin_to_gray(uint_i);
    uint_q = bin_to_gray(uint_q);

    mat_ret(sz_sym) = mat_s(uint_i, uint_q);

  }

  return mat_ret;
}

bitset<> qam::demod_to_bits(const cmatrix <double>& cmat_arg)
{
  size_t    sz_num_syms = cmat_arg.size();
  bitset<>  ret(sz_num_syms * uint_bps);

  for (size_t index = 0; index < sz_num_syms; index++)
  {
    unsigned uint_i(0), uint_q(0);

    std::tuple<unsigned, unsigned> sym_indx = demod_to_tuple(cmat_arg(index));
    uint_i = gray_to_bin(std::get<0>(sym_indx));
    uint_q = gray_to_bin(std::get<1>(sym_indx));

    for (unsigned uint_bit = 0; uint_bit < (uint_bps / 2); uint_bit++)
    {
      unsigned mask = 1 << uint_bit;
      if (uint_i & mask) ret.set(index * uint_bps + uint_bit);
      if (uint_q & mask) ret.set(index * uint_bps + (uint_bps / 2) + uint_bit);
    }
  }

  return ret;
}

qam::~qam() {}

} // NAMESPACE SUSA
