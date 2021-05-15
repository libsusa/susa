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
 * @file modulation.h
 * @brief base-band modulation and demodulation schemes (declaration).
 *
 * <i>Modulation</i> in the communications theory is the technique of mapping
 * digital and discrete information to a measurable quantity of an analog
 * signal such as amplitude, frequency and phase.
 *
 * @author Behrooz Kamary
 */

#ifndef SUSA_MODULATION_H
#define SUSA_MODULATION_H

namespace susa {

/**
 * @brief BPSK Modulation
 * converts the input matrix with 0,1 elements to -1,1 elements.
 * @param mat_arg input binary matrix
 * @ingroup Communications
 */
matrix <int8_t> bpsk(const matrix <uint8_t>& mat_arg);

/**
 * @brief BPSK Modulation
 * converts the input matrix with 0,1 elements to -1,1 elements.
 * @param bset_arg input binary bitset
 * @ingroup Communications
 */
matrix <int8_t> bpsk(const bitset<>& bset_arg);

/**
 * @brief QAM Modulation
 * This is a wrapper class for QAM modulation and demodulation methods.
 * @ingroup Communications
 */
class qam
{

  public:

    /**
     * @brief The constructor of the QAM modulation class
     * @param uint_m The modulation order, e.g., 4, 8, 16, 64, 128, 256
     */
    qam(unsigned int uint_m);

    qam() = delete;
    ~qam();

    //! returns the constructed QAM constellation matrix
    cmatrix <double> get_constellation()
    {
      return mat_s;
    }

    /**
     * @brief compute the AWGN noise deviation
     * @param dbl_arg Eb / N_0 in dB.
     */
    double get_noise_deviation(double dbl_arg);

    /**
     * @brief detects a symbol using Eucleadian distance metrics
     */
    std::complex <double>           demod_sym(std::complex <double> complex_arg);

    /**
     * @brief detects a symbol two diemensioanl indecies using Eucleadian distance metrics
     */
    std::tuple<unsigned, unsigned>  demod_to_tuple(std::complex <double> complex_arg);

    /**
     * @brief detects a symbol linear index using Eucleadian distance metrics
     */
    unsigned int                    demod_to_index(std::complex <double> complex_arg);

    /**
     * @brief detects bits set using Eucleadian distance metrics
     */
    bitset<>                        demod_to_bits(const cmatrix <double>& cmat_arg);

    /**
     * @brief maps the input bits set to symbols matrix
     */
    cmatrix <double>                mod_bits(const bitset<>& bset_arg);

  private:

    double          dbl_es;          // Energy per symbol
    double          dbl_eb;          // Energy per bit
    unsigned int    uint_bps;        // Number of bits per symbol
    unsigned int    uint_m;          // Number of symbols in constellation

    cmatrix <double> mat_s;
    matrix  <double> mat_axis;

    unsigned bin_to_gray(unsigned uint_arg)
    {
      return (uint_arg ^ (uint_arg / 2));
    }

    unsigned gray_to_bin(unsigned uint_arg)
    {
      unsigned mask = uint_arg;
      while (mask)
      {
        mask /= 2;
        uint_arg ^= mask;
      }
      return uint_arg;
    }
};

/**
 * @brief add Additive White Gaussian Noise (AWGN) to the input signal
 *
 * AWGN is a basic noise model used in information theory to mimic
 * the effect of many random processes that occur in nature.
 *
 * @param mat_signal the input signal can be a real or a complex matrix
 * @param uint_ovs the input signal oversamplaing rate
 *
 * @ingroup Communications
 */
template <typename T, template <typename> typename Allocator>
auto awgn(matrix<T, Allocator<T>> mat_signal, float flt_snr_db, unsigned uint_ovs = 1)
{
  // Es = k * Eb
  // k = log2(M)
  // gamma_s = Es / N0 = k * gamma_b
  // N0 = P / gamma
  // Noise_var = N0 / 2

  double P;
  if constexpr (is_complex<T>::value)
  {
    P = uint_ovs * sum(sum(mag(mat_signal)))(0) / mat_signal.size();
  }
  else
  {
    P = uint_ovs * sum(sum(mat_signal * mat_signal))(0) / mat_signal.size();
  }

  double gamma = std::pow(10, flt_snr_db / 10);
  double N0 = P / gamma;
  double dbl_std_dev = std::sqrt(0.5f * N0);

  rng random;

  std::conditional_t<is_complex<T>::value,
                     matrix<std::complex<double>, Allocator<std::complex<double>>>,
                     matrix<double, Allocator<double>>>
      mat_noise(mat_signal.shape());

  if constexpr (is_complex<T>::value)
  {
    for (size_t indx = 0; indx < mat_noise.size(); indx++)
      mat_noise(indx) = dbl_std_dev * std::complex<double>(random.randn(), random.randn());
  }
  else
  {
    for (size_t indx = 0; indx < mat_noise.size(); indx++)
      mat_noise(indx) = dbl_std_dev * random.randn();
  }

  return mat_noise;
}

} // NAMESPACE SUSA
#endif // SUSA_MODULATION_H
