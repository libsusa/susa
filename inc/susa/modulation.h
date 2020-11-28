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

#ifndef MODULATION_H
#define MODULATION_H

namespace susa {

/**
 * @brief BPSK Modulation
 * Converts the input matrix with 0,1 elements to -1,1 elements.
 * @ingroup Communications
 */
matrix <int8_t> bpsk(const matrix <uint8_t> &mat_arg);


/**
 * @brief QAM Modulation
 * This is a wrapper class for QAM modulation and demodulation methods.
 * @ingroup Communications
 */
class qam {

  public:
    /**
     * @brief The constructor of the QAM modulation class
     * @param uint_m The modulation order, e.g., 4, 8, 16, 64, 128, 256
     */
    qam(unsigned int uint_m);

    ~qam();

    matrix < std::complex <double> > get_constellation()
    {
        return mat_s;
    }

    /**
     * @brief compute the AWGN noise deviation
     * @param dbl_arg Eb / N_0 in dB.
     */
    double get_noise_deviation(double dbl_arg);

    std::complex <double> demodulate_symbol(std::complex <double> complex_arg);

  private:
    double          dbl_es;          // Energy per symbol
    double          dbl_eb;          // Energy per bit
    unsigned int    uint_bps;        // Number of bits per symbol
    unsigned int    uint_m;          // Number of symbols in constellation

    matrix < std::complex <double> > mat_s;
    matrix < std::complex <double> > mat_axis;

};

/**
 * @brief add Additive white Gaussian noise (AWGN) to the given signal
 *
 * AWGN is a basic noise model used in information theory to mimic
 * the effect of many random processes that occur in nature.
 *
 * @param mat_signal the input signal can be a real or a complex matrix
 * @param ovs oversamplaing rate
 *
 * @ingroup Communications
 */
template <typename T, template <typename> typename Allocator>
auto awgn(matrix <T, Allocator<T>> mat_signal, float flt_snr_db, unsigned ovs = 1)
{
    // Es = k * Eb
    // k = log2(M)
    // gamma_s = Es / N0 = k * gamma_b
    // N0 = P / gamma
    // Noise_var = N0 / 2

    double P;
    if constexpr (is_complex<T>::value)
    {
        P = ovs * sum(sum(mag(mat_signal)))(0)/mat_signal.size();
    }
    else
    {
        P = ovs * sum(sum(mat_signal * mat_signal))(0)/mat_signal.size();
    }

    double gamma        = std::pow(10, flt_snr_db/10);
    double N0           = P / gamma;
    double dbl_std_dev  = std::sqrt(0.5f * N0);

    rng  random;

    std::conditional_t<is_complex<T>::value,
                      matrix <std::complex<double>, Allocator<std::complex<double>>>,
                      matrix <double, Allocator<double>>> mat_noise(mat_signal.shape());


    if constexpr (is_complex<T>::value)
    {
        for (size_t indx = 0; indx < mat_noise.size(); indx++)
            mat_noise(indx) = dbl_std_dev * std::complex<double> (random.randn(), random.randn());
    }
    else
    {
        for (size_t indx = 0; indx < mat_noise.size(); indx++)
            mat_noise(indx) = dbl_std_dev * random.randn();
    }

    return mat_noise;
}

} // NAMESPACE SUSA
#endif // MODULATION_H
