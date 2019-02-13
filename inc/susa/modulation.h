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
 * @brief The base-band modulation and demodulation
 *
 * <i>Modulation</i> in the communications theory is the technique of mapping
 * digital and discrete information to a measurable quantity of an analog
 * signal such as amplitude, frequency and phase.
 *
 * @author Behrooz Kamary Aliabadi
 * @version 1.0.0
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
     * @param uint_m The modulation order, e.g., 4,8,16,64,256
     */
    qam(unsigned int uint_m);

    ~qam();

    matrix < std::complex <double> > get_constellation()
    {
        return mat_s;
    }


    /**
     * @brief This method computes the AWGN noise deviation.
     * @param dbl_arg Eb / N_0 in dB.
     */
    double get_noise_deviation(double dbl_arg);

    matrix < std::complex <double> > modulate_bits(matrix <char> mat_bits);

    matrix <char>   demodulate_bits(const matrix < std::complex <double> >& mat_symbols);

    matrix <int>    demodulate_symbols(const matrix < std::complex <double> >& mat_symbols);

    std::complex <double> demodulate_symbol(std::complex <double> complex_arg);

  private:
    double          dbl_es;          // Energy per symbol
    double          dbl_eb;          // Energy per bit
    unsigned int    uint_bps;        // Number of bits per symbol
    unsigned int    uint_m;          // Number of symbols in constellation

    matrix < std::complex <double> > mat_s;
    matrix < std::complex <double> > mat_axis;

    unsigned int log2(unsigned int uint_x);

};

} // NAMESPACE SUSA
#endif // MODULATION_H
