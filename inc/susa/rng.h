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
 * @file rng.h
 * @brief Random Number Generator (declaration).
 *
 *
 * @author Behrooz Kamary
 */

#ifndef SUSA_RNG_H
#define SUSA_RNG_H

#include <cstdint>

namespace susa {

/**
* @brief Random Number Generator class.
*
*
* This implementation of "Mersenne Twister" generates uniformly
* 32 bit integers in the range [0, 2^32 - 1] using MT19937 algorithm.<br>
* Mersenne Twister - MT19937<br>
* (w, n, m, r) = (32, 624, 397, 31)<br>
* a = 0x9908B0DF<br>
* u = 11<br>
* (s, b) = (7, 0x9D2C5680)<br>
* (t, c) = (15, 0xEFC60000)<br>
* l = 18<br>
*
*/
class rng {

  public:

    //! Constructor
    rng();

    /**
     * @brief Constructor
     *
     * @param uint_seed RNG initial seed
     */
    rng(unsigned int uint_seed);

    /**
     * @brief Initializes the RNG
     *
     * @param uint_seed RNG initial seed
     */
    void init(unsigned int uint_seed);

    /**
     * @brief gaussian distributed random double
     * with mean value equal to zero and unit variance.
     *
     */
    double randn();

    /**
     * @brief uniformly distributed random double
     *
     */
    double rand();

    /**
     * @brief uniformly distributed random unsigned integers
     *
     * The output numbers will be masked to get uniformly distributed numbers
     * in a limited range. As an example masking with 0x0F would produce random
     * numbers between 0 to 15.
     *
     * @param uint_mask used to mask output numbers
     */
    unsigned int rand_mask(unsigned int uint_mask);

    /**
     * @brief uniformly distributed random unsigned integers
     *
     * The output numbers will be masked to get uniformly distributed numbers
     * in a limited range. As an example masking with 0x0F would produce random
     * numbers between 0 to 15.
     *
     * @param uint_mask used to mask output numbers
     * @param uint_num number of output samples
     */
    matrix <unsigned int> rand_mask(unsigned int uint_mask, unsigned int uint_num);

    /**
     * @brief gaussian distributed random double
     *
     * @param uint_num Number of random numbers
     * @return a column vector of type, susa::matrix<double>
     */
    matrix <double> randn(unsigned int uint_num);

    /**
     * @brief uniformly distributed random double
     *
     * @param uint_num Number of random numbers
     * @return a column vector of type susa::matrix<double>
     */
    matrix <double> rand(unsigned int uint_num);

    /**
     * @brief uniformly distributed random double
     *
     * @param uint_rows Number of rows
     * @param uint_cols Number of columns
     * @return a random susa::matrix<double> instance
     */
    matrix <double> rand(unsigned int uint_rows, unsigned int uint_cols);

    /**
     * @brief Bernoulli random samples
     * 
     * @param size_num number of samples i.e. vector size
     * @return a column vector of type susa::matrix<unsigned char>
     */
    matrix <unsigned char> bernoulli(size_t size_num);

    enum
    {
        W = 32,
        N = 624,
        M = 397,
        R = 31,
        A = 0x9908B0DF,

        F = 1812433253,

        U = 11,
        D = 0xFFFFFFFF,

        S = 7,
        B = 0x9D2C5680,

        T = 15,
        C = 0xEFC60000,

        L = 18,

        MASK_LOWER = (1ull << R) - 1,
        MASK_UPPER = (1ull << R)
    };

  private:
    // Create a length 624 array to store the state of the generator
    uint32_t MT[N];
    uint16_t uint_index;

    void generate_numbers();
    unsigned int extract_number();


    unsigned int GetUInt();
    unsigned int GetNonUniform(double* pr, unsigned int n);
    unsigned int nonUniform(std::vector <float>);
    double       GetDouble();
};

}      // NAMESPACE SUSA
#endif // SUSA_RNG_H
