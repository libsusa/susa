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
 * @brief Random Number Generator
 *
 *
 * @author Behrooz, Kamary Aliabadi
 * @version 1.0.0
 */

#ifndef RNG_H
#define RNG_H

namespace susa {

/**
* @brief Random Number Generator wrapper class.
*
* @author Behrooz, Kamary Aliabadi
*
* This implementation of "Mersenne Twister" generates uniformly 32 bit integers in the range [0, 2^32 - 1] using MT19937 algorithm.<br>
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
     * @param intSeed RNG initial seed
     */
    rng(unsigned int intSeed);

    /**
     * @brief Initializes the RNG
     *
     * @param intSeed RNG initial seed
     */
    void init(unsigned int intSeed);

    /**
     * @brief Returns uniformly distributed random double
     *
     * @see rand()
     */
    double GetDouble();

    /**
     * @brief Returns gaussian distributed random double
     *
     * With mean value equal to zero and unit variance.
     *
     */
    double randn();

    /**
     * @brief Returns uniformly distributed random double
     *
     * @see GetDouble()
     */
    double rand();

    /**
     * @brief Returns uniformly distributed random unsigned integers
     *
     * The output numbers will be masked to get uniformly distributed numbers
     * in a limited range. As an example masking with 0x0F would produce random
     * numbers between 0 to 15.
     *
     * @param uint_mask used to mask output numbers
     */
    unsigned int rand_mask(unsigned int uint_mask);

    /**
     * @brief Returns uniformly distributed random unsigned integers
     *
     * The output numbers will be masked to get uniformly distributed numbers
     * in a limited range. As an example masking with 0x0F would produce random
     * numbers between 0 to 15.
     *
     * @param uint_mask used to mask output numbers
     * @param uint_N number of output samples
     */
    matrix <unsigned int> rand_mask(unsigned int uint_mask, unsigned int uint_N);

    /**
     * @brief Returns gaussian distributed random double
     *
     * @param uint_N Number of random numbers
     * @return Susa matrix
     */
    matrix <double> randn(unsigned int uint_N);

    /**
     * @brief Returns uniformly distributed random double
     *
     * @param uint_N Number of random numbers
     * @return Susa matrix
     */
    matrix <double> rand(unsigned int uint_N);

    unsigned int GetUInt();

    unsigned int GetNonUniform(double* pr, unsigned int n);

    unsigned int nonUniform(std::vector <float>);

  private:
    // Create a length 624 array to store the state of the generator
    unsigned int MT[624];
    unsigned int y;
    int intIndex;

    void generateNumbers();
    unsigned int extractNumber(int);
};
}

#endif // RNG_H


