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
 * @file mt.h
 * @brief Mersenne Twister Random Number Generator (declaration).
 *
 *
 * @author Behrooz Kamary
 * @version 1.0.0
 *
 * @defgroup RNG Random Number Generator
 */

#ifndef SUSA_MT_H
#define SUSA_MT_H

#include <cstdint>

#define _MT_N 624
#define _MT_M 397
#define _MT_MATRIX_A   0x9908b0dfUL /* constant vector a */
#define _MT_UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define _MT_LOWER_MASK 0x7fffffffUL /* least significant r bits */


namespace susa {

/**
 * @brief  The Mersenne Twister pseudorandom number generator class.
 *
 * This class uses the same code as the creators of MT (The Mersenne Twister) have been published.
 *
 * @ingroup RNG
*/

class mt {

  private:
    matrix <double>       mat_probabilities;
    matrix <unsigned int> mat_sort_indices;
    matrix <double>       mat_prob_sorted;
    uint32_t*             uint_mt; /* the array for the state vector  */
    int                   mti;   /* mti==N+1 means mt[N] is not initialized */

    void init_genrand(uint32_t s);
    void init_by_array(uint32_t init_key[], int key_length);

    uint32_t genrand_int32(void);      // generates unsigned 32-bit integers.
    unsigned int genrand_int31(void);  // generates unsigned 31-bit integers.

    double genrand_real1(void);        // generates uniform real in [0,1] (32-bit resolution).
    double genrand_real2(void);        // generates uniform real in [0,1) (32-bit resolution).
    double genrand_real3(void);        // generates uniform real in (0,1) (32-bit resolution).
    double genrand_res53(void);        // generates uniform real in [0,1) with 53-bit resolution.

  public:
    //! Constructor
    mt(void);

    /**
     * @brief Constructor
     *
     * @param uint_seed the random number generator initial seed
     */
    mt(uint32_t uint_seed);

    //! Destructor
    ~mt();

    /**
     * @brief Bernoulli
     *
     * @param dbl_p the head probability i.e. 1's probability
     */
    unsigned int bernoulli(double dbl_p);

    /**
     * @brief Random integer
     *
     * @param uint_max the maximum random integer
     */
    unsigned int randint(unsigned int uint_max);
    
    unsigned int rand_mask(unsigned int uint_mask);

    matrix <unsigned int> rand_mask(unsigned int uint_mask, unsigned int uint_num);

    matrix <double> randn(unsigned int uint_num);

    matrix <double> rand(unsigned int uint_num);

    void set_probabilities(matrix <double> mat_probabilities);

    unsigned int nonuniform(matrix <double> mat_probabilities);

    unsigned int nonuniform();
    
    matrix <unsigned int> nonuniform(matrix <double> mat_probabilities, unsigned int uint_length);
    
    matrix <unsigned int> nonuniform(unsigned int uint_length);
};

}        // NAMESPACE SUSA
#endif   // SUSA_MT_H
