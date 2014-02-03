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
 * @brief Random Number Generator
 * This class uses the same code as creators of MT (The Mersenne Twister)
 * made available on the web.
 *
 * @author Behrooz, Kamary Aliabadi
 * @version 1.0.0
 *
 * @defgroup RNG Random Number Generator
 */

#ifndef MT_H
#define MT_H

#define _MT_N 624
#define _MT_M 397
#define _MT_MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define _MT_UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define _MT_LOWER_MASK 0x7fffffffUL /* least significant r bits */


namespace susa {

/**
 * @brief  The Mersenne Twister pseudorandom number generator wrapper class.
 *
 *
 * @ingroup RNG
*/

class mt {
  private:

    unsigned long *ul_mt; /* the array for the state vector  */
    int mti;              /* mti==N+1 means mt[N] is not initialized */

  public:

    mt(void);
    mt(unsigned long ul_seed);
    ~mt();

    void init_genrand(unsigned long s);
    void init_by_array(unsigned long init_key[], int key_length);

    unsigned long genrand_int32(void); // generates unsigned 32-bit integers.
    long genrand_int31(void);          // generates unsigned 31-bit integers.

    double genrand_real1(void);        // generates uniform real in [0,1] (32-bit resolution).
    double genrand_real2(void);        // generates uniform real in [0,1) (32-bit resolution).
    double genrand_real3(void);        // generates uniform real in (0,1) (32-bit resolution).
    double genrand_res53(void);        // generates uniform real in [0,1) with 53-bit resolution.

    unsigned int rand_mask(unsigned int uint_mask);
    matrix <unsigned int> rand_mask(unsigned int uint_mask, unsigned int uint_N);
    matrix <double> randn(unsigned int uint_N);
    matrix <double> rand(unsigned int uint_N);

};
} // NAMESPACE SUSA
#endif
