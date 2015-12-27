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
 * @file mt.cpp
 * @brief Random Number Generator
 * @author Behrooz, Kamary Aliabadi
 * @version 1.0.0
 */

#include "../inc/susa.h"

namespace susa {

mt::mt() {
    ul_mt = new unsigned long[_MT_N];
    mti = _MT_N+1;
}

mt::mt(unsigned long ul_seed) {
    ul_mt = new unsigned long[_MT_N];
    mti = _MT_N+1;

    init_genrand(ul_seed);
}

mt::~mt() {
    delete [] ul_mt;
}

/* initializes ul_mt[_MT_N] with a seed */
void mt::init_genrand(unsigned long s) {
    ul_mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<_MT_N; mti++) {
        ul_mt[mti] = (1812433253UL * (ul_mt[mti-1] ^ (ul_mt[mti-1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array ul_mt[].                     */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        ul_mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void mt::init_by_array(unsigned long init_key[], int key_length) {
    int i, j, k;
    init_genrand(19650218UL);
    i=1;
    j=0;
    k = (_MT_N>key_length ? _MT_N : key_length);
    for (; k; k--) {
        ul_mt[i] = (ul_mt[i] ^ ((ul_mt[i-1] ^ (ul_mt[i-1] >> 30)) * 1664525UL))
                   + init_key[j] + j; /* non linear */
        ul_mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        j++;
        if (i>=_MT_N) {
            ul_mt[0] = ul_mt[_MT_N-1];
            i=1;
        }
        if (j>=key_length) j=0;
    }
    for (k=_MT_N-1; k; k--) {
        ul_mt[i] = (ul_mt[i] ^ ((ul_mt[i-1] ^ (ul_mt[i-1] >> 30)) * 1566083941UL))
                   - i; /* non linear */
        ul_mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=_MT_N) {
            ul_mt[0] = ul_mt[_MT_N-1];
            i=1;
        }
    }

    ul_mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long mt::genrand_int32(void) {
    unsigned long y;
    static unsigned long mag01[2]= {0x0UL, _MT_MATRIX_A};
    /* mag01[x] = x * _MT_MATRIX_A  for x=0,1 */

    if (mti >= _MT_N) { /* generate _MT_N words at one time */
        int kk;

        if (mti == _MT_N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0; kk<_MT_N-_MT_M; kk++) {
            y = (ul_mt[kk]&_MT_UPPER_MASK)|(ul_mt[kk+1]&_MT_LOWER_MASK);
            ul_mt[kk] = ul_mt[kk+_MT_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (; kk<_MT_N-1; kk++) {
            y = (ul_mt[kk]&_MT_UPPER_MASK)|(ul_mt[kk+1]&_MT_LOWER_MASK);
            ul_mt[kk] = ul_mt[kk+(_MT_M-_MT_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (ul_mt[_MT_N-1]&_MT_UPPER_MASK)|(ul_mt[0]&_MT_LOWER_MASK);
        ul_mt[_MT_N-1] = ul_mt[_MT_M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = ul_mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long mt::genrand_int31(void) {
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double mt::genrand_real1(void) {
    return genrand_int32()*(1.0/4294967295.0);
    /* divided by 2^32-1 */
}

/* generates a random number on [0,1)-real-interval */
double mt::genrand_real2(void) {
    return genrand_int32()*(1.0/4294967296.0);
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double mt::genrand_real3(void) {
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0);
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double mt::genrand_res53(void) {
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6;
    return(a*67108864.0+b)*(1.0/9007199254740992.0);
}

/*************************************/

matrix <double> mt::rand(unsigned int uint_N) {
    matrix <double> mat_rnd(uint_N,1);

    for (unsigned int uint_i=0; uint_i<uint_N; uint_i++) mat_rnd(uint_i) = genrand_real1();

    return mat_rnd;
}

matrix <double> mt::randn(unsigned int uint_N) {
    //float pi = 3.141592653589793;
    matrix <double> mat_rndn(uint_N,1);

    //for (unsigned int uint_i = 0; uint_i < uint_N; uint_i ++) {
    //  mat_rndn(uint_i) = std::sqrt(-2.0 * std::log(genrand_real1())) * std::cos ( 2.0 * pi * genrand_real1());
    //}

    // Box-Muller Transform
    double x1, x2, w, y1, y2;

    for (unsigned int i=0; i<uint_N; i+=2) {
        do {
            x1 = 2.0 * genrand_real1() - 1.0;
            x2 = 2.0 * genrand_real1() - 1.0;
            w = x1 * x1 + x2 * x2;
        } while ( w >= 1.0 );

        w = std::sqrt( (-2.0 * std::log( w ) ) / w );
        y1 = x1 * w;
        y2 = x2 * w;

        mat_rndn(i) = y1;
        if ((i+1) < uint_N) mat_rndn(i+1) = y2;
    }

    return mat_rndn;
}

unsigned int mt::rand_mask(unsigned int uint_mask) {
    return (genrand_int32() & uint_mask);
}
matrix <unsigned int> mt::rand_mask(unsigned int uint_mask, unsigned int uint_N) {
    matrix <unsigned int> mat_rnd(uint_N,1);

    for (unsigned int uint_i=0; uint_i<uint_N; uint_i++) mat_rnd(uint_i) = genrand_int32() & uint_mask;

    return mat_rnd;
}
} //NAMESPACE SUSA
