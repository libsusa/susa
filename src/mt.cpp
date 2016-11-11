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
 * along with Susa. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file mt.cpp
 * @brief The Mersenne Twister pseudorandom number generator class (definition).
 *
 * @author Behrooz Kamary Aliabadi
 * @version 1.0.0
 */

#include <susa.h>

namespace susa {

mt::mt()
{
    uint_mt = new uint32_t[_MT_N];
    mti = _MT_N + 1;
}

mt::mt(uint32_t uint_seed)
{
    uint_mt = new uint32_t[_MT_N];
    mti = _MT_N + 1;

    init_genrand(uint_seed);
}

mt::~mt()
{
    delete [] uint_mt;
}

/* initializes uint_mt[_MT_N] with a seed */
void mt::init_genrand(uint32_t uint_seed)
{
    uint_mt[0]= uint_seed & 0xffffffffUL;
    for (mti=1; mti<_MT_N; mti++) {
        uint_mt[mti] = (1812433253UL * (uint_mt[mti - 1] ^ (uint_mt[mti - 1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array uint_mt[].                     */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        uint_mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void mt::init_by_array(uint32_t init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1;
    j=0;
    k = (_MT_N>key_length ? _MT_N : key_length);
    for (; k; k--)
    {
        uint_mt[i] = (uint_mt[i] ^ ((uint_mt[i - 1] ^ (uint_mt[i - 1] >> 30)) * 1664525UL)) + init_key[j] + j; /* non linear */
        uint_mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        j++;
        if (i >= _MT_N)
        {
            uint_mt[0] = uint_mt[_MT_N - 1];
            i = 1;
        }
        if (j >= key_length) j = 0;
    }
    for (k = _MT_N - 1; k; k--)
    {
        uint_mt[i] = (uint_mt[i] ^ ((uint_mt[i - 1] ^ (uint_mt[i - 1] >> 30)) * 1566083941UL))
                   - i; /* non linear */
        uint_mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i >= _MT_N)
        {
            uint_mt[0] = uint_mt[_MT_N - 1];
            i = 1;
        }
    }

    uint_mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0,0xffffffff]-interval */
uint32_t mt::genrand_int32(void)
{
    uint32_t y;
    static uint32_t mag01[2] = {0x0UL, _MT_MATRIX_A};
    /* mag01[x] = x * _MT_MATRIX_A  for x=0,1 */

    if (mti >= _MT_N)
    { /* generate _MT_N words at one time */
        int kk;

        if (mti == _MT_N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk = 0; kk < _MT_N - _MT_M; kk++)
        {
            y = (uint_mt[kk]&_MT_UPPER_MASK)|(uint_mt[kk+1]&_MT_LOWER_MASK);
            uint_mt[kk] = uint_mt[kk+_MT_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (; kk<_MT_N-1; kk++)
        {
            y = (uint_mt[kk]&_MT_UPPER_MASK)|(uint_mt[kk+1]&_MT_LOWER_MASK);
            uint_mt[kk] = uint_mt[kk+(_MT_M-_MT_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (uint_mt[_MT_N-1]&_MT_UPPER_MASK)|(uint_mt[0]&_MT_LOWER_MASK);
        uint_mt[_MT_N - 1] = uint_mt[_MT_M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = uint_mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7)  & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
unsigned int mt::genrand_int31(void)
{
    return (unsigned int)(genrand_int32() >> 1);
}

/* generates a random number on [0,1]-real-interval */
double mt::genrand_real1(void)
{
    /* divided by 2^32-1 */
    return genrand_int32() * (1.0/4294967295.0);
}

/* generates a random number on [0,1)-real-interval */
double mt::genrand_real2(void)
{
    /* divided by 2^32 */
    return genrand_int32() * (1.0/4294967296.0);
}

/* generates a random number on (0,1)-real-interval */
double mt::genrand_real3(void)
{
    /* divided by 2^32 */
    return (((double)genrand_int32()) + 0.5)*(1.0 / 4294967296.0);
}

/* generates a random number on [0,1) with 53-bit resolution*/
double mt::genrand_res53(void)
{
    uint32_t a = genrand_int32() >> 5;
    uint32_t b = genrand_int32() >> 6;
    return (a * 67108864.0 + b) * (1.0/9007199254740992.0);
}

matrix <double> mt::rand(unsigned int uint_num)
{
    matrix <double> mat_ret(uint_num, 1);

    for (unsigned int uint_i = 0; uint_i < uint_num; uint_i++)
        mat_ret(uint_i) = genrand_real1();

    return mat_ret;
}

matrix <double> mt::randn(unsigned int uint_num)
{
    matrix <double> mat_ret(uint_num, 1);

    // Box-Muller Transform
    double x1, x2, w, y1, y2;

    for (unsigned int indx = 0; indx < uint_num; indx += 2)
    {
        do
        {
            x1 = 2.0 * genrand_real1() - 1.0;
            x2 = 2.0 * genrand_real1() - 1.0;
            w = x1 * x1 + x2 * x2;
        } while ( w >= 1.0 );

        w = std::sqrt( (-2.0 * std::log( w ) ) / w );
        y1 = x1 * w;
        y2 = x2 * w;

        mat_ret(indx) = y1;
        if ((indx + 1) < uint_num) mat_ret(indx + 1) = y2;
    }

    return mat_ret;
}

unsigned int mt::rand_mask(unsigned int uint_mask)
{
    return (genrand_int32() & uint_mask);
}

matrix <unsigned int> mt::rand_mask(unsigned int uint_mask, unsigned int uint_num)
{
    matrix <unsigned int> mat_ret(uint_num, 1);

    for (unsigned int uint_i = 0; uint_i < uint_num; uint_i++)
         mat_ret(uint_i) = genrand_int32() & uint_mask;

    return mat_ret;
}

unsigned int mt::bernoulli(double dbl_p)
{
  if (genrand_real3() <= dbl_p) return 1;
  return 0;
}

unsigned int mt::randint(unsigned int uint_max)
{
  return (genrand_int32() % uint_max + 1);
}

unsigned int mt::nonuniform(matrix <double> mat_probabilities)
{

  unsigned int           uint_prob_size       = mat_probabilities.size();
  matrix <unsigned int>  mat_sort_indices     = sort_indices(mat_probabilities);
  matrix <double>        mat_prob_sorted      = sort(mat_probabilities);
  double                 dbl_th               = 0;
  double                 dbl_rnd              = genrand_real3();

  dbl_th = 0;
  for ( unsigned int uint_indx = 0; uint_indx < uint_prob_size; uint_indx++ )
  {
    dbl_th += mat_prob_sorted(uint_indx);
    if ( dbl_th > dbl_rnd)
    {
      return mat_sort_indices(uint_indx);
    }
  }

  return 0;
}


unsigned int mt::nonuniform()
{

  unsigned int uint_prob_size = this->mat_probabilities.size();
  double       dbl_th         = 0;
  double       dbl_rnd        = genrand_real3();

  for ( unsigned int uint_indx = 0; uint_indx < uint_prob_size; uint_indx++ )
  {
    dbl_th += this->mat_prob_sorted(uint_indx);
    if ( dbl_th > dbl_rnd)
    {
      return this->mat_sort_indices(uint_indx);
    }
  }

  return 0;
}


matrix <unsigned int> mt::nonuniform(matrix <double> mat_probabilities, unsigned int uint_length)
{
  matrix <unsigned int> mat_rand(uint_length, 1);
  unsigned int uint_prob_size = mat_probabilities.size();
  double       dbl_th         = 0;
  double       dbl_rnd        = genrand_real3();

  for (unsigned int uint_mat_indx = 0; uint_mat_indx < uint_length; uint_mat_indx++)
  {
    for ( unsigned int uint_indx = 0; uint_indx < uint_prob_size; uint_indx++ )
    {
      dbl_th += mat_probabilities(uint_indx);
      if ( dbl_th > dbl_rnd)
      {
        mat_rand(uint_mat_indx) = uint_indx;
        break;
      }
    }
    dbl_th  = 0;
    dbl_rnd = genrand_real3();
  }

  return mat_rand;
}

void mt::set_probabilities(matrix <double> mat_probabilities)
{
  this->mat_probabilities = mat_probabilities;
  this->mat_sort_indices  = sort_indices(mat_probabilities);
  this->mat_prob_sorted   = sort(mat_probabilities);
}

} //NAMESPACE SUSA
