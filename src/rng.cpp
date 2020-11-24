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
 * @file rng.cpp
 * @brief Random Number Generator (definition).
 *
 * @author Behrooz Kamary
 */



#include <susa.h>

namespace susa {

// Constructors and Destructor
rng::rng()
{
    uint_index = 0;
}

rng::rng(unsigned int uint_seed)
{
    uint_index = 0;
    init(uint_seed);
}


unsigned int rng::rand_mask(unsigned int uint_mask)
{
    return (extract_number() & uint_mask);
}

matrix <unsigned int> rng::rand_mask(unsigned int uint_mask, unsigned int uint_num)
{
    matrix <unsigned int> mat_ret(uint_num, 1);

    for (unsigned int uint_i = 0; uint_i < uint_num; uint_i++)
    {
        mat_ret(uint_i) = extract_number() & uint_mask;
    }

    return mat_ret;
}

double rng::get_double()
{

    return ((double)extract_number() / std::numeric_limits<uint32_t>::max());

}

double rng::rand()
{

    return ((double)extract_number() / std::numeric_limits<uint32_t>::max());
}

double rng::randn()
{
    // Returns a single Gaussian random number
    // Box-Muller transform
    double x1, x2, w, y1;

    do
    {
        x1 = 2.0 * rand() - 1.0;
        x2 = 2.0 * rand() - 1.0;
        w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );

    w = std::sqrt( (-2.0 * std::log( w ) ) / w );
    y1 = x1 * w;

    return (y1);
}

matrix <double> rng::rand(unsigned int uint_num)
{
    matrix <double> mat_ret(uint_num, 1);

    for (unsigned int i = 0; i < uint_num; i++)
    {
        mat_ret(i) = (double)extract_number() / std::numeric_limits<uint32_t>::max();
    }

    return mat_ret;
}

matrix <double> rng::rand(unsigned int uint_rows, unsigned int uint_cols)
{
    matrix <double> mat_ret(uint_rows, uint_cols);
    unsigned int uint_size = mat_ret.size();

    for (unsigned int i = 0; i < uint_size; i++)
    {
        mat_ret(i) = (double)extract_number() / std::numeric_limits<uint32_t>::max();
    }

    return mat_ret;
}

matrix <double> rng::randn(unsigned int uint_num)
{
    // Box-Muller Transform
    double x1, x2, w, y1, y2;
    matrix <double> mat_ret(uint_num, 1);

    for (unsigned int indx = 0; indx < uint_num; indx += 2)
    {
        do
        {
            x1 = 2.0 * rand() - 1.0;
            x2 = 2.0 * rand() - 1.0;
            w  = x1 * x1 + x2 * x2;
        } while ( w >= 1.0 );

        w  = std::sqrt( (-2.0 * std::log( w ) ) / w );
        y1 = x1 * w;
        y2 = x2 * w;

        mat_ret(indx) = y1;
        if ((indx + 1) < uint_num) mat_ret(indx + 1) = y2;
    }

    return mat_ret;
}

void rng::init(unsigned int uint_seed)
{
    // Initialise the generator from a seed
    MT[0] = uint_seed;
    for (unsigned int indx = 1; indx < N; indx++)
    {
        MT[indx] = F * (MT[indx - 1] ^ (MT[indx - 1]) >> 30) + indx;
    }
}

size_t rng::get_nonuniform(double* prob, size_t n)
{

    double  t = 0;
    double  x = get_double();
    for (size_t indx = 0; indx < n; indx++)
    {
        t += prob[indx];
        if ( t > x ) return indx;
    }

    return 0;
}

size_t rng::get_nonuniform(const std::vector <float>& vec_prob)
{

    size_t n   = vec_prob.size();
    double t   = 0;
    double x   = get_double();

    for (size_t i = 0; i < n; i++)
    {
        t += vec_prob[i];
        if ( t > x ) return i;
    }

    return 0;
}

matrix <int8_t> rng::bernoulli(size_t size_num)
{
    // TODO: optimize the generation time
    matrix <int8_t> mat_ret(size_num, 1);

    for (unsigned int uint_i = 0; uint_i < size_num; uint_i++)
    {
        mat_ret(uint_i) = extract_number() & 0x01;
    }

    return mat_ret;
 }


// Private methods
void rng::generate_numbers()
{
    uint32_t x;
    for (uint32_t indx = 0; indx < N; indx++)
    {
        x = (MT[indx] & MASK_UPPER) + ((MT[(indx + 1) % N]) & MASK_LOWER);

        if ((x & 0x01) == 0)  MT[indx] = MT[(indx + M) % N] ^ (x >> 1);
        else  MT[indx] = MT[(indx + M) % N] ^ (x >> 1) ^ A;
    }
}


unsigned int rng::extract_number()
{
    if (uint_index >= N)
    {
        generate_numbers();
        uint_index = 0;
    }

    uint32_t x = MT[uint_index++];
    x = x ^ (x  >> U);
    x = x ^ ((x << S)  & B);
    x = x ^ ((x << T)  ^ C);
    x = x ^ (x  >> L);
    return x;
}

}
