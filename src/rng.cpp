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
 * @file rng.cpp
 * @brief Random Number Generator
 * @author Behrooz, Kamary Aliabadi
 * @version 1.0.0
 */



#include "../inc/susa.h"

namespace susa {

// Constructors and Destructor

rng::rng() {
  intIndex = -1;
}

rng::rng(unsigned int intSeed) {
  intIndex = -1;
  init(intSeed);
}

// Public methods
unsigned int rng::GetUInt() {

  intIndex++;

  if (intIndex >623) {
    generateNumbers();
    intIndex=0;
    }

  return extractNumber(intIndex);

}

unsigned int rng::rand_mask(unsigned int uint_mask) {
  intIndex++;

  if (intIndex >623) {
    generateNumbers();
    intIndex=0;
    }

  return (extractNumber(intIndex) & uint_mask);
}

matrix <unsigned int> rng::rand_mask(unsigned int uint_mask, unsigned int uint_N) {
  matrix <unsigned int> mat_ret(uint_N,1);

  for (unsigned int uint_i = 0; uint_i < uint_N; uint_i++) {
    intIndex++;

    if (intIndex >623) {
      generateNumbers();
      intIndex=0;
    }

    mat_ret(uint_i) = extractNumber(intIndex) & uint_mask;
  }

  return mat_ret;
}

double rng::GetDouble() {

  unsigned int UIntMax = (unsigned int)(0xFFFFFFFF);

  intIndex++;

  if (intIndex >623) {
    generateNumbers();
    intIndex=0;
    }

  return ((double)extractNumber(intIndex)/UIntMax);

}

double rng::rand() {

  unsigned int UIntMax = (unsigned int)(0xFFFFFFFF);

  intIndex++;

  if (intIndex >623) {
    generateNumbers();
    intIndex=0;
    }

  return ((double)extractNumber(intIndex)/UIntMax);

}

double rng::randn() { 
  // Returns a single Gaussian random number 
  // Box-Muller transform
  double x1, x2, w, y1, y2;

  do {
    x1 = 2.0 * rand() - 1.0;
    x2 = 2.0 * rand() - 1.0;
    w = x1 * x1 + x2 * x2;
  } while ( w >= 1.0 );

  w = std::sqrt( (-2.0 * std::log( w ) ) / w );
  y1 = x1 * w;
  y2 = x2 * w;

  return (y1);
}

matrix <double> rng::rand(unsigned int uint_N) { // Returns a single uniform random number
  matrix <double> mat_rnd(uint_N,1);
  unsigned int UIntMax = (unsigned int)(0xFFFFFFFF);

  for (unsigned int i=0;i<uint_N;i++) {
    intIndex++;

    if (intIndex >623) {
      generateNumbers();
      intIndex=0;
    }
    mat_rnd(i) = ((double)extractNumber(intIndex)/UIntMax);
  }

  return mat_rnd;
}

matrix <double> rng::randn(unsigned int uint_N) {
  // Box-Muller Transform
  double x1, x2, w, y1, y2;
  matrix <double> mat_rndn(uint_N,1);
  /*unsigned int UIntMax = (unsigned int)(0xFFFFFFFF);*/

  for (unsigned int i=0; i<uint_N; i+=2) {
    do{
      x1 = 2.0 * rand() - 1.0;
      x2 = 2.0 * rand() - 1.0;
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

void rng::init(unsigned int intSeed) {
  // Initialise the generator from a seed
  MT[0] = intSeed;
  for (int i=1;i<624;i++) {
    MT[i] = 1812433253 * (MT[i-1] ^(MT[i-1])>>30) + i;
  }
}

unsigned int rng::GetNonUniform( double* pr, unsigned int n ) {

  double  t = 0;
  double x = GetDouble();
  for ( unsigned int i = 0; i < n; i++ ) {
    t += pr[i];
    if ( t > x ) return i;
  }

  return 0; // Avoids GCC warning
}

unsigned int rng::nonUniform(std::vector <float> pr) {

  int n = pr.size();
  double  t = 0;
  double x = GetDouble();
  for ( int i = 0; i < n; i++ ) {
    t += pr[i];
    if ( t > x ) return i;
  }

  return 0; // Avoids GCC warning
}

// Private methods
void rng::generateNumbers() {
  // Generate an array of 624 untempered numbers
  for (int i=0;i<624;i++) {
    y = (MT[i]&0x80000000) + ((MT[(i+1)%624])&0x7FFFFFFF);
    if (y%2==0)  MT[i] = MT[(i + 397) % 624] ^ (y>>1);
    else  MT[i] = MT[(i + 397) % 624] ^ (y>>1) ^ 0x9908b0df; // (a) Parameter

  }
}


unsigned int rng::extractNumber(int i) {
  // Extract a tempered pseudorandom number based on the i-th value
  // generateNumbers() will have to be called again once the array of 624 numbers is exhausted
     y = MT[i];
     y = y ^ (y>>11);
     y = y ^ ((y<<7) & 0x9d2c5680);
     y = y ^ ((y<<15) ^ 0xefc60000);
     y = y ^ (y>>18);
     return y;
 }

}

