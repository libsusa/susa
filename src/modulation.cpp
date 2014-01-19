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
 * @file modulation.cpp
 * @brief This file contains methods that implement digital base-band communication schemes.
 * <i>Modulation</i> is the term that is being used in the communications theory.
 * @author Behrooz, Kamary Aliabadi
 * @version 1.0.0
 */


#include "../inc/susa.h"

namespace susa {

unsigned int qam::log2(unsigned int uint_x) {
  unsigned int r = 0;
  
  while( (uint_x >> r) != 0) {
    r++;
  }
  return (r-1);
}


qam::qam(unsigned int uint_m) {

  uint_K = log2(uint_m);
  double dP = sqrt(uint_m);
  unsigned int  uintL = (unsigned int)(uint_m/dP);

  // Constellation generation
  mat_s =  matrix < std::complex <double> > (uintL, uintL, std::complex <double>(0,0));
  mat_axis = matrix < std::complex <double> > ((unsigned int)dP, 1, std::complex <double>(0,0));
  
  for (unsigned int i=0; i<dP; i++) mat_axis(i) = std::complex <double> (i - (dP - 1)/2,i - (dP - 1)/2);
  
  for (unsigned int i=0; i<uintL; i++) {
    for (unsigned int j=0;j<uintL;j++) {
      mat_s(i,j) = std::complex <double> (mat_axis(i).real(),mat_axis(j).imag());
    }
  }

  // Energy calculation
  matrix <double> mat_mag = mag(mat_s);
  dEs = (sum(mat_mag) * (1.00/uint_m))(0);
  dEb = dEs/uint_K;
}


qam::~qam() {}

} // NAMESPACE SUSA
