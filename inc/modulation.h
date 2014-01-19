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
 * @brief The matrix operation wrapper class for non-complex data sets
 * @author Behrooz, Kamary Aliabadi
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
template <class T> matrix <T> bpsk(const matrix <T> &mat_arg);

/**
 * @brief QAM Modulation
 * This is a wrapper class for QAM modulation and demodulation.
 * @ingroup Communications
*/
class qam {

  public :
  qam(unsigned int uint_m);
  ~qam();

  matrix < std::complex <double> > modulate_bits(matrix <char> mat_bits);

  
  matrix <char> demodulate_bits(matrix < std::complex <double> > mat_symbols);
  matrix <int> demodulate_symbols(matrix < std::complex <double> > mat_symbols);

  private :
  double dEs; // Energy per symbol
  double dEb; // Energy per bit
  unsigned int uint_K; // Number of bits per symbol

  matrix < std::complex <double> > mat_s;
  matrix < std::complex <double> > mat_axis;
  
  unsigned int log2(unsigned int uint_x);

};

template <class T> matrix <T> bpsk(const matrix <T> &mat_arg) {
  matrix <char> mat_ret(mat_arg.no_rows(), mat_arg.no_cols());
  unsigned int uint_length = mat_arg.no_cols() * mat_arg.no_rows();
  for (unsigned int uint_index = 0; uint_index < uint_length; uint_index++) mat_ret(uint_index) = (mat_arg(uint_index) == 1) ? 1 : -1;
  return mat_ret;
}



/*
std::vector < std::complex <double> > qam::modulate_bits(std::vector <char> vec_bits) {
  unsigned int uint_length = vec_bits.size();
  
  uint_length -= uint_length % uint_K; // Trim the bit stream if the size is not a factor of K

  unsigned int uint_symbol_length = uint_length / uint_K;

  std::vector < std::complex <double> >  vec_modulated(uint_symbol_length, std::complex <double> (0,0));

  return vec_modulated;
    
}

std::vector < std::complex <double> > qam::modulate_symbols(std::vector <short> vec_symbols) {
  //unsigned int uint_length = vec_symbols.size();
  

  std::vector < std::complex <double> >  vec_modulated(1, std::complex <double> (0,0));

  return vec_modulated;
    
}
*/
/*
 * std::vector <char> demodulate_symbols(std::vector <double> vec_symbols) {
  unsigned int uint_length = vec-symbols.size();

  for (unsigned int i=0; i < uint_length; i++) {
    cmat_t_symbols(i,0) = vec_s[];
  }
}*/


} // NAMESPACE SUSA
#endif // MODULATION_H
