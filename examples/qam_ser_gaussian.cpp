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


#include "susa.h"

using namespace std;
using namespace susa;

int main(void){
	
	fstream fs_matlab;
	fs_matlab.open("result.txt",fstream::out);
	rng _rng(2987549);



	unsigned int uintM = 16;
	unsigned int uintN = 1000000;
  double dbl_min_noise_db = -2;
  double dbl_max_noise_db = 10;
  unsigned int uint_num_steps = 12;

  matrix <unsigned int> mat_num_errors(uint_num_steps,1);

	qam _qam(uintM);

	/* Uniform symbol generation */
	matrix < std::complex <double> > cmat_symbols =
        matrix < std::complex <double> > (uintN, 1, std::complex <double>(0,0));
  for (unsigned int uint_i = 0; uint_i < uintN; uint_i++) cmat_symbols(uint_i) = _qam.get_constellation()(_rng.rand_mask(0xF));
  
  /* AWGN channel generation */
  matrix < std::complex <double> > cmat_noise(uintN,1);
  matrix < std::complex <double> > cmat_noisy_symbols;

  for (unsigned int uint_noise_step = 0; uint_noise_step < uint_num_steps; uint_noise_step++) {

    double dbl_noise_db = dbl_min_noise_db + uint_noise_step * (dbl_max_noise_db - dbl_min_noise_db)/(uint_num_steps - 1);
    double dbl_noise_dev = _qam.get_noise_deviation(dbl_noise_db);

    for (unsigned int uint_i = 0; uint_i < uintN; uint_i++)
      cmat_noise(uint_i) = dbl_noise_dev * std::complex <double> (_rng.randn(), _rng.randn());

    /* The QAM symbols pass AWGN channel */
    cmat_noisy_symbols = cmat_symbols + cmat_noise;


    /* Demodulate the noisy symbols */
    for (unsigned int uint_i = 0; uint_i < uintN; uint_i++)
      if (_qam.demodulate_symbol(cmat_noisy_symbols(uint_i)) != cmat_symbols(uint_i)) mat_num_errors(uint_noise_step)++;


    cout << "Noise : " << dbl_noise_db << " \t \t" << "SER : " << (double) mat_num_errors(uint_noise_step)/uintN << endl;
    fs_matlab << dbl_noise_db << "   " <<  (double) mat_num_errors(uint_noise_step)/uintN << endl;
  }

	fs_matlab.close();
	return 0;
}
