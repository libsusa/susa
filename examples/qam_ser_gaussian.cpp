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

/*
 This program simulates a 16-QAM modem over AWGN channel using Susa.
 The output is a list of Eb/N0 vs. Symbol Error Rate (SER). 
 */



#include "susa.h"

using namespace std;
using namespace susa;

int main(void) {

    fstream fs_result;
    fs_result.open("result.txt",fstream::out);
    rng _rng(2987549);



    unsigned int uint_m = 16;         // Modulation order
    unsigned int uint_n = 1000000;    // Number of transmitted symbols
    double dbl_min_noise_db = -2;     // Minimum Eb/N0 in dB
    double dbl_max_noise_db = 10;     // Maximum Eb/N0 in dB
    unsigned int uint_num_steps = 12; // Number of simulation points

    matrix <unsigned int> mat_num_errors(uint_num_steps,1);

    qam _qam(uint_m);

    /* Uniform symbol generation */
    matrix < std::complex <double> > cmat_symbols =
        matrix < std::complex <double> > (uint_n, 1, std::complex <double>(0,0));
    for (unsigned int uint_i = 0; uint_i < uint_n; uint_i++) cmat_symbols(uint_i) = _qam.get_constellation()(_rng.rand_mask(0xF));

    /* AWGN channel generation */
    matrix < std::complex <double> > cmat_noise(uint_n,1);
    matrix < std::complex <double> > cmat_noisy_symbols;

    for (unsigned int uint_noise_step = 0; uint_noise_step < uint_num_steps; uint_noise_step++) {

        double dbl_noise_db = dbl_min_noise_db + uint_noise_step * (dbl_max_noise_db - dbl_min_noise_db)/(uint_num_steps - 1);
        double dbl_noise_dev = _qam.get_noise_deviation(dbl_noise_db);

        for (unsigned int uint_i = 0; uint_i < uint_n; uint_i++)
            cmat_noise(uint_i) = dbl_noise_dev * std::complex <double> (_rng.randn(), _rng.randn());

        /* The QAM symbols pass AWGN channel */
        cmat_noisy_symbols = cmat_symbols + cmat_noise;


        /* Demodulate the noisy symbols */
        for (unsigned int uint_i = 0; uint_i < uint_n; uint_i++)
            if (_qam.demodulate_symbol(cmat_noisy_symbols(uint_i)) != cmat_symbols(uint_i)) mat_num_errors(uint_noise_step)++;


        cout << "Eb/N0 = " << dbl_noise_db << " \t \t" << "SER = " << (double) mat_num_errors(uint_noise_step)/uint_n << endl;
        fs_result << dbl_noise_db << "   " <<  (double) mat_num_errors(uint_noise_step)/uint_n << endl;
    }

    fs_result.close();
    return 0;
}
