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
 * This program simulates a communication system with BPSK modulation
 * with a R=1/2 convolutional code G(1 + D^2, 1 + D + D^2) over an AWGN noisy channel,
 */

#include "susa.h"


#define _N 1000000

using namespace std;
using namespace susa;


int main(void) {

    // Random Number Generator (RNG)
    rng _rng(23098);

    // Encoder and State machine
    susa::ccode state(2,1,2);
    state.set_generator(7,0);
    state.set_generator(5,1);


    matrix <char> mat_bits = _rng.rand_mask(1, _N);
    matrix <char> mat_coded;

    mat_bits(0) = 0;
    mat_bits(1) = 0;
    mat_bits(_N - 1) = 0;
    mat_bits(_N - 2) = 0;

    state.build_trellis();
    mat_coded = state.encode(mat_bits);
    matrix <double> mod_bpsk = bpsk(mat_coded);


    matrix <unsigned int>   mat_err(10,1);
    matrix <double>         awgn_mod_bpsk;
    matrix <double>         awgn;
    matrix <double>         mat_l;
    matrix <double>         mat_l_ln;
    matrix <double>         mat_l_sign;
    double                  EbN0;


    for (int EbN0db = 0; EbN0db < 10; EbN0db++)
    {

        // AWGN channel
        EbN0 = std::pow(10,(double)EbN0db/10);
        awgn = sqrt(1 / EbN0) * _rng.randn(mat_coded.size());
        awgn_mod_bpsk = mod_bpsk + awgn;


        // BCJR decoder
        mat_l = state.decode_bcjr(awgn_mod_bpsk, EbN0);


        unsigned int uint_num_stages = awgn_mod_bpsk.size() / 2;

        mat_l_ln = susa::log(mat_l);
        mat_l_sign = susa::sign(mat_l_ln);

        for (unsigned int inti = 0; inti < uint_num_stages; inti++) if (mat_l_sign(inti) == -1) mat_l_sign(inti) = 0;

        for (unsigned int inti = 0; inti < uint_num_stages; inti++) if (mat_l_sign(inti) != mat_bits(inti)) mat_err(EbN0db)++;


        cout << "Eb/N0 = " << EbN0db << " \t \t" << "BER  =  " << ((double)mat_err(EbN0db)/_N) << endl;

    }


    fstream fs_result;
    fs_result.open("result.txt",fstream::out);
    fs_result.precision(5);

    double dbl_tmp;

    for (int inti = 0; inti < 10; inti++)
    {
        dbl_tmp = (double)mat_err(inti)/_N;
        fs_result << inti << " " << scientific << dbl_tmp << endl;
    }

    fs_result.close();

    return 0;
}
