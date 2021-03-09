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

#include <cmath>
#include <iostream>
#include <complex>

#include <susa.h>

using namespace std;
using namespace susa;

int main(void)
{
	rng _rng(29009);
	rrcosine rrc(11, 0.65, 7);

	cout << "Root-Raised Cosine (RRC) : " << endl;
	cout << "\t Alpha = 65 %" << endl;
	cout << "\t Order = 7" << endl;
	cout << "\t Upsampling ratio = 8" << endl;
	cout << "\t Fs/Fd = L = 11" << endl;

	unsigned int uintM = 16;
	unsigned int uintN = 100000;

	qam mapper(uintM);

	/* Uniformly distributed symbol generation */
	cmatrix<double> cmat_t_symbols(uintN, 1);

	for (unsigned int i = 0; i < uintN; i++)
	{
		cmat_t_symbols(i, 0) = mapper.get_constellation()(_rng.rand_mask(0xF));
	}

	/*Root-Raised Cosine filtering on symbols*/
	cmatrix<double> cmat_t_rrc_symbols;
	cmat_t_rrc_symbols = filter(rrc.get(), susa::matrix<double>(1, 1, 1), upsample(cmat_t_symbols, 8));

	/* AWGN channel generation */
	cmatrix<double> cmat_noise(cmat_t_rrc_symbols.no_rows(), 1);
	double dDev;
	vector<double> dSER(15, 0);
	int int_cnt = 0;

	for (int i = -2; i < 11; i++)
	{
		dDev = mapper.get_noise_deviation(i);
		for (size_t j = 0; j < cmat_t_rrc_symbols.no_rows(); j++)
		{
			cmat_noise(j, 0) = std::complex<double>(dDev * _rng.randn(), dDev * _rng.randn());
		}
		cmatrix<double> cmat_noisy_symbols = cmat_t_rrc_symbols + cmat_noise;

		/*Receiver*/

		/*Root-Raised Cosine filtering on noisy symbols*/
		cmatrix<double> cmat_noisy_rrc_symbols = filter(rrc.get(), susa::matrix<double>(1, 1, 1), cmat_noisy_symbols);
		cmatrix<double> cmat_noisy_d_symbols(uintN, 1);

		for (size_t j = 0; j < uintN; j++)
		{
			cmat_noisy_d_symbols(j, 0) = cmat_noisy_rrc_symbols(7 + j * 8 - 1, 0);
		}

		/* Estimate Symbol Error Rate (SER) */
		for (size_t j = 0; j < uintN; j++)
			if (mapper.demodulate_symbol(cmat_noisy_d_symbols(j)) != cmat_t_symbols(j))
				dSER[int_cnt]++;

		cout << " i = " << int_cnt << "\t   N0/2 = " << dDev << "   -    SER = " << dSER[int_cnt] / uintN << "\t  t = " << (clock() / CLOCKS_PER_SEC) << " Sec" << endl;
		int_cnt++;
	}

	return 0;
}
