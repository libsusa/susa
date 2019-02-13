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
 * @file convolutional.cpp
 * @brief The convolutional ecoder/decoder class.
 * @author Behrooz Kamary Aliabadi
 * @version 1.0.0
 */


#include <susa.h>

namespace susa {

// Constructors and Destructors

ccode::ccode()
{
    uint_k = 1;
    uint_n = 0;
    uint_m = 0;
}

ccode::~ccode()
{
    if (uint_gen != NULL)
    {
        delete [] uint_gen;
        uint_gen = nullptr;
    }
}

ccode::ccode(uint32_t uint_n, uint32_t uint_k, uint32_t uint_m)
{
    SUSA_ASSERT_MESSAGE(uint_k == 1, "the number of inputs must be one.");
    this->uint_k     = uint_k;
    this->uint_n     = uint_n;
    this->uint_m     = uint_m;
    this->uint_mmask = susa::pow(2, uint_m) - 1;
    this->uint_gen   = new uint32_t[uint_n];
}


// Public domain methods

void ccode::set_generator(uint32_t uint_gen, uint32_t uint_gen_id)
{
    SUSA_ASSERT_MESSAGE(uint_gen_id < uint_n, "id exceeded the number of generators.");
    this->uint_gen[uint_gen_id] = OctToDec(uint_gen);
}

uint32_t ccode::next_state(uint32_t uint_state, bool b_input)
{
    if (b_input)
        return ((uint_state>>1) | (1<<(uint_m - 1)));
    else
        return ((uint_state>>1) & (~(1<<(uint_m - 1))));
}

uint32_t ccode::next_state(bool b_input)
{
    uint_last_state = uint_current_state;
    uint_current_state = next_state(uint_current_state, b_input);
    return uint_current_state;
}

std::vector <uint32_t> ccode::prev_states(uint32_t uint_state)
{
    std::vector <uint32_t> prev_states(2);
    prev_states[0]  = (uint_state<<1);
    prev_states[0] &= uint_mmask;
    prev_states[1]  = (uint_state<<1) + 1;
    prev_states[1] &= uint_mmask;

    return prev_states;
}

void ccode::zero_state()
{
    uint_current_state = 0;
    uint_last_state    = 0;
}

uint8_t ccode::next_output(uint32_t uint_state, bool b_input)
{
    uint8_t uintTmp = 0;
    if (b_input)
    {
        for (uint32_t j=0; j<uint_n; j++) uintTmp |=  ((count_1bits((uint_state | (1<<uint_m)) & uint_gen[j]) % 2) << j);
    }
    else
    {
        for (uint32_t j=0; j<uint_n; j++) uintTmp |=  ((count_1bits((uint_state & (~(1<<uint_m))) & uint_gen[j]) % 2) << j);
    }

    return uintTmp;
}

uint8_t ccode::next_output(bool b_input)
{
    uint_last_state = uint_current_state;
    uint8_t uintTmp = 0;
    if (b_input)
    {
        for (uint32_t j=0; j<uint_n; j++) uintTmp |=  ((count_1bits((uint_current_state | (1<<uint_m)) & uint_gen[j]) % 2) << j);
    }
    else
    {
        for (uint32_t j=0; j<uint_n; j++) uintTmp |=  ((count_1bits((uint_current_state & (~(1<<uint_m))) & uint_gen[j]) % 2) << j);
    }
    
    next_state(b_input);
    
    return uintTmp;
}


matrix <uint8_t> ccode::encode(const matrix <uint8_t>& mat_arg)
{
    uint32_t uint_rate      = uint_n/uint_k;
    uint32_t uint_in_size   = mat_arg.size();
    uint32_t uint_out_size  = uint_in_size * uint_rate;
    uint32_t uint_counter   = 0;

    matrix <uint8_t> mat_out(uint_out_size, 1);

    //init
    zero_state();

    for (uint32_t uinti = 0; uinti < uint_in_size; uinti++)
    {

        if (mat_arg(uinti) != 0)
        {
            for (uint32_t j = 0; j < uint_n; j++) mat_out(uint_counter + j) =  (count_1bits((uint_current_state | (1<<uint_m)) & uint_gen[j]) % 2);
            next_state(true);
        }
        else
        {
            for (uint32_t j = 0; j < uint_n; j++) mat_out(uint_counter + j) =  (count_1bits((uint_current_state & (~(1<<uint_m))) & uint_gen[j]) % 2);
            next_state(false);
        }
        uint_counter += uint_n;
    }

    return mat_out;
}

matrix <double> ccode::decode_bcjr(const matrix <double> &mat_arg, double dbl_ebn0)
{
    double a       = 1;
    double l_c     = 4 * a * dbl_ebn0;
    double c_k     = 0.5; // for equiprobable binary signal (bernolli process)
    double dbl_sum = 0;

    zero_state();

    uint32_t uint_num_stages = mat_arg.size() / uint_n;
    uint32_t uint_num_states = (1 << uint_m);

    // Gamma
    matrix <double> mat_gamma(uint_num_states, uint_num_states, 0);
    std::vector <matrix <double> > vec_gamma(uint_num_stages);

    // Alpha
    matrix <double> mat_alpha(uint_num_states, 1, 0);
    std::vector <matrix <double> > vec_alpha(uint_num_stages + 1, matrix <double> (uint_num_states,1));
    vec_alpha[0](0,0) = 1;

    // Beta
    matrix <double> mat_beta(uint_num_states, 1, 0);
    std::vector <matrix <double> > vec_beta(uint_num_stages + 1, matrix <double> (uint_num_states,1));
    vec_beta[uint_num_stages](0,0) = 1;


    // Forward calculation of Gamma/Alpha
    std::vector <uint32_t> vec_uint_prev_states;

    uint32_t uint_next_zero, uint_next_one;
    double dbl_arg = 0;

    for (uint32_t uint_stage = 0; uint_stage < uint_num_stages; uint_stage++)
    {
        // Gamma Calculation
        for (uint32_t uint_state = 0; uint_state < uint_num_states; uint_state++)
        {

            uint_next_zero  =  this->next_state(uint_state, false);
            uint_next_one   =  this->next_state(uint_state, true);

            dbl_arg = 0;
            for (uint32_t uint_i = 0; uint_i < uint_n; uint_i++)
            {
                dbl_arg += mat_arg(uint_n * uint_stage + uint_i) * (((this->next_output(uint_state,false) >> uint_i) & 0x1) == 0 ? -1 : 1);
            }

            mat_gamma(uint_state, uint_next_zero) = c_k * exp(0.5 * l_c * dbl_arg);

            dbl_arg = 0;
            for (uint32_t uint_i = 0; uint_i < uint_n; uint_i++)
            {
                dbl_arg += mat_arg(uint_n * uint_stage + uint_i) * (((this->next_output(uint_state,true) >> uint_i) & 0x1) == 0 ? -1 : 1);
            }

            mat_gamma(uint_state, uint_next_one) = c_k * exp(0.5 * l_c * dbl_arg);
        }

        vec_gamma[uint_stage] = mat_gamma;


        // Alpha Calculation
        for (uint32_t uint_state = 0; uint_state < uint_num_states; uint_state++)
        {
            vec_uint_prev_states = this->prev_states(uint_state);

            for (uint32_t uinti = 0; uinti < vec_uint_prev_states.size(); uinti++)
                mat_alpha(uint_state) += vec_alpha[uint_stage](vec_uint_prev_states[uinti]) * mat_gamma(vec_uint_prev_states[uinti], uint_state);

        }

        // Normalize the alpha
        dbl_sum = 0;
        for (uint32_t uinti = 0; uinti < mat_alpha.size(); uinti++) dbl_sum += mat_alpha(uinti);
        mat_alpha = mat_alpha / dbl_sum;

        vec_alpha[uint_stage + 1] = mat_alpha;
        mat_alpha.set_all(0);
    }


    // Beta Calculation
    for (uint32_t uint_stage = uint_num_stages; uint_stage > 0; uint_stage--)
    {
        for (uint32_t uint_state = 0; uint_state < uint_num_states; uint_state++)
        {
            mat_gamma = vec_gamma[uint_stage - 1];

            uint_next_zero =  this->next_state(uint_state, false);
            uint_next_one =  this->next_state(uint_state, true);

            mat_beta(uint_state) += vec_beta[uint_stage](uint_next_zero) * mat_gamma(uint_state, uint_next_zero);
            mat_beta(uint_state) += vec_beta[uint_stage](uint_next_one) * mat_gamma(uint_state, uint_next_one);

        }


        // Normalize the beta
        dbl_sum = 0;
        for (uint32_t uinti = 0; uinti < mat_beta.size(); uinti++) dbl_sum += mat_beta(uinti);
        mat_beta = mat_beta / dbl_sum;

        vec_beta[uint_stage - 1] = mat_beta;

        mat_beta.set_all(0);
    }


    matrix <double> mat_p_norm(uint_num_stages, 1, 0);
    matrix <double> mat_p_one(uint_num_stages, 1, 0);
    matrix <double> mat_p_zero(uint_num_stages, 1, 0);
    matrix <double> mat_lr(mat_arg.shape(), 0);

    for (uint32_t uint_stage = 0; uint_stage < uint_num_stages; uint_stage++)
    {

        mat_alpha   = vec_alpha[uint_stage];
        mat_beta    = vec_beta[uint_stage + 1];
        mat_gamma   = vec_gamma[uint_stage]; // Note ! code : Gamma_1 = Gamma(0)


        for (uint32_t uint_state = 0; uint_state < uint_num_states; uint_state++)
        {

            uint_next_zero =  this->next_state(uint_state,false);
            uint_next_one =  this->next_state(uint_state,true);

            mat_p_zero(uint_stage) += mat_alpha(uint_state) * mat_gamma(uint_state,uint_next_zero) * mat_beta(uint_next_zero);
            mat_p_one(uint_stage) += mat_alpha(uint_state) * mat_gamma(uint_state,uint_next_one) * mat_beta(uint_next_one);
        }

        mat_p_norm(uint_stage) = mat_p_one(uint_stage) + mat_p_zero(uint_stage);

        mat_p_one(uint_stage) = mat_p_one(uint_stage) / mat_p_norm(uint_stage);
        mat_p_zero(uint_stage) = mat_p_zero(uint_stage) / mat_p_norm(uint_stage);

        mat_lr(uint_stage) = mat_p_one(uint_stage) / mat_p_zero(uint_stage);
    }


    return mat_lr;
}



// Private methods

uint32_t ccode::OctToDec(uint32_t uintMask)
{

    uint32_t uintTmp = 0;

    for (uint32_t j=1; uintMask!=0; j*=8)
    {
        uintTmp += (uintMask%10)*j;
        uintMask = uintMask/10;
    }
    return uintTmp;
}

uint8_t ccode::count_1bits(uint32_t x)
{
    x = x - ((x >> 1) & 0x55555555);
    x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
    x = (x & 0x0F0F0F0F) + ((x >> 4) & 0x0F0F0F0F);
    x = (x & 0x00FF00FF) + ((x >> 8) & 0x00FF00FF);
    x = (x & 0x0000FFFF) + ((x >> 16) & 0x0000FFFF);
    return (x & 0x3F);
}

uint8_t ccode::count_0bits(uint32_t x)
{
    return (32 - count_1bits(x));
}

}      // NAMESPACE SUSA
