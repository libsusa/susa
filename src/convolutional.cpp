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
 * @author Behrooz, Kamary Aliabadi
 * @version 1.0.0
 */


#include "../inc/susa.h"

namespace susa {
// Constructors and Destructors

convolutional_codec::convolutional_codec() {
    uint_k = 1;
    uint_n = 0;
    uint_m = 0;
}

convolutional_codec::~convolutional_codec() {
    if (uint_gen != NULL) {
        //std::cout << "\n del \n";
        delete [] uint_gen;
        uint_gen = 0;
    }
}

convolutional_codec::convolutional_codec(unsigned int uint_n, unsigned int uint_k, unsigned int uint_m) {
    this->uint_k = 1;
    this->uint_n = uint_n;
    this->uint_m = uint_m;
    this->uint_gen = new unsigned int[uint_n];
}


// Public domain methods

void convolutional_codec::set_generator(unsigned int uint_gen,unsigned int uint_gen_id) {
    this->uint_gen[uint_gen_id] = OctToDec(uint_gen);
}

void convolutional_codec::out_outputs() {

    unsigned int uintTmp=0;

    for (unsigned int i=0; i< (unsigned int)(1<<uint_m); i++) {

        // Input = 0
        for (unsigned int j=0; j<uint_n; j++) uintTmp |=  (EvenParity((i & (~(1<<uint_m)))&uint_gen[j]) << j);


        std::cout << UIntToBinChar(i) << " - 0 - " << UIntToBinChar(uintTmp) << std::endl;

        uintTmp = 0;

        // Input = 1
        for (unsigned int j=0; j<uint_n; j++) uintTmp |=  (EvenParity((i | (1<<uint_m))&uint_gen[j]) << j);


        std::cout << UIntToBinChar(i) << " - 1 - " << UIntToBinChar(uintTmp)<<std::endl<<std::endl;

        uintTmp = 0;
    }
}


void convolutional_codec::out_states() {

    unsigned int uint_next_state;

    for (unsigned int uinti = 0; uinti < (unsigned int)(1<<uint_m); uinti++) {

        // Input = 0
        uint_next_state = (uinti>>1) & (~(1<<(uint_m - 1)));

        std::cout << UIntToBinChar(uinti) << " - 0 - " << UIntToBinChar(uint_next_state)<<std::endl;

        // Input = 1
        uint_next_state = (uinti>>1) | (1<<(uint_m - 1));

        std::cout << UIntToBinChar(uinti) << " - 1 - " << UIntToBinChar(uint_next_state)<<std::endl<<std::endl;

    }

}

void convolutional_codec::build_trellis() {

    vec_uint_prev_states = std::vector <std::vector <unsigned int> > (1<<uint_m);
    vec_uint_next_states = std::vector <std::vector <unsigned int> > (1<<uint_m);

    unsigned int uint_next_state = 0;

    for (unsigned int uint_prev_state = 0; uint_prev_state < (unsigned int)(1<<uint_m); uint_prev_state++) {

        // Input = 0
        uint_next_state = (uint_prev_state>>1) & (~(1<<(uint_m - 1)));
        vec_uint_prev_states[uint_next_state].push_back(uint_prev_state);
        vec_uint_next_states[uint_prev_state].push_back(uint_next_state);

        // Input = 1
        uint_next_state = (uint_prev_state>>1) | (1<<(uint_m - 1));
        vec_uint_prev_states[uint_next_state].push_back(uint_prev_state);
        vec_uint_next_states[uint_prev_state].push_back(uint_next_state);

    }
}

unsigned int convolutional_codec::next_state(unsigned int uint_state, bool b_input) {

    if (b_input)
        return ((uint_state>>1) | (1<<(uint_m - 1)));
    else
        return ((uint_state>>1) & (~(1<<(uint_m - 1))));
}

unsigned int convolutional_codec::next_state(bool b_input) {
    uint_last_state = uint_current_state;
    uint_current_state = next_state(uint_current_state,b_input);
    return uint_current_state;
}


unsigned int convolutional_codec::next_output(unsigned int uint_state, bool b_input) {

    unsigned int uintTmp = 0;

    if (!b_input) {
        // Input = 0
        for (unsigned int j=0; j<uint_n; j++) uintTmp |=  (EvenParity((uint_state & (~(1<<uint_m)))&uint_gen[j]) << j);
        return uintTmp;
    } else {

        // Input = 1
        for (unsigned int j=0; j<uint_n; j++) uintTmp |=  (EvenParity((uint_state | (1<<uint_m))&uint_gen[j]) << j);
        return uintTmp;
    }
}

unsigned int convolutional_codec::next_output(bool b_input) {
    uint_last_state = uint_current_state;
    unsigned int uintTmp = 0;
    if (!b_input) {
        // Input = 0
        for (unsigned int j=0; j<uint_n; j++) uintTmp |=  (EvenParity((uint_current_state & (~(1<<uint_m)))&uint_gen[j]) << j);
        uint_current_state = next_state(uint_current_state,b_input);
        return uintTmp;
    } else {

        // Input = 1
        for (unsigned int j=0; j<uint_n; j++) uintTmp |=  (EvenParity((uint_current_state | (1<<uint_m))&uint_gen[j]) << j);

        uint_current_state = next_state(uint_current_state,b_input);
        return uintTmp;
    }

}

unsigned int convolutional_codec::prev_state(unsigned int uint_state,bool b_input) {
    return 0;
}

std::vector <unsigned int> convolutional_codec::prev_states(unsigned int uint_state) {
    return vec_uint_prev_states[uint_state];
}
unsigned int convolutional_codec::prev_internal_state() {
    return uint_last_state;
}

void convolutional_codec::zero_state() {
    uint_current_state = 0;
}

matrix <char> convolutional_codec::encode(const matrix <char> &mat_arg) {
    unsigned int uint_rate = uint_n/uint_k;
    unsigned int uint_out_size = mat_arg.no_rows() * mat_arg.no_cols() * uint_rate;
    unsigned int uint_counter = 0;

    matrix <char> mat_out(uint_out_size,1);

    //init
    zero_state();

    for (unsigned int uinti = 0; uinti < (mat_arg.no_rows() * mat_arg.no_cols()); uinti++) {

        if (mat_arg(uinti) == 1) {
            // Input = 1
            for (unsigned int j=0; j<uint_n; j++) mat_out(uint_counter + j) =  EvenParity((uint_current_state | (1<<uint_m)) & uint_gen[j]);
            uint_current_state = next_state(uint_current_state,true);

        } else {
            // Input = 0
            for (unsigned int j=0; j<uint_n; j++) mat_out(uint_counter + j) =  EvenParity((uint_current_state & (~(1<<uint_m))) & uint_gen[j]) ;
            uint_current_state = next_state(uint_current_state,false);
        }
        uint_counter += 2;
    }

    return mat_out;
}

matrix <double> convolutional_codec::decode_bcjr(const matrix <double> &mat_arg, double dbl_ebn0) {
    double a =1;
    double l_c = 4 * a * dbl_ebn0;
    double c_k = 0.5; // For equiprobable binary signal (bernolli process)
    double dbl_sum = 0;

    unsigned int uint_num_stages = mat_arg.size() / uint_n;
    unsigned int uint_num_states = (1 << uint_m);

    // Gamma
    matrix <double> mat_gamma(uint_num_states, uint_num_states);
    std::vector <matrix <double> > vec_gamma(uint_num_stages);

    // Alpha
    matrix <double> mat_alpha(uint_num_states, 1);
    std::vector <matrix <double> > vec_alpha(uint_num_stages + 1, matrix <double> (uint_num_states,1));
    vec_alpha[0](0,0) = 1;

    // Beta
    matrix <double> mat_beta(uint_num_states, 1);
    std::vector <matrix <double> > vec_beta(uint_num_stages + 1, matrix <double> (uint_num_states,1));
    vec_beta[uint_num_stages](0,0) = 1;


    // Forward calculation of Gamma/Alpha
    std::vector <unsigned int> vec_uint_prev_states;

    unsigned int uint_next_zero, uint_next_one;
    double dbl_arg = 0;

    for (unsigned int uint_stage = 0; uint_stage < uint_num_stages; uint_stage++) {

        // Gamma Calculation
        for (unsigned int uint_state = 0; uint_state < uint_num_states; uint_state++) {

            uint_next_zero =  this->next_state(uint_state,false);
            uint_next_one =  this->next_state(uint_state,true);

            dbl_arg = 0;
            for (unsigned int uint_i = 0; uint_i < uint_n; uint_i++) {
                dbl_arg += mat_arg(uint_n * uint_stage + uint_i) * (((this->next_output(uint_state,false) >> uint_i) & 0x1) == 0 ? -1 : 1);
            }

            mat_gamma(uint_state,uint_next_zero) = c_k * exp(0.5 * l_c * dbl_arg);

            dbl_arg = 0;
            for (unsigned int uint_i = 0; uint_i < uint_n; uint_i++) {
                dbl_arg += mat_arg(uint_n * uint_stage + uint_i) * (((this->next_output(uint_state,true) >> uint_i) & 0x1) == 0 ? -1 : 1);
            }

            mat_gamma(uint_state,uint_next_one) = c_k * exp(0.5 * l_c * dbl_arg);

        }

        vec_gamma[uint_stage] = mat_gamma;


        // Alpha Calculation
        for (unsigned int uint_state = 0; uint_state < uint_num_states; uint_state++) {

            vec_uint_prev_states = this->prev_states(uint_state);

            for (unsigned int uinti = 0; uinti < vec_uint_prev_states.size(); uinti++)
                mat_alpha(uint_state) += vec_alpha[uint_stage](vec_uint_prev_states[uinti]) * mat_gamma(vec_uint_prev_states[uinti], uint_state);

        }

        // Normalize the alpha
        dbl_sum = 0;
        for (unsigned int uinti = 0; uinti < mat_alpha.size(); uinti++) dbl_sum += mat_alpha(uinti);
        mat_alpha = mat_alpha / dbl_sum;

        vec_alpha[uint_stage + 1] = mat_alpha;
        mat_alpha.set_all(0);
    }


    // Beta Calculation
    for (unsigned int uint_stage = uint_num_stages; uint_stage > 0; uint_stage--) {

        for (unsigned int uint_state = 0; uint_state < uint_num_states; uint_state++) {

            mat_gamma = vec_gamma[uint_stage - 1];

            uint_next_zero =  this->next_state(uint_state,false);
            uint_next_one =  this->next_state(uint_state,true);

            mat_beta(uint_state) += vec_beta[uint_stage](uint_next_zero) * mat_gamma(uint_state, uint_next_zero);
            mat_beta(uint_state) += vec_beta[uint_stage](uint_next_one) * mat_gamma(uint_state, uint_next_one);

        }


        // Normalize the beta
        dbl_sum = 0;
        for (unsigned int uinti = 0; uinti < mat_beta.size(); uinti++) dbl_sum += mat_beta(uinti);
        mat_beta = mat_beta / dbl_sum;

        vec_beta[uint_stage - 1] = mat_beta;

        mat_beta.set_all(0);
    }


    matrix <double> mat_p_norm(uint_num_stages,1);
    matrix <double> mat_p_one(uint_num_stages,1);
    matrix <double> mat_p_zero(uint_num_stages,1);
    matrix <double> mat_lr(uint_num_stages,1);

    for (unsigned int uint_stage = 0; uint_stage < uint_num_stages; uint_stage++) {

        mat_alpha = vec_alpha[uint_stage];
        mat_beta = vec_beta[uint_stage + 1];
        mat_gamma = vec_gamma[uint_stage]; // Note ! code : Gamma_1 = Gamma(0)


        for (unsigned int uint_state = 0; uint_state < uint_num_states; uint_state++) {

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

unsigned int convolutional_codec::OctToDec(unsigned int uintMask) {

    unsigned int uintTmp = 0;

    for (unsigned int j=1; uintMask!=0; j*=8) {
        uintTmp += (uintMask%10)*j;
        uintMask = uintMask/10;
    }
    return uintTmp;
}

bool convolutional_codec::EvenParity(unsigned int uint_register) {
    bool boolTmp =false;
    for (int i=0; i<31; i++) {
        boolTmp ^= (bool)(uint_register%2);
        uint_register /= 2;
    }
    return boolTmp;
}

std::string convolutional_codec::UIntToBinChar(unsigned int uint_register) {
    std::string strBuffer;

    int intMax = (uint_m > uint_n) ? uint_m:uint_n;
    for (int i=0; i<intMax; i++) {
        strBuffer = (char)(0x30 + uint_register%2) + strBuffer;
        uint_register /= 2;
    }

    return strBuffer;
}

}      // NAMESPACE SUSA

