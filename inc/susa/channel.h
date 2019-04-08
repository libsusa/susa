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
 * @file channel.h
 * @brief Channel Equalizers for ISI channels
 * @author Behrooz Kamary
 * @version 1.0.0
 *
 * @defgroup Communications
 */

#ifndef CHANNEL_H
#define CHANNEL_H

namespace susa {

/**
 * @brief The Inter-Symbole Interference (ISI) class.
 *
 * This class contains a set of methods for signal generation
 * and equalization of ISI channels that currently supports only
 * BPSK modulation.
 *
 * @ingroup Communications
 */
template <class T> class channel
{
  private:

    void init(const matrix <T> &mat_taps, const matrix <T> &mat_pam);

    unsigned int    uint_num_taps;
    unsigned int    uint_num_pam;
    unsigned int    uint_num_states;
    unsigned int    uint_num_states_mem;

    matrix <T>      mat_pam;
    matrix <T>      mat_taps;

    matrix <T>      mat_main_taps;
    matrix <T>      mat_offset_taps;

    matrix <T>      mat_outputs;

    matrix <unsigned int> mat_source_states;
    matrix <unsigned int> mat_decoded_states;

    //! next state
    unsigned int next_state(unsigned int uint_state, unsigned int uint_pam_index);

    //! previous state
    matrix <unsigned int> prev_states(unsigned int uint_state);

    T next_output(unsigned int uint_state, unsigned int uint_pam_index);

    matrix <T> get_full_state(unsigned int uint_state);

    matrix <T> get_mem_state(unsigned int uint_state);

  public:

    /**
     * @brief Constructor
     *
     * Each symbol in the constellation shall be transduced to
     * a signal amplitude. The series of transduced symbols is
     * the signal that passes the channel (which is represented by the channel impulse response).
     *
     * @param mat_taps the Channel Impulse Response (CIR)
     * @param mat_pam the vector representing the corresponding signal for each constellation symbol
     */
    channel(const matrix <T> &mat_taps,const matrix <T> &mat_pam);

    //! Destructor
    ~channel();

    /**
     * @brief Inter-Symbol Interference (ISI) channel encoder
     *
     * The method delicately convolves the modulated input signal
     * with the Channel Impulse Response (CIR).
     *
     * @param mat_arg the channel input signal
     * @param uint_init_state the initial state
     */
    matrix <T> encode_isi(const matrix <T> &mat_arg, unsigned int uint_init_state);

    /**
     * @brief Bahl, Cocke, Jelinek and Raviv (BCJR) algorithm
     * is a maximum a posteriori probability (MAP) based algorithm for channel equalization.
     *
     * @param mat_arg the encoded signal vector
     * @param dbl_ebn0 the signal to noise ratio Eb/N0
     */
    matrix <T> decode_bcjr(const matrix <T> &mat_arg, T dbl_ebn0);

    /**
     * @brief Viterbi i.e. maximum likelihood sequence estimation (MLSE) algorithm
     * is a maximum likelihood (ML) based algorithm for channel equalization.
     *
     * @param mat_arg the encoded signal vector
     * @param uint_init_state the initial state of the trellis
     */
    matrix <T> decode_mlse(const matrix <T> &mat_arg, unsigned int uint_init_state);

    /**
     * @brief Viterbi algorithm
     * Viterbi is a Maximum Likelihood (ML) based algorithm for channel equalization
     * with an unknown initial state. Viterbi is also known as Maximum likelihood sequence estimation (MLSE).
     *
     * @param mat_arg the encoded signal vector
     */
    matrix <T> decode_mlse(const matrix <T> &mat_arg);


    /**
     * @brief ZF-DFE with BPSK quantizer
     *
     * Zero-Forcing (ZF) Decision–Feedback Equalization (DFE)
     * for BPSK modulated signals.
     *
     * @param mat_arg the input signal
     */
    matrix <T> decode_bpsk_dfe(const matrix <T> &mat_arg);

};



// Constructor

template <class T> channel <T>::channel(const matrix <T> &mat_taps, const matrix <T> &mat_pam)
{
    init(mat_taps, mat_pam);
}

// Destructor

template <class T> channel <T>::~channel() {}

// Private methods

template <class T> void channel <T>::init(const matrix <T> &mat_taps, const matrix <T> &mat_pam)
{
    uint_num_taps       = mat_taps.size();
    uint_num_pam        = mat_pam.size();
    uint_num_states     = susa::pow(uint_num_pam, uint_num_taps);
    uint_num_states_mem = susa::pow(uint_num_pam, uint_num_taps - 1);

    this->mat_pam       = mat_pam;
    this->mat_taps      = mat_taps;

    matrix <T> mat_trellis = matrix <T> (uint_num_taps, uint_num_states, 0.0f);

    unsigned int uint_tmp_state;

    for (unsigned int uint_state = 0; uint_state < uint_num_states; uint_state++)
    {
        uint_tmp_state = uint_state;
        for (unsigned int uint_tap = 0; uint_tap < uint_num_taps; uint_tap++)
        {
            mat_trellis(uint_tap, uint_state) = mat_pam(uint_tmp_state % uint_num_pam);
            uint_tmp_state /= uint_num_pam;
        }
    }

    mat_outputs = matmul(mat_taps, mat_trellis);
} // INIT



template <class T>  matrix <T> channel <T>::get_full_state(unsigned int uint_state)
{
    matrix <T> mat_ret(uint_num_taps, 1, 0);
    unsigned int uint_tmp_state = uint_state;
    for (unsigned int uint_tap = 0; uint_tap < uint_num_taps; uint_tap++)
    {
        mat_ret(uint_tap) = mat_pam(uint_tmp_state % uint_num_pam);
        uint_tmp_state /= uint_num_pam;
    }
    return mat_ret;
}

template <class T>  matrix <T> channel <T>::get_mem_state(unsigned int uint_state)
{
    matrix <T> mat_ret(uint_num_taps - 1, 1, 0);
    unsigned int uint_tmp_state = uint_state;
    for (unsigned int uint_tap = 1; uint_tap < uint_num_taps; uint_tap++)
    {
        mat_ret(uint_tap - 1) = mat_pam(uint_tmp_state % uint_num_pam);
        uint_tmp_state /= uint_num_pam;
    }
    return mat_ret;
}

template <class T> unsigned int channel <T>::next_state(unsigned int uint_state_mem_arg, unsigned int uint_pam_index)
{
    unsigned int uint_base = 1;
    for (unsigned int uint_i = 0; uint_i < (uint_num_taps - 1); uint_i++)
    {
        uint_base *= uint_num_pam;
    }

    unsigned int uint_state = (uint_state_mem_arg * uint_num_pam + uint_pam_index) % uint_base;

    return uint_state;
}

template <class T> T channel <T>::next_output(unsigned int uint_state_mem_arg, unsigned int uint_pam_index)
{
    return mat_outputs(uint_state_mem_arg * uint_num_pam + uint_pam_index);
}

template <class T> matrix <unsigned int> channel <T>::prev_states(unsigned int uint_state_mem_arg)
{

    matrix <unsigned int> mat_ret(uint_num_pam, 1, 0);
    unsigned int uint_shift = uint_state_mem_arg / uint_num_pam;

    unsigned int uint_base = 1;
    for (unsigned int uint_i = 0; uint_i < (uint_num_taps - 2); uint_i++) uint_base *= uint_num_pam;
    for (unsigned int uint_i = 0; uint_i < uint_num_pam; uint_i++) mat_ret(uint_i) = uint_shift + uint_i * uint_base;

    return mat_ret;

} // PREV_STATES


// Public methods

template <class T> matrix <T> channel <T>::encode_isi(const matrix <T> &mat_arg, unsigned int uint_init_state)
{
    // encode the signal stream
    // This convolve the signal with a CIR and start channel from a known state
    matrix <T> mat_ret;
    matrix <T> mat_init_state = flipud(get_mem_state(uint_init_state));

    mat_ret = conv(concat(concat(mat_init_state, mat_arg),mat_init_state), mat_taps).mid(mat_taps.size() - 1, mat_arg.size() + 2 * mat_taps.size() - 3);

    return mat_ret;
} // ENCODE


template <class T> matrix <T> channel <T>::decode_mlse(const matrix <T> &mat_arg, unsigned int uint_init_state) { // DECODE_MLSE Initial state

    T dbl_metric = 0;
    unsigned int uint_l = mat_taps.size() - 1;
    unsigned int uint_num_stages = mat_arg.size();
    unsigned int uint_num_states = uint_num_states_mem;


    matrix <T> mat_metric_past(uint_num_states, 1);
    matrix <T> mat_metric_next(uint_num_states, 1,std::numeric_limits<T>::max());


    matrix <unsigned int> mat_survivor_path(uint_num_states, uint_num_stages);
    matrix <unsigned int> mat_previous_state(uint_num_states, uint_num_stages);

    matrix <unsigned int> mat_visited(uint_num_states,1);
    matrix <unsigned int> mat_visited_forward(uint_num_states,1);
    mat_visited(uint_init_state) = 1;

    T dbl_next_output;
    unsigned int uint_next_state;

    for (unsigned int uint_stage = 0; uint_stage < uint_num_stages; uint_stage++)
    {
        for (unsigned int uint_state = 0; uint_state < uint_num_states; uint_state++)
        {
            if (mat_visited(uint_state) == 1)
            {
                for (unsigned int uint_pam_index = 0; uint_pam_index < uint_num_pam; uint_pam_index++)
                {

                    dbl_next_output = this->next_output(uint_state,uint_pam_index);
                    uint_next_state = this->next_state(uint_state,uint_pam_index);

                    mat_visited_forward(uint_next_state) = 1;

                    dbl_metric = (dbl_next_output - mat_arg(uint_stage)) * (dbl_next_output - mat_arg(uint_stage));


                    if ((dbl_metric + mat_metric_past(uint_state)) < mat_metric_next(uint_next_state))
                    {
                        mat_metric_next(uint_next_state) = dbl_metric + mat_metric_past(uint_state);
                        mat_survivor_path(uint_next_state,uint_stage) = uint_pam_index;
                        mat_previous_state(uint_next_state,uint_stage) = uint_state;
                    }
                }
            }
        }
        mat_metric_past = mat_metric_next;
        mat_metric_next = matrix <T> (uint_num_states, 1,std::numeric_limits<T>::max());
        mat_visited = mat_visited_forward;
        mat_visited_forward = matrix <unsigned int> (uint_num_states,1);
    }



    // Trace Back the Trellis (decoding)
    unsigned int uint_curr_state = uint_init_state;
    unsigned int uint_curr_input;
    matrix <T> mat_ret(uint_num_stages,1);


    for (unsigned int uint_stage = uint_num_stages; uint_stage > 0; uint_stage--)
    {
        uint_curr_input = mat_survivor_path(uint_curr_state,uint_stage - 1);
        mat_ret(uint_stage - 1) = mat_pam(uint_curr_input);
        uint_curr_state = mat_previous_state(uint_curr_state,uint_stage - 1);
    }


    mat_ret = mat_ret.left(uint_num_stages - uint_l);

    return mat_ret;
}// DECODE_MLSE Initial state

template <class T> matrix <T> channel <T>::decode_mlse(const matrix <T> &mat_arg)
{ // DECODE_MLSE

    // NOTES
    // This decoder must be used with CONV as the channel encoder.
    // It detects the initial and final states.

    T dbl_metric = 0;
    unsigned int uint_l = mat_taps.size() - 1;
    unsigned int uint_num_stages = mat_arg.size() - 2 * uint_l;
    unsigned int uint_num_states = uint_num_states_mem;

    matrix <T> mat_new_arg = mat_arg.mid(uint_l, mat_arg.size() - uint_l - 1);



    ////matrix <T> mat_metric(uint_num_states, uint_num_stages + 1,std::numeric_limits<T>::max());
    ////for (unsigned int uint_i = 0; uint_i < uint_num_states; uint_i++) mat_metric(uint_i, 0) = 0;
    matrix <T> mat_metric_past(uint_num_states, 1);
    matrix <T> mat_metric_next(uint_num_states, 1,std::numeric_limits<T>::max());

    matrix <unsigned int> mat_survivor_path(uint_num_states, uint_num_stages);
    matrix <unsigned int> mat_previous_state(uint_num_states, uint_num_stages);


    // First

    matrix <T> mat_first(uint_num_states,1);

    for (unsigned int uint_i = 0; uint_i < uint_num_states; uint_i++)
    {
        mat_first(uint_i) = norm(conv(flipud(get_mem_state(uint_i)),mat_taps).left(uint_l) - mat_arg.left(uint_l))(0);
        mat_first(uint_i) *= mat_first(uint_i);
        mat_metric_past(uint_i) = mat_first(uint_i)/uint_l;
    }


    T dbl_next_output;
    unsigned int uint_next_state;

    for (unsigned int uint_stage = 0; uint_stage < uint_num_stages; uint_stage++) {
        for (unsigned int uint_state = 0; uint_state < uint_num_states; uint_state++) {
            for (unsigned int uint_pam_index = 0; uint_pam_index < uint_num_pam; uint_pam_index++) {

                dbl_next_output = this->next_output(uint_state,uint_pam_index);
                uint_next_state = this->next_state(uint_state,uint_pam_index);

                dbl_metric = (dbl_next_output - mat_new_arg(uint_stage)) * (dbl_next_output - mat_new_arg(uint_stage));

                if ((dbl_metric + mat_metric_past(uint_state)) < mat_metric_next(uint_next_state)) {
                    mat_metric_next(uint_next_state) = dbl_metric + mat_metric_past(uint_state);
                    mat_survivor_path(uint_next_state,uint_stage) = uint_pam_index;
                    mat_previous_state(uint_next_state,uint_stage) = uint_state;
                }
            }
        }
        mat_metric_past = mat_metric_next;
        mat_metric_next = matrix <T> (uint_num_states, 1,std::numeric_limits<T>::max());
    }

    // Last
    matrix <T> mat_last(uint_num_states,1);

    for (unsigned int uint_i = 0; uint_i < uint_num_states; uint_i++)
    {
        mat_last(uint_i) = norm(conv(flipud(get_mem_state(uint_i)),mat_taps).right(uint_l) - mat_arg.right(uint_l))(0);
        mat_last(uint_i) *= mat_last(uint_i);
        mat_metric_next(uint_i) += mat_last(uint_i)/uint_l;
    }


    // To find the state which has the minimum metric at the very last stage

    unsigned int uint_min_state = min(mat_metric_past)(0);


    // Trace Back the Trellis (decoding)
    unsigned int uint_curr_state = uint_min_state;
    unsigned int uint_curr_input;
    matrix <T> mat_ret(uint_num_stages,1);


    for (unsigned int uint_stage = uint_num_stages; uint_stage > 0; uint_stage--)
    {
        uint_curr_input = mat_survivor_path(uint_curr_state,uint_stage - 1);
        mat_ret(uint_stage - 1) = mat_pam(uint_curr_input);
        uint_curr_state = mat_previous_state(uint_curr_state,uint_stage - 1);
    }
    mat_ret = concat(flipud(get_mem_state(min(mat_first)(0))),mat_ret);

    return mat_ret;
} // DECODE_MLSE

template <class T> matrix <T> channel <T>::decode_bcjr(const matrix <T> &mat_arg, T dbl_ebn0)
{
// The data sequence must begin/end to zero state.

    T dbl_sum = 0;

    unsigned int uint_num_stages = mat_arg.size();
    unsigned int uint_num_states = uint_num_states_mem;

    // Gamma
    matrix <T> mat_gamma(uint_num_states, uint_num_states);
    std::vector <matrix <T> > vec_gamma(uint_num_stages, matrix <T> (uint_num_states, uint_num_states));


    // Alpha
    matrix <T> mat_alpha(uint_num_states, 1, 0);
    std::vector <matrix <T> > vec_alpha(uint_num_stages + 1, matrix <T> (uint_num_states, 1, 0));
    vec_alpha[0](0,0) = 1;

    // Beta
    matrix <T> mat_beta(uint_num_states, 1, 0);
    std::vector <matrix <T> > vec_beta(uint_num_stages + 1, matrix <T> (uint_num_states, 1, 0));
    vec_beta[uint_num_stages](0,0) = 1;


    // Forward calculation of Gamma/Alpha

    matrix <unsigned int> mat_uint_prev_states;

    unsigned int uint_next_zero, uint_next_one;

    for (unsigned int uint_stage = 0; uint_stage < uint_num_stages; uint_stage++)
    {
        // Gamma Calculation

        for (unsigned int uint_state = 0; uint_state < uint_num_states; uint_state++)
        {

            uint_next_zero  =  this->next_state(uint_state, 0);
            uint_next_one   =  this->next_state(uint_state, 1);

            mat_gamma(uint_state,uint_next_zero) = exp((-dbl_ebn0) * std::pow(mat_arg(uint_stage) - this->next_output(uint_state, 0),2));

            mat_gamma(uint_state,uint_next_one) = exp((-dbl_ebn0) * std::pow(mat_arg(uint_stage) - this->next_output(uint_state, 1),2));

        }
        vec_gamma[uint_stage] = mat_gamma;



        // Alpha Calculation


        for (unsigned int uint_state = 0; uint_state < uint_num_states; uint_state++)
        {

            mat_uint_prev_states = this->prev_states(uint_state);

            for (unsigned int uinti = 0; uinti < mat_uint_prev_states.size(); uinti++)
                mat_alpha(uint_state) += vec_alpha[uint_stage](mat_uint_prev_states(uinti)) * mat_gamma(mat_uint_prev_states(uinti), uint_state);

        }

        // Normalize the alpha
        dbl_sum = 0;
        for (unsigned int uinti = 0; uinti < mat_alpha.size(); uinti++) dbl_sum += mat_alpha(uinti);
        mat_alpha = mat_alpha / dbl_sum;

        vec_alpha[uint_stage + 1] = mat_alpha;
        mat_alpha.set_all(0);
    }


    // Beta Calculation

    for (unsigned int uint_stage = uint_num_stages; uint_stage > 0; uint_stage--)
    {

        for (unsigned int uint_state = 0; uint_state < uint_num_states; uint_state++)
        {

            mat_gamma = vec_gamma[uint_stage - 1];

            uint_next_zero  =  this->next_state(uint_state, 0);
            uint_next_one   =  this->next_state(uint_state, 1);

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


    matrix <T> mat_p_norm(uint_num_stages, 1, 0);
    matrix <T> mat_p_one(uint_num_stages, 1, 0);
    matrix <T> mat_p_zero(uint_num_stages, 1, 0);

    // likelihood ratios
    matrix <T> mat_lr(mat_arg.shape(), 0);

    for (unsigned int uint_stage = 0; uint_stage < uint_num_stages; uint_stage++)
    {

        mat_alpha   = vec_alpha[uint_stage];
        mat_beta    = vec_beta[uint_stage + 1];
        mat_gamma   = vec_gamma[uint_stage]; // Note ! code : Gamma_1 = Gamma(0)


        for (unsigned int uint_state = 0; uint_state < uint_num_states; uint_state++)
        {

            uint_next_zero          =  this->next_state(uint_state,0);
            uint_next_one           =  this->next_state(uint_state,1);

            mat_p_zero(uint_stage)  += mat_alpha(uint_state) * mat_gamma(uint_state,uint_next_zero) * mat_beta(uint_next_zero);
            mat_p_one(uint_stage)   += mat_alpha(uint_state) * mat_gamma(uint_state,uint_next_one) * mat_beta(uint_next_one);
        }

        mat_p_norm(uint_stage)      = mat_p_one(uint_stage) + mat_p_zero(uint_stage);

        mat_p_one(uint_stage)       = mat_p_one(uint_stage) / mat_p_norm(uint_stage);
        mat_p_zero(uint_stage)      = mat_p_zero(uint_stage) / mat_p_norm(uint_stage);

        mat_lr(uint_stage)          = mat_p_one(uint_stage) / mat_p_zero(uint_stage);
    }


    return mat_lr;
} // DECODE_BCJR

} // NAMESPACE SUSA
#endif // CHANNEL_H
