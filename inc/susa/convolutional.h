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
 * @file convolutional.h
 * @brief The convolutional encoder/decoder
 * @author Behrooz Kamary Aliabadi
 * @version 1.0.0
 */

#ifndef CONVOLUTIONAL_H
#define CONVOLUTIONAL_H


namespace susa {

/**
 * @brief Convolutional encoder/decoder.
 *
 * This class implements encoding and decoding of the convolutional error correction codes.
 * The constructor sets up the internal wiring of the state-machine. This is defined through
 * the constructor arguments which are a set of numbers in octal format.
 * To know more about it refer to <i>Fundamentals of convolutional coding, Rolf Johannesson and Kamil Zigangirov, IEEE Press, 1999</i>.
 * Note that input/output of the class methods are all in binary format (boolean to represent a single bit).
 *
 * @ingroup Communications
 */
class ccode
    {
  public:
    ccode(unsigned int uint_n, unsigned int uint_k, unsigned int uint_m);

    ccode();

    ~ccode();


    /**
     * @brief Set the generator polynomail for each output
     *
     * @param uintGen Generator in octal format
     * @param uintGenID  Generator polynomial output identifier
     **/
    void set_generator(unsigned int uint_gen, unsigned int uint_gen_id);  // Set generators

    /**
     * @brief Set the internal memory directly through this method
     * @param uint_state The internal state to be set
     **/
    void set_internal_state(unsigned int uint_state) {
        uint_current_state = uint_state;   // Set Internal State
    }

    /**
     * @brief next state
     *
     * param uint_state the current state
     * param b_input the input
     */
    unsigned int next_state(unsigned int uint_state, bool b_input);
    unsigned int next_state(bool b_input);                           // Next state using internal state
    unsigned int next_output(unsigned int uint_state, bool b_input);
    unsigned int next_output(bool b_input);

    unsigned int prev_state(unsigned int uint_state, bool b_input);  // Previous state (Current state,input)
    unsigned int prev_state(bool b_input);                           // Previous state using internal state
    std::vector <unsigned int> prev_states(unsigned int uint_state); // Previous states
    unsigned int prev_internal_state();                              // Previous state using internal state
    unsigned int prev_output(unsigned int uint_state, bool b_input);
    unsigned int prev_output(bool b_input);

    unsigned int get_current_state()
    {
        return uint_current_state;
    }

    unsigned int get_last_state()
    {
        return uint_last_state;
    }

    void zero_state();         // Zero internal state
    void out_states();         // Print out states
    void out_outputs();        // Print out outputs

    /**
     * @brief MUST BE CALLED after initialization of the generator polynomials are set
     *
     * This methos builds Previous and Next vectors
     **/
    void build_trellis();

    /**
     * @brief Returns the rate of the convolutional code
     **/
    double rate();

    /**
     * @brief 1/n Convolutional encoder
     * @param mat_arg hard bit matrix to be encoded
     **/
    matrix <char> encode(const matrix <char> &mat_arg);

    /**
     * @brief BCJR decoder
     * @param mat_arg Input matrix to be decoded
     * @param dbl_ebn0 Eb/N0 in linear scale with AWGN assumption
     **/
    matrix <double> decode_bcjr(const matrix <double> &mat_arg, double dbl_ebn0);

  private:
    unsigned int uint_k;     // Number of inputs (currently implemented for ONE INPUT ONLY)
    unsigned int uint_n;     // Number of outputs
    unsigned int uint_m;     // Number of memories
    unsigned int* uint_gen;  // Dynamic array of generators

    std::vector <std::vector <unsigned int> > vec_uint_prev_states;
    std::vector <std::vector <unsigned int> > vec_uint_next_states;

    std::vector <std::vector <unsigned int> > vec_uint_outputs;

    unsigned int OctToDec(unsigned int);        // Converts Octal to Decimal
    std::string UIntToBinChar(unsigned int);    // Converts unsigned integers to string
    bool EvenParity(unsigned int);              // Checks even parity of unsigned integers

    unsigned int uint_current_state;
    unsigned int uint_last_state;
};

}
#endif // CONVOLUTIONAL_H
