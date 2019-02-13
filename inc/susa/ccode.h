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
 * @brief Convolutional Codes
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

    /**
     * @brief Constructor
     * 
     * A convolutional code is defined as (n,k,m) where is <i>n</i>
     * is the number of outputs, <i>k</i> is the number of inputs
     * and <i>m</i> is the memory size.
     *
     * @param uint_n number of outputs
     * @param uint_k number of inputs (has to be one)
     * @param uint_m number of memories
     */
    ccode(uint32_t uint_n, uint32_t uint_k, uint32_t uint_m);

    /**
     * @brief Constructor
     */
    ccode();

    //! Destructor
    ~ccode();


    /**
     * @brief Set the generator polynomail for each output
     *
     * @param uint_gen Generator in octal format
     * @param uint_gen_id Generator polynomial output identifier
     **/
    void set_generator(uint32_t uint_gen, uint32_t uint_gen_id);

    /**
     * @brief Set the internal memory directly through this method
     * @param uint_state The internal state to be set
     **/
    void set_internal_state(uint32_t uint_state)
    {
        uint_current_state = uint_state;
    }

    /**
     * @brief the rate of the convolutional code
     **/
    float rate()
    {
      return (uint_k / uint_n);
    }

    /**
     * @brief 1/n Convolutional encoder
     * @param mat_arg hard bit matrix to be encoded
     **/
    matrix <uint8_t> encode(const matrix <uint8_t>& mat_arg);

    /**
     * @brief BCJR decoder
     * @param mat_arg Input matrix to be decoded
     * @param dbl_ebn0 Eb/N0 in linear scale with AWGN assumption
     **/
    matrix <double> decode_bcjr(const matrix <double> &mat_arg, double dbl_ebn0);

  private:

    /**
     * @brief get next state
     *
     * @param uint_state the current state
     * @param b_input the input
     */
    uint32_t next_state(uint32_t uint_state, bool b_input);

    /**
     * @brief get next state using internal state
     *
     */
    uint32_t next_state(bool b_input);

    uint8_t next_output(uint32_t uint_state, bool b_input);

    uint8_t next_output(bool b_input);
    
    /**
     * @brief get all possible previous states
     *
     * @param uint_state the current state
     */
    std::vector <uint32_t> prev_states(uint32_t uint_state);

    uint32_t prev_output(uint32_t uint_state, bool b_input);

    uint32_t prev_output(bool b_input);
    
    uint32_t get_current_state()
    {
        return uint_current_state;
    }
    
    uint32_t get_last_state()
    {
        return uint_last_state;
    }
    
    uint8_t count_1bits(uint32_t x);

    uint8_t count_0bits(uint32_t x);

    void zero_state();

    uint32_t  uint_k;     // Number of inputs (currently implemented for A SINGLE INPUT ONLY)
    uint32_t  uint_n;     // Number of outputs
    uint32_t  uint_m;     // Number of memories
    uint32_t  uint_mmask; // Memory mask
    uint32_t* uint_gen;   // Dynamic array of generators

    uint32_t OctToDec(uint32_t);        // Converts Octal to Decimal

    uint32_t uint_current_state;
    uint32_t uint_last_state;

};

}
#endif // CONVOLUTIONAL_H
