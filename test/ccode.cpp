#include <susa.h>

#include "test.h"

using namespace std;
using namespace susa;


int main(void)
{

    {
        susa::ccode state(2, 1, 2);
        state.set_generator(7, 0);
        state.set_generator(5, 1);

        matrix<uint8_t> mat_bits("1 0 1 1 0");

        matrix<uint8_t> mat_expected("1 1 1 0 0 0 0 1 0 1");

        matrix<uint8_t> mat_coded = state.encode(mat_bits);

        SUSA_TEST_EQ(mat_expected, transpose(mat_coded), "Convolutional Code Encoder");
    }

    {
        susa::ccode state(2, 1, 2);
        state.set_generator(7, 0);
        state.set_generator(5, 1);

        matrix<uint8_t> mat_bits("1 0 1 0 0 0");

        matrix<uint8_t> mat_expected("1 1 1 0 0 0 1 0 1 1 0 0");

        matrix<uint8_t> mat_coded = state.encode(mat_bits);

        SUSA_TEST_EQ(mat_expected, transpose(mat_coded), "Convolutional Code Encoder");
    }

    {
        susa::ccode state(2, 1, 2);
        state.set_generator(7, 0);
        state.set_generator(5, 1);

        matrix<uint8_t> mat_bits("0 1 0 1 1 1 0 0 1 0 1 0 0 0 1 0 0");

        matrix<uint8_t> mat_expected("0 0 1 1 1 0 0 0 0 1 1 0 0 1 1 1 1 1 1 0 0 0 1 0 1 1 0 0 1 1 1 0 1 1");

        matrix<uint8_t> mat_coded = state.encode(mat_bits);

        SUSA_TEST_EQ(mat_expected, transpose(mat_coded), "Convolutional Code Encoder");
    }

    SUSA_TEST_PRINT_STATS();

    return (uint_failed);
}