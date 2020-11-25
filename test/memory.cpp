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
 * along with Susa. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file memory.cpp
 * @brief Unit Test Suit
 * @author Behrooz Kamary
 */


#include <vector>
#include <list>
#include <set>

#ifndef SUSA_NDEBUG
#define SUSA_NDEBUG
#endif

#include "test.h"


template <template <typename T, typename AT> typename TYPE> void test()
{
    SUSA_LOG_DBG("[before] total allocated memory : " << susa::memory_tacker::instance().read() << " bytes.");
    TYPE <unsigned, susa::allocator_log<unsigned>> container;
    for (unsigned cnt = 0; cnt < 1e3; cnt++)
    {
        container.insert(std::end(container), cnt);
    }

    SUSA_LOG_DBG("[after] total allocated memory : " << susa::memory_tacker::instance().read() << " bytes.");
}

template <typename T, typename AT> using set_default_comparator = std::set<T, std::less<>, AT>;

int main()
{
    SUSA_TEST_EQ(susa::memory_tacker::instance().read(), 0, "initial memory tracker must be zero");

    {
    std::vector<int, susa::allocator_log<int> > int_vec(1024, 0, susa::allocator_log<int>());

    test <std::vector>();
    test <std::list>();
    test <set_default_comparator>();
    }

    SUSA_TEST_EQ(susa::memory_tacker::instance().read(), 0, "std types memory leak with susa allocator");


    {
    susa::array<int, susa::allocator_log<int> > int_array_a({2,3,4});
    susa::array<int, susa::allocator_log<int> > int_array_b({20,313,473,5});
    susa::array<int, susa::allocator_log<int> > int_array_c;
    int_array_b = int_array_a;
    susa::array<int, susa::allocator_log<int> > int_array_d(std::move(int_array_a));

    susa::array <int, susa::allocator_log<int> > arr_a({21,6,5,15,43});
    arr_a(2,4,3,0,1) = 55;
    arr_a(12,4,3,5,1) = 32;
    susa::array <int, susa::allocator_log<int> > arr_b = arr_a;
    }

    {
        susa::matrix <float, susa::allocator_log<float> > mat_m("[0 1 1; 2 3 2; 1 3 2; 4 2 2]");
        susa::matrix <float> result = susa::mean(mat_m);
    }

    SUSA_TEST_EQ(susa::memory_tacker::instance().read(), 0, "susa::array memory leak with susa allocator");

    SUSA_TEST_PRINT_STATS();

    return (uint_failed);
}