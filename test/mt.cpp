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
 * @file mt_test.cpp
 * @brief Mersenne Twister Random Number Generator Unit Tests
 * @author Test Suite
 */

#include "test.h"
#include <cmath>

int main(void)
{
    // ========== Constructor Tests ==========
    
    // Test default constructor
    susa::mt rng_default;
    SUSA_TEST_EQ(true, true, "mt default constructor");

    // Test constructor with seed
    susa::mt rng_seeded(12345);
    SUSA_TEST_EQ(true, true, "mt seeded constructor");

    // ========== Determinism Test ==========
    // Same seed should produce same sequence
    susa::mt rng_a(42);
    susa::mt rng_b(42);
    
    unsigned int val_a = rng_a.randint(1000);
    unsigned int val_b = rng_b.randint(1000);
    SUSA_TEST_EQ(val_a, val_b, "mt: same seed produces same sequence");

    // ========== randint Tests ==========
    susa::mt rng_int(54321);
    
    // Test randint range - max
    unsigned int rand_val_max = rng_int.randint(100);
    SUSA_TEST_EQ((rand_val_max < 100), true, "mt: randint produces value less than max");

    // Test randint produces various values
    susa::mt rng_int2(12341);
    unsigned int val1 = rng_int2.randint(1000000);
    unsigned int val2 = rng_int2.randint(1000000);
    unsigned int val3 = rng_int2.randint(1000000);
    SUSA_TEST_EQ((val1 != val2 || val2 != val3), true, "mt: randint produces different values");

    // ========== rand Tests ==========
    susa::mt rng_real(99999);
    
    // Test rand returns matrix of correct size
    susa::matrix<double> rand_vals = rng_real.rand(10);
    SUSA_TEST_EQ(rand_vals.size(), 10, "mt: rand produces correct size");

    // Test rand values are in [0, 1)
    bool all_in_range = true;
    for (size_t i = 0; i < rand_vals.size(); i++)
    {
        if (rand_vals(i) < 0.0 || rand_vals(i) >= 1.0)
        {
            all_in_range = false;
            break;
        }
    }
    SUSA_TEST_EQ(all_in_range, true, "mt: rand values in [0, 1) range");

    // ========== randn Tests ==========
    susa::mt rng_normal(555);
    
    // Test randn returns matrix of correct size
    susa::matrix<double> normal_vals = rng_normal.randn(50);
    SUSA_TEST_EQ(normal_vals.size(), 50, "mt: randn produces correct size");

    // Test randn values - should have some variation
    double min_val = normal_vals(0);
    double max_val = normal_vals(0);
    for (size_t i = 0; i < normal_vals.size(); i++)
    {
        if (normal_vals(i) < min_val) min_val = normal_vals(i);
        if (normal_vals(i) > max_val) max_val = normal_vals(i);
    }
    SUSA_TEST_EQ((max_val > min_val), true, "mt: randn produces varied values");

    // ========== rand_mask Tests ==========
    susa::mt rng_mask(777);
    
    // Test rand_mask with mask value
    unsigned int mask_val = rng_mask.rand_mask(0xFF); // 255 mask
    SUSA_TEST_EQ((mask_val <= 0xFF), true, "mt: rand_mask respects mask");

    // Test rand_mask matrix output
    susa::matrix<unsigned int> mask_vals = rng_mask.rand_mask(0x0F, 20); // 4-bit mask, 20 values
    SUSA_TEST_EQ(mask_vals.size(), 20, "mt: rand_mask matrix produces correct size");
    
    bool all_masked = true;
    for (size_t i = 0; i < mask_vals.size(); i++)
    {
        if (mask_vals(i) > 0x0F)
        {
            all_masked = false;
            break;
        }
    }
    SUSA_TEST_EQ(all_masked, true, "mt: rand_mask matrix values respect mask");

    // ========== bernoulli Tests ==========
    susa::mt rng_bernoulli(888);
    
    // Test bernoulli with p=0.5 - should produce mix of 0s and 1s
    unsigned int heads = 0;
    unsigned int num_trials = 200;
    for (unsigned int i = 0; i < num_trials; i++)
    {
        heads += rng_bernoulli.bernoulli(0.5);
    }
    // With p=0.5 and 200 trials, expect roughly 100, but allow 30-170 range
    SUSA_TEST_EQ((heads > 30 && heads < 170), true, "mt: bernoulli(0.5) produces reasonable distribution");

    // Test bernoulli with p=1.0 - should always return 1
    unsigned int ones = 0;
    for (unsigned int i = 0; i < 100; i++)
    {
        ones += rng_bernoulli.bernoulli(1.0);
    }
    SUSA_TEST_EQ(ones, 100, "mt: bernoulli(1.0) always returns 1");

    // Test bernoulli with p=0.0 - should always return 0
    susa::mt rng_bernoulli_zero(999);
    unsigned int zeros = 0;
    for (unsigned int i = 0; i < 100; i++)
    {
        if (rng_bernoulli_zero.bernoulli(0.0) == 0)
            zeros++;
    }
    SUSA_TEST_EQ(zeros, 100, "mt: bernoulli(0.0) always returns 0");

    // ========== Reproducibility Test ==========
    susa::mt rng_repro1(12321);
    susa::matrix<double> sequence1 = rng_repro1.rand(5);
    
    susa::mt rng_repro2(12321);
    susa::matrix<double> sequence2 = rng_repro2.rand(5);
    
    bool sequences_match = true;
    for (size_t i = 0; i < sequence1.size(); i++)
    {
        if (sequence1(i) != sequence2(i))
        {
            sequences_match = false;
            break;
        }
    }
    SUSA_TEST_EQ(sequences_match, true, "mt: reproducible with same seed");

    // ========== Different Seeds Test ==========
    susa::mt rng_diff1(11111);
    susa::mt rng_diff2(22222);
    
    susa::matrix<double> diff_seq1 = rng_diff1.rand(5);
    susa::matrix<double> diff_seq2 = rng_diff2.rand(5);
    
    bool sequences_differ = false;
    for (size_t i = 0; i < diff_seq1.size(); i++)
    {
        if (diff_seq1(i) != diff_seq2(i))
        {
            sequences_differ = true;
            break;
        }
    }
    SUSA_TEST_EQ(sequences_differ, true, "mt: different seeds produce different sequences");

    SUSA_TEST_PRINT_STATS();

    return (uint_failed);
}
