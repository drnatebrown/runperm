// Unit tests for split_by_length_capping and split_by_balancing.
// These are simple assert-based tests, no external framework.

#include "internal/move/move_splitting.hpp"
#include "internal/ds/packed_vector.hpp"
#include "internal/permutation.hpp"

#include <cassert>
#include <iostream>
#include <vector>

using std::size_t;
using std::vector;

// Use the packed IntVector as IntVectorType for tests
using TestIntVector = IntVector;
using TestSplitResult = SplitResult<TestIntVector>;

void test_split_by_length_capping_no_splitting() {
    // Three equal-length runs that should not be split.
    const vector<ulint> lengths_vec = {4, 4, 4};
    const vector<ulint> perm_vec = {0, 4, 8};
    const vector<ulint> tau_inv_vec = {0, 1, 2};
    const ulint domain = 12; // Sum of lengths.
    const double factor = 1.0;

    TestIntVector lengths(lengths_vec);
    TestIntVector tau_inv(tau_inv_vec);

    TestSplitResult result;
    split_by_length_capping(lengths, tau_inv, domain, factor, result);

    assert(result.lengths.size() == lengths_vec.size());
    assert(result.tau_inv.size() == tau_inv_vec.size());

    for (size_t i = 0; i < lengths_vec.size(); ++i) {
        assert(result.lengths[i] == lengths_vec[i]);
        assert(result.tau_inv[i] == tau_inv_vec[i]);
    }

    // With these parameters, max_allowed_length is 4, so max_length should be 4.
    assert(result.max_length == 4);
}

void test_split_by_length_capping_with_splitting() {
    // One run is longer than the allowed maximum and should be split.
    const vector<ulint> lengths_vec = {2, 3, 1, 2, 2, 1, 1, 1, 3};
    const vector<ulint> perm_vec = {1, 9, 3, 12, 4, 14, 0, 15, 6};
    const vector<ulint> tau_inv_vec = {6, 0, 2, 4, 8, 1, 3, 5, 7};
    const ulint domain = 15; // Sum of lengths.
    size_t max_length = 2;

    TestIntVector lengths(lengths_vec);
    TestIntVector tau_inv(tau_inv_vec);

    TestSplitResult result;
    split_by_max_allowed_length(lengths, tau_inv, domain, max_length, result);

    const vector<ulint> expected_lengths = {2, 2, 1, 1, 2, 2, 1, 1, 1, 2, 1};
    const vector<ulint> expected_perm = {1, 9, 11, 3, 12, 4, 14, 0, 15, 6, 8};
    const vector<ulint> expected_tau_inv = {7, 0, 3, 5, 9, 10, 1, 2, 4, 6, 8};


    std::cout << "result.lengths: " << std::endl;
    for (size_t i = 0; i < result.lengths.size(); ++i) {
        std::cout << result.lengths[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "result.tau_inv: " << std::endl;
    for (size_t i = 0; i < result.tau_inv.size(); ++i) {
        std::cout << result.tau_inv[i] << " ";
    }
    std::cout << std::endl;
    
    assert(result.lengths.size() == expected_lengths.size());
    assert(result.tau_inv.size() == expected_tau_inv.size());

    for (size_t i = 0; i < expected_lengths.size(); ++i) {
        assert(result.lengths[i] == expected_lengths[i]);
        assert(result.tau_inv[i] == expected_tau_inv[i]);
    }
    assert(result.max_length == max_length);
}

void test_split_by_length_capping_mixed_lengths() {
    // Mixture of small and large runs; ensure max_length tracks the maximum.
    const vector<ulint> lengths_vec = {1, 16, 2};
    const vector<ulint> perm_vec    = {0, 1, 17};
    const vector<ulint> tau_inv_vec = {0, 1, 2};
    const ulint domain = 19;
    const double factor = 0.5;

    TestIntVector lengths(lengths_vec);
    TestIntVector tau_inv(tau_inv_vec);

    TestSplitResult result;
    split_by_length_capping(lengths, tau_inv, domain, factor, result);

    // avg_run_length = 19 / 3 ≈ 6.333...
    // ceil(avg_run_length * 0.5) = ceil(3.166...) = 4
    // bit_width(4) = 3, so max_allowed_length = 2^3 - 1 = 7
    //
    // So the run of length 16 is split into [7, 7, 2]:
    // [1, 16, 2] -> [1, 7, 7, 2, 2]
    const vector<ulint> expected_lengths = {1, 7, 7, 2, 2};
    const vector<ulint> expected_perm    = {0, 1, 8, 15, 17};
    const vector<ulint> expected_tau_inv = {0, 1, 2, 3, 4};

    assert(result.lengths.size() == expected_lengths.size());
    assert(result.tau_inv.size()   == expected_tau_inv.size());

    for (size_t i = 0; i < expected_lengths.size(); ++i) {
        assert(result.lengths[i] == expected_lengths[i]);
        assert(result.tau_inv[i] == expected_tau_inv[i]);
    }
    assert(result.max_length == 7);
}

// Dummy test for balancing: just ensure the function is callable.
void test_split_by_balancing_dummy() {
    const vector<ulint> lengths_vec = {3, 5, 2};
    const vector<ulint> perm_vec    = {0, 3, 8};
    const vector<ulint> tau_inv_vec = {0, 1, 2};
    const ulint domain = 10;
    const ulint factor = 4;

    TestIntVector lengths(lengths_vec);
    TestIntVector tau_inv(tau_inv_vec);

    TestSplitResult result;
    split_by_balancing(lengths, tau_inv, domain, factor, result);
}

int main() {
    test_split_by_length_capping_no_splitting();
    test_split_by_length_capping_with_splitting();
    test_split_by_length_capping_mixed_lengths();
    test_split_by_balancing_dummy();

    std::cout << "move_splitting tests passed" << std::endl;
    return 0;
}
