// Unit tests for split_by_length_capping and split_by_balancing.
// These are simple assert-based tests, no external framework.

#include "internal/move/move_splitting.hpp"

#include <cassert>
#include <iostream>
#include <vector>

using std::size_t;
using std::vector;

void test_split_by_length_capping_no_splitting() {
    // Three equal-length runs that should not be split.
    const vector<ulint> lengths = {4, 4, 4};
    const vector<ulint> perm = {0, 4, 8};
    const ulint domain = 12; // Sum of lengths.
    const double factor = 1.0;

    SplitResult result;
    split_by_length_capping(lengths, perm, domain, factor, result);

    assert(result.lengths.size() == lengths.size());
    assert(result.interval_permutations.size() == perm.size());

    for (size_t i = 0; i < lengths.size(); ++i) {
        assert(result.lengths[i] == lengths[i]);
        assert(result.interval_permutations[i] == perm[i]);
    }

    // With these parameters, max_allowed_length is 4, so max_length should be 4.
    assert(result.max_length == 4);
}

void test_split_by_length_capping_with_splitting() {
    // One run is longer than the allowed maximum and should be split.
    const vector<ulint> lengths = {2, 10, 3};
    const vector<ulint> perm = {0, 2, 12};
    const ulint domain = 15; // Sum of lengths.
    const double factor = 1.0;

    SplitResult result;
    split_by_length_capping(lengths, perm, domain, factor, result);

    // avg_run_length = 15 / 3 = 5
    // max_allowed_length = next_power_of_two(ceil(5 * 1.0)) = next_power_of_two(5) = 8
    // The middle run (length 10, start 2) is split into [8, 2].
    const vector<ulint> expected_lengths = {2, 8, 2, 3};
    const vector<ulint> expected_perm = {0, 2, 10, 12};

    assert(result.lengths == expected_lengths);
    assert(result.interval_permutations == expected_perm);
    assert(result.max_length == 8);
}

void test_split_by_length_capping_mixed_lengths() {
    // Mixture of small and large runs; ensure max_length tracks the maximum.
    const vector<ulint> lengths = {1, 16, 2};
    const vector<ulint> perm = {0, 1, 17};
    const ulint domain = 19;
    const double factor = 0.5;

    SplitResult result;
    split_by_length_capping(lengths, perm, domain, factor, result);

    // avg_run_length = 19 / 3 ≈ 6.333...
    // ceil(avg_run_length * 0.5) = ceil(3.166...) = 4
    // max_allowed_length = next_power_of_two(4) = 4
    //
    // So the run of length 16 is split into four chunks of length 4:
    // [1, 16, 2] -> [1, 4, 4, 4, 4, 2]
    const vector<ulint> expected_lengths = {1, 4, 4, 4, 4, 2};
    const vector<ulint> expected_perm = {0, 1, 5, 9, 13, 17};

    assert(result.lengths == expected_lengths);
    assert(result.interval_permutations == expected_perm);
    assert(result.max_length == 4);
}

// Dummy test for balancing: just ensure the function is callable.
void test_split_by_balancing_dummy() {
    const vector<ulint> lengths = {3, 5, 2};
    const vector<ulint> perm = {0, 3, 8};
    const ulint domain = 10;
    const ulint factor = 4;

    SplitResult result;
    split_by_balancing(lengths, perm, domain, factor, result);
}

int main() {
    test_split_by_length_capping_no_splitting();
    test_split_by_length_capping_with_splitting();
    test_split_by_length_capping_mixed_lengths();
    test_split_by_balancing_dummy();

    std::cout << "move_splitting tests passed" << std::endl;
    return 0;
}
