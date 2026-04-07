// Unit tests for split_by_length_capping and split_by_balancing.
// These are simple assert-based tests, no external framework.

#include "orbit/internal/move/move_splitting.hpp"
#include "orbit/internal/ds/packed_vector.hpp"
#include "orbit/interval_encoding.hpp"

#include <cassert>
#include <iostream>
#include <vector>

using std::size_t;
using std::vector;

using namespace orbit;

// Use the packed int_vector as IntVectorType for tests
using test_int_vector = int_vector;
using test_split_result = split_result<test_int_vector>;

void test_split_by_length_capping_no_splitting() {
    // Three equal-length runs that should not be split.
    const vector<ulint> lengths_vec = {4, 4, 4};
    const vector<ulint> perm_vec = {0, 4, 8};
    const vector<ulint> img_rank_inv_vec = {0, 1, 2};
    const ulint domain = 12; // Sum of lengths.
    const double factor = 1.0;

    test_int_vector lengths(lengths_vec);
    test_int_vector img_rank_inv(img_rank_inv_vec);

    test_split_result result;
    move_splitting::split_by_length_capping(lengths, img_rank_inv, domain, factor, result);

    assert(result.lengths.size() == lengths_vec.size());
    assert(result.img_rank_inv.size() == img_rank_inv_vec.size());

    for (size_t i = 0; i < lengths_vec.size(); ++i) {
        assert(result.lengths[i] == lengths_vec[i]);
        assert(result.img_rank_inv[i] == img_rank_inv_vec[i]);
    }

    // With these parameters, max_allowed_length is 4, so max_length should be 4.
    assert(result.max_length == 4);
}

void test_split_by_length_capping_with_splitting() {
    // One run is longer than the allowed maximum and should be split.
    const vector<ulint> lengths_vec = {2, 3, 1, 2, 2, 1, 1, 1, 3};
    const vector<ulint> perm_vec = {1, 9, 3, 12, 4, 14, 0, 15, 6};
    const vector<ulint> img_rank_inv_vec = {6, 0, 2, 4, 8, 1, 3, 5, 7};
    const ulint domain = 15; // Sum of lengths.
    size_t max_length = 2;

    test_int_vector lengths(lengths_vec);
    test_int_vector img_rank_inv(img_rank_inv_vec);

    test_split_result result;
    move_splitting::split_by_max_allowed_length(lengths, img_rank_inv, domain, max_length, result);

    const vector<ulint> expected_lengths = {2, 2, 1, 1, 2, 2, 1, 1, 1, 2, 1};
    const vector<ulint> expected_perm = {1, 9, 11, 3, 12, 4, 14, 0, 15, 6, 8};
    const vector<ulint> expected_img_rank_inv = {7, 0, 3, 5, 9, 10, 1, 2, 4, 6, 8};

    assert(result.lengths.size() == expected_lengths.size());
    assert(result.img_rank_inv.size() == expected_img_rank_inv.size());

    for (size_t i = 0; i < expected_lengths.size(); ++i) {
        assert(result.lengths[i] == expected_lengths[i]);
        assert(result.img_rank_inv[i] == expected_img_rank_inv[i]);
    }
    assert(result.max_length == max_length);
}

void test_split_by_length_capping_mixed_lengths() {
    // Mixture of small and large runs; ensure max_length tracks the maximum.
    const vector<ulint> lengths_vec = {1, 16, 2};
    const vector<ulint> perm_vec    = {0, 1, 17};
    const vector<ulint> img_rank_inv_vec = {0, 1, 2};
    const ulint domain = 19;
    const double factor = 0.5;

    test_int_vector lengths(lengths_vec);
    test_int_vector img_rank_inv(img_rank_inv_vec);

    test_split_result result;
    move_splitting::split_by_length_capping(lengths, img_rank_inv, domain, factor, result);

    // avg_run_length = 19 / 3 ≈ 6.333...
    // ceil(avg_run_length * 0.5) = ceil(3.166...) = 4
    // bit_width(4) = 3, so max_allowed_length = 2^3 - 1 = 7
    //
    // So the run of length 16 is split into [7, 7, 2]:
    // [1, 16, 2] -> [1, 7, 7, 2, 2]
    const vector<ulint> expected_lengths = {1, 7, 7, 2, 2};
    const vector<ulint> expected_perm    = {0, 1, 8, 15, 17};
    const vector<ulint> expected_img_rank_inv = {0, 1, 2, 3, 4};

    assert(result.lengths.size() == expected_lengths.size());
    assert(result.img_rank_inv.size()   == expected_img_rank_inv.size());

    for (size_t i = 0; i < expected_lengths.size(); ++i) {
        assert(result.lengths[i] == expected_lengths[i]);
        assert(result.img_rank_inv[i] == expected_img_rank_inv[i]);
    }
    assert(result.max_length == 7);
}

// Dummy test for balancing: just ensure the function is callable.
void test_split_by_balancing_no_change() {
    const vector<ulint> lengths_vec = {2, 3, 1, 2, 2, 1, 1, 1, 3};
    const vector<ulint> perm_vec = {1, 9, 3, 12, 4, 14, 0, 15, 6};
    const vector<ulint> img_rank_inv_vec = {6, 0, 2, 4, 8, 1, 3, 5, 7};
    const ulint domain = 16;
    const ulint factor = 4;

    test_int_vector lengths(lengths_vec);
    test_int_vector img_rank_inv(img_rank_inv_vec);

    test_split_result result;
    move_splitting::split_by_balancing(lengths, img_rank_inv, domain, factor, result);

    for (size_t i = 0; i < lengths_vec.size(); ++i) {
        assert(result.lengths[i] == lengths_vec[i]);
        assert(result.img_rank_inv[i] == img_rank_inv_vec[i]);
    }
    assert(result.max_length == 3);
}

int main() {
    test_split_by_length_capping_no_splitting();
    test_split_by_length_capping_with_splitting();
    test_split_by_length_capping_mixed_lengths();
    test_split_by_balancing_no_change();

    std::cout << "move_splitting tests passed" << std::endl;
    return 0;
}
