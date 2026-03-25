// Unit tests for PermutationImpl and permutation helpers.
// Simple assert-based tests, no external framework.

#include "orbit/permutation.hpp"
#include "orbit/common.hpp"
#include "orbit/internal/ds/packed_vector.hpp"
#include "orbit/internal/move/move_splitting.hpp"

#include <cassert>
#include <iostream>
#include <vector>

using std::size_t;
using std::vector;

using namespace orbit;

// Use the packed int_vector as IntVectorType for tests.
using test_int_vector = int_vector;
using test_permutation = permutation_impl<test_int_vector>;

// Basic helper to check that img_rank_inv is a permutation of [0, n).
static void assert_img_rank_inv_is_permutation(test_permutation &perm) {
    const size_t r = perm.runs();
    vector<bool> seen(r, false);
    for (size_t i = 0; i < r; ++i) {
        ulint v = perm.get_img_rank_inv(i);
        assert(v < r);
        assert(!seen[v]);
        seen[v] = true;
    }
}

void test_inverse_permutation() {
    {
        const vector<ulint> perm = {0, 1, 2, 3, 4};
        const vector<ulint> expected_inv_perm = {0, 1, 2, 3, 4};
        assert(get_inverse_permutation(perm) == expected_inv_perm);
    }
}

void test_permutation_intervals_trivial_and_runs() {
    {
        // Empty permutation -> empty intervals.
        vector<ulint> perm;
        auto [lengths, intervals] = get_permutation_intervals(perm);
        assert(lengths.empty());
        assert(intervals.empty());
    }
    {
        // Single element permutation.
        vector<ulint> perm = {5};
        auto [lengths, intervals] = get_permutation_intervals(perm);
        assert(lengths.size() == 1);
        assert(intervals.size() == 1);
        assert(lengths[0] == 1);
        assert(intervals[0] == 5);
    }
    {
        // Increasing consecutive block followed by a jump.
        // permutation: [3,4,5, 10,11]
        vector<ulint> perm = {3, 4, 5, 10, 11};
        auto [lengths, intervals] = get_permutation_intervals(perm);

        assert(lengths.size() == 2);
        assert(intervals.size() == 2);

        // First interval [3,4,5] of length 3 starting at 3.
        assert(lengths[0] == 3);
        assert(intervals[0] == 3);

        // Second interval [10,11] of length 2 starting at 10.
        assert(lengths[1] == 2);
        assert(intervals[1] == 10);
    }
}

static void test_permutation_helpers() {
    // starts_to_lengths
    {
        const vector<ulint> starts = {0, 3, 8};
        const ulint domain = 10;
        auto [lengths, max_len] = starts_to_lengths(starts, domain);
        const vector<ulint> expected_lengths = {3, 5, 2};
        for (size_t i = 0; i < lengths.size(); ++i) {
            assert(lengths[i] == expected_lengths[i]);
        }
        assert(max_len == 5);
    }

    // sum_and_max
    {
        const vector<ulint> data = {1, 5, 3};
        auto [sum, max_val] = sum_and_max(data);
        assert(sum == 9);
        assert(max_val == 5);
    }

    // get_img_rank_inv
    {
        const vector<ulint> interval_starts = {5, 1, 10};
        auto img_rank_inv = compute_img_rank_inv(interval_starts);
        // Sorted by interval_starts: indices [1, 0, 2].
        const vector<ulint> expected_img_rank_inv = {1, 0, 2};
        for (size_t i = 0; i < img_rank_inv.size(); ++i) {
            assert(img_rank_inv[i] == expected_img_rank_inv[i]);
        }
    }

    // get_permutation_intervals (from common.hpp)
    {
        const vector<ulint> perm = {2, 3, 4, 10, 11};
        auto [lengths, interval_perm] = get_permutation_intervals(perm);
        // Runs: [2,3,4] and [10,11]
        const vector<ulint> expected_lengths = {3, 2};
        const vector<ulint> expected_interval_perm = {2, 10};
        assert(lengths == expected_lengths);
        assert(interval_perm == expected_interval_perm);
    }
}

static void test_all_construction_paths_no_splitting() {
    // Use the same example permutation as other tests, and derive all
    // run-level structures from it to avoid invalid combinations.
    const vector<ulint> permutation = {1, 2, 9, 10, 11, 3, 12, 13, 4, 5, 14, 0, 15, 6, 7, 8};
    const ulint domain = static_cast<ulint>(permutation.size());

    ulint max_length = 0;
    auto [lengths, interval_perm] = get_permutation_intervals(permutation, &max_length);
    assert(!interval_perm.empty());
    // 0 must appear in interval_perm since 0 cannot be consecutive with a predecessor.
    assert(std::find(interval_perm.begin(), interval_perm.end(), 0) != interval_perm.end());

    // img_rank_inv as defined by the library (ranks intervals by their output start).
    auto img_rank_inv_packed_vector = compute_img_rank_inv(interval_perm);
    std::vector<ulint> img_rank_inv(img_rank_inv_packed_vector.size());
    for (size_t i = 0; i < img_rank_inv_packed_vector.size(); ++i) {
        img_rank_inv[i] = img_rank_inv_packed_vector.get(i);
    }

    // img_rank is the inverse of img_rank_inv: img_rank[rank] = original_interval_index.
    vector<ulint> img_rank = get_inverse_permutation(img_rank_inv);

    // starts derived from lengths.
    vector<ulint> starts;
    starts.reserve(lengths.size());
    ulint pref = 0;
    for (size_t i = 0; i < lengths.size(); ++i) {
        starts.push_back(pref);
        pref += lengths[i];
    }
    assert(pref == domain);

    // 1) lengths + img_rank_inv (with and without explicit domain/max_length)
    {
        test_permutation p1 = test_permutation::from_lengths_and_img_rank_inv(
            lengths, img_rank_inv, NO_SPLITTING
        );
        assert(p1.domain() == domain);
        assert(p1.runs() == lengths.size());
        assert(p1.intervals() == lengths.size());
        assert(p1.max_length() == max_length);
        for (size_t i = 0; i < lengths.size(); ++i) {
            assert(p1.get_length(i) == lengths[i]);
            assert(p1.get_img_rank_inv(i) == img_rank_inv[i]);
        }
        assert_img_rank_inv_is_permutation(p1);

        test_permutation p2 = test_permutation::from_lengths_and_img_rank_inv(
            lengths, img_rank_inv, domain, max_length, NO_SPLITTING
        );
        assert(p2.domain() == domain);
        assert(p2.runs() == lengths.size());
        assert(p2.intervals() == lengths.size());
        assert(p2.max_length() == max_length);
        for (size_t i = 0; i < lengths.size(); ++i) {
            assert(p2.get_length(i) == lengths[i]);
            assert(p2.get_img_rank_inv(i) == img_rank_inv[i]);
        }
        assert_img_rank_inv_is_permutation(p2);
    }

    // 2) lengths + img_rank (with and without explicit domain/max_length)
    {
        // Expected img_rank_inv is inverse of img_rank (and must match the one we computed above).
        vector<ulint> expected_img_rank_inv(img_rank.size());
        for (size_t rank = 0; rank < img_rank.size(); ++rank) {
            expected_img_rank_inv[img_rank[rank]] = static_cast<ulint>(rank);
        }
        assert(expected_img_rank_inv == img_rank_inv);

        test_permutation p1 = test_permutation::from_lengths_and_img_rank(
            lengths, img_rank, NO_SPLITTING
        );
        assert(p1.domain() == domain);
        assert(p1.runs() == lengths.size());
        assert(p1.intervals() == lengths.size());
        assert(p1.max_length() == max_length);
        for (size_t i = 0; i < lengths.size(); ++i) {
            assert(p1.get_length(i) == lengths[i]);
            assert(p1.get_img_rank_inv(i) == expected_img_rank_inv[i]);
        }
        assert_img_rank_inv_is_permutation(p1);

        test_permutation p2 = test_permutation::from_lengths_and_img_rank(
            lengths, img_rank, domain, max_length, NO_SPLITTING
        );
        assert(p2.domain() == domain);
        assert(p2.runs() == lengths.size());
        assert(p2.intervals() == lengths.size());
        assert(p2.max_length() == max_length);
        for (size_t i = 0; i < lengths.size(); ++i) {
            assert(p2.get_length(i) == lengths[i]);
            assert(p2.get_img_rank_inv(i) == expected_img_rank_inv[i]);
        }
        assert_img_rank_inv_is_permutation(p2);
    }

    // 3) lengths + image
    {
        test_permutation p1 = test_permutation::from_lengths_and_images(
            lengths, interval_perm, NO_SPLITTING
        );
        assert(p1.domain() == domain);
        assert(p1.runs() == lengths.size());
        assert(p1.intervals() == lengths.size());
        assert(p1.max_length() == max_length);
        assert_img_rank_inv_is_permutation(p1);

        test_permutation p2 = test_permutation::from_lengths_and_images(
            lengths, interval_perm, domain, max_length, NO_SPLITTING
        );
        assert(p2.domain() == domain);
        assert(p2.runs() == lengths.size());
        assert(p2.intervals() == lengths.size());
        assert(p2.max_length() == max_length);
        assert_img_rank_inv_is_permutation(p2);
    }

    // 4) starts + img_rank_inv
    {
        test_permutation p1 = test_permutation::from_starts_and_img_rank_inv(
            starts, img_rank_inv, domain, NO_SPLITTING
        );
        assert(p1.domain() == domain);
        assert(p1.runs() == lengths.size());
        assert(p1.intervals() == lengths.size());
        assert(p1.max_length() == max_length);
        for (size_t i = 0; i < lengths.size(); ++i) {
            assert(p1.get_length(i) == lengths[i]);
            assert(p1.get_img_rank_inv(i) == img_rank_inv[i]);
        }
        assert_img_rank_inv_is_permutation(p1);

        test_permutation p2 = test_permutation::from_starts_and_img_rank_inv(
            starts, img_rank_inv, domain, max_length, NO_SPLITTING
        );
        assert(p2.domain() == domain);
        assert(p2.runs() == lengths.size());
        assert(p2.intervals() == lengths.size());
        assert(p2.max_length() == max_length);
        for (size_t i = 0; i < lengths.size(); ++i) {
            assert(p2.get_length(i) == lengths[i]);
            assert(p2.get_img_rank_inv(i) == img_rank_inv[i]);
        }
        assert_img_rank_inv_is_permutation(p2);
    }

    // 5) starts + img_rank
    {
        vector<ulint> expected_img_rank_inv(img_rank.size());
        for (size_t rank = 0; rank < img_rank.size(); ++rank) {
            expected_img_rank_inv[img_rank[rank]] = static_cast<ulint>(rank);
        }
        assert(expected_img_rank_inv == img_rank_inv);

        test_permutation p1 = test_permutation::from_starts_and_img_rank(
            starts, img_rank, domain, NO_SPLITTING
        );
        assert(p1.domain() == domain);
        assert(p1.runs() == lengths.size());
        assert(p1.intervals() == lengths.size());
        assert(p1.max_length() == max_length);
        for (size_t i = 0; i < lengths.size(); ++i) {
            assert(p1.get_length(i) == lengths[i]);
            assert(p1.get_img_rank_inv(i) == expected_img_rank_inv[i]);
        }
        assert_img_rank_inv_is_permutation(p1);

        test_permutation p2 = test_permutation::from_starts_and_img_rank(
            starts, img_rank, domain, max_length, NO_SPLITTING
        );
        assert(p2.domain() == domain);
        assert(p2.runs() == lengths.size());
        assert(p2.intervals() == lengths.size());
        assert(p2.max_length() == max_length);
        for (size_t i = 0; i < lengths.size(); ++i) {
            assert(p2.get_length(i) == lengths[i]);
            assert(p2.get_img_rank_inv(i) == expected_img_rank_inv[i]);
        }
        assert_img_rank_inv_is_permutation(p2);
    }

    // 6) starts + image
    {
        test_permutation p1 = test_permutation::from_starts_and_images(
            starts, interval_perm, domain, NO_SPLITTING
        );
        assert(p1.domain() == domain);
        assert(p1.runs() == lengths.size());
        assert(p1.intervals() == lengths.size());
        assert(p1.max_length() == max_length);
        assert_img_rank_inv_is_permutation(p1);

        test_permutation p2 = test_permutation::from_starts_and_images(
            starts, interval_perm, domain, max_length, NO_SPLITTING
        );
        assert(p2.domain() == domain);
        assert(p2.runs() == lengths.size());
        assert(p2.intervals() == lengths.size());
        assert(p2.max_length() == max_length);
        assert_img_rank_inv_is_permutation(p2);
    }
}

static void test_from_permutation_and_split_run_data() {
    // Permutation with three runs of consecutive values.
    const vector<ulint> perm_vec = {0, 1, 2, 5, 6, 7, 8, 10};
    const ulint domain = static_cast<ulint>(perm_vec.size());

    test_permutation perm_no_split = test_permutation::from_permutation(
        perm_vec, NO_SPLITTING
    );

    assert(perm_no_split.domain() == domain);
    assert(perm_no_split.runs() == 3);
    assert(perm_no_split.intervals() == 3);

    const vector<ulint> expected_lengths = {3, 4, 1};
    const vector<ulint> expected_img_rank_inv = {0, 1, 2};
    assert(perm_no_split.max_length() == 4);
    for (size_t i = 0; i < expected_lengths.size(); ++i) {
        assert(perm_no_split.get_length(i) == expected_lengths[i]);
        assert(perm_no_split.get_img_rank_inv(i) == expected_img_rank_inv[i]);
    }
    assert_img_rank_inv_is_permutation(perm_no_split);

    // Now build a permutation where splitting should occur and test split_run_data_with_copy.
    const vector<ulint> original_lengths = {2, 10, 3};
    const vector<ulint> img_rank_inv = {0, 1, 2};
    const ulint domain2 = 15;
    const ulint max_length2 = 10;

    test_permutation perm_split = test_permutation::from_lengths_and_img_rank_inv(
        original_lengths, img_rank_inv, domain2, max_length2, DEFAULT_SPLITTING
    );

    assert(perm_split.domain() == domain2);
    assert(perm_split.runs() == original_lengths.size());
    assert(perm_split.intervals() >= perm_split.runs());
    assert_img_rank_inv_is_permutation(perm_split);

    // Run data: just store the original run index.
    vector<ulint> run_data(original_lengths.size());
    for (size_t i = 0; i < original_lengths.size(); ++i) {
        run_data[i] = static_cast<ulint>(i);
    }

    auto new_run_data = perm_split.split_run_data_with_copy(original_lengths, run_data);
    assert(new_run_data.size() == perm_split.intervals());

    // Each split interval derived from original run i should carry run_data[i].
    size_t new_idx = 0;
    for (size_t i = 0; i < original_lengths.size(); ++i) {
        ulint remaining = original_lengths[i];
        while (remaining > 0) {
            ulint chunk = perm_split.get_length(new_idx);
            assert(new_run_data[new_idx] == run_data[i]);
            remaining -= chunk;
            ++new_idx;
        }
    }
    assert(new_idx == perm_split.intervals());
}

int main() {
    test_inverse_permutation();
    test_permutation_intervals_trivial_and_runs();
    test_permutation_helpers();
    test_all_construction_paths_no_splitting();
    test_from_permutation_and_split_run_data();

    std::cout << "permutation tests passed" << std::endl;
    return 0;
}

