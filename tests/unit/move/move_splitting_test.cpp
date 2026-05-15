// Unit tests for move_splitting (length capping, balancing, union).
// These are simple assert-based tests, no external framework.

#include "orbit/internal/move/move_splitting.hpp"
#include "orbit/internal/ds/packed_vector.hpp"
#include "orbit/interval_encoding.hpp"
#include "orbit/permutation.hpp"

#include <cassert>
#include <iostream>
#include <vector>

using std::size_t;
using std::vector;

using namespace orbit;

// Use the packed int_vector as IntVectorType for tests
using test_int_vector = int_vector;
using test_split_result = split_result<test_int_vector>;
using test_split_result_union = split_result_union<test_int_vector>;
using test_invertible_encoding = interval_encoding_impl<true, test_int_vector>;

static void assert_img_rank_inv_is_permutation(const test_int_vector &img, size_t n) {
    std::vector<bool> seen(n, false);
    assert(img.size() == n);
    for (size_t i = 0; i < n; ++i) {
        ulint v = img[i];
        assert(v < n);
        assert(!seen[v]);
        seen[v] = true;
    }
}

static void assert_union_result_invariants(const test_split_result_union &result, ulint domain) {
    assert(result.lengths.size() == result.img_rank_inv.size());
    assert(result.lengths.size() == result.is_fwd_interval.size());
    assert(result.lengths.size() == result.is_inv_interval.size());
    ulint sum = 0;
    for (size_t i = 0; i < result.lengths.size(); ++i) {
        ulint len = static_cast<ulint>(result.lengths[i]);
        assert(len > 0);
        assert(len <= result.max_length);
        sum += len;
    }
    assert(sum == domain);
    assert_img_rank_inv_is_permutation(result.img_rank_inv, result.lengths.size());
}

static void assert_union_matches_invertible_encoding(
    const vector<ulint> &lengths_vec,
    const vector<ulint> &img_rank_inv_vec,
    const vector<ulint> &permutation,
    ulint domain,
    ulint max_length
) {
    test_int_vector lengths(lengths_vec);
    test_int_vector img_rank_inv(img_rank_inv_vec);

    test_split_result_union result;
    move_splitting::split_by_union(lengths, img_rank_inv, domain, max_length, result);
    assert_union_result_invariants(result, domain);

    test_invertible_encoding inv = test_invertible_encoding::from_permutation(permutation, NO_SPLITTING);
    assert(inv.intervals() == result.lengths.size());
    assert(inv.max_length() == result.max_length);
    for (size_t i = 0; i < inv.intervals(); ++i) {
        assert(inv.get_length(i) == static_cast<ulint>(result.lengths[i]));
        assert(inv.get_img_rank_inv(i) == static_cast<ulint>(result.img_rank_inv[i]));
        assert(inv.get_is_fwd_interval(i) == static_cast<bool>(result.is_fwd_interval[i]));
        assert(inv.get_is_inv_interval(i) == static_cast<bool>(result.is_inv_interval[i]));
    }
}

static void assert_union_preserves_permutation(
    const vector<ulint> &lengths_vec,
    const vector<ulint> &images_vec,
    ulint domain,
    ulint max_length
) {
    test_int_vector lengths(lengths_vec);
    test_int_vector img_rank_inv = compute_img_rank_inv<test_int_vector>(images_vec);

    test_split_result_union union_result;
    move_splitting::split_by_union(lengths, img_rank_inv, domain, max_length, union_result);
    assert_union_result_invariants(union_result, domain);

    permutation_absolute<> before_perm(lengths_vec, images_vec, NO_SPLITTING);
    auto enc_after = interval_encoding_impl<>::from_lengths_and_img_rank_inv(
        union_result.lengths, union_result.img_rank_inv, domain, union_result.max_length, NO_SPLITTING);
    permutation_absolute<> after_perm(enc_after);

    auto before_pos = before_perm.first();
    auto after_pos = after_perm.first();
    assert(before_pos == after_pos);
    for (size_t i = 0; i < domain; ++i) {
        before_pos = before_perm.next(before_pos);
        after_pos = after_perm.next(after_pos);
        assert(before_pos.idx == after_pos.idx);
    }
    assert(before_pos == before_perm.first());
    assert(after_pos == after_perm.first());
}

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
    move_splitting::split_by_max_allowed_length(lengths, img_rank_inv, max_length, result);

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

// Dummy test for balancing: just ensure the function is callable.
void test_split_by_balancing_splitting() {
    const vector<ulint> lengths_vec = {6, 1, 1, 1, 1, 1, 1};
    const vector<ulint> img_rank_inv_vec = {3, 1, 5, 6, 2, 4, 0};
    const vector<ulint> expected_lengths = {3, 3, 1, 1, 1, 1, 1, 1};
    const vector<ulint> expected_img_rank_inv = {4, 2, 6, 7, 3, 5, 0, 1};
    const ulint domain = 12;
    const ulint factor = 2;

    test_int_vector lengths(lengths_vec);
    test_int_vector img_rank_inv(img_rank_inv_vec);

    test_split_result result;
    move_splitting::split_by_balancing(lengths, img_rank_inv, domain, factor, result);

    for (size_t i = 0; i < expected_lengths.size(); ++i) {
        assert(result.lengths[i] == expected_lengths[i]);
        assert(result.img_rank_inv[i] == expected_img_rank_inv[i]);
    }

    assert(result.max_length == 3);
}

template<typename int_vector_t>
ulint compute_max_input_weight(const int_vector_t &lengths, const int_vector_t &tau_inv, size_t domain) {
    ulint max_weight = 0;

    ulint curr_input_start = 0;
    ulint curr_output_start = 0;
    ulint curr_output_idx = 0;
    for (size_t i = 0; i < lengths.size(); ++i) {
        ulint curr_length = lengths[i];
        ulint curr_weight = 0;
        while (curr_output_start < curr_input_start + curr_length && curr_output_idx < tau_inv.size()) {
            if (curr_output_start > curr_input_start) {
                ++curr_weight;
            }
            curr_output_start += lengths[tau_inv[curr_output_idx]];
            ++curr_output_idx;
        }
        max_weight = std::max(max_weight, curr_weight);
        curr_input_start += curr_length;
    }
    assert(curr_input_start == domain);
    assert(curr_output_start == domain);
    assert(curr_output_idx == tau_inv.size());

    return max_weight;
}

template<typename int_vector_t>
ulint compute_max_output_weight(const int_vector_t &lengths, const int_vector_t &tau_inv, size_t domain) {
    ulint max_weight = 0;

    ulint curr_output_start = 0;
    ulint curr_input_start = 0;
    ulint curr_input_idx = 0;
    for (size_t i = 0; i < lengths.size(); ++i) {
        ulint curr_length = lengths[tau_inv[i]];
        ulint curr_weight = 0;
        while (curr_input_start < curr_output_start + curr_length && curr_input_idx < tau_inv.size()) {
            if (curr_input_start > curr_output_start) {
                ++curr_weight;
            }
            curr_input_start += lengths[curr_input_idx];
            ++curr_input_idx;
        }
        max_weight = std::max(max_weight, curr_weight);
        curr_output_start += curr_length;
    }
    assert(curr_output_start == domain);
    assert(curr_input_start == domain);
    assert(curr_input_idx == tau_inv.size());

    return max_weight;
}

void test_split_by_union_two_equal_runs() {
    // Two runs of length 2 in input order; output order matches input (img_rank_inv identity).
    const vector<ulint> lengths_vec = {2, 2};
    const vector<ulint> img_rank_inv_vec = {0, 1};
    const ulint domain = 4;
    const ulint max_length = 2;

    test_int_vector lengths(lengths_vec);
    test_int_vector img_rank_inv(img_rank_inv_vec);

    test_split_result_union result;
    move_splitting::split_by_union(lengths, img_rank_inv, domain, max_length, result);

    assert(result.lengths.size() == 2);
    assert(result.max_length == 2);
    ulint sum = 0;
    for (size_t i = 0; i < result.lengths.size(); ++i) {
        sum += static_cast<ulint>(result.lengths[i]);
    }
    assert(sum == domain);
    assert_img_rank_inv_is_permutation(result.img_rank_inv, result.lengths.size());
    assert(static_cast<ulint>(result.lengths[0]) == 2);
    assert(static_cast<ulint>(result.lengths[1]) == 2);
    assert(static_cast<ulint>(result.img_rank_inv[0]) == 0);
    assert(static_cast<ulint>(result.img_rank_inv[1]) == 1);
    for (size_t i = 0; i < result.is_fwd_interval.size(); ++i) {
        assert(static_cast<bool>(result.is_fwd_interval[i]));
        assert(static_cast<bool>(result.is_inv_interval[i]));
    }
}

void test_split_by_union_nine_run_regression() {
    // Run lengths / img_rank_inv from the 16-element permutation used in interval_encoding tests
    // (get_permutation_intervals on {1,2,9,...,8}); domain is sum(lengths) == 16.
    const vector<ulint> lengths_vec = {2, 3, 1, 2, 2, 1, 1, 1, 3};
    const vector<ulint> img_rank_inv_vec = {6, 0, 2, 4, 8, 1, 3, 5, 7};
    const ulint domain = 16;
    const ulint max_length = 3;

    test_int_vector lengths(lengths_vec);
    test_int_vector img_rank_inv(img_rank_inv_vec);

    test_split_result_union result;
    move_splitting::split_by_union(lengths, img_rank_inv, domain, max_length, result);

    const vector<ulint> expected_lengths = {1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1};
    const vector<ulint> expected_img_rank_inv = {10, 0, 1, 5, 7, 8, 12, 13, 14, 2, 3, 4, 6, 9, 11};
    const vector<ulint> expected_fwd = {1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0};
    const vector<ulint> expected_inv = {1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1};

    assert(result.lengths.size() == expected_lengths.size());
    assert(result.max_length == 2);
    ulint sum = 0;
    for (size_t i = 0; i < result.lengths.size(); ++i) {
        sum += static_cast<ulint>(result.lengths[i]);
        assert(static_cast<ulint>(result.lengths[i]) == expected_lengths[i]);
        assert(static_cast<ulint>(result.img_rank_inv[i]) == expected_img_rank_inv[i]);
        assert(static_cast<ulint>(result.is_fwd_interval[i]) == expected_fwd[i]);
        assert(static_cast<ulint>(result.is_inv_interval[i]) == expected_inv[i]);
    }
    assert(sum == domain);
    assert_img_rank_inv_is_permutation(result.img_rank_inv, result.lengths.size());
}

void test_split_by_union_single_run_no_split() {
    const vector<ulint> lengths_vec = {4};
    const vector<ulint> img_rank_inv_vec = {0};
    const ulint domain = 4;
    const ulint max_length = 4;

    test_int_vector lengths(lengths_vec);
    test_int_vector img_rank_inv(img_rank_inv_vec);

    test_split_result_union result;
    move_splitting::split_by_union(lengths, img_rank_inv, domain, max_length, result);

    assert(result.lengths.size() == 1);
    assert(static_cast<ulint>(result.lengths[0]) == 4);
    assert(static_cast<ulint>(result.img_rank_inv[0]) == 0);
    assert(static_cast<bool>(result.is_fwd_interval[0]));
    assert(static_cast<bool>(result.is_inv_interval[0]));
    assert(result.max_length == 4);
    assert_union_result_invariants(result, domain);
}

void test_split_by_union_disjoint_runs_no_split() {
    // Input and output interval orders align; no refinement splits.
    const vector<ulint> lengths_vec = {2, 2, 2};
    const vector<ulint> img_rank_inv_vec = {0, 1, 2};
    const ulint domain = 6;
    const ulint max_length = 2;

    test_int_vector lengths(lengths_vec);
    test_int_vector img_rank_inv(img_rank_inv_vec);

    test_split_result_union result;
    move_splitting::split_by_union(lengths, img_rank_inv, domain, max_length, result);

    assert(result.lengths.size() == lengths_vec.size());
    assert(result.max_length == max_length);
    for (size_t i = 0; i < lengths_vec.size(); ++i) {
        assert(static_cast<ulint>(result.lengths[i]) == lengths_vec[i]);
        assert(static_cast<ulint>(result.img_rank_inv[i]) == img_rank_inv_vec[i]);
        assert(static_cast<bool>(result.is_fwd_interval[i]));
        assert(static_cast<bool>(result.is_inv_interval[i]));
    }
    assert_union_result_invariants(result, domain);
}

void test_split_by_union_staggered_two_runs() {
    // Output run 0 (len 2) nests inside input run 0 (len 3); exercises intersecting_offset > 0.
    const vector<ulint> lengths_vec = {3, 2};
    const vector<ulint> img_rank_inv_vec = {1, 0};
    const ulint domain = 5;
    const ulint max_length = 3;

    test_int_vector lengths(lengths_vec);
    test_int_vector img_rank_inv(img_rank_inv_vec);

    test_split_result_union result;
    move_splitting::split_by_union(lengths, img_rank_inv, domain, max_length, result);

    const vector<ulint> expected_lengths = {2, 1, 2};
    const vector<ulint> expected_img_rank_inv = {2, 0, 1};
    const vector<ulint> expected_fwd = {1, 0, 1};
    const vector<ulint> expected_inv = {1, 1, 0};

    assert(result.lengths.size() == expected_lengths.size());
    assert(result.max_length == 2);
    for (size_t i = 0; i < expected_lengths.size(); ++i) {
        assert(static_cast<ulint>(result.lengths[i]) == expected_lengths[i]);
        assert(static_cast<ulint>(result.img_rank_inv[i]) == expected_img_rank_inv[i]);
        assert(static_cast<ulint>(result.is_fwd_interval[i]) == expected_fwd[i]);
        assert(static_cast<ulint>(result.is_inv_interval[i]) == expected_inv[i]);
    }
    assert_union_result_invariants(result, domain);
}

void test_split_by_union_three_element_permutation() {
    // Permutation {2,0,1}: one run stays whole with fwd=true, inv=false on the tail piece.
    const vector<ulint> permutation = {2, 0, 1};
    ulint max_length = 0;
    auto [lengths_vec, images_vec] = get_permutation_intervals(permutation, &max_length);
    const ulint domain = static_cast<ulint>(permutation.size());
    test_int_vector img_rank_inv = compute_img_rank_inv<test_int_vector>(images_vec);

    test_int_vector lengths(lengths_vec);
    test_split_result_union result;
    move_splitting::split_by_union(lengths, img_rank_inv, domain, max_length, result);

    const vector<ulint> expected_lengths = {1, 1, 1};
    const vector<ulint> expected_img_rank_inv = {1, 2, 0};
    const vector<ulint> expected_fwd = {1, 1, 0};
    const vector<ulint> expected_inv = {1, 1, 1};

    assert(result.lengths.size() == expected_lengths.size());
    assert(result.max_length == 1);
    for (size_t i = 0; i < expected_lengths.size(); ++i) {
        assert(static_cast<ulint>(result.lengths[i]) == expected_lengths[i]);
        assert(static_cast<ulint>(result.img_rank_inv[i]) == expected_img_rank_inv[i]);
        assert(static_cast<ulint>(result.is_fwd_interval[i]) == expected_fwd[i]);
        assert(static_cast<ulint>(result.is_inv_interval[i]) == expected_inv[i]);
    }
    assert_union_result_invariants(result, domain);
    assert_union_matches_invertible_encoding(lengths_vec, vector<ulint>(img_rank_inv.begin(), img_rank_inv.end()), permutation, domain, max_length);
    assert_union_preserves_permutation(lengths_vec, images_vec, domain, max_length);
}

void test_split_by_union_matches_invertible_encoding_nine_run() {
    const vector<ulint> permutation = {1, 2, 9, 10, 11, 3, 12, 13, 4, 5, 14, 0, 15, 6, 7, 8};
    ulint max_length = 0;
    auto [lengths_vec, images_vec] = get_permutation_intervals(permutation, &max_length);
    const ulint domain = static_cast<ulint>(permutation.size());
    test_int_vector img_rank_inv = compute_img_rank_inv<test_int_vector>(images_vec);

    assert_union_matches_invertible_encoding(
        lengths_vec, vector<ulint>(img_rank_inv.begin(), img_rank_inv.end()), permutation, domain, max_length);
    assert_union_preserves_permutation(lengths_vec, images_vec, domain, max_length);
}

void test_split_by_balancing_invariant() {
    const vector<ulint> lengths_vec = {4, 3, 1, 2, 2, 1, 4, 1, 1, 1, 2, 1, 1, 2, 1, 2, 1, 24, 1, 5};
    const vector<ulint> images_vec = {30, 47, 34, 50, 35, 52, 37, 53, 41, 54, 42, 55, 44, 56, 45, 58, 46, 1, 0, 29};
    const ulint domain = 60;
    const ulint factor = 2;

    test_int_vector lengths(lengths_vec);
    test_int_vector tau_inv = compute_img_rank_inv<test_int_vector>(images_vec);
    ulint max_output_weight_before = compute_max_output_weight(lengths, tau_inv, domain);
    ulint max_input_weight_before = compute_max_input_weight(lengths, tau_inv, domain);
    assert(max_output_weight_before == 13);
    assert(max_input_weight_before == 12);

    test_split_result result;
    move_splitting::split_by_balancing(lengths, tau_inv, domain, factor, result);
    ulint max_output_weight_after = compute_max_output_weight(result.lengths, result.img_rank_inv, domain);
    ulint max_input_weight_after = compute_max_input_weight(result.lengths, result.img_rank_inv, domain);
    assert(max_output_weight_after < 2 * factor);
    assert(max_input_weight_after < 2 * factor);

    permutation_absolute<> before_perm(lengths, images_vec, NO_SPLITTING);
    auto enc_after = interval_encoding_impl<>::from_lengths_and_img_rank_inv(result.lengths, result.img_rank_inv, domain, result.max_length, NO_SPLITTING);
    permutation_absolute<> after_perm(enc_after);
    
    auto before_pos = before_perm.first();
    auto after_pos = after_perm.first();
    assert(before_pos == after_pos);
    for (size_t i = 0; i < domain; ++i) {
        before_pos = before_perm.next(before_pos);
        after_pos = after_perm.next(after_pos);
        assert(before_pos.idx == after_pos.idx);
    }
    assert(before_pos == before_perm.first());
    assert(after_pos == after_perm.first());
}

int main() {
    test_split_by_length_capping_no_splitting();
    test_split_by_length_capping_with_splitting();
    test_split_by_length_capping_mixed_lengths();
    test_split_by_union_two_equal_runs();
    test_split_by_union_nine_run_regression();
    test_split_by_union_single_run_no_split();
    test_split_by_union_disjoint_runs_no_split();
    test_split_by_union_staggered_two_runs();
    test_split_by_union_three_element_permutation();
    test_split_by_union_matches_invertible_encoding_nine_run();
    test_split_by_balancing_no_change();
    test_split_by_balancing_splitting();
    test_split_by_balancing_invariant();

    std::cout << "move_splitting tests passed" << std::endl;
    return 0;
}
