// Unit tests for RLBWT-specific permutation builder (`RLBWTPermutationImpl`).
// Focus: construction invariants and head/alphabet mapping (with/without splitting).

#include "internal/rlbwt/specializations/rlbwt_permutation.hpp"

#include <cassert>
#include <iostream>
#include <vector>

using std::size_t;
using std::vector;

template<class Perm>
static vector<uchar> expected_heads_for_permutation(
    const Perm& perm,
    const vector<uchar>& orig_heads,
    const vector<ulint>& orig_lengths
) {
    vector<uchar> expected;
    expected.reserve(perm.intervals());
    size_t new_interval_idx = 0;
    for (size_t i = 0; i < orig_lengths.size(); ++i) {
        ulint remaining = orig_lengths[i];
        while (remaining > 0) {
            ulint chunk = perm.get_length(new_interval_idx);
            expected.push_back(orig_heads[i]);
            remaining -= chunk;
            ++new_interval_idx;
        }
    }
    assert(new_interval_idx == perm.intervals());
    return expected;
}

void test_rlbwt_lf_permutation_heads_and_alphabet() {
    // Known-valid RLBWT example (see integration tests).
    vector<uchar> heads =       {'T','C','G','A','T', 1 ,'A','T','A'};
    vector<ulint> lengths  =    { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };

    // No splitting: intervals == runs == original runs.
    auto perm = RLBWTPermutationImpl<>::lf_permutation(heads, lengths, NO_SPLITTING);
    assert(perm.domain() == 27);
    assert(perm.runs() == heads.size());
    assert(perm.intervals() == heads.size());
    assert(perm.sigma() == Nucleotide::size());

    const auto& alpha = perm.get_alphabet();
    const auto& perm_heads = perm.get_heads();
    assert(perm_heads.size() == perm.intervals());
    for (size_t i = 0; i < perm.intervals(); ++i) {
        assert(alpha.unmap_char(static_cast<uchar>(perm_heads[i])) == heads[i]);
    }

    // With splitting: heads must be expanded consistently with new intervals.
    auto perm_split = RLBWTPermutationImpl<>::lf_permutation(heads, lengths, DEFAULT_SPLITTING);
    assert(perm_split.domain() == 27);
    assert(perm_split.runs() == heads.size());
    assert(perm_split.intervals() >= heads.size());
    assert(perm_split.get_split_params() == DEFAULT_SPLITTING);

    const auto& alpha_split = perm_split.get_alphabet();
    const auto& perm_heads_split = perm_split.get_heads();
    assert(perm_heads_split.size() == perm_split.intervals());

    auto expected = expected_heads_for_permutation(perm_split, heads, lengths);
    for (size_t i = 0; i < expected.size(); ++i) {
        assert(alpha_split.unmap_char(static_cast<uchar>(perm_heads_split[i])) == expected[i]);
    }
}

void test_rlbwt_fl_permutation_heads_and_alphabet() {
    // Known-valid RLBWT example (see integration tests).
    vector<uchar> heads =       {'T','C','G','A','T', 1 ,'A','T','A'};
    vector<ulint> lengths  =    { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };

    // For FL, the internal builder derives "F runs" (different heads/lengths).
    auto perm = RLBWTPermutationImpl<>::fl_permutation(heads, lengths, NO_SPLITTING);
    assert(perm.domain() == 27);
    assert(perm.runs() == heads.size());
    assert(perm.sigma() == Nucleotide::size());

    auto [head_counts, F_lens_and_origin_run, n, max_length] = fl::get_FL_head_counts(heads, lengths);
    (void)head_counts;
    (void)max_length;
    assert(n == perm.domain());

    auto [F_heads, F_lens, F_tau_inv] = fl::get_FL_runs_and_tau_inv(heads.size(), F_lens_and_origin_run);
    (void)F_tau_inv;
    assert(F_heads.size() == heads.size());
    assert(F_lens.size() == lengths.size());

    const auto& alpha = perm.get_alphabet();
    const auto& perm_heads = perm.get_heads();
    assert(perm_heads.size() == perm.intervals());

    auto expected = expected_heads_for_permutation(perm, F_heads, F_lens);
    for (size_t i = 0; i < expected.size(); ++i) {
        assert(alpha.unmap_char(static_cast<uchar>(perm_heads[i])) == expected[i]);
    }
}

int main() {
    test_rlbwt_lf_permutation_heads_and_alphabet();
    test_rlbwt_fl_permutation_heads_and_alphabet();
    std::cout << "rlbwt_permutation unit tests passed" << std::endl;
    return 0;
}

