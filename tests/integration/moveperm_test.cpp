// Integration-style tests for MovePerm.
// These are simple assert-based tests, no external framework.

#include "runperm.hpp"

#include <cassert>
#include <iostream>
#include <vector>

using std::size_t;
using std::vector;

// Helper to compute global index for a relative MovePerm position.
static ulint compute_idx_relative(const MovePermRelative &mp,
                                  typename MovePermRelative::Position pos) {
    ulint idx = 0;
    for (ulint i = 0; i < pos.interval; ++i) {
        idx += mp.get_length(i);
    }
    idx += pos.offset;
    return idx;
}

static MovePermRelative::Position make_pos_relative(const MovePermRelative &mp,
                                                    ulint idx) {
    MovePermRelative::Position pos{};
    ulint prefix = 0;
    for (ulint interval = 0; interval < mp.move_runs(); ++interval) {
        ulint len = mp.get_length(interval);
        if (idx < prefix + len) {
            pos.interval = interval;
            pos.offset = idx - prefix;
            return pos;
        }
        prefix += len;
    }
    assert(false && "index out of range");
    return pos;
}

static MovePermAbsolute::Position make_pos_absolute(const MovePermAbsolute &mp,
                                                    ulint idx) {
    MovePermAbsolute::Position pos{};
    pos.idx = idx;
    ulint prefix = 0;
    for (ulint interval = 0; interval < mp.move_runs(); ++interval) {
        ulint len = mp.get_length(interval);
        if (idx < prefix + len) {
            pos.interval = interval;
            pos.offset = idx - prefix;
            return pos;
        }
        prefix += len;
    }
    assert(false && "index out of range");
    return pos;
}

static void integration_move_perm_relative_and_absolute() {
    // Same permutation used in examples.cpp example2.
    vector<ulint> perm = {1, 2, 9, 10, 11, 3, 12, 13, 4, 5, 14, 0, 15, 6, 7, 8};
    const ulint domain = static_cast<ulint>(perm.size());

    // Relative constructor from full permutation.
    MovePermRelative mp_rel(perm);
    assert(mp_rel.domain() == domain);

    for (ulint idx = 0; idx < domain; ++idx) {
        auto pos = make_pos_relative(mp_rel, idx);
        auto next_pos = mp_rel.next(pos);
        ulint idx_next = compute_idx_relative(mp_rel, next_pos);
        assert(idx_next == perm[idx]);
    }

    // Absolute constructor from lengths + interval permutation.
    auto [lengths, interval_perm] = get_permutation_intervals(perm);
    MovePermAbsolute mp_abs(lengths, interval_perm, domain);
    assert(mp_abs.domain() == domain);

    for (ulint idx = 0; idx < domain; ++idx) {
        auto pos = make_pos_absolute(mp_abs, idx);
        auto next_pos = mp_abs.next(pos);
        assert(next_pos.idx == perm[idx]);
    }
}

static void integration_move_perm_splitting_equivalence() {
    // Example3-like configuration: long run that can be split.
    const vector<ulint> lengths = {2, 1, 8};
    const vector<ulint> interval_perm = {9, 0, 1};
    const ulint domain = 11;

    // Build equivalent full permutation to validate semantics.
    vector<ulint> perm(domain);
    // Reconstruct perm[i] = interval_perm[j] + offset within interval j.
    ulint src_idx = 0;
    for (size_t j = 0; j < lengths.size(); ++j) {
        for (ulint o = 0; o < lengths[j]; ++o) {
            perm[src_idx++] = interval_perm[j] + o;
        }
    }

    MovePermRelative mp_no_split(lengths, interval_perm, domain);

    SplitParams split;
    split.length_capping_factor = 1.0;
    MovePermRelative mp_split(lengths, interval_perm, domain, split);

    assert(mp_no_split.domain() == domain);
    assert(mp_split.domain() == domain);

    for (ulint idx = 0; idx < domain; ++idx) {
        auto pos0 = make_pos_relative(mp_no_split, idx);
        auto pos1 = make_pos_relative(mp_split, idx);

        auto next0 = mp_no_split.next(pos0);
        auto next1 = mp_split.next(pos1);

        ulint idx0 = compute_idx_relative(mp_no_split, next0);
        ulint idx1 = compute_idx_relative(mp_split, next1);

        assert(idx0 == perm[idx]);
        assert(idx1 == perm[idx]);
    }
}

int main() {
    integration_move_perm_relative_and_absolute();
    integration_move_perm_splitting_equivalence();

    std::cout << "moveperm integration tests passed" << std::endl;
    return 0;
}

