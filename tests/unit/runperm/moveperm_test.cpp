// Unit tests for MovePerm (MovePermImpl via public alias).
// Simple assert-based tests, no external framework.

#include "runperm.hpp"

#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>

using std::size_t;
using std::vector;

// Helper to compute the global index for a relative MovePerm position.
static ulint compute_idx_relative(const MovePermRelative &mp,
                                  typename MovePermRelative::Position pos) {
    ulint idx = 0;
    for (ulint i = 0; i < pos.interval; ++i) {
        idx += mp.get_length(i);
    }
    idx += pos.offset;
    return idx;
}

// Helper to build a relative position from a global index.
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

// Helper to build an absolute position from a global index.
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

static void test_move_perm_relative_from_permutation() {
    // Simple permutation with several runs of contiguous values.
    vector<ulint> perm = {1, 2, 9, 10, 11, 3, 12, 13, 4, 5, 14, 0, 15, 6, 7, 8};

    MovePermRelative mp(perm);
    assert(mp.domain() == perm.size());

    for (ulint idx = 0; idx < mp.domain(); ++idx) {
        auto pos = make_pos_relative(mp, idx);
        auto next_pos = mp.next(pos);
        ulint idx_next = compute_idx_relative(mp, next_pos);
        assert(idx_next == perm[idx]);
    }
}

static void test_move_perm_absolute_from_lengths() {
    const vector<ulint> perm = {1, 2, 9, 10, 11, 3, 12, 13, 4, 5, 14, 0, 15, 6, 7, 8};
    const ulint domain = static_cast<ulint>(perm.size());

    auto [lengths, interval_perm] = get_permutation_intervals(perm);

    MovePermAbsolute mp(lengths, interval_perm, domain);
    assert(mp.domain() == domain);
    assert(mp.move_runs() == lengths.size());

    // Check that each absolute position maps according to perm.
    for (ulint idx = 0; idx < domain; ++idx) {
        auto pos = make_pos_absolute(mp, idx);
        auto next_pos = mp.next(pos);
        assert(next_pos.idx == perm[idx]);
    }
}

static void test_move_perm_relative_splitting_preserves_permutation() {
    // Permutation with a long run to trigger length capping.
    vector<ulint> perm = {4, 5, 6, 7, 0, 1, 2, 3};
    const ulint domain = static_cast<ulint>(perm.size());

    MovePermRelative mp_no_split(perm);

    SplitParams split;
    split.length_capping_factor = 1.0; // Strong capping to force splitting if beneficial.
    MovePermRelative mp_split(perm, split);

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

static void test_move_perm_serialize_roundtrip_absolute() {
    const vector<ulint> perm = {3, 0, 1, 2};
    const ulint domain = static_cast<ulint>(perm.size());
    auto [lengths, interval_perm] = get_permutation_intervals(perm);

    MovePermAbsolute mp(lengths, interval_perm, domain);

    std::stringstream ss;
    size_t bytes = mp.serialize(ss);
    assert(bytes > 0);

    MovePermAbsolute loaded;
    loaded.load(ss);

    assert(loaded.domain() == mp.domain());
    assert(loaded.move_runs() == mp.move_runs());

    for (ulint idx = 0; idx < domain; ++idx) {
        auto pos = make_pos_absolute(loaded, idx);
        auto next_pos = loaded.next(pos);
        assert(next_pos.idx == perm[idx]);
    }
}

static void test_move_perm_next_with_steps() {
    // Small permutation to exercise multi-step next.
    vector<ulint> perm = {1, 2, 0, 4, 3};
    const ulint domain = static_cast<ulint>(perm.size());

    MovePermRelative mp_rel(perm);

    for (ulint idx = 0; idx < domain; ++idx) {
        auto pos = make_pos_relative(mp_rel, idx);

        // Apply next with steps = 3.
        auto pos_steps = mp_rel.next(pos, 3);

        // Apply single-step next three times.
        auto pos_iter = pos;
        for (int k = 0; k < 3; ++k) {
            pos_iter = mp_rel.next(pos_iter);
        }

        ulint idx_steps = compute_idx_relative(mp_rel, pos_steps);
        ulint idx_iter = compute_idx_relative(mp_rel, pos_iter);
        assert(idx_steps == idx_iter);
    }
}
int main() {
    test_move_perm_relative_from_permutation();
    test_move_perm_absolute_from_lengths();
    test_move_perm_relative_splitting_preserves_permutation();
    test_move_perm_serialize_roundtrip_absolute();
    test_move_perm_next_with_steps();

    std::cout << "moveperm unit tests passed" << std::endl;
    return 0;
}
