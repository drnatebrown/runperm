// Integration-style tests for MoveStructure navigation semantics.
// These are simple assert-based tests, no external framework.

#include "internal/move/move_structure.hpp"

#include <cassert>
#include <iostream>
#include <vector>

using std::size_t;
using std::vector;

template <typename ColumnsT>
static ulint compute_global_index(
    const MoveStructure<ColumnsT>& ms,
    typename MoveStructure<ColumnsT>::Position pos
) {
    if constexpr (MoveColsTraits<ColumnsT>::RELATIVE) {
        // Convert (interval, offset) to absolute index via prefix sums of lengths.
        ulint idx = 0;
        for (size_t i = 0; i < pos.interval; ++i) {
            idx += ms.get_length(i);
        }
        idx += pos.offset;
        return idx;
    } else {
        return pos.idx;
    }
}

static void test_move_structure_move_matches_interval_permutation_relative() {
    const vector<ulint> lengths = {3, 2, 1, 2, 2};
    const vector<ulint> perm = {4, 0, 9, 2, 7};
    const ulint domain = 10;

    using MS = MoveStructure<MoveCols>;
    MS ms(lengths, perm, NO_SPLITTING);
    using Position = typename MS::Position;

    // For each interval and offset, move once and check absolute index.
    ulint start_prefix = 0;
    for (size_t j = 0; j < lengths.size(); ++j) {
        for (ulint o = 0; o < lengths[j]; ++o) {
            Position pos{static_cast<ulint>(j), o};
            pos = ms.move(pos);
            const ulint idx = compute_global_index(ms, pos);
            const ulint expected = perm[j] + o;
            assert(idx == expected);
        }
        start_prefix += lengths[j];
    }
}

static void test_move_structure_move_matches_interval_permutation_absolute() {
    const vector<ulint> lengths = {3, 2, 1, 2, 2};
    const vector<ulint> perm = {4, 0, 9, 2, 7};
    const ulint domain = 10;

    using MS = MoveStructure<MoveColsIdx>;
    MS ms(lengths, perm, NO_SPLITTING);
    using Position = typename MS::Position;

    for (size_t j = 0; j < lengths.size(); ++j) {
        for (ulint o = 0; o < lengths[j]; ++o) {
            // Absolute position: idx is the source index.
            Position pos{static_cast<ulint>(j), o, 0};
            pos.idx = 0;
            // Reconstruct idx from (interval, offset).
            for (size_t i = 0; i < j; ++i) {
                pos.idx += ms.get_length(i);
            }
            pos.idx += o;

            pos = ms.move(pos);
            const ulint idx = compute_global_index(ms, pos);
            const ulint expected = perm[j] + o;
            assert(idx == expected);
        }
    }
}

static void test_move_structure_move_exponential_agrees_with_move_absolute() {
    const vector<ulint> lengths = {3, 2, 1, 2, 2};
    const vector<ulint> perm = {4, 0, 9, 2, 7};
    const ulint domain = 10;

    using MS = MoveStructure<MoveColsIdx>;
    MS ms(lengths, perm, NO_SPLITTING);
    using Position = typename MS::Position;

    for (size_t j = 0; j < lengths.size(); ++j) {
        for (ulint o = 0; o < lengths[j]; ++o) {
            Position pos{static_cast<ulint>(j), o, 0};
            for (size_t i = 0; i < j; ++i) {
                pos.idx += ms.get_length(i);
            }
            pos.idx += o;

            Position p1 = ms.move(pos);
            Position p2 = ms.move_exponential(pos);

            assert(p1.interval == p2.interval);
            assert(p1.offset == p2.offset);
            assert(p1.idx == p2.idx);
        }
    }
}

static void test_move_structure_splitting_preserves_semantics_relative() {
    // Same example, but force length capping so some runs are split.
    const vector<ulint> lengths = {3, 2, 1, 2, 2};
    const vector<ulint> perm = {4, 0, 9, 2, 7};
    const ulint domain = 10;

    using MS = MoveStructure<MoveCols>;
    MS base_ms(lengths, perm, NO_SPLITTING);

    // Use only length capping with a small factor to trigger splitting.
    SplitParams params = ONLY_LENGTH_CAPPING;
    MS capped_ms(lengths, perm, params);

    using Position = typename MS::Position;

    // Check semantics for every source index in [0, domain).
    for (ulint idx = 0; idx < domain; ++idx) {
        // Find (interval, offset) in base_ms.
        Position pos_base;
        {
            ulint prefix = 0;
            for (size_t j = 0; j < base_ms.runs(); ++j) {
                const ulint len = base_ms.get_length(j);
                if (idx < prefix + len) {
                    pos_base.interval = j;
                    pos_base.offset = idx - prefix;
                    break;
                }
                prefix += len;
            }
        }

        // Find (interval, offset) in capped_ms.
        Position pos_cap;
        {
            ulint prefix = 0;
            for (size_t j = 0; j < capped_ms.runs(); ++j) {
                const ulint len = capped_ms.get_length(j);
                if (idx < prefix + len) {
                    pos_cap.interval = j;
                    pos_cap.offset = idx - prefix;
                    break;
                }
                prefix += len;
            }
        }

        pos_base = base_ms.move(pos_base);
        pos_cap = capped_ms.move(pos_cap);

        const ulint idx_base = compute_global_index(base_ms, pos_base);
        const ulint idx_cap = compute_global_index(capped_ms, pos_cap);
        assert(idx_base == idx_cap);
    }
}

int main() {
    test_move_structure_move_matches_interval_permutation_relative();
    test_move_structure_move_matches_interval_permutation_absolute();
    test_move_structure_move_exponential_agrees_with_move_absolute();
    test_move_structure_splitting_preserves_semantics_relative();

    std::cout << "move_structure integration tests passed" << std::endl;
    return 0;
}
