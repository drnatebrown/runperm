// Unit tests for MoveStructure build-time invariants.
// These are simple assert-based tests, no external framework.

#include "internal/move/move_structure.hpp"

#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>

using std::size_t;
using std::vector;

static vector<ulint> compute_starts(const vector<ulint>& lengths) {
    vector<ulint> starts(lengths.size() + 1, 0);
    for (size_t i = 0; i < lengths.size(); ++i) {
        starts[i + 1] = starts[i] + lengths[i];
    }
    return starts;
}

static size_t interval_of(const vector<ulint>& starts, ulint idx) {
    // returns i such that starts[i] <= idx < starts[i+1]
    for (size_t i = 0; i + 1 < starts.size(); ++i) {
        if (starts[i] <= idx && idx < starts[i + 1]) return i;
    }
    assert(false && "idx not in any interval");
    return 0;
}

template <typename ColumnsT>
static void assert_pointer_offset_correct(
    const MoveStructure<ColumnsT>& ms,
    const vector<ulint>& lengths,
    const vector<ulint>& interval_permutation
) {
    const auto starts = compute_starts(lengths);
    assert(starts.back() == ms.domain());
    assert(ms.runs() == lengths.size());

    for (size_t j = 0; j < lengths.size(); ++j) {
        const ulint mapped_start = interval_permutation[j];
        const size_t k = interval_of(starts, mapped_start);
        const ulint expected_offset = mapped_start - starts[k];

        assert(ms.get_pointer(j) == k);
        assert(ms.get_offset(j) == expected_offset);
    }
}

static void test_move_structure_relative_build_invariants() {
    // A small valid run-decomposition of a permutation of [0..9]:
    // π = [4,5,6, 0,1, 9, 2,3, 7,8]
    const vector<ulint> lengths = {3, 2, 1, 2, 2};
    const vector<ulint> perm = {4, 0, 9, 2, 7};
    const ulint domain = 10;

    MoveStructure<MoveCols> ms(lengths, perm, NO_SPLITTING);
    assert(ms.domain() == domain);
    assert(ms.runs() == lengths.size());

    for (size_t i = 0; i < lengths.size(); ++i) {
        assert(ms.get_length(i) == lengths[i]);
    }

    assert_pointer_offset_correct(ms, lengths, perm);
}

static void test_move_structure_absolute_build_invariants() {
    // Same run-decomposition as above, but using absolute columns.
    const vector<ulint> lengths = {3, 2, 1, 2, 2};
    const vector<ulint> perm = {4, 0, 9, 2, 7};
    const ulint domain = 10;

    MoveStructure<MoveColsIdx> ms(lengths, perm, NO_SPLITTING);
    assert(ms.domain() == domain);
    assert(ms.runs() == lengths.size());

    const auto starts = compute_starts(lengths);
    for (size_t i = 0; i < lengths.size(); ++i) {
        assert(ms.get_start(i) == starts[i]);
        assert(ms.get_length(i) == lengths[i]);
    }
    assert(ms.get_start(lengths.size()) == domain);

    assert_pointer_offset_correct(ms, lengths, perm);
}

static void test_move_structure_serialize_roundtrip() {
    const vector<ulint> lengths = {3, 2, 1, 2, 2};
    const vector<ulint> perm = {4, 0, 9, 2, 7};
    const ulint domain = 10;

    MoveStructure<MoveColsIdx> ms(lengths, perm, NO_SPLITTING);

    std::stringstream ss;
    const size_t bytes = ms.serialize(ss);
    assert(bytes > 0);

    MoveStructure<MoveColsIdx> loaded;
    loaded.load(ss);

    assert(loaded.domain() == ms.domain());
    assert(loaded.runs() == ms.runs());

    for (size_t i = 0; i < lengths.size(); ++i) {
        assert(loaded.get_start(i) == ms.get_start(i));
        assert(loaded.get_pointer(i) == ms.get_pointer(i));
        assert(loaded.get_offset(i) == ms.get_offset(i));
    }
}

static void test_move_structure_widths_relative_and_absolute_with_and_without_splitting() {
    // Base example from earlier tests.
    const vector<ulint> lengths = {3, 2, 1, 2, 2};
    const vector<ulint> perm = {4, 0, 9, 2, 7};
    const ulint domain = 10;

    // Relative, no splitting.
    MoveStructure<MoveCols> ms_rel(lengths, perm, NO_SPLITTING);
    auto widths_rel = ms_rel.get_widths();

    // PRIMARY and OFFSET widths should match for relative representation.
    uchar w_primary_rel = widths_rel[static_cast<size_t>(MoveColsTraits<MoveCols>::PRIMARY)];
    uchar w_pointer_rel = widths_rel[static_cast<size_t>(MoveColsTraits<MoveCols>::POINTER)];
    uchar w_offset_rel  = widths_rel[static_cast<size_t>(MoveColsTraits<MoveCols>::OFFSET)];
    assert(w_primary_rel == w_offset_rel);
    // Pointer width must be enough to index all runs.
    assert(w_pointer_rel >= bit_width(ms_rel.runs()));

    // Relative, with length-capping splitting.
    SplitParams split = ONLY_LENGTH_CAPPING;
    MoveStructure<MoveCols> ms_rel_split(lengths, perm, split);
    auto widths_rel_split = ms_rel_split.get_widths();

    uchar w_primary_rel_split =
        widths_rel_split[static_cast<size_t>(MoveColsTraits<MoveCols>::PRIMARY)];
    uchar w_pointer_rel_split =
        widths_rel_split[static_cast<size_t>(MoveColsTraits<MoveCols>::POINTER)];
    uchar w_offset_rel_split =
        widths_rel_split[static_cast<size_t>(MoveColsTraits<MoveCols>::OFFSET)];

    // Still must have PRIMARY == OFFSET after splitting.
    assert(w_primary_rel_split == w_offset_rel_split);
    // Pointer width must match the (possibly increased) number of runs.
    assert(w_pointer_rel_split >= bit_width(ms_rel_split.runs()));

    // Absolute, no splitting.
    MoveStructure<MoveColsIdx> ms_abs(lengths, perm, NO_SPLITTING);
    auto widths_abs = ms_abs.get_widths();

    uchar w_primary_abs =
        widths_abs[static_cast<size_t>(MoveColsTraits<MoveColsIdx>::PRIMARY)];
    uchar w_pointer_abs =
        widths_abs[static_cast<size_t>(MoveColsTraits<MoveColsIdx>::POINTER)];
    uchar w_offset_abs =
        widths_abs[static_cast<size_t>(MoveColsTraits<MoveColsIdx>::OFFSET)];

    // PRIMARY width should be enough for domain indices.
    assert(w_primary_abs >= bit_width(domain));
    // Pointer width enough for runs, offset width reused for relative length.
    assert(w_pointer_abs >= bit_width(ms_abs.runs()));
    assert(w_offset_abs >= bit_width(ms_abs.get_length(0)));
}

int main() {
    test_move_structure_relative_build_invariants();
    test_move_structure_absolute_build_invariants();
    test_move_structure_serialize_roundtrip();
    test_move_structure_widths_relative_and_absolute_with_and_without_splitting();

    std::cout << "move_structure unit tests passed" << std::endl;
    return 0;
}
