// Unit tests for RLBWTMoveStructure (RLBWT-specific MoveStructure wrapper).
// Simple assert-based tests, no external framework.

#include "internal/rlbwt/specializations/rlbwt_structure.hpp"

#include <cassert>
#include <iostream>
#include <vector>

using std::size_t;
using std::vector;

static void test_rlbwt_move_structure_relative_chars_and_widths() {
    // Small RLBWT example with 4 runs.
    // Characters are already mapped into [0, sigma) as RLBWTMoveStructure expects.
    const vector<uchar> rlbwt_chars = {0, 1, 2, 3};
    const vector<ulint> lengths     = {1, 2, 1, 2};
    const vector<ulint> interval_perm = {0, 1, 3, 5};
    const ulint domain = 6;
    const uchar sigma = 4; // Alphabet size (A,C,G,T).

    using MS = RLBWTMoveStructure<RLBWTCols>;
    MS ms(rlbwt_chars, lengths, interval_perm, domain, sigma, NO_SPLITTING);

    assert(ms.domain() == domain);
    assert(ms.runs() == lengths.size());

    auto widths = ms.get_widths();
    using ColsTraits = MoveColsTraits<RLBWTCols>;

    uchar w_char = widths[static_cast<size_t>(ColsTraits::CHARACTER)];
    // Character width should be enough to encode sigma-1.
    assert(w_char == bit_width(static_cast<ulint>(sigma - 1)));

    // Row-level characters should match run heads.
    for (size_t i = 0; i < lengths.size(); ++i) {
        assert(ms.get_character(i) == rlbwt_chars[i]);
    }

    // Position-based accessor should agree with interval-based accessor.
    typename ColsTraits::Position pos{};
    pos.interval = 1;
    pos.offset = 0;
    assert(ms.get_character(pos) == ms.get_character(pos.interval));
}

static void test_rlbwt_move_structure_relative_splitting_preserves_chars() {
    // Two runs: a long run (symbol 0) followed by a short run (symbol 1).
    const vector<uchar> rlbwt_chars = {0, 1};
    const vector<ulint> lengths     = {5, 1};
    const vector<ulint> interval_perm = {0, 5};
    const ulint domain = 6;
    const uchar sigma = 4;

    using MS = RLBWTMoveStructure<RLBWTCols>;

    SplitParams split;
    split.length_capping = 1.0; // force aggressive splitting if beneficial
    MS ms(rlbwt_chars, lengths, interval_perm, domain, sigma, split);

    assert(ms.domain() == domain);
    assert(ms.runs() >= lengths.size()); // first run may be split into multiple

    // Sum up lengths per character to ensure splitting preserved run heads.
    ulint sym0_len = 0;
    ulint sym1_len = 0;
    for (ulint i = 0; i < ms.intervals(); ++i) {
        uchar ch = ms.get_character(i);
        if (ch == 0) {
            sym0_len += ms.get_length(i);
        } else if (ch == 1) {
            sym1_len += ms.get_length(i);
        } else {
            assert(false && "Unexpected character in RLBWTMoveStructure");
        }
    }
    assert(sym0_len == 5);
    assert(sym1_len == 1);
}

static void test_rlbwt_move_structure_absolute_chars_and_widths() {
    // Same example as relative, but using absolute columns.
    const vector<uchar> rlbwt_chars = {0, 1, 2, 3};
    const vector<ulint> lengths     = {1, 2, 1, 2};
    const vector<ulint> interval_perm = {0, 1, 3, 5};
    const ulint domain = 6;
    const uchar sigma = 4;

    using MS = RLBWTMoveStructure<RLBWTColsIdx>;
    MS ms(rlbwt_chars, lengths, interval_perm, domain, sigma, NO_SPLITTING);

    assert(ms.domain() == domain);
    assert(ms.runs() == lengths.size());

    auto widths = ms.get_widths();
    using ColsTraits = MoveColsTraits<RLBWTColsIdx>;

    uchar w_char = widths[static_cast<size_t>(ColsTraits::CHARACTER)];
    assert(w_char == bit_width(static_cast<ulint>(sigma - 1)));

    for (size_t i = 0; i < lengths.size(); ++i) {
        assert(ms.get_character(i) == rlbwt_chars[i]);
    }
}

int main() {
    test_rlbwt_move_structure_relative_chars_and_widths();
    test_rlbwt_move_structure_relative_splitting_preserves_chars();
    test_rlbwt_move_structure_absolute_chars_and_widths();

    std::cout << "rlbwt_structure unit tests passed" << std::endl;
    return 0;
}
