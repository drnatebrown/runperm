// Unit tests for RLBWTRow and its traits.
// Simple assert-based tests, no external framework.

#include "internal/rlbwt/specializations/rlbwt_row.hpp"

#include <cassert>
#include <iostream>
#include <vector>

using std::size_t;

static void test_rlbwt_row_set_get_default() {
    using Row = RLBWTRow<RLBWTCols>;
    using RowTraits = Row::RowTraits;
    using ColsTraits = MoveColsTraits<RLBWTCols>;

    constexpr size_t NumCols = static_cast<size_t>(RLBWTCols::COUNT);

    Row row;

    // Values chosen within bit-width limits (use ulint to avoid overflow).
    ulint length    = (ulint(1) << (RowTraits::PRIMARY_BITS   - 1)) - 1;
    ulint pointer   = (ulint(1) << (RowTraits::POINTER_BITS   - 1)) - 1;
    ulint offset    = (ulint(1) << (RowTraits::OFFSET_BITS    - 1)) - 1;
    ulint character = (ulint(1) << (RowTraits::CHARACTER_BITS - 1)) - 1;

    row.set<ColsTraits::PRIMARY>(length);
    row.set<ColsTraits::POINTER>(pointer);
    row.set<ColsTraits::OFFSET>(offset);
    row.set<ColsTraits::CHARACTER>(character);

    assert(row.get<ColsTraits::PRIMARY>() == length);
    assert(row.get<ColsTraits::POINTER>() == pointer);
    assert(row.get<ColsTraits::OFFSET>() == offset);
    assert(row.get<ColsTraits::CHARACTER>() == character);

    // Round-trip via array-based API.
    std::array<ulint, NumCols> values{};
    values[static_cast<size_t>(ColsTraits::PRIMARY)]   = length;
    values[static_cast<size_t>(ColsTraits::POINTER)]   = pointer;
    values[static_cast<size_t>(ColsTraits::OFFSET)]    = offset;
    values[static_cast<size_t>(ColsTraits::CHARACTER)] = character;

    Row row2(values);
    auto roundtrip = row2.get();

    assert(roundtrip[static_cast<size_t>(ColsTraits::PRIMARY)]   == length);
    assert(roundtrip[static_cast<size_t>(ColsTraits::POINTER)]   == pointer);
    assert(roundtrip[static_cast<size_t>(ColsTraits::OFFSET)]    == offset);
    assert(roundtrip[static_cast<size_t>(ColsTraits::CHARACTER)] == character);
}

static void test_rlbwt_row_assert_widths() {
    using Row = RLBWTRow<RLBWTCols>;
    using RowTraits = Row::RowTraits;
    using ColsTraits = MoveColsTraits<RLBWTCols>;

    constexpr size_t NumCols = static_cast<size_t>(RLBWTCols::COUNT);
    std::array<uchar, NumCols> widths{};

    widths[static_cast<size_t>(ColsTraits::PRIMARY)]   = static_cast<uchar>(RowTraits::PRIMARY_BITS);
    widths[static_cast<size_t>(ColsTraits::POINTER)]   = static_cast<uchar>(RowTraits::POINTER_BITS);
    widths[static_cast<size_t>(ColsTraits::OFFSET)]    = static_cast<uchar>(RowTraits::OFFSET_BITS);
    widths[static_cast<size_t>(ColsTraits::CHARACTER)] = static_cast<uchar>(RowTraits::CHARACTER_BITS);

    Row::assert_widths(widths);
}

static void test_rlbwt_row_traits_aliases() {
    // Default traits for RLBWTCols should alias the configured default bit-layout.
    using DefaultBits = RLBWTColsDefault;
    using TraitsFromBits = MoveRowTraits<DefaultBits>;
    using TraitsFromCols = MoveRowTraits<RLBWTCols>;

    static_assert(TraitsFromBits::CHARACTER_BITS == TraitsFromCols::CHARACTER_BITS, "CHARACTER_BITS mismatch");
    static_assert(TraitsFromBits::PRIMARY_BITS   == TraitsFromCols::PRIMARY_BITS,   "PRIMARY_BITS mismatch");
    static_assert(TraitsFromBits::POINTER_BITS   == TraitsFromCols::POINTER_BITS,   "POINTER_BITS mismatch");
    static_assert(TraitsFromBits::OFFSET_BITS    == TraitsFromCols::OFFSET_BITS,    "OFFSET_BITS mismatch");
}

int main() {
    test_rlbwt_row_set_get_default();
    test_rlbwt_row_assert_widths();
    test_rlbwt_row_traits_aliases();

    std::cout << "rlbwt_row unit tests passed" << std::endl;
    return 0;
}

