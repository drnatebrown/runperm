// Unit tests for RLBWTMoveStructure (RLBWT-specific MoveStructure wrapper).
// Simple assert-based tests, no external framework.

#include "orbit/internal/move/move_structure_impl.hpp"
#include "orbit/internal/rlbwt/specializations/rlbwt_structure.hpp"
#include "orbit/rlbwt.hpp"

#include <cassert>
#include <iostream>
#include <type_traits>
#include <vector>

using std::size_t;
using std::vector;

using namespace orbit;
using namespace orbit::rlbwt;

static void test_rlbwt_move_structure_relative_chars_and_widths() {
    // Small RLBWT example with 4 runs.
    // Characters are already mapped into [0, sigma) as RLBWTMoveStructure expects.
    const vector<uchar> rlbwt_chars = {0, 1, 2, 3};
    const vector<ulint> lengths     = {1, 2, 1, 2};
    const vector<ulint> interval_perm = {0, 1, 3, 5};
    const ulint domain = 6;
    const uchar sigma = 4; // alphabet size (A,C,G,T).

    using MS = rlbwt_move_structure<rlbwt_columns>;
    MS ms(rlbwt_chars, lengths, interval_perm, domain, sigma, NO_SPLITTING);

    assert(ms.domain() == domain);
    assert(ms.runs() == lengths.size());

    auto widths = ms.get_widths();
    using ColsTraits = move_cols_traits<rlbwt_columns>;

    uchar w_char = widths[static_cast<size_t>(ColsTraits::CHARACTER)];
    // Character width should be enough to encode sigma-1.
    assert(w_char == bit_width(static_cast<ulint>(sigma - 1)));

    // Row-level characters should match run heads.
    for (size_t i = 0; i < lengths.size(); ++i) {
        assert(ms.get_character(i) == rlbwt_chars[i]);
    }

    // position-based accessor should agree with interval-based accessor.
    typename ColsTraits::position pos{};
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

    using MS = rlbwt_move_structure<rlbwt_columns>;

    split_params split;
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

// Known-valid RLBWT example (see integration tests).
static const vector<uchar> kBwtHeads =       {'T','C','G','A','T', 1 ,'A','T','A'};
static const vector<ulint> kBwtRunLengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };

using rlbwt_invertible_structure = rlbwt_move_structure<rlbwt_invertible_columns>;
using rlbwt_invertible_structure_idx = rlbwt_move_structure<rlbwt_invertible_columns_idx>;

static_assert(std::is_same_v<
    move_structure<rlbwt_invertible_columns, move_vector>,
    invertible_structure_impl<rlbwt_invertible_columns, move_vector>>);

template <typename MS>
static ulint global_index_relative(const MS& ms, typename MS::position pos) {
    ulint idx = 0;
    for (size_t i = 0; i < pos.interval; ++i) {
        idx += ms.get_length(i);
    }
    return idx + pos.offset;
}

template <typename MS>
static typename MS::position position_from_index_relative(const MS& ms, ulint idx) {
    ulint prefix = 0;
    for (ulint interval = 0; interval < ms.intervals(); ++interval) {
        const ulint len = ms.get_length(interval);
        if (idx < prefix + len) {
            return {interval, idx - prefix};
        }
        prefix += len;
    }
    assert(false && "index out of range");
    return {};
}

template <typename MS>
static typename MS::position position_from_index_absolute(const MS& ms, ulint idx) {
    typename MS::position pos{};
    pos.idx = idx;
    ulint prefix = 0;
    for (ulint interval = 0; interval < ms.intervals(); ++interval) {
        const ulint len = ms.get_length(interval);
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

static void test_rlbwt_invertible_structure_build_from_lf_encoding() {
    const auto enc = invertible_rlbwt_interval_encoding<>::lf_interval_encoding(
        kBwtHeads, kBwtRunLengths, NO_SPLITTING);
    rlbwt_invertible_structure ms(enc);

    assert(ms.domain() == enc.domain());
    assert(ms.intervals() == enc.intervals());
    assert(ms.get_file_extension() == ".imove");
    assert(ms.get_fwd_interval(0));

    const auto& alpha = enc.get_alphabet();
    const auto& heads = enc.get_heads();
    for (size_t i = 0; i < enc.intervals(); ++i) {
        assert(ms.get_length(i) == enc.get_length(i));
        assert(ms.get_fwd_interval(i) == enc.get_is_fwd_interval(i));
        assert(ms.get_inv_interval(i) == enc.get_is_inv_interval(i));
        assert(alpha.unmap_char(static_cast<uchar>(ms.get_character(i)))
               == alpha.unmap_char(static_cast<uchar>(heads[i])));
    }
}

static void test_rlbwt_invertible_structure_move_fwd_matches_lf_move() {
    const auto enc = invertible_rlbwt_interval_encoding<>::lf_interval_encoding(
        kBwtHeads, kBwtRunLengths, NO_SPLITTING);
    rlbwt_invertible_structure ms(enc);
    lf_move<> ref(kBwtHeads, kBwtRunLengths, NO_SPLITTING);

    assert(ms.domain() == ref.domain());
    for (ulint idx = 0; idx < ms.domain(); ++idx) {
        auto pos = position_from_index_relative(ms, idx);
        auto ref_pos = position_from_index_relative(ref, idx);
        pos = ms.move_fwd(pos);
        ref_pos = ref.next(ref_pos);
        assert(global_index_relative(ms, pos) == global_index_relative(ref, ref_pos));
    }
}

static void test_rlbwt_invertible_structure_fwd_inv_roundtrip() {
    const auto enc = invertible_rlbwt_interval_encoding<>::lf_interval_encoding(
        kBwtHeads, kBwtRunLengths, NO_SPLITTING);
    rlbwt_invertible_structure ms(enc);

    for (ulint idx = 0; idx < ms.domain(); ++idx) {
        auto pos = position_from_index_relative(ms, idx);
        pos = ms.move_fwd(pos);
        pos = ms.move_inv(pos);
        assert(global_index_relative(ms, pos) == idx);
    }
}

static void test_rlbwt_invertible_structure_absolute_move_fwd() {
    const auto enc = invertible_rlbwt_interval_encoding<>::lf_interval_encoding(
        kBwtHeads, kBwtRunLengths, NO_SPLITTING);
    rlbwt_invertible_structure_idx ms(enc);
    lf_move<true> ref(kBwtHeads, kBwtRunLengths, NO_SPLITTING);

    for (ulint idx = 0; idx < ms.domain(); ++idx) {
        auto pos = position_from_index_absolute(ms, idx);
        auto ref_pos = position_from_index_absolute(ref, idx);
        pos = ms.move_fwd(pos);
        ref_pos = ref.next(ref_pos);
        assert(pos.idx == ref_pos.idx);
    }
}

static void test_rlbwt_invertible_structure_widths() {
    const auto enc = invertible_rlbwt_interval_encoding<>::lf_interval_encoding(
        kBwtHeads, kBwtRunLengths, NO_SPLITTING);
    rlbwt_invertible_structure ms(enc);

    const auto widths = ms.get_widths();
    using ColsTraits = move_cols_traits<rlbwt_invertible_columns>;
    assert(widths[static_cast<size_t>(ColsTraits::POINTER_FWD)] >= bit_width(ms.intervals()));
    assert(widths[static_cast<size_t>(ColsTraits::POINTER_INV)] >= bit_width(ms.intervals()));
    assert(widths[static_cast<size_t>(ColsTraits::FWD_INTERVAL)] == 1);
    assert(widths[static_cast<size_t>(ColsTraits::INV_INTERVAL)] == 1);
    assert(widths[static_cast<size_t>(ColsTraits::CHARACTER)] == bit_width(enc.sigma() - 1));
}

static void test_rlbwt_move_structure_absolute_chars_and_widths() {
    // Same example as relative, but using absolute columns.
    const vector<uchar> rlbwt_chars = {0, 1, 2, 3};
    const vector<ulint> lengths     = {1, 2, 1, 2};
    const vector<ulint> interval_perm = {0, 1, 3, 5};
    const ulint domain = 6;
    const uchar sigma = 4;

    using MS = rlbwt_move_structure<rlbwt_columns_idx>;
    MS ms(rlbwt_chars, lengths, interval_perm, domain, sigma, NO_SPLITTING);

    assert(ms.domain() == domain);
    assert(ms.runs() == lengths.size());

    auto widths = ms.get_widths();
    using ColsTraits = move_cols_traits<rlbwt_columns_idx>;

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
    test_rlbwt_invertible_structure_build_from_lf_encoding();
    test_rlbwt_invertible_structure_move_fwd_matches_lf_move();
    test_rlbwt_invertible_structure_fwd_inv_roundtrip();
    test_rlbwt_invertible_structure_absolute_move_fwd();
    test_rlbwt_invertible_structure_widths();

    std::cout << "rlbwt_structure unit tests passed" << std::endl;
    return 0;
}
