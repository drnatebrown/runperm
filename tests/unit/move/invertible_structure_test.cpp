// Unit tests for invertible move structures (invertible_structure_impl).
// Simple assert-based tests, no external framework.

#include "orbit/internal/move/move_structure_impl.hpp"
#include "orbit/internal/perm/permutation_impl.hpp"
#include "orbit/interval_encoding.hpp"
#include "orbit/move_structure.hpp"

#include <cassert>
#include <iostream>
#include <sstream>
#include <type_traits>
#include <vector>

using std::size_t;
using std::vector;

using namespace orbit;

static const vector<ulint> kRunnyPerm = {1, 2, 9, 10, 11, 3, 12, 13, 4, 5, 14, 0, 15, 6, 7, 8};
static const vector<ulint> kSmallPerm = {1, 0, 3, 2};

static_assert(std::is_same_v<
    move_structure<invertible_columns, move_vector>,
    invertible_structure_impl<invertible_columns, move_vector>>);
static_assert(std::is_same_v<
    move_structure<move_columns, move_vector>,
    move_structure_impl<move_columns, move_vector>>);

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

static invertible_interval_encoding make_invertible_encoding(const vector<ulint>& permutation,
                                                               const split_params& sp = NO_SPLITTING) {
    return invertible_interval_encoding::from_permutation(permutation, sp);
}

static void assert_build_matches_encoding(
    const invertible_structure_vec& ms,
    const invertible_interval_encoding& enc
) {
    assert(ms.domain() == enc.domain());
    assert(ms.intervals() == enc.intervals());
    for (size_t i = 0; i < enc.intervals(); ++i) {
        assert(ms.get_length(i) == enc.get_length(i));
        assert(ms.get_fwd_interval(i) == enc.get_is_fwd_interval(i));
        assert(ms.get_inv_interval(i) == enc.get_is_inv_interval(i));
    }
    assert(ms.get_fwd_interval(0));
}

// Mirrors populate_structure pointer assignments in invertible_structure_impl.
static void assert_pointers_match_encoding(
    const invertible_structure_vec& ms,
    const invertible_interval_encoding& enc
) {
    size_t input_start_val = 0;
    size_t output_start_val = 0;
    size_t img_rank_inv_idx = 0;
    for (size_t i = 0; i < enc.intervals(); ++i) {
        const size_t length = enc.get_length(i);
        while (img_rank_inv_idx < enc.intervals() && output_start_val < input_start_val + length) {
            const ulint img_rank_inv_val = enc.get_img_rank_inv(img_rank_inv_idx);
            if (enc.get_is_fwd_interval(img_rank_inv_val) && output_start_val == input_start_val) {
                assert(ms.get_pointer_fwd(img_rank_inv_val) == i);
            }
            if (enc.get_is_inv_interval(i) && output_start_val == input_start_val) {
                assert(ms.get_pointer_inv(i) == img_rank_inv_val);
            }
            output_start_val += enc.get_length(img_rank_inv_val);
            ++img_rank_inv_idx;
        }
        input_start_val += length;
    }
}

template <typename MS>
static void assert_move_fwd_matches_perm(const MS& ms, const vector<ulint>& permutation) {
    for (ulint idx = 0; idx < ms.domain(); ++idx) {
        auto pos = position_from_index_relative(ms, idx);
        pos = ms.move_fwd(pos);
        assert(global_index_relative(ms, pos) == permutation[idx]);
    }
}

template <typename MS>
static void assert_move_inv_matches_inverse(const MS& ms, const vector<ulint>& permutation) {
    const vector<ulint> inverse = get_inverse_permutation(permutation);
    for (ulint idx = 0; idx < ms.domain(); ++idx) {
        auto pos = position_from_index_relative(ms, idx);
        pos = ms.move_inv(pos);
        assert(global_index_relative(ms, pos) == inverse[idx]);
    }
}

template <typename MS>
static void assert_move_fwd_matches_perm_absolute(const MS& ms, const vector<ulint>& permutation) {
    for (ulint idx = 0; idx < ms.domain(); ++idx) {
        auto pos = position_from_index_absolute(ms, idx);
        pos = ms.move_fwd(pos);
        assert(pos.idx == permutation[idx]);
    }
}

template <typename MS>
static void assert_move_inv_matches_inverse_absolute(const MS& ms, const vector<ulint>& permutation) {
    const vector<ulint> inverse = get_inverse_permutation(permutation);
    for (ulint idx = 0; idx < ms.domain(); ++idx) {
        auto pos = position_from_index_absolute(ms, idx);
        pos = ms.move_inv(pos);
        assert(pos.idx == inverse[idx]);
    }
}

static void test_invertible_structure_build_from_encoding() {
    const auto enc = make_invertible_encoding(kRunnyPerm);
    invertible_structure_vec ms(enc);
    assert_build_matches_encoding(ms, enc);
    assert_pointers_match_encoding(ms, enc);
    assert(ms.get_file_extension() == ".imove");
}

static void test_invertible_structure_build_from_lengths_and_images() {
    const vector<ulint> lengths = {3, 2, 1, 2, 2};
    const vector<ulint> images = {4, 0, 9, 2, 7};
    invertible_structure_vec ms(lengths, images, NO_SPLITTING);
    assert(ms.domain() == 10);
    assert_move_fwd_matches_perm(ms, {4, 5, 6, 0, 1, 9, 2, 3, 7, 8});
}

static void test_invertible_structure_move_fwd_matches_permutation() {
    const auto enc = make_invertible_encoding(kRunnyPerm);
    invertible_structure_vec ms(enc);
    assert_move_fwd_matches_perm(ms, kRunnyPerm);
}

static void test_invertible_structure_move_inv_matches_inverse() {
    const auto enc = make_invertible_encoding(kRunnyPerm);
    invertible_structure_vec ms(enc);
    assert_move_inv_matches_inverse(ms, kRunnyPerm);
}

static void test_invertible_structure_fwd_inv_roundtrip() {
    const auto enc = make_invertible_encoding(kRunnyPerm);
    invertible_structure_vec ms(enc);

    for (ulint idx = 0; idx < ms.domain(); ++idx) {
        auto pos = position_from_index_relative(ms, idx);
        pos = ms.move_fwd(pos);
        pos = ms.move_inv(pos);
        assert(global_index_relative(ms, pos) == idx);
    }
}

static void test_invertible_structure_inv_fwd_roundtrip() {
    const auto enc = make_invertible_encoding(kRunnyPerm);
    invertible_structure_vec ms(enc);

    for (ulint idx = 0; idx < ms.domain(); ++idx) {
        auto pos = position_from_index_relative(ms, idx);
        pos = ms.move_inv(pos);
        pos = ms.move_fwd(pos);
        assert(global_index_relative(ms, pos) == idx);
    }
}

static void test_invertible_structure_move_from_non_fwd_interval() {
    const auto enc = make_invertible_encoding(kRunnyPerm);
    invertible_structure_vec ms(enc);

    // interval 1 is inv-only for the runny permutation (see interval_encoding_test).
    assert(!ms.get_fwd_interval(1));
    assert(ms.get_inv_interval(1));

    ulint prefix = 0;
    for (ulint interval = 0; interval < ms.intervals(); ++interval) {
        if (interval == 1) {
            for (ulint o = 0; o < ms.get_length(interval); ++o) {
                const ulint idx = prefix + o;
                auto pos = position_from_index_relative(ms, idx);
                assert(pos.interval == 1);
                pos = ms.move_fwd(pos);
                assert(global_index_relative(ms, pos) == kRunnyPerm[idx]);
            }
            break;
        }
        prefix += ms.get_length(interval);
    }
}

static void test_invertible_structure_absolute_semantics() {
    const auto enc = make_invertible_encoding(kRunnyPerm);
    invertible_structure_vec_idx ms(enc);
    assert_move_fwd_matches_perm_absolute(ms, kRunnyPerm);
    assert_move_inv_matches_inverse_absolute(ms, kRunnyPerm);
}

static void test_invertible_structure_absolute_move_exponential() {
    const auto enc = make_invertible_encoding(kRunnyPerm);
    const vector<ulint> inverse = get_inverse_permutation(kRunnyPerm);
    invertible_structure_vec_idx ms(enc);

    for (ulint idx = 0; idx < ms.domain(); ++idx) {
        auto pos_fwd = position_from_index_absolute(ms, idx);
        auto pos_lin = ms.move_fwd(pos_fwd);
        auto pos_exp = ms.move_fwd_exponential(pos_fwd);
        assert(pos_lin.interval == pos_exp.interval);
        assert(pos_lin.offset == pos_exp.offset);
        assert(pos_lin.idx == pos_exp.idx);
        assert(pos_lin.idx == kRunnyPerm[idx]);

        auto pos_inv = position_from_index_absolute(ms, idx);
        pos_lin = ms.move_inv(pos_inv);
        pos_exp = ms.move_inv_exponential(pos_inv);
        assert(pos_lin.interval == pos_exp.interval);
        assert(pos_lin.offset == pos_exp.offset);
        assert(pos_lin.idx == pos_exp.idx);
        assert(pos_lin.idx == inverse[idx]);
    }
}

static void test_invertible_structure_small_and_identity_permutations() {
    {
        const auto enc = make_invertible_encoding(kSmallPerm);
        invertible_structure_vec ms(enc);
        assert_pointers_match_encoding(ms, enc);
        assert_move_fwd_matches_perm(ms, kSmallPerm);
        assert_move_inv_matches_inverse(ms, kSmallPerm);
    }
    {
        vector<ulint> identity(8);
        for (ulint i = 0; i < identity.size(); ++i) {
            identity[i] = i;
        }
        const auto enc = make_invertible_encoding(identity);
        invertible_structure_vec ms(enc);
        assert_move_fwd_matches_perm(ms, identity);
        assert_move_inv_matches_inverse(ms, identity);
    }
}

static void test_invertible_structure_serialize_roundtrip() {
    const auto enc = make_invertible_encoding(kRunnyPerm);
    invertible_structure_vec ms(enc);

    std::stringstream ss;
    const size_t bytes = ms.serialize(ss);
    assert(bytes > 0);

    invertible_structure_vec loaded;
    loaded.load(ss);

    assert(loaded.domain() == ms.domain());
    assert(loaded.runs() == ms.runs());
    assert(loaded.intervals() == ms.intervals());
    for (size_t i = 0; i < ms.intervals(); ++i) {
        assert(loaded.get_length(i) == ms.get_length(i));
        assert(loaded.get_pointer_fwd(i) == ms.get_pointer_fwd(i));
        assert(loaded.get_pointer_inv(i) == ms.get_pointer_inv(i));
        assert(loaded.get_fwd_interval(i) == ms.get_fwd_interval(i));
        assert(loaded.get_inv_interval(i) == ms.get_inv_interval(i));
    }

    assert_move_fwd_matches_perm(loaded, kRunnyPerm);
    assert_move_inv_matches_inverse(loaded, kRunnyPerm);
}

static void test_invertible_structure_widths() {
    const auto enc = make_invertible_encoding(kRunnyPerm);
    invertible_structure_vec ms(enc);

    const auto widths = ms.get_widths();
    const uchar w_pointer_fwd =
        widths[static_cast<size_t>(move_cols_traits<invertible_columns>::POINTER_FWD)];
    const uchar w_pointer_inv =
        widths[static_cast<size_t>(move_cols_traits<invertible_columns>::POINTER_INV)];
    const uchar w_fwd_interval =
        widths[static_cast<size_t>(move_cols_traits<invertible_columns>::FWD_INTERVAL)];
    const uchar w_inv_interval =
        widths[static_cast<size_t>(move_cols_traits<invertible_columns>::INV_INTERVAL)];

    assert(w_pointer_fwd >= bit_width(ms.intervals()));
    assert(w_pointer_inv >= bit_width(ms.intervals()));
    assert(w_fwd_interval == 1);
    assert(w_inv_interval == 1);
}

// Permutation-level tests for prev() / next() with invertible move columns.
using invertible_move_perm = permutation_impl<
    empty_data_columns, false, false, false,
    invertible_columns, move_structure, move_vector>;

using invertible_move_perm_absolute = permutation_impl<
    empty_data_columns, false, true, false,
    invertible_columns, move_structure, move_vector>;

static ulint compute_idx_relative_perm(const invertible_move_perm& p,
                                       invertible_move_perm::position pos) {
    ulint idx = 0;
    for (ulint i = 0; i < pos.interval; ++i) {
        idx += p.get_length(i);
    }
    return idx + pos.offset;
}

static invertible_move_perm::position make_pos_relative_perm(const invertible_move_perm& p, ulint idx) {
    ulint prefix = 0;
    for (ulint interval = 0; interval < p.intervals(); ++interval) {
        const ulint len = p.get_length(interval);
        if (idx < prefix + len) {
            return {interval, idx - prefix};
        }
        prefix += len;
    }
    assert(false && "index out of range");
    return {};
}

static invertible_move_perm_absolute::position make_pos_absolute_perm(
    const invertible_move_perm_absolute& p, ulint idx
) {
    invertible_move_perm_absolute::position pos{};
    pos.idx = idx;
    ulint prefix = 0;
    for (ulint interval = 0; interval < p.intervals(); ++interval) {
        const ulint len = p.get_length(interval);
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

static void test_invertible_permutation_next_prev() {
    const auto enc = make_invertible_encoding(kRunnyPerm);
    invertible_move_perm p(enc);

    assert(p.domain() == kRunnyPerm.size());
    assert(p.intervals() == enc.intervals());

    for (ulint idx = 0; idx < p.domain(); ++idx) {
        auto pos = make_pos_relative_perm(p, idx);
        auto next_pos = p.next(pos);
        assert(compute_idx_relative_perm(p, next_pos) == kRunnyPerm[idx]);

        auto back_pos = p.prev(next_pos);
        assert(compute_idx_relative_perm(p, back_pos) == idx);
    }
}

static void test_invertible_permutation_next_prev_multi_step() {
    const auto enc = make_invertible_encoding(kRunnyPerm);
    invertible_move_perm p(enc);

    for (ulint idx = 0; idx < p.domain(); ++idx) {
        auto pos = make_pos_relative_perm(p, idx);
        auto stepped = p.next(pos, 3);
        auto manual = pos;
        for (ulint s = 0; s < 3; ++s) {
            manual = p.next(manual);
        }
        assert(compute_idx_relative_perm(p, stepped) == compute_idx_relative_perm(p, manual));

        stepped = p.prev(stepped, 3);
        assert(compute_idx_relative_perm(p, stepped) == idx);
    }
}

static void test_invertible_permutation_absolute_exponential_prev() {
    const auto enc = make_invertible_encoding(kRunnyPerm);
    invertible_move_perm_absolute p(enc);

    for (ulint idx = 0; idx < p.domain(); ++idx) {
        auto pos = make_pos_absolute_perm(p, idx);

        auto next_lin = p.next_linear(pos);
        auto next_exp = p.next_exponential(pos);
        assert(next_lin.interval == next_exp.interval);
        assert(next_lin.offset == next_exp.offset);
        assert(next_lin.idx == next_exp.idx);
        assert(next_lin.idx == kRunnyPerm[idx]);

        auto next_default = p.next(pos);
        assert(next_default.idx == kRunnyPerm[idx]);

        auto prev_lin = p.prev_linear(next_lin);
        auto prev_exp = p.prev_exponential(next_lin);
        assert(prev_lin.interval == prev_exp.interval);
        assert(prev_lin.offset == prev_exp.offset);
        assert(prev_lin.idx == prev_exp.idx);
        assert(prev_lin.idx == idx);

        auto prev_default = p.prev(next_lin);
        assert(prev_default.idx == idx);
    }
}

static void test_invertible_permutation_up_down() {
    const auto enc = make_invertible_encoding(kSmallPerm);
    invertible_move_perm p(enc);

    auto pos = p.first();
    for (ulint i = 0; i < p.intervals(); ++i) {
        pos = p.down(pos);
    }
    auto first_pos = p.first();
    assert(pos.interval == first_pos.interval);
    assert(pos.offset == first_pos.offset);

    pos = p.first();
    for (ulint i = 0; i < p.intervals(); ++i) {
        pos = p.up(pos);
    }
    first_pos = p.first();
    assert(pos.interval == first_pos.interval);
    assert(pos.offset == first_pos.offset);
}

int main() {
    test_invertible_structure_build_from_encoding();
    test_invertible_structure_build_from_lengths_and_images();
    test_invertible_structure_move_fwd_matches_permutation();
    test_invertible_structure_move_inv_matches_inverse();
    test_invertible_structure_fwd_inv_roundtrip();
    test_invertible_structure_inv_fwd_roundtrip();
    test_invertible_structure_move_from_non_fwd_interval();
    test_invertible_structure_absolute_semantics();
    test_invertible_structure_absolute_move_exponential();
    test_invertible_structure_small_and_identity_permutations();
    test_invertible_structure_serialize_roundtrip();
    test_invertible_structure_widths();
    test_invertible_permutation_next_prev();
    test_invertible_permutation_next_prev_multi_step();
    test_invertible_permutation_absolute_exponential_prev();
    test_invertible_permutation_up_down();

    std::cout << "invertible_structure unit tests passed" << std::endl;
    return 0;
}
