// Small API-level tests for public headers `move_structure.hpp`, `interval_encoding.hpp`,
// `rlbwt.hpp`, and `permutation.hpp`.
// Goal: ensure convenience aliases and basic constructors compile and link.
// These are intentionally light on logic/algorithmic checks.

#include "orbit/internal/rlbwt/specializations/rlbwt_structure.hpp"
#include "orbit/interval_encoding.hpp"
#include "orbit/move_structure.hpp"
#include "orbit/rlbwt.hpp"
#include "orbit/permutation.hpp"

#include <cassert>
#include <iostream>
#include <type_traits>
#include <vector>

using std::size_t;
using std::vector;

using namespace orbit;
using namespace orbit::rlbwt;

void test_move_api() {
    // Simple permutation of size 4 split into two intervals.
    const vector<ulint> lengths = {2, 2};
    const vector<ulint> interval_perm = {2, 0};
    const ulint domain = 4;

    // Exercise all public move_structure aliases (4 column types x 2 table backends).
    move_structure_vec ms_vec(lengths, interval_perm, NO_SPLITTING);
    move_structure_vec_idx ms_vec_idx(lengths, interval_perm, NO_SPLITTING);
    move_structure_tbl ms_tbl(lengths, interval_perm, NO_SPLITTING);
    move_structure_tbl_idx ms_tbl_idx(lengths, interval_perm, NO_SPLITTING);

    static_assert(std::is_same_v<move_structure_vec,
        move_structure<move_columns, move_vector>>);
    static_assert(std::is_same_v<move_structure_tbl,
        move_structure<move_columns, move_table>>);
    static_assert(std::is_same_v<move_structure_vec_idx,
        move_structure<move_columns_idx, move_vector>>);
    static_assert(std::is_same_v<move_structure_tbl_idx,
        move_structure<move_columns_idx, move_table>>);

    assert(ms_vec.domain() == domain);
    assert(ms_vec_idx.runs() == lengths.size());
    assert(ms_tbl.domain() == domain);
    assert(ms_tbl_idx.domain() == domain);
    (void)ms_tbl;
    (void)ms_tbl_idx;
}

void test_move_invertible_api() {
    const vector<ulint> permutation = {1, 0, 3, 2};
    const auto enc = invertible_interval_encoding::from_permutation(permutation, NO_SPLITTING);

    invertible_structure_vec ms_vec(enc);
    invertible_structure_vec_idx ms_vec_idx(enc);
    invertible_structure_tbl ms_tbl(enc);
    invertible_structure_tbl_idx ms_tbl_idx(enc);

    static_assert(std::is_same_v<
        move_structure<invertible_columns, move_vector>,
        invertible_structure_impl<invertible_columns, move_vector>>);
    static_assert(std::is_same_v<invertible_structure_vec, move_structure<invertible_columns, move_vector>>);
    static_assert(std::is_same_v<invertible_structure_vec_idx, move_structure<invertible_columns_idx, move_vector>>);
    static_assert(std::is_same_v<invertible_structure_tbl, move_structure<invertible_columns, move_table>>);
    static_assert(std::is_same_v<invertible_structure_tbl_idx, move_structure<invertible_columns_idx, move_table>>);

    assert(ms_vec.domain() == permutation.size());
    assert(ms_vec.get_file_extension() == ".imove");
    assert(ms_vec.get_fwd_interval(0));

    auto pos = ms_vec.first();
    pos = ms_vec.move_fwd(pos);
    (void)ms_vec.move_inv(pos);
    assert(ms_vec_idx.domain() == permutation.size());
    assert(ms_tbl.domain() == permutation.size());
    assert(ms_tbl_idx.domain() == permutation.size());
    (void)ms_vec_idx;
    (void)ms_tbl;
    (void)ms_tbl_idx;
}

void test_interval_encoding_api() {
    const vector<ulint> permutation = {1, 0, 3, 2};
    interval_encoding enc = interval_encoding::from_permutation(permutation, NO_SPLITTING);
    invertible_interval_encoding inv_enc =
        invertible_interval_encoding::from_permutation(permutation, NO_SPLITTING);

    static_assert(!interval_encoding::invertible_tag);
    static_assert(invertible_interval_encoding::invertible_tag);

    assert(enc.domain() == permutation.size());
    assert(inv_enc.domain() == permutation.size());
    (void)enc.get_split_params();
    (void)inv_enc.get_split_params();
}

void test_rlbwt_api() {
    // Tiny artificial RLBWT: three runs whose lengths sum to 5.
    vector<uchar> heads = {static_cast<uchar>('A'),
                           static_cast<uchar>('C'),
                           static_cast<uchar>('G')};
    vector<ulint> run_lengths = {2, 1, 2};

    // LF/FL move aliases (relative vs absolute).
    lf_move_relative move_lf_rel(heads, run_lengths);
    lf_move_absolute move_lf_abs(heads, run_lengths);
    fl_move_relative move_fl_rel(heads, run_lengths);
    fl_move_absolute move_fl_abs(heads, run_lengths);

    static_assert(std::is_same_v<lf_move<>, lf_move_relative>);
    static_assert(std::is_same_v<fl_move<>, fl_move_relative>);

    assert(move_lf_rel.domain() == 5);
    assert(move_lf_abs.domain() == 5);
    assert(move_fl_rel.domain() == 5);
    assert(move_fl_abs.domain() == 5);

    const auto lf_enc = rlbwt_interval_encoding<>::lf_interval_encoding(heads, run_lengths, NO_SPLITTING);

    // All rlbwt_move_structure column x table combinations (non-invertible).
    rlbwt_move_structure<rlbwt_columns> ms_lf_vec(lf_enc);
    rlbwt_move_structure<rlbwt_columns, move_table> ms_lf_tbl(lf_enc);
    rlbwt_move_structure<rlbwt_columns_idx> ms_lf_vec_idx(lf_enc);
    rlbwt_move_structure<rlbwt_columns_idx, move_table> ms_lf_tbl_idx(lf_enc);

    assert(ms_lf_vec.domain() == 5);
    assert(ms_lf_tbl.domain() == 5);
    assert(ms_lf_vec_idx.domain() == 5);
    assert(ms_lf_tbl_idx.domain() == 5);
    (void)ms_lf_vec_idx;
    (void)ms_lf_tbl_idx;

    // Exercise one RunPermLF / RunPermFL instantiation with trivial run data.
    enum class run_cols {
        VAL,
        COUNT
    };
    static constexpr size_t NUM_FIELDS = static_cast<size_t>(run_cols::COUNT);
    vector<std::array<ulint, NUM_FIELDS>> run_data(heads.size());
    for (size_t i = 0; i < run_data.size(); ++i) {
        run_data[i][0] = static_cast<ulint>(i);
    }

    // All public LF/FL permutation layout aliases (integrated x absolute).
    lf_permutation_separated<run_cols> rp_lf_sep(heads, run_lengths, run_data);
    lf_permutation_integrated<run_cols> rp_lf_int(heads, run_lengths, run_data);
    lf_permutation_absolute<run_cols> rp_lf_abs(heads, run_lengths, run_data);
    lf_permutation_integrated_absolute<run_cols> rp_lf_int_abs(heads, run_lengths, run_data);

    fl_permutation_separated<run_cols> rp_fl_sep(heads, run_lengths, run_data);
    fl_permutation_integrated<run_cols> rp_fl_int(heads, run_lengths, run_data);
    fl_permutation_absolute<run_cols> rp_fl_abs(heads, run_lengths, run_data);
    fl_permutation_integrated_absolute<run_cols> rp_fl_int_abs(heads, run_lengths, run_data);

    static_assert(std::is_same_v<lf_permutation<run_cols>, lf_permutation_separated<run_cols>>);
    static_assert(std::is_same_v<fl_permutation<run_cols>, fl_permutation_separated<run_cols>>);

    using PosLF = typename lf_permutation<run_cols>::position;
    using PosFL = typename fl_permutation<run_cols>::position;

    PosLF pos_lf = rp_lf_sep.first();
    PosFL pos_fl = rp_fl_sep.first();
    (void)rp_lf_sep.LF(pos_lf);
    (void)rp_fl_sep.FL(pos_fl);

    assert(rp_lf_sep.domain() == 5);
    assert(rp_lf_int.domain() == 5);
    assert(rp_lf_abs.domain() == 5);
    assert(rp_lf_int_abs.domain() == 5);
    assert(rp_fl_sep.domain() == 5);
    assert(rp_fl_int.domain() == 5);
    assert(rp_fl_abs.domain() == 5);
    assert(rp_fl_int_abs.domain() == 5);
    (void)rp_lf_int;
    (void)rp_lf_abs;
    (void)rp_lf_int_abs;
    (void)rp_fl_int;
    (void)rp_fl_abs;
    (void)rp_fl_int_abs;

    // Exercise the public Phi / phi_inv convenience wrappers on a known-valid example.
    // TEXT: GATTACATGATTACATAGATTACATT$
    // BWT:  TTTTTCCCGGGAAAT$ATTTTAAAAAA
    // RLBWT: TCGAT$ATA
    vector<uchar> bwt_heads =       {'T','C','G','A','T', 1 ,'A','T','A'};
    vector<ulint> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };

    auto phi_perm = rlbwt_to_phi(bwt_heads, bwt_run_lengths, NO_SPLITTING);
    auto phi_inv_perm = rlbwt_to_phi_inv(bwt_heads, bwt_run_lengths, NO_SPLITTING);

    assert(phi_perm.domain() == phi_inv_perm.domain());
    assert(phi_perm.domain() == 27);
    assert(phi_perm.get_split_params() == NO_SPLITTING);
    assert(phi_inv_perm.get_split_params() == NO_SPLITTING);

    // Phi / phi_inv wrappers need a valid BWT (use the GATTACA example above).
    phi_move move_phi(bwt_heads, bwt_run_lengths);
    phi_inv_move move_phi_inv(bwt_heads, bwt_run_lengths);
    phi_permutation<> rp_phi(bwt_heads, bwt_run_lengths);
    phi_permutation_integrated<> rp_phi_int(bwt_heads, bwt_run_lengths);
    phi_inv_permutation<> rp_phi_inv(bwt_heads, bwt_run_lengths);
    phi_inv_permutation_integrated<> rp_phi_inv_int(bwt_heads, bwt_run_lengths);

    assert(move_phi.domain() == 27);
    assert(move_phi_inv.domain() == 27);
    assert(rp_phi.domain() == 27);
    assert(rp_phi_int.domain() == 27);
    assert(rp_phi_inv.domain() == 27);
    assert(rp_phi_inv_int.domain() == 27);

    auto pos_phi = rp_phi.first();
    (void)rp_phi.phi(pos_phi);
    (void)move_phi.phi(move_phi.first());
    (void)move_phi_inv.phi_inv(move_phi_inv.first());

    // Phi permutations with run data (from derived interval images).
    size_t phi_domain = 0;
    ulint max_length = 0;
    auto [phi_lengths, phi_images] =
        rlbwt_to_phi_images(bwt_heads, bwt_run_lengths, &phi_domain, &max_length);
    size_t inv_domain = 0;
    ulint max_length_inv = 0;
    auto [inv_lengths, inv_images] =
        rlbwt_to_phi_inv_images(bwt_heads, bwt_run_lengths, &inv_domain, &max_length_inv);

    vector<std::array<ulint, NUM_FIELDS>> run_data_phi(phi_lengths.size());
    vector<std::array<ulint, NUM_FIELDS>> run_data_inv(inv_lengths.size());
    phi_permutation<run_cols> rp_phi_data(phi_lengths, phi_images, run_data_phi);
    phi_inv_permutation<run_cols> rp_phi_inv_data(inv_lengths, inv_images, run_data_inv);

    assert(rp_phi_data.domain() == phi_domain);
    assert(rp_phi_inv_data.domain() == inv_domain);
    (void)rp_phi_data.first();
    (void)rp_phi_inv_data.first();
    (void)move_phi;
    (void)move_phi_inv;
    (void)rp_phi;
    (void)rp_phi_int;
    (void)rp_phi_inv;
    (void)rp_phi_inv_int;
}

void test_rlbwt_invertible_api() {
    vector<uchar> heads = {static_cast<uchar>('A'),
                           static_cast<uchar>('C'),
                           static_cast<uchar>('G')};
    vector<ulint> run_lengths = {2, 1, 2};

    const auto enc = invertible_rlbwt_interval_encoding<>::lf_interval_encoding(
        heads, run_lengths, NO_SPLITTING);
    // All rlbwt_move_structure column x table combinations (invertible).
    rlbwt_move_structure<rlbwt_invertible_columns> ms(enc);
    rlbwt_move_structure<rlbwt_invertible_columns_idx> ms_idx(enc);
    rlbwt_move_structure<rlbwt_invertible_columns, move_table> ms_tbl(enc);
    rlbwt_move_structure<rlbwt_invertible_columns_idx, move_table> ms_tbl_idx(enc);

    static_assert(invertible_rlbwt_interval_encoding<>::invertible_tag == true);
    static_assert(std::is_same_v<
        typename table_row_for<rlbwt_invertible_columns>::type,
        invertible_rlbwt_row<rlbwt_invertible_columns>>);
    static_assert(std::is_same_v<
        typename table_row_for<rlbwt_invertible_columns_idx>::type,
        invertible_rlbwt_row<rlbwt_invertible_columns_idx>>);

    assert(ms.domain() == 5);
    assert(ms.get_file_extension() == ".imove");
    assert(ms.get_fwd_interval(0));

    auto pos = ms.first();
    pos = ms.move_fwd(pos);
    (void)ms.move_inv(pos);
    assert(ms_idx.domain() == 5);
    assert(ms_tbl.domain() == 5);
    assert(ms_tbl_idx.domain() == 5);
    (void)ms_idx;
    (void)ms_tbl;
    (void)ms_tbl_idx;
}

void test_permutation_api() {
    // All public permutation layout aliases (integrated x absolute).
    enum class test_cols {
        X,
        COUNT
    };

    const vector<ulint> lengths = {2, 2};
    const vector<ulint> interval_perm = {2, 0};
    vector<ulint> perm = {1, 0, 3, 2};

    using run_data = typename permutation_separated<test_cols>::data_tuple;
    vector<run_data> data(lengths.size());

    permutation_separated<test_cols> rp_sep(lengths, interval_perm, data);
    permutation_integrated<test_cols> rp_int(lengths, interval_perm, data);
    permutation_absolute<test_cols> rp_abs(lengths, interval_perm, data);
    permutation_integrated_absolute<test_cols> rp_int_abs(lengths, interval_perm, data);

    static_assert(std::is_same_v<permutation<test_cols>, permutation_separated<test_cols>>);
    static_assert(std::is_same_v<permutation_absolute<test_cols>, permutation_separated_absolute<test_cols>>);

    assert(rp_sep.domain() == perm.size());
    assert(rp_int.domain() == perm.size());
    assert(rp_abs.domain() == perm.size());
    assert(rp_int_abs.domain() == perm.size());
    (void)rp_sep.first();
    (void)rp_int.first();
    (void)rp_abs.first();
    (void)rp_int_abs.first();

    // move_permutation aliases (no run data).
    move_permutation_relative mp_rel(perm);
    move_permutation_absolute mp_abs(lengths, interval_perm);

    static_assert(std::is_same_v<move_permutation<>, move_permutation_relative>);

    assert(mp_rel.domain() == perm.size());
    assert(mp_abs.domain() == perm.size());
    (void)mp_rel.first();
    (void)mp_abs.first();
}

int main() {
    test_move_api();
    test_move_invertible_api();
    test_interval_encoding_api();
    test_permutation_api();
    test_rlbwt_api();
    test_rlbwt_invertible_api();

    std::cout << "api tests passed" << std::endl;
    return 0;
}
