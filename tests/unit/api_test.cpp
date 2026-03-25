// Small API-level tests for public headers `move.hpp` and `rlbwt.hpp` and minor `runperm.hpp` tests.
// Goal: ensure the convenience aliases and basic constructors are usable.
// These are intentionally light on logic/algorithmic checks.

#include "orbit/move_structure.hpp"
#include "orbit/rlbwt.hpp"
#include "orbit/runperm.hpp"

#include <cassert>
#include <iostream>
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

    // Exercise the public aliases from `move.hpp`.
    move_structure_vec ms_vec(lengths, interval_perm, NO_SPLITTING);
    move_structure_vec_idx ms_vec_idx(lengths, interval_perm, NO_SPLITTING);

    // Just touch a couple of basic APIs.
    (void)ms_vec;
    (void)ms_vec_idx;
    assert(ms_vec.domain() == domain);
    assert(ms_vec_idx.runs() == lengths.size());
}

void test_rlbwt_api() {
    // Tiny artificial RLBWT: three runs whose lengths sum to 5.
    vector<uchar> heads = {static_cast<uchar>('A'),
                           static_cast<uchar>('C'),
                           static_cast<uchar>('G')};
    vector<ulint> run_lengths = {2, 1, 2};

    // Exercise MoveLF / MoveFL default aliases.
    move_lf<> move_lf(heads, run_lengths);
    move_fl<> move_fl(heads, run_lengths);

    assert(move_lf.domain() == 5);
    assert(move_fl.domain() == 5);

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

    runperm_lf<run_cols> rp_lf(heads, run_lengths, run_data);
    runperm_fl<run_cols> rp_fl(heads, run_lengths, run_data);

    // Touch a couple of simple methods to ensure the API is usable.
    using PosLF = typename runperm_lf<run_cols>::position;
    using PosFL = typename runperm_fl<run_cols>::position;

    PosLF pos_lf = rp_lf.first();
    PosFL pos_fl = rp_fl.first();

    (void)pos_lf;
    (void)pos_fl;
    (void)rp_lf.domain();
    (void)rp_fl.domain();

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
}

void test_runperm_header_is_available() {
    // Simple compile/link test for `runperm.hpp` top-level aliases.
    enum class test_cols {
        X,
        COUNT
    };

    using rp_separated = runperm_separated<test_cols>;
    using rp_integrated_absolute = runperm_integrated_absolute<test_cols>;

    const vector<ulint> lengths = {3};
    const vector<ulint> interval_perm = {0};
    const ulint domain = 3;

    using run_data = typename rp_separated::data_tuple;
    vector<run_data> data_sep(1);
    vector<run_data> data_int(1);

    rp_separated rp_sep(lengths, interval_perm, data_sep);
    rp_integrated_absolute rp_int(lengths, interval_perm, data_int);

    (void)rp_sep;
    (void)rp_int;
}

int main() {
    test_move_api();
    test_rlbwt_api();
    test_runperm_header_is_available();

    std::cout << "api tests passed" << std::endl;
    return 0;
}
