// Small API-level tests for public headers `move.hpp` and `rlbwt.hpp` and minor `runperm.hpp` tests.
// Goal: ensure the convenience aliases and basic constructors are usable.
// These are intentionally light on logic/algorithmic checks.

#include "move.hpp"
#include "rlbwt.hpp"
#include "runperm.hpp"

#include <cassert>
#include <iostream>
#include <vector>

using std::size_t;
using std::vector;

void test_move_api() {
    // Simple permutation of size 4 split into two intervals.
    const vector<ulint> lengths = {2, 2};
    const vector<ulint> interval_perm = {2, 0};
    const ulint domain = 4;

    // Exercise the public aliases from `move.hpp`.
    MoveStructureVec ms_vec(lengths, interval_perm, NO_SPLITTING);
    MoveStructureVecIdx ms_vec_idx(lengths, interval_perm, NO_SPLITTING);

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
    MoveLF<> move_lf(heads, run_lengths);
    MoveFL<> move_fl(heads, run_lengths);

    assert(move_lf.domain() == 5);
    assert(move_fl.domain() == 5);

    // Exercise one RunPermLF / RunPermFL instantiation with trivial run data.
    enum class RunCols {
        VAL,
        COUNT
    };
    static constexpr size_t NUM_FIELDS = static_cast<size_t>(RunCols::COUNT);
    vector<std::array<ulint, NUM_FIELDS>> run_data(heads.size());
    for (size_t i = 0; i < run_data.size(); ++i) {
        run_data[i][0] = static_cast<ulint>(i);
    }

    RunPermLF<RunCols> rp_lf(heads, run_lengths, run_data);
    RunPermFL<RunCols> rp_fl(heads, run_lengths, run_data);

    // Touch a couple of simple methods to ensure the API is usable.
    using PosLF = typename RunPermLF<RunCols>::Position;
    using PosFL = typename RunPermFL<RunCols>::Position;

    PosLF pos_lf = rp_lf.first();
    PosFL pos_fl = rp_fl.first();

    (void)pos_lf;
    (void)pos_fl;
    (void)rp_lf.domain();
    (void)rp_fl.domain();

    // Exercise the public Phi / InvPhi convenience wrappers on a known-valid example.
    // TEXT: GATTACATGATTACATAGATTACATT$
    // BWT:  TTTTTCCCGGGAAAT$ATTTTAAAAAA
    // RLBWT: TCGAT$ATA
    vector<uchar> bwt_heads =       {'T','C','G','A','T', 1 ,'A','T','A'};
    vector<ulint> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };

    auto phi_perm = phi::phi(bwt_heads, bwt_run_lengths, NO_SPLITTING);
    auto invphi_perm = invphi::invphi(bwt_heads, bwt_run_lengths, NO_SPLITTING);

    assert(phi_perm.domain() == invphi_perm.domain());
    assert(phi_perm.domain() == 27);
    assert(phi_perm.get_split_params() == NO_SPLITTING);
    assert(invphi_perm.get_split_params() == NO_SPLITTING);
}

void test_runperm_header_is_available() {
    // Simple compile/link test for `runperm.hpp` top-level aliases.
    enum class TestCols {
        X,
        COUNT
    };

    using RPSeparated = RunPermSeparated<TestCols>;
    using RPIntegratedAbs = RunPermIntegratedAbsolute<TestCols>;

    const vector<ulint> lengths = {3};
    const vector<ulint> interval_perm = {0};
    const ulint domain = 3;

    using RunData = typename RPSeparated::RunData;
    vector<RunData> data_sep(1);
    vector<RunData> data_int(1);

    RPSeparated rp_sep(lengths, interval_perm, data_sep);
    RPIntegratedAbs rp_int(lengths, interval_perm, data_int);

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
