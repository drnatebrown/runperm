// Tiny unit tests for LF/FL-based RLBWT components.
// Focus: wrapper methods (LF/FL) behave like next(), not full text reconstruction.

#include "rlbwt.hpp"

#include <cassert>
#include <iostream>
#include <vector>

using std::size_t;
using std::vector;

void test_move_lf_wrapper_equivalence() {
    // Same small RLBWT example as in the integration test, but we only
    // check a few LF steps.
    vector<uchar> bwt_heads =       {'T','C','G','A','T', 1 ,'A','T','A'};
    vector<ulint> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };

    MoveLF<> move_lf(bwt_heads, bwt_run_lengths);

    using Position = typename MoveLF<>::Position;
    Position start = move_lf.first();

    // Single-step LF must match next().
    Position lf1 = move_lf.LF(start);
    Position next1 = move_lf.next(start);
    assert(lf1.interval == next1.interval);
    assert(lf1.offset == next1.offset);

    // Multi-step LF(pos, k) must match k repeated LF calls.
    Position lf3 = move_lf.LF(start, 3);
    Position iter = start;
    iter = move_lf.LF(iter);
    iter = move_lf.LF(iter);
    iter = move_lf.LF(iter);
    assert(lf3.interval == iter.interval);
    assert(lf3.offset == iter.offset);
}

void test_runperm_lf_wrapper_equivalence() {
    vector<uchar> bwt_heads =       {'T','C','G','A','T', 1 ,'A','T','A'};
    vector<ulint> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };

    enum class RunCols {
        V,
        COUNT
    };
    static constexpr size_t NUM_FIELDS = static_cast<size_t>(RunCols::COUNT);
    vector<std::array<ulint, NUM_FIELDS>> run_data(bwt_heads.size());
    for (size_t i = 0; i < run_data.size(); ++i) {
        run_data[i][0] = static_cast<ulint>(i);
    }

    RunPermLF<RunCols> rp_lf(bwt_heads, bwt_run_lengths, run_data);

    using Position = typename RunPermLF<RunCols>::Position;
    Position start = rp_lf.first();

    Position lf1 = rp_lf.LF(start);
    Position next1 = rp_lf.next(start);
    assert(lf1.interval == next1.interval);
    assert(lf1.offset == next1.offset);

    Position lf2 = rp_lf.LF(start, 2);
    Position iter = start;
    iter = rp_lf.LF(iter);
    iter = rp_lf.LF(iter);
    assert(lf2.interval == iter.interval);
    assert(lf2.offset == iter.offset);
}

void test_move_fl_wrapper_equivalence() {
    vector<uchar> bwt_heads =       {'T','C','G','A','T', 1 ,'A','T','A'};
    vector<ulint> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };

    MoveFL<> move_fl(bwt_heads, bwt_run_lengths);

    using Position = typename MoveFL<>::Position;
    Position start = move_fl.first();

    Position fl1 = move_fl.FL(start);
    Position next1 = move_fl.next(start);
    assert(fl1.interval == next1.interval);
    assert(fl1.offset == next1.offset);

    Position fl4 = move_fl.FL(start, 4);
    Position iter = start;
    for (int i = 0; i < 4; ++i) {
        iter = move_fl.FL(iter);
    }
    assert(fl4.interval == iter.interval);
    assert(fl4.offset == iter.offset);
}

void test_runperm_fl_wrapper_equivalence() {
    vector<uchar> bwt_heads =       {'T','C','G','A','T', 1 ,'A','T','A'};
    vector<ulint> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };

    enum class RunCols {
        V,
        COUNT
    };
    static constexpr size_t NUM_FIELDS = static_cast<size_t>(RunCols::COUNT);
    vector<std::array<ulint, NUM_FIELDS>> run_data(bwt_heads.size());
    for (size_t i = 0; i < run_data.size(); ++i) {
        run_data[i][0] = static_cast<ulint>(i * 2);
    }

    RunPermFL<RunCols> rp_fl(bwt_heads, bwt_run_lengths, run_data);

    using Position = typename RunPermFL<RunCols>::Position;
    Position start = rp_fl.first();

    Position fl1 = rp_fl.FL(start);
    Position next1 = rp_fl.next(start);
    assert(fl1.interval == next1.interval);
    assert(fl1.offset == next1.offset);

    Position fl2 = rp_fl.FL(start, 2);
    Position iter = start;
    iter = rp_fl.FL(iter);
    iter = rp_fl.FL(iter);
    assert(fl2.interval == iter.interval);
    assert(fl2.offset == iter.offset);
}

void test_runpermlf_construct_from_precomputed_permutation_no_splitting() {
    vector<uchar> bwt_heads =       {'T','C','G','A','T', 1 ,'A','T','A'};
    vector<ulint> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };

    auto perm = RLBWTPermutation::lf_permutation(bwt_heads, bwt_run_lengths, NO_SPLITTING);

    enum class RunCols {
        V,
        COUNT
    };
    using RP = RunPermLF<RunCols>;
    using RunData = typename RP::RunDataTuple;

    vector<RunData> run_data(perm.intervals());
    for (size_t i = 0; i < run_data.size(); ++i) {
        run_data[i][0] = static_cast<ulint>(1000 + i);
    }

    RP rp_from_perm(perm, run_data);
    RP rp_direct(bwt_heads, bwt_run_lengths, NO_SPLITTING, run_data);

    assert(rp_from_perm.domain() == rp_direct.domain());
    assert(rp_from_perm.runs() == rp_direct.runs());
    assert(rp_from_perm.intervals() == rp_direct.intervals());
    assert(rp_from_perm.get_split_params() == NO_SPLITTING);

    auto p1 = rp_from_perm.first();
    auto p2 = rp_direct.first();
    for (size_t i = 0; i < rp_direct.domain(); ++i) {
        assert(rp_from_perm.get_character(p1) == rp_direct.get_character(p2));
        p1 = rp_from_perm.LF(p1);
        p2 = rp_direct.LF(p2);
    }
    assert(p1.interval == 0 && p1.offset == 0);
    assert(p2.interval == 0 && p2.offset == 0);

    for (size_t i = 0; i < rp_from_perm.intervals(); ++i) {
        assert(rp_from_perm.template get<RunCols::V>(static_cast<ulint>(i)) == static_cast<ulint>(1000 + i));
        assert(rp_direct .template get<RunCols::V>(static_cast<ulint>(i)) == static_cast<ulint>(1000 + i));
    }
}

void test_runpermlf_construct_from_precomputed_permutation_with_splitting() {
    vector<uchar> bwt_heads =       {'T','C','G','A','T', 1 ,'A','T','A'};
    vector<ulint> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };

    auto perm = RLBWTPermutation::lf_permutation(bwt_heads, bwt_run_lengths, DEFAULT_SPLITTING);
    assert(perm.intervals() >= perm.runs());

    enum class RunCols {
        V,
        COUNT
    };
    using RP = RunPermLF<RunCols>;
    using RunData = typename RP::RunDataTuple;

    vector<RunData> run_data_per_run(perm.runs());
    for (size_t i = 0; i < run_data_per_run.size(); ++i) {
        run_data_per_run[i][0] = static_cast<ulint>(2000 + i);
    }
    auto run_data_split = perm.split_run_data_with_copy(bwt_run_lengths, run_data_per_run);
    assert(run_data_split.size() == perm.intervals());

    RP rp_from_perm(perm, run_data_split);
    assert(rp_from_perm.get_split_params() == DEFAULT_SPLITTING);
    assert(rp_from_perm.intervals() == perm.intervals());

    for (size_t i = 0; i < rp_from_perm.intervals(); ++i) {
        assert(rp_from_perm.template get<RunCols::V>(static_cast<ulint>(i)) == run_data_split[i][0]);
    }
}

int main() {
    test_move_lf_wrapper_equivalence();
    test_runperm_lf_wrapper_equivalence();
    test_move_fl_wrapper_equivalence();
    test_runperm_fl_wrapper_equivalence();
    test_runpermlf_construct_from_precomputed_permutation_no_splitting();
    test_runpermlf_construct_from_precomputed_permutation_with_splitting();

    std::cout << "runperm_lf_fl unit tests passed" << std::endl;
    return 0;
}
