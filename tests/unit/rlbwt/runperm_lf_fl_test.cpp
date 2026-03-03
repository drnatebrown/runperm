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

int main() {
    test_move_lf_wrapper_equivalence();
    test_runperm_lf_wrapper_equivalence();
    test_move_fl_wrapper_equivalence();
    test_runperm_fl_wrapper_equivalence();

    std::cout << "runperm_lf_fl unit tests passed" << std::endl;
    return 0;
}

