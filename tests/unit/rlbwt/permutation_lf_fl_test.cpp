// Tiny unit tests for LF/FL-based RLBWT components.
// Focus: wrapper methods (LF/FL) behave like next(), not full text reconstruction.

#include "orbit/rlbwt.hpp"

#include <cassert>
#include <iostream>
#include <vector>

using namespace orbit;
using namespace orbit::rlbwt;

using std::size_t;
using std::vector;

void test_move_lf_wrapper_equivalence() {
    // Same small RLBWT example as in the integration test, but we only
    // check a few LF steps.
    vector<uchar> bwt_heads =       {'T','C','G','A','T', 1 ,'A','T','A'};
    vector<ulint> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };

    move_lf<> move_lf(bwt_heads, bwt_run_lengths);

    using position = typename move_lf<>::position;
    position start = move_lf.first();

    // Single-step LF must match next().
    position lf1 = move_lf.LF(start);
    position next1 = move_lf.next(start);
    assert(lf1.interval == next1.interval);
    assert(lf1.offset == next1.offset);

    // Multi-step LF(pos, k) must match k repeated LF calls.
    position lf3 = move_lf.LF(start, 3);
    position iter = start;
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

    runperm_lf<RunCols> rp_lf(bwt_heads, bwt_run_lengths, run_data);

    using position = typename runperm_lf<RunCols>::position;
    position start = rp_lf.first();

    position lf1 = rp_lf.LF(start);
    position next1 = rp_lf.next(start);
    assert(lf1.interval == next1.interval);
    assert(lf1.offset == next1.offset);

    position lf2 = rp_lf.LF(start, 2);
    position iter = start;
    iter = rp_lf.LF(iter);
    iter = rp_lf.LF(iter);
    assert(lf2.interval == iter.interval);
    assert(lf2.offset == iter.offset);
}

void test_move_fl_wrapper_equivalence() {
    vector<uchar> bwt_heads =       {'T','C','G','A','T', 1 ,'A','T','A'};
    vector<ulint> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };

    move_fl<> move_fl(bwt_heads, bwt_run_lengths);

    using position = typename move_fl<>::position;
    position start = move_fl.first();

    position fl1 = move_fl.FL(start);
    position next1 = move_fl.next(start);
    assert(fl1.interval == next1.interval);
    assert(fl1.offset == next1.offset);

    position fl4 = move_fl.FL(start, 4);
    position iter = start;
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

    runperm_fl<RunCols> rp_fl(bwt_heads, bwt_run_lengths, run_data);

    using position = typename runperm_fl<RunCols>::position;
    position start = rp_fl.first();

    position fl1 = rp_fl.FL(start);
    position next1 = rp_fl.next(start);
    assert(fl1.interval == next1.interval);
    assert(fl1.offset == next1.offset);

    position fl2 = rp_fl.FL(start, 2);
    position iter = start;
    iter = rp_fl.FL(iter);
    iter = rp_fl.FL(iter);
    assert(fl2.interval == iter.interval);
    assert(fl2.offset == iter.offset);
}

void test_runpermlf_construct_from_precomputed_permutation_no_splitting() {
    vector<uchar> bwt_heads =       {'T','C','G','A','T', 1 ,'A','T','A'};
    vector<ulint> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };

    auto perm = rlbwt_interval_encoding<>::lf_interval_encoding(bwt_heads, bwt_run_lengths, NO_SPLITTING);

    enum class RunCols {
        V,
        COUNT
    };
    using RP = runperm_lf<RunCols>;
    using RunData = typename RP::data_tuple;

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

    auto perm = rlbwt_interval_encoding<>::lf_interval_encoding(bwt_heads, bwt_run_lengths, DEFAULT_SPLITTING);
    assert(perm.intervals() >= perm.runs());

    enum class RunCols {
        V,
        COUNT
    };
    using RP = runperm_lf<RunCols>;
    using RunData = typename RP::data_tuple;

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

void test_bwt_to_rlbwt_basic() {
    {
        // Empty input -> empty output.
        vector<uchar> bwt;
        auto [heads, lengths] = bwt_to_rlbwt(bwt);
        assert(heads.empty());
        assert(lengths.empty());
    }
    {
        // Single run.
        vector<uchar> bwt = {'A', 'A', 'A'};
        auto [heads, lengths] = bwt_to_rlbwt(bwt);
        assert(heads.size() == 1);
        assert(lengths.size() == 1);
        assert(heads[0] == static_cast<uchar>('A'));
        assert(lengths[0] == 3);
    }
    {
        // Multiple runs.
        vector<uchar> bwt = {'A', 'A', 'C', 'C', 'C', 'G', 'G'};
        auto [heads, lengths] = bwt_to_rlbwt(bwt);

        assert(heads.size() == 3);
        assert(lengths.size() == 3);

        assert(heads[0] == static_cast<uchar>('A'));
        assert(lengths[0] == 2);

        assert(heads[1] == static_cast<uchar>('C'));
        assert(lengths[1] == 3);

        assert(heads[2] == static_cast<uchar>('G'));
        assert(lengths[2] == 2);

        // Sum of run lengths equals original length.
        ulint total = 0;
        for (ulint v : lengths) {
            total += v;
        }
        assert(total == bwt.size());
    }
}

int main() {
    test_move_lf_wrapper_equivalence();
    test_runperm_lf_wrapper_equivalence();
    test_move_fl_wrapper_equivalence();
    test_runperm_fl_wrapper_equivalence();
    test_runpermlf_construct_from_precomputed_permutation_no_splitting();
    test_runpermlf_construct_from_precomputed_permutation_with_splitting();
    test_bwt_to_rlbwt_basic();

    std::cout << "runperm_lf_fl unit tests passed" << std::endl;
    return 0;
}
