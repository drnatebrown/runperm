// Unit tests for RunPerm (RunPermImpl via public alias).
// Simple assert-based tests, no external framework.

#include "runperm.hpp"

#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>

using std::size_t;
using std::vector;

// Simple run-data schema with two fields.
enum class TestRunCols {
    VAL1,
    VAL2,
    COUNT
};

using TestRunData = ColumnsTuple<TestRunCols>;

// Helper to build an absolute RunPerm position from a global index.
template <typename RP>
static typename RP::Position make_pos_absolute(const RP &rp, ulint idx) {
    using Position = typename RP::Position;
    Position pos{};
    pos.idx = idx;
    ulint prefix = 0;
    for (ulint interval = 0; interval < rp.intervals(); ++interval) {
        ulint len = rp.get_length(interval);
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

static void test_runperm_separated_absolute_basic_mapping_and_run_data() {
    const vector<ulint> perm = {1, 2, 9, 10, 11, 3, 12, 13, 4, 5, 14, 0, 15, 6, 7, 8};
    const ulint domain = static_cast<ulint>(perm.size());
    auto [lengths, interval_perm] = get_permutation_intervals(perm);

    // Set up simple run data: VAL1 = interval index, VAL2 = interval index + 100.
    vector<TestRunData> run_data(lengths.size());
    for (size_t i = 0; i < lengths.size(); ++i) {
        run_data[i] = {static_cast<ulint>(i), static_cast<ulint>(i + 100)};
    }

    using RP = RunPermSeparatedAbsolute<TestRunCols>;
    RP rp(lengths, interval_perm, run_data);

    assert(rp.domain() == domain);
    assert(rp.runs() == lengths.size());
    assert(rp.intervals() == lengths.size());

    // Interval-level checks.
    for (ulint i = 0; i < rp.intervals(); ++i) {
        assert(rp.get_length(i) == lengths[i]);
        assert(rp.get<TestRunCols::VAL1>(i) == static_cast<ulint>(i));
        assert(rp.get<TestRunCols::VAL2>(i) == static_cast<ulint>(i + 100));
        // get_row must match per-column get
        auto row = rp.get_row(i);
        assert(row[0] == rp.get<TestRunCols::VAL1>(i));
        assert(row[1] == rp.get<TestRunCols::VAL2>(i));
        assert(row == run_data[i]);
    }

    // Position-level mapping: next() must follow perm, and run data must agree.
    for (ulint idx = 0; idx < domain; ++idx) {
        auto pos = make_pos_absolute<RP>(rp, idx);
        auto next_pos = rp.next(pos);
        assert(next_pos.idx == perm[idx]);

        // Run-data for a position should equal the data for its interval.
        ulint interval = next_pos.interval;
        assert(rp.get<TestRunCols::VAL1>(next_pos) == rp.get<TestRunCols::VAL1>(interval));
        assert(rp.get<TestRunCols::VAL2>(next_pos) == rp.get<TestRunCols::VAL2>(interval));
        // get_row(Position) must equal get_row(interval)
        assert(rp.get_row(next_pos) == rp.get_row(interval));
    }
}

static void test_runperm_up_down_navigation() {
    const vector<ulint> perm = {1, 0, 3, 2};
    const ulint domain = static_cast<ulint>(perm.size());
    auto [lengths, interval_perm] = get_permutation_intervals(perm);

    vector<TestRunData> run_data(lengths.size());
    for (size_t i = 0; i < lengths.size(); ++i) {
        run_data[i] = {static_cast<ulint>(i), 0};
    }

    using RP = RunPermSeparatedAbsolute<TestRunCols>;
    RP rp(lengths, interval_perm, run_data);

    auto pos = rp.first();
    // Going down move_runs() times should wrap back to first().
    for (ulint i = 0; i < rp.intervals(); ++i) {
        pos = rp.down(pos);
    }
    auto first_pos = rp.first();
    assert(pos.interval == first_pos.interval);
    assert(pos.offset == first_pos.offset);

    // Going up move_runs() times from first() should also wrap back.
    pos = rp.first();
    for (ulint i = 0; i < rp.intervals(); ++i) {
        pos = rp.up(pos);
    }
    first_pos = rp.first();
    assert(pos.interval == first_pos.interval);
    assert(pos.offset == first_pos.offset);
}

static void test_runperm_serialize_roundtrip_separated_absolute() {
    const vector<ulint> perm = {2, 0, 3, 1};
    const ulint domain = static_cast<ulint>(perm.size());
    auto [lengths, interval_perm] = get_permutation_intervals(perm);

    vector<TestRunData> run_data(lengths.size());
    for (size_t i = 0; i < lengths.size(); ++i) {
        run_data[i] = {static_cast<ulint>(i * 10), static_cast<ulint>(i * 10 + 1)};
    }

    using RP = RunPermSeparatedAbsolute<TestRunCols>;
    RP rp(lengths, interval_perm, run_data);

    std::stringstream ss;
    size_t bytes = rp.serialize(ss);
    assert(bytes > 0);

    RP loaded;
    loaded.load(ss);

    assert(loaded.domain() == rp.domain());
    assert(loaded.intervals() == rp.intervals());
    assert(loaded.runs() == rp.runs());
    assert(loaded.get_split_params() == NO_SPLITTING);
    assert(loaded.get_split_params() == rp.get_split_params());

    for (ulint i = 0; i < loaded.intervals(); ++i) {
        assert(loaded.get_length(i) == rp.get_length(i));
        assert(loaded.get<TestRunCols::VAL1>(i) == rp.get<TestRunCols::VAL1>(i));
        assert(loaded.get<TestRunCols::VAL2>(i) == rp.get<TestRunCols::VAL2>(i));
        assert(loaded.get_row(i) == rp.get_row(i));
    }

    for (ulint idx = 0; idx < domain; ++idx) {
        auto pos = make_pos_absolute<RP>(loaded, idx);
        auto next_pos = loaded.next(pos);
        assert(next_pos.idx == perm[idx]);
    }
}

static void test_runperm_next_with_steps_and_pred_succ() {
    // Two intervals, each with distinct VAL1 values to search.
    const vector<ulint> perm = {1, 0, 3, 2};
    const ulint domain = static_cast<ulint>(perm.size());
    auto [lengths, interval_perm] = get_permutation_intervals(perm);

    vector<TestRunData> run_data(lengths.size());
    for (size_t i = 0; i < lengths.size(); ++i) {
        run_data[i] = {static_cast<ulint>(i), static_cast<ulint>(10 + i)};
    }

    using RP = RunPermSeparatedAbsolute<TestRunCols>;
    RP rp(lengths, interval_perm, run_data);

    // next with steps: compare to repeated single-step.
    for (ulint idx = 0; idx < domain; ++idx) {
        auto pos = make_pos_absolute<RP>(rp, idx);

        auto pos_steps = rp.next(pos, 2);
        auto pos_iter = pos;
        pos_iter = rp.next(pos_iter);
        pos_iter = rp.next(pos_iter);

        assert(pos_steps.idx == pos_iter.idx);
    }

    // pred/succ on VAL1 column.
    using Col = TestRunCols;
    auto start = rp.first();

    // succ from first, looking for VAL1 == 1 should find interval 1.
    auto succ1 = rp.succ<Col::VAL1>(start, 1);
    assert(succ1.has_value());
    assert(succ1->interval == 1);

    // pred from last, looking for VAL1 == 0 should find interval 0.
    auto last = rp.last();
    auto pred0 = rp.pred<Col::VAL1>(last, 0);
    assert(pred0.has_value());
    assert(pred0->interval == 0);

    // succ from last with value not present should return nullopt.
    auto succ_missing = rp.succ<Col::VAL1>(last, 42);
    assert(!succ_missing.has_value());

    // pred from first with value not present should return nullopt.
    auto pred_missing = rp.pred<Col::VAL1>(start, 42);
    assert(!pred_missing.has_value());
}

int main() {
    test_runperm_separated_absolute_basic_mapping_and_run_data();
    test_runperm_up_down_navigation();
    test_runperm_serialize_roundtrip_separated_absolute();
    test_runperm_next_with_steps_and_pred_succ();

    std::cout << "runperm unit tests passed" << std::endl;
    return 0;
}
