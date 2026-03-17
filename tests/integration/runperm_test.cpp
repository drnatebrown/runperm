// Integration-style tests for RunPerm (generic, non-RLBWT).
// These are simple assert-based tests, no external framework.

#include "runperm.hpp"

#include <cassert>
#include <iostream>
#include <vector>

using std::size_t;
using std::vector;

// Simple run-data schema with two fields.
enum class TestRunColsInt {
    VAL1,
    VAL2,
    COUNT
};

using TestRunDataInt = ColumnsTuple<TestRunColsInt>;

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

static void integration_runperm_separated_and_integrated_absolute() {
    // Use the same run-level representation as in examples.cpp example1.
    const vector<ulint> lengths = {2, 3, 1, 2, 2, 1, 1, 1, 3};
    const vector<ulint> interval_perm = {1, 9, 3, 12, 4, 14, 0, 15, 6};
    const ulint domain = 16;

    // Reconstruct the underlying position-level permutation induced by
    // (lengths, interval_perm).
    vector<ulint> full_perm(domain);
    ulint src_idx = 0;
    for (size_t j = 0; j < lengths.size(); ++j) {
        for (ulint o = 0; o < lengths[j]; ++o) {
            full_perm[src_idx++] = interval_perm[j] + o;
        }
    }

    vector<TestRunDataInt> run_data(lengths.size());
    for (size_t i = 0; i < lengths.size(); ++i) {
        run_data[i] = {static_cast<ulint>(i), static_cast<ulint>(i + 10)};
    }

    using RPSeparatedAbs = RunPermSeparatedAbsolute<TestRunColsInt>;
    using RPIntegratedAbs = RunPermIntegratedAbsolute<TestRunColsInt>;

    RPSeparatedAbs rp_sep(lengths, interval_perm, run_data);
    RPIntegratedAbs rp_int(lengths, interval_perm, run_data);

    assert(rp_sep.domain() == domain);
    assert(rp_int.domain() == domain);
    assert(rp_sep.intervals() == lengths.size());
    assert(rp_int.intervals() == lengths.size());

    // Both configurations must represent the same permutation and run data.
    for (ulint idx = 0; idx < domain; ++idx) {
        auto pos_sep = make_pos_absolute<RPSeparatedAbs>(rp_sep, idx);
        auto pos_int = make_pos_absolute<RPIntegratedAbs>(rp_int, idx);

        auto next_sep = rp_sep.next(pos_sep);
        auto next_int = rp_int.next(pos_int);

        assert(next_sep.idx == full_perm[idx]);
        assert(next_int.idx == full_perm[idx]);

        ulint ival_sep = next_sep.interval;
        ulint ival_int = next_int.interval;

        assert(rp_sep.get<TestRunColsInt::VAL1>(ival_sep) ==
               rp_int.get<TestRunColsInt::VAL1>(ival_int));
        assert(rp_sep.get<TestRunColsInt::VAL2>(ival_sep) ==
               rp_int.get<TestRunColsInt::VAL2>(ival_int));
    }
}

int main() {
    integration_runperm_separated_and_integrated_absolute();

    std::cout << "runperm integration tests passed" << std::endl;
    return 0;
}
