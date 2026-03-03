// Small randomized tests for RunPerm invariants.
// These are simple assert-based tests, no external framework.
//
// Goal:
// - For small, random permutations, verify that different RunPerm
//   configurations (separated vs integrated, absolute) all represent
//   the same underlying permutation and expose consistent run data.

#include "runperm.hpp"

#include <cassert>
#include <iostream>
#include <random>
#include <vector>

using std::size_t;
using std::vector;

// Simple run-data schema with two fields.
enum class TestRunCols {
    VAL1,
    VAL2,
    COUNT
};

using TestRunData = DataTuple<TestRunCols>;

// Helper to build an absolute RunPerm position from a global index.
template <typename RP>
static typename RP::Position make_pos_absolute(const RP &rp, ulint idx) {
    using Position = typename RP::Position;
    Position pos{};
    pos.idx = idx;
    ulint prefix = 0;
    for (ulint interval = 0; interval < rp.move_runs(); ++interval) {
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

static vector<ulint> make_random_permutation(size_t n, std::mt19937 &rng) {
    vector<ulint> perm(n);
    for (size_t i = 0; i < n; ++i) {
        perm[i] = static_cast<ulint>(i);
    }
    std::shuffle(perm.begin(), perm.end(), rng);
    return perm;
}

static void test_runperm_random_small_permutations() {
    // Fixed seed for determinism across runs.
    std::mt19937 rng(42);

    // A few small sizes and a handful of permutations for each.
    const vector<size_t> sizes = {10, 31, 50};
    const size_t num_perms_per_size = 5;

    for (size_t n : sizes) {
        const ulint domain = static_cast<ulint>(n);

        for (size_t iter = 0; iter < num_perms_per_size; ++iter) {
            vector<ulint> perm = make_random_permutation(n, rng);

            // Build run decomposition of the permutation.
            auto [lengths, interval_perm] = get_permutation_intervals(perm);

            // Simple, deterministic run data: VAL1 = i, VAL2 = i + 100.
            vector<TestRunData> run_data(lengths.size());
            for (size_t i = 0; i < lengths.size(); ++i) {
                run_data[i] = {static_cast<ulint>(i), static_cast<ulint>(i + 100)};
            }

            using RPSeparatedAbs = RunPermSeparatedAbsolute<TestRunCols>;
            using RPIntegratedAbs = RunPermIntegratedAbsolute<TestRunCols>;

            RPSeparatedAbs rp_sep(lengths, interval_perm, domain, run_data);
            RPIntegratedAbs rp_int(lengths, interval_perm, domain, run_data);

            assert(rp_sep.domain() == domain);
            assert(rp_int.domain() == domain);
            assert(rp_sep.move_runs() == lengths.size());
            assert(rp_int.move_runs() == lengths.size());

            // For every starting index, next() must follow the base permutation,
            // and both configurations must agree on the mapping and run data.
            for (ulint idx = 0; idx < domain; ++idx) {
                auto pos_sep = make_pos_absolute<RPSeparatedAbs>(rp_sep, idx);
                auto pos_int = make_pos_absolute<RPIntegratedAbs>(rp_int, idx);

                auto next_sep = rp_sep.next(pos_sep);
                auto next_int = rp_int.next(pos_int);

                assert(next_sep.idx == perm[static_cast<size_t>(idx)]);
                assert(next_int.idx == perm[static_cast<size_t>(idx)]);

                ulint ival_sep = next_sep.interval;
                ulint ival_int = next_int.interval;

                // Run data seen through a position must correspond to its interval,
                // and both configurations must expose the same values.
                auto val1_sep_pos = rp_sep.template get<TestRunCols::VAL1>(next_sep);
                auto val2_sep_pos = rp_sep.template get<TestRunCols::VAL2>(next_sep);
                auto val1_int_pos = rp_int.template get<TestRunCols::VAL1>(next_int);
                auto val2_int_pos = rp_int.template get<TestRunCols::VAL2>(next_int);

                auto val1_sep_run = rp_sep.template get<TestRunCols::VAL1>(ival_sep);
                auto val2_sep_run = rp_sep.template get<TestRunCols::VAL2>(ival_sep);
                auto val1_int_run = rp_int.template get<TestRunCols::VAL1>(ival_int);
                auto val2_int_run = rp_int.template get<TestRunCols::VAL2>(ival_int);

                assert(val1_sep_pos == val1_sep_run);
                assert(val2_sep_pos == val2_sep_run);
                assert(val1_int_pos == val1_int_run);
                assert(val2_int_pos == val2_int_run);

                assert(val1_sep_pos == val1_int_pos);
                assert(val2_sep_pos == val2_int_pos);
            }
        }
    }
}

int main() {
    test_runperm_random_small_permutations();

    std::cout << "runperm random tests passed" << std::endl;
    return 0;
}
