/* These benchmarks were made using generative AI, they need to be improved upon
    I wrote some of the functions myself, but need to manually write better
    benchmarks. Leaving them here for now since they somewhat 
    verify performance. Better benchmarks should be written. */

#include "orbit/runperm.hpp"

#include <iostream>
#include <set>
#include <cassert>
#include <sstream>
#include <optional>
#include <vector>
#include <chrono>
#include <functional>
#include <random>

using namespace orbit;

using namespace std;
using namespace std::chrono;

// Deterministic RNG for reproducible benchmarks
static std::mt19937 rng(42);
    
std::vector<ulint> random_permutation(size_t n) {
    std::vector<ulint> permutation(n);
    for (size_t i = 0; i < n; ++i) {
        permutation[i] = i;
    }
    std::shuffle(permutation.begin(), permutation.end(), rng);
    return permutation;
}

std::vector<ulint> random_runny_permutation(size_t n, size_t r) {
    if (r == 0 || r > n) {
        return random_permutation(n);
    }
    
    // Step 1: Choose r-1 random break points in range [1, n-1] to create r intervals
    std::vector<size_t> break_points;
    if (r > 1) {
        std::vector<size_t> possible_breaks;
        for (size_t i = 1; i < n; ++i) {
            possible_breaks.push_back(i);
        }
        std::shuffle(possible_breaks.begin(), possible_breaks.end(), rng);
        
        for (size_t i = 0; i < r - 1; ++i) {
            break_points.push_back(possible_breaks[i]);
        }
        std::sort(break_points.begin(), break_points.end());
    }
    
    // Step 2: Define the intervals [start, end) for each of the r intervals
    std::vector<std::pair<size_t, size_t>> intervals;
    size_t start = 0;
    for (size_t bp : break_points) {
        intervals.push_back({start, bp});
        start = bp;
    }
    intervals.push_back({start, n}); // Last interval goes to n
    
    // Step 3: Create random order for processing the intervals
    std::vector<size_t> interval_order(r);
    for (size_t i = 0; i < r; ++i) {
        interval_order[i] = i;
    }
    std::shuffle(interval_order.begin(), interval_order.end(), rng);
    
    // Step 4: Fill intervals with consecutive numbers
    std::vector<ulint> result(n);
    ulint current_value = 0;
    
    for (size_t idx = 0; idx < r; ++idx) {
        size_t interval_idx = interval_order[idx];
        size_t start_pos = intervals[interval_idx].first;
        size_t end_pos = intervals[interval_idx].second;
        
        // Fill this interval with consecutive values
        for (size_t pos = start_pos; pos < end_pos; ++pos) {
            result[pos] = current_value++;
        }
    }
    
    return result;
}

bool verify_permutation(const std::vector<ulint>& perm) {
    std::vector<ulint> sorted_perm(perm);
    std::sort(sorted_perm.begin(), sorted_perm.end());
    for (size_t i = 0; i < sorted_perm.size(); ++i) {
        if (sorted_perm[i] != i) {
            return false;
        }
    }
    return true;
}

// std::vector<ulint> test_n = {1048576, 2097152, 4194304, 8388608};
std::vector<ulint> test_n = {18388608};  
std::vector<size_t> percentage_runs = {1, 2, 5, 10};
split_params sp = NO_SPLITTING;
std::optional<double> length_capping_factor = 4.0;

template<typename RunData, typename RunPermType>
std::string get_runperm_type_name() {
    if constexpr (std::is_same_v<RunPermType, runperm<RunData, true, true>>) {
        return "RunPermIntegratedAbsolute";
    } else if constexpr (std::is_same_v<RunPermType, runperm<RunData, true, false>>) {
        return "RunPermIntegratedRelative";
    } else if constexpr (std::is_same_v<RunPermType, runperm<RunData, false, true>>) {
        return "RunPermSeperatedAbsolute";
    } else if constexpr (std::is_same_v<RunPermType, runperm<RunData, false, false>>) {
        return "RunPermSeperatedRelative";
    } else {
        return "Unknown";
    }
}

template<typename RunData, typename RunPermType>
void bench_runperm(const std::vector<ulint>& lengths, 
                   const std::vector<ulint>& images, 
                   const std::vector<std::array<ulint, 2>>& run_data, 
                   size_t n) {

    auto start_time = high_resolution_clock::now();
    auto runperm = RunPermType(lengths, images, sp, run_data);
    auto creation_time = high_resolution_clock::now();

    auto pos = runperm.first();
    for (size_t i = 0; i < n; ++i) {
        pos = runperm.next(pos);
    }
    auto move_time = high_resolution_clock::now();

    pos = runperm.first();
    for (size_t i = 0; i < n; ++i) {
        pos = runperm.next(pos);
        assert(runperm.template get<RunData::VAL_1>(pos) == run_data[pos.interval][0]);
        assert(runperm.template get<RunData::VAL_2>(pos) == run_data[pos.interval][1]);
    }
    auto getter_time = high_resolution_clock::now();

    auto creation_duration = duration_cast<microseconds>(creation_time - start_time);
    auto move_duration = duration_cast<microseconds>(move_time - creation_time);
    auto getter_duration = duration_cast<microseconds>(getter_time - move_time);
    auto total_duration = duration_cast<microseconds>(getter_time - start_time);

    std::cout << "  " << get_runperm_type_name<RunData, RunPermType>() << ":" << std::endl;
    std::cout << "    Creation: " << creation_duration.count() << "μs" << std::endl;
    std::cout << "    Moves: " << move_duration.count() << "μs" << std::endl;
    std::cout << "    Getters: " << getter_duration.count() << "μs" << std::endl;
    std::cout << "    Total: " << total_duration.count() << "μs" << std::endl;
    std::stringstream ss;
    std::cout << "    Size: " << runperm.serialize(ss) << std::endl;
}

// Benchmarks for different RunPerm configurations
void run_benchmarks() {
    std::cout << "=== RunPerm Benchmarks ===" << endl << endl;

    for (size_t n : test_n) {
        for (size_t percentage_run : percentage_runs) {
            size_t r = (n * percentage_run) / 100;
            std::vector<ulint> test_perm = random_runny_permutation(n, r);
            if (!verify_permutation(test_perm)) {
                std::cout << "Test permutation is not correct" << std::endl;
                continue;
            }
                
            auto [lengths, images] = get_permutation_intervals(test_perm);

            size_t val1_max = max_val(32);
            size_t val2_max = max_val(55);
            
            enum class RunData {
                VAL_1,
                VAL_2,
                COUNT
            };
            static constexpr size_t NUM_FIELDS = static_cast<size_t>(RunData::COUNT);
            std::vector<std::array<ulint, NUM_FIELDS>> run_data(lengths.size());
            for (size_t i = 0; i < lengths.size(); ++i) {
                run_data[i][0] = rand() % val1_max;
                run_data[i][1] = rand() % val2_max;
            }

            std::cout << "Testing n=" << n << ", r=" << r << " (n/r=" << n/r << "):" << std::endl;

            bench_runperm<RunData, runperm<RunData, true, true>>(lengths, images, run_data, n);
            bench_runperm<RunData, runperm<RunData, true, false>>(lengths, images, run_data, n);
            bench_runperm<RunData, runperm<RunData, false, true>>(lengths, images, run_data, n);
            bench_runperm<RunData, runperm<RunData, false, false>>(lengths, images, run_data, n);

            std::cout << std::endl;
        }
    }
}

int main() {
    run_benchmarks();
    std::cout << "Benchmarks complete" << std::endl;
    return 0;
}
