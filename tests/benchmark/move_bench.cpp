/* These benchmarks were made using generative AI, they need to be improved upon
    I wrote some of the functions myself, but need to manually write better
    benchmarks. Leaving them here for now since they somewhat 
    verify performance. Better benchmarks should be written. */

#include "move.hpp"

#include <iostream>
#include <cassert>
#include <sstream>
#include <optional>
#include <vector>
#include <chrono>
#include <functional>
#include <algorithm>
#include <random>

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
std::optional<double> length_capping_factor = 4.0;
SplitParams split_params = NO_SPLITTING; 

// Helper function to get readable type name
template<typename MoveStructType>
std::string get_type_name() {
    if constexpr (std::is_same_v<MoveStructType, MoveStructureTbl>) {
        return "MoveStructureTbl";
    } else if constexpr (std::is_same_v<MoveStructType, MoveStructureVec>) {
        return "MoveStructureVec";
    } else if constexpr (std::is_same_v<MoveStructType, MoveStructureTblIdx>) {
        return "MoveStructureTblIdx";
    } else if constexpr (std::is_same_v<MoveStructType, MoveStructureVecIdx>) {
        return "MoveStructureVecIdx";
    } else {
        return "Unknown";
    }
}

// Benchmark function template (no splitting: reference data matches unsplit layout)
template<typename MoveStructType>
void bench_move_structure(const std::vector<ulint>& lengths, 
                          const std::vector<ulint>& interval_permutation,
                          const std::vector<ulint>& starts,
                          const std::vector<ulint>& pointers,
                          const std::vector<ulint>& offsets,
                          const std::vector<ulint>& test_perm,
                          size_t n) {
    // Creation
    auto t0 = high_resolution_clock::now();
    auto move_structure = MoveStructType(lengths, interval_permutation, split_params);
    auto t1 = high_resolution_clock::now();

    // Getter phase
    for (size_t i = 0; i < lengths.size(); ++i) {
        assert(move_structure.get_pointer(i) == pointers[i]);
        assert(move_structure.get_offset(i) == offsets[i]);
        
        if constexpr (std::is_same_v<MoveStructType, MoveStructureTblIdx> || 
                      std::is_same_v<MoveStructType, MoveStructureVecIdx>) {
            assert(move_structure.get_start(i) == starts[i]);
        }
        else {
            assert(move_structure.get_length(i) == lengths[i]);
        }
    }
    auto t2 = high_resolution_clock::now();

    // Benchmark move() operations
    typename MoveStructType::Position pos;
    if constexpr (std::is_same_v<MoveStructType, MoveStructureTblIdx> || 
                  std::is_same_v<MoveStructType, MoveStructureVecIdx>) {
        pos = {0, 0, 0};  // idx types have 3 fields
    } else {
        pos = {0, 0};     // normal types have 2 fields
    }

    auto t3 = high_resolution_clock::now();
    for (size_t i = 0; i < n; ++i) {
        pos = move_structure.move(pos);
    }
    auto t4 = high_resolution_clock::now();

    // Benchmark move_exponential() operations (no splitting)
    typename MoveStructType::Position pos_exp;
    if constexpr (std::is_same_v<MoveStructType, MoveStructureTblIdx> || 
                  std::is_same_v<MoveStructType, MoveStructureVecIdx>) {
        pos_exp = {0, 0, 0};
    } else {
        pos_exp = {0, 0};
    }

    auto t5 = high_resolution_clock::now();
    for (size_t i = 0; i < n; ++i) {
        pos_exp = move_structure.move_exponential(pos_exp);
    }
    auto t6 = high_resolution_clock::now();

    // Print timing results
    auto creation_duration   = duration_cast<microseconds>(t1 - t0);
    auto getter_duration     = duration_cast<microseconds>(t2 - t1);
    auto move_duration       = duration_cast<microseconds>(t4 - t3);
    auto move_exp_duration   = duration_cast<microseconds>(t6 - t5);
    auto total_duration      = duration_cast<microseconds>(t6 - t0);
    
    std::cout << "  " << get_type_name<MoveStructType>() << " (no splitting):" << std::endl;
    std::cout << "    Creation:        " << creation_duration.count() << "μs" << std::endl;
    std::cout << "    Getters:         " << getter_duration.count() << "μs" << std::endl;
    std::cout << "    Moves (move):    " << move_duration.count() << "μs" << std::endl;
    std::cout << "    Moves (exp):     " << move_exp_duration.count() << "μs" << std::endl;
    std::cout << "    Total:           " << total_duration.count() << "μs" << std::endl;
    std::stringstream ss;
    std::cout << "    Size: " << move_structure.serialize(ss) << std::endl;
}

// Splitting-aware benchmark: measures both move() and move_exponential() under splitting.
template<typename MoveStructType>
void bench_move_structure_with_splitting(const std::vector<ulint>& lengths, 
                                         const std::vector<ulint>& interval_permutation,
                                         const std::vector<ulint>& test_perm,
                                         size_t n) {
    auto t0 = high_resolution_clock::now();

    SplitParams split_params_split(length_capping_factor, std::nullopt);
    auto move_structure = MoveStructType(lengths, interval_permutation, split_params_split);

    auto t1 = high_resolution_clock::now();

    assert(move_structure.intervals() >= lengths.size() && "splitting can only add runs");

    // Benchmark move() in splitting mode
    typename MoveStructType::Position pos;
    if constexpr (std::is_same_v<MoveStructType, MoveStructureTblIdx> || 
                  std::is_same_v<MoveStructType, MoveStructureVecIdx>) {
        pos = {0, 0, 0};
    } else {
        pos = {0, 0};
    }

    auto t2 = high_resolution_clock::now();
    for (size_t i = 0; i < n; ++i) {
        pos = move_structure.move(pos);
    }
    auto t3 = high_resolution_clock::now();

    // Benchmark move_exponential() in splitting mode
    typename MoveStructType::Position pos_exp;
    if constexpr (std::is_same_v<MoveStructType, MoveStructureTblIdx> || 
                  std::is_same_v<MoveStructType, MoveStructureVecIdx>) {
        pos_exp = {0, 0, 0};
    } else {
        pos_exp = {0, 0};
    }

    auto t4 = high_resolution_clock::now();
    for (size_t i = 0; i < n; ++i) {
        pos_exp = move_structure.move_exponential(pos_exp);
    }
    auto t5 = high_resolution_clock::now();

    auto creation_duration   = duration_cast<microseconds>(t1 - t0);
    auto move_duration       = duration_cast<microseconds>(t3 - t2);
    auto move_exp_duration   = duration_cast<microseconds>(t5 - t4);
    auto total_duration      = duration_cast<microseconds>(t5 - t0);

    std::cout << "  " << get_type_name<MoveStructType>() << " (splitting):" << std::endl;
    std::cout << "    Creation:        " << creation_duration.count() << "μs" << std::endl;
    std::cout << "    Moves (move):    " << move_duration.count() << "μs" << std::endl;
    std::cout << "    Moves (exp):     " << move_exp_duration.count() << "μs" << std::endl;
    std::cout << "    Total:           " << total_duration.count() << "μs" << std::endl;
    std::stringstream ss;
    std::cout << "    Size: " << move_structure.serialize(ss) << std::endl;
}

void run_benchmarks() {
    cout << "=== MoveStructure Benchmarks ===" << endl << endl;

    for (size_t n : test_n) {
        for (size_t percentage_run : percentage_runs) {
            size_t r = (n * percentage_run) / 100;
            std::vector<ulint> test_perm = random_runny_permutation(n, r);
            if (!verify_permutation(test_perm)) {
                std::cout << "Test permutation is not correct" << std::endl;
                continue;
            }
                
            auto [lengths, interval_permutation] = get_permutation_intervals(test_perm);
            std::vector<ulint> starts(lengths.size());
            size_t start = 0;
            for (size_t i = 0; i < lengths.size(); ++i) {
                starts[i] = start;
                start += lengths[i];
            }

            vector<ulint> pointers(lengths.size());
            vector<ulint> offsets(lengths.size());
            for (size_t i = 0; i < lengths.size(); ++i) {
                assert(interval_permutation[i] == test_perm[starts[i]]);
                auto it = std::upper_bound(starts.begin(), starts.end(), interval_permutation[i]);
                if (it != starts.begin()) {
                    size_t index = std::distance(starts.begin(), it) - 1;
                    size_t element = starts[index];
                    pointers[i] = index;
                    offsets[i] = interval_permutation[i] - element;
                }
            }

            std::cout << "Testing n=" << n << ", r=" << r << " (n/r=" << n/r << "):" << std::endl;
            
            // Benchmark all table types (no splitting + splitting)
            bench_move_structure<MoveStructureTbl>(lengths, interval_permutation, starts, pointers, offsets, test_perm, n);
            bench_move_structure_with_splitting<MoveStructureTbl>(lengths, interval_permutation, test_perm, n);
            bench_move_structure<MoveStructureVec>(lengths, interval_permutation, starts, pointers, offsets, test_perm, n);
            bench_move_structure_with_splitting<MoveStructureVec>(lengths, interval_permutation, test_perm, n);
            bench_move_structure<MoveStructureTblIdx>(lengths, interval_permutation, starts, pointers, offsets, test_perm, n);
            bench_move_structure_with_splitting<MoveStructureTblIdx>(lengths, interval_permutation, test_perm, n);
            bench_move_structure<MoveStructureVecIdx>(lengths, interval_permutation, starts, pointers, offsets, test_perm, n);
            bench_move_structure_with_splitting<MoveStructureVecIdx>(lengths, interval_permutation, test_perm, n);

            std::cout << std::endl;
        }
    }
}

int main() {
    run_benchmarks();
    return 0;
}
