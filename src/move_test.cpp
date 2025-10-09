/* These tests were made using generative AI, they need to be improved upon
    I wrote some of the functions myself, but need to manually write better
    tests and benchmarks. Leaving them here for now since they somewhat 
    verify correctness. Better unit tests and benchmarks should be written. */

#include "move/move_structure.hpp"
#include "move/move_table.hpp"
#include "runperm/runperm.hpp"

#include <iostream>
#include <cassert>
#include <sstream>
#include <optional>
#include <vector>
#include <chrono>
#include <functional>

using namespace std;
using namespace std::chrono;

std::vector<ulint> random_permutation(size_t n) {
    std::vector<ulint> permutation(n);
    for (size_t i = 0; i < n; ++i) {
        permutation[i] = i;
    }
    std::random_shuffle(permutation.begin(), permutation.end());
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
        std::random_shuffle(possible_breaks.begin(), possible_breaks.end());
        
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
    std::random_shuffle(interval_order.begin(), interval_order.end());
    
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
std::optional<ulint> max_allowed_length = 255;
// std::optional<ulint> max_allowed_length = std::nullopt;

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

template<typename RunData, typename RunPermType>
std::string get_runperm_type_name() {
    if constexpr (std::is_same_v<RunPermType, RunPerm<RunData, true, true>>) {
        return "RunPermIntegratedAbsolute";
    } else if constexpr (std::is_same_v<RunPermType, RunPerm<RunData, true, false>>) {
        return "RunPermIntegratedRelative";
    } else if constexpr (std::is_same_v<RunPermType, RunPerm<RunData, false, true>>) {
        return "RunPermSeperatedAbsolute";
    } else if constexpr (std::is_same_v<RunPermType, RunPerm<RunData, false, false>>) {
        return "RunPermSeperatedRelative";
    } else {
        return "Unknown";
    }
}

// Test function template
template<typename MoveStructType>
void test_move_structure(const std::vector<ulint>& lengths, 
                        const std::vector<ulint>& interval_permutation,
                        const std::vector<ulint>& starts,
                        const std::vector<ulint>& pointers,
                        const std::vector<ulint>& offsets,
                        const std::vector<ulint>& test_perm,
                        const std::vector<std::pair<size_t, size_t>>& full_cycle_pos,
                        size_t n) {
    auto start_time = high_resolution_clock::now();
    
    SplitParams split_params;
    split_params.max_allowed_length = max_allowed_length;

    // Create move structure
    auto move_structure = MoveStructType(lengths, interval_permutation, n);
    
    auto creation_time = high_resolution_clock::now();
    
    // Test getters
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
        // ulint pointer = move_structure.get_pointer(i);
        // ulint offset = move_structure.get_offset(i);
        // if constexpr (std::is_same_v<MoveStructType, MoveStructureTblIdx> || 
        //               std::is_same_v<MoveStructType, MoveStructureVecIdx>) {
        //     ulint start = move_structure.get_start(i);
        // }
        // else {
        //     ulint length = move_structure.get_length(i);
        // }
    }
    
    auto getter_time = high_resolution_clock::now();
    
    // Test move operations
    typename MoveStructType::Position pos;
    if constexpr (std::is_same_v<MoveStructType, MoveStructureTblIdx> || 
                  std::is_same_v<MoveStructType, MoveStructureVecIdx>) {
        pos = {0, 0, 0};  // idx types have 3 fields
    } else {
        pos = {0, 0};     // normal types have 2 fields
    }

    ulint real_pos = 0;
    
    for (size_t i = 0; i < n; ++i) {
        size_t last_real_pos = real_pos;
        real_pos = test_perm[real_pos];
        pos = move_structure.move(pos);
        // if constexpr (std::is_same_v<MoveStructType, MoveStructureTblIdx> || 
        //               std::is_same_v<MoveStructType, MoveStructureVecIdx>) {
        //     assert(pos.idx == real_pos);
        // }
        // else {
        //     assert(pos.interval == full_cycle_pos[last_real_pos].first);
        //     assert(pos.offset == full_cycle_pos[last_real_pos].second);
        // }
    }
    
    auto move_time = high_resolution_clock::now();
    
    // Print timing results
    auto creation_duration = duration_cast<microseconds>(creation_time - start_time);
    auto getter_duration = duration_cast<microseconds>(getter_time - creation_time);
    auto move_duration = duration_cast<microseconds>(move_time - getter_time);
    auto total_duration = duration_cast<microseconds>(move_time - start_time);
    
    std::cout << "  " << get_type_name<MoveStructType>() << ":" << std::endl;
    std::cout << "    Creation: " << creation_duration.count() << "μs" << std::endl;
    std::cout << "    Getters: " << getter_duration.count() << "μs" << std::endl;
    std::cout << "    Moves: " << move_duration.count() << "μs" << std::endl;
    std::cout << "    Total: " << total_duration.count() << "μs" << std::endl;
    std::stringstream ss;
    std::cout << "    Size: " << move_structure.serialize(ss) << std::endl;
}

template<typename RunData, typename RunPermType>
void test_runperm(const std::vector<ulint>& lengths, 
                  const std::vector<ulint>& interval_permutation, 
                  const std::vector<std::array<ulint, 2>>& run_data, 
                  size_t n) {
    
    SplitParams split_params;
    split_params.max_allowed_length = max_allowed_length;

    auto start_time = high_resolution_clock::now();
    auto runperm = RunPermType(lengths, interval_permutation, n, run_data);
    auto creation_time = high_resolution_clock::now();

    runperm.first();
    for (size_t i = 0; i < n; ++i) {
        typename RunPermType::Position pos = runperm.get_position();
        // assert(runperm.template get<RunData::VAL_1>() == run_data[pos.interval][0]);
        // assert(runperm.template get<RunData::VAL_2>() == run_data[pos.interval][1]);
        runperm.next();
    }

    auto move_time = high_resolution_clock::now();

    auto creation_duration = duration_cast<microseconds>(creation_time - start_time);
    auto move_duration = duration_cast<microseconds>(move_time - creation_time);
    auto total_duration = duration_cast<microseconds>(move_time - start_time);

    std::cout << "  " << get_runperm_type_name<RunData, RunPermType>() << ":" << std::endl;
    std::cout << "    Creation: " << creation_duration.count() << "μs" << std::endl;
    std::cout << "    Moves: " << move_duration.count() << "μs" << std::endl;
    std::cout << "    Total: " << total_duration.count() << "μs" << std::endl;
    std::stringstream ss;
    std::cout << "    Size: " << runperm.serialize(ss) << std::endl;
}

void tests() {
    cout << "=== MoveStructure Table Tests ===" << endl << endl;

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

            vector<std::pair<size_t, size_t>> full_cycle_pos(n);
            for (size_t i = 0; i < n; ++i) {
                ulint result = test_perm[i];
                auto it = std::upper_bound(starts.begin(), starts.end(), result);
                if (it != starts.begin()) {
                    size_t index = std::distance(starts.begin(), it) - 1;
                    size_t element = starts[index];
                    full_cycle_pos[i] = {index, result - element};
                }
            }

            std::cout << "Testing n=" << n << ", r=" << r << " (n/r=" << n/r << "):" << std::endl;
            
            // Test all table types
            test_move_structure<MoveStructureTbl>(lengths, interval_permutation, starts, pointers, offsets, test_perm, full_cycle_pos, n);
            test_move_structure<MoveStructureVec>(lengths, interval_permutation, starts, pointers, offsets, test_perm, full_cycle_pos, n);
            test_move_structure<MoveStructureTblIdx>(lengths, interval_permutation, starts, pointers, offsets, test_perm, full_cycle_pos, n);
            test_move_structure<MoveStructureVecIdx>(lengths, interval_permutation, starts, pointers, offsets, test_perm, full_cycle_pos, n);
            
            std::cout << std::endl;
        }
    }

    std::cout << "=== RunPerm Tests ===" << endl << endl;

    for (size_t n : test_n) {
        for (size_t percentage_run : percentage_runs) {
            size_t r = (n * percentage_run) / 100;
            std::vector<ulint> test_perm = random_runny_permutation(n, r);
            if (!verify_permutation(test_perm)) {
                std::cout << "Test permutation is not correct" << std::endl;
                continue;
            }
                
            auto [lengths, interval_permutation] = get_permutation_intervals(test_perm);

            size_t val1_max = MAX_VAL(32);
            size_t val2_max = MAX_VAL(55);
            
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

            test_runperm<RunData, RunPerm<RunData, true, true>>(lengths, interval_permutation, run_data, n);
            test_runperm<RunData, RunPerm<RunData, true, false>>(lengths, interval_permutation, run_data, n);
            test_runperm<RunData, RunPerm<RunData, false, true>>(lengths, interval_permutation, run_data, n);
            test_runperm<RunData, RunPerm<RunData, false, false>>(lengths, interval_permutation, run_data, n);
            std::cout << std::endl;
        }
    }
}

int main() {
    tests();
    std::cout << "Tests complete" << std::endl;

    // std::cout << "=== RunPerm Tests ===" << endl << endl;

    // enum class RunData {
    //     EASY,
    //     HARD,
    //     COUNT
    // };
    // std::vector<std::array<ulint, static_cast<size_t>(RunData::COUNT)>> run_data = {{
    //     {1, 2},
    //     {3, 4}
    // }};
    // std::vector<ulint> permutation = {6, 7, 8, 9, 10, 0, 1, 2, 3, 4, 5};
    // std::cout << "Permutation: ";
    // for (size_t i = 0; i < permutation.size(); ++i) {
    //     std::cout << permutation[i] << " ";
    // }
    // std::cout << std::endl;
    // auto [lengths, interval_permutation] = get_permutation_intervals(permutation);
    // RunPerm<RunData> run_perm(lengths, interval_permutation, run_data);
    // run_perm.first();
    // std::cout << "Size: " << run_perm.size() << std::endl;
    // std::cout << "Move Runs: " << run_perm.move_runs() << std::endl;
    // std::cout << "Permutation Runs: " << run_perm.permutation_runs() << std::endl;
    // using Position = typename RunPerm<RunData>::Position;
    // Position pos = run_perm.get_position();
    // for (size_t i = 0; i <= run_perm.size(); ++i) {
    //     std::cout << "Position: " << pos.interval << ", " << pos.offset << ", " << pos.idx << " --> ";
    //     run_perm.next();
    //     pos = run_perm.get_position();

    //     std::cout << "Run Data: " << run_perm.template get<RunData::EASY>() << ", " << run_perm.template get<RunData::HARD>() << std::endl;
    // }
    
    // std::cout << "=== MovePerm Tests ===" << endl << endl;
    // std::cout << "Permutation: ";
    // for (size_t i = 0; i < permutation.size(); ++i) {
    //     std::cout << permutation[i] << " ";
    // }
    // using PositionMove = typename MovePerm<>::Position;
    // std::cout << std::endl;
    // MovePerm<> move_perm(permutation);
    // move_perm.first();
    // std::cout << "Size: " << move_perm.size() << std::endl;
    // std::cout << "Move Runs: " << move_perm.move_runs() << std::endl;
    // std::cout << "Permutation Runs: " << move_perm.permutation_runs() << std::endl;
    // PositionMove pos2 = move_perm.get_position();
    // for (size_t i = 0; i <= move_perm.size(); ++i) {
    //     std::cout << "Position: " << pos2.interval << ", " << pos2.offset << std::endl;
    //     move_perm.next();
    //     pos2 = move_perm.get_position();
    // }
    return 0;
}