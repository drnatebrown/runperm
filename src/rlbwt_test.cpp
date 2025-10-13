#include <iostream>
#include <cassert>
#include <sstream>
#include <optional>
#include <vector>
#include <chrono>
#include <functional>

#include "rlbwt.hpp"

void test_move_lf(std::vector<uchar> bwt_heads, std::vector<ulint> bwt_run_lengths) {
    std::cout << "Testing MoveLF" << std::endl;
    MoveLF<> move_lf(bwt_heads, bwt_run_lengths);

    std::cout << "Size: " << move_lf.size() << std::endl;
    std::cout << "Move Runs: " << move_lf.move_runs() << std::endl;
    std::cout << "Permutation Runs: " << move_lf.permutation_runs() << std::endl;
    using Position = typename MoveLF<>::Position;
    Position pos = move_lf.get_position();
    for (size_t i = 0; i <= move_lf.size(); ++i) {
        std::cout << "Position: " << pos.interval << ", " << pos.offset << " --> ";
        std::cout << "Run Data: " << move_lf.get_character() << std::endl;
        move_lf.next();
        pos = move_lf.get_position();
    }
}

void test_move_fl(std::vector<uchar> bwt_heads, std::vector<ulint> bwt_run_lengths) {
    std::cout << "Testing MoveFL" << std::endl;
    MoveFL<> move_fl(bwt_heads, bwt_run_lengths);

    std::cout << "Size: " << move_fl.size() << std::endl;
    std::cout << "Move Runs: " << move_fl.move_runs() << std::endl;
    std::cout << "Permutation Runs: " << move_fl.permutation_runs() << std::endl;
    using Position = typename MoveFL<>::Position;
    Position pos = move_fl.get_position();
    for (size_t i = 0; i <= move_fl.size(); ++i) {
        std::cout << "Position: " << pos.interval << ", " << pos.offset << " --> ";
        std::cout << "Run Data: " << move_fl.get_character() << std::endl;
        move_fl.next();
        pos = move_fl.get_position();
    }
}

void test_runperm_lf(std::vector<uchar> bwt_heads, std::vector<ulint> bwt_run_lengths) {
    std::cout << "Testing RunPermLF" << std::endl;
    
    // Define dummy run data similar to RunPerm tests
    enum class RunData {
        VAL_1,
        VAL_2,
        COUNT
    };
    static constexpr size_t NUM_FIELDS = static_cast<size_t>(RunData::COUNT);
    std::vector<std::array<ulint, NUM_FIELDS>> run_data(bwt_heads.size());
    
    // Fill with dummy data
    for (size_t i = 0; i < bwt_heads.size(); ++i) {
        run_data[i][0] = i * 10;  // Some dummy value
        run_data[i][1] = i * 5;   // Another dummy value
    }
    
    // Try the constructor without SplitParams
    RunPermLF<RunData> runperm_lf(bwt_heads, bwt_run_lengths, run_data);
    
    std::cout << "Size: " << runperm_lf.size() << std::endl;
    std::cout << "Move Runs: " << runperm_lf.move_runs() << std::endl;
    std::cout << "Permutation Runs: " << runperm_lf.permutation_runs() << std::endl;
    
    runperm_lf.first();
    using Position = typename RunPermLF<RunData>::Position;
    Position pos = runperm_lf.get_position();
    
    for (size_t i = 0; i <= runperm_lf.size(); ++i) {
        std::cout << "Position: " << pos.interval << ", " << pos.offset << " --> ";
        std::cout << "Character: " << (char)runperm_lf.get_character() << " ";
        std::cout << "Run Data: " << runperm_lf.template get<RunData::VAL_1>() << ", " 
                  << runperm_lf.template get<RunData::VAL_2>() << std::endl;
        runperm_lf.LF();
        pos = runperm_lf.get_position();
    }
}

void test_runperm_fl(std::vector<uchar> bwt_heads, std::vector<ulint> bwt_run_lengths) {
    std::cout << "Testing RunPermFL" << std::endl;
    
    // Define dummy run data similar to RunPerm tests
    enum class RunData {
        VAL_1,
        VAL_2,
        COUNT
    };
    static constexpr size_t NUM_FIELDS = static_cast<size_t>(RunData::COUNT);
    std::vector<std::array<ulint, NUM_FIELDS>> run_data(bwt_heads.size());
    
    // Fill with dummy data
    for (size_t i = 0; i < bwt_heads.size(); ++i) {
        run_data[i][0] = i * 7;   // Some dummy value
        run_data[i][1] = i * 3;   // Another dummy value
    }
    
    // Try the constructor without SplitParams
    RunPermFL<RunData> runperm_fl(bwt_heads, bwt_run_lengths, run_data);
    
    std::cout << "Size: " << runperm_fl.size() << std::endl;
    std::cout << "Move Runs: " << runperm_fl.move_runs() << std::endl;
    std::cout << "Permutation Runs: " << runperm_fl.permutation_runs() << std::endl;
    
    runperm_fl.first();
    using Position = typename RunPermFL<RunData>::Position;
    Position pos = runperm_fl.get_position();
    
    for (size_t i = 0; i <= runperm_fl.size(); ++i) {
        std::cout << "Position: " << pos.interval << ", " << pos.offset << " --> ";
        std::cout << "Character: " << (char)runperm_fl.get_character() << " ";
        std::cout << "Run Data: " << runperm_fl.template get<RunData::VAL_1>() << ", " 
                  << runperm_fl.template get<RunData::VAL_2>() << std::endl;
        runperm_fl.FL();
        pos = runperm_fl.get_position();
    }
}

int main() {
    // TEXT: GATTACATGATTACATAGATTACATT$
    // BWT:  TTTTTCCCGGGAAAT$ATTTTAAAAAA
    // RLBWT: TCGAT$ATA
    std::vector<uchar> bwt_heads =       {'T','C','G','A','T', 1 ,'A','T','A'};
    std::vector<ulint> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };
    test_move_lf(bwt_heads, bwt_run_lengths);\
    std::cout << std::endl;
    test_move_fl(bwt_heads, bwt_run_lengths);
    std::cout << std::endl;
    test_runperm_lf(bwt_heads, bwt_run_lengths);
    std::cout << std::endl;
    test_runperm_fl(bwt_heads, bwt_run_lengths);
    std::cout << std::endl;
    return 0;
}