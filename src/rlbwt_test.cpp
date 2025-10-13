#include <iostream>
#include <cassert>
#include <sstream>
#include <optional>
#include <vector>
#include <chrono>
#include <functional>

#include "rlbwt.hpp"

void print_vector(std::vector<ulint> vec) {
    std::cout << "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        std::cout << vec[i];
        if (i < vec.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;
}

void test_move_lf(std::vector<uchar> bwt_heads, std::vector<ulint> bwt_run_lengths, std::string text) {
    std::cout << "Testing MoveLF" << std::endl;
    MoveLF<> move_lf(bwt_heads, bwt_run_lengths);

    std::cout << "Size: " << move_lf.size() << std::endl;
    std::cout << "Move Runs: " << move_lf.move_runs() << std::endl;
    std::cout << "Permutation Runs: " << move_lf.permutation_runs() << std::endl;
    using Position = typename MoveLF<>::Position;
    Position pos = move_lf.get_position();
    for (size_t i = 0; i < move_lf.size(); ++i) {
        move_lf.next();
        pos = move_lf.get_position();
    }
    assert(pos.interval == 0);
    assert(pos.offset == 0);
    std::cout << "Successfully moved to start" << std::endl;

    std::string recovered_text(text.size(), '\0');
    move_lf.first();
    for (size_t i = 1; i < move_lf.size(); ++i) {
        recovered_text[text.size() - i] = (char)move_lf.get_character();
        move_lf.LF();
    }
    std::cout << "Recovered text: " << recovered_text << std::endl;
    std::cout << "Original text:  " << text << std::endl;
    assert(recovered_text.compare(text) == 0);
    std::cout << "Successfully recovered text" << std::endl;
}

void test_move_fl(std::vector<uchar> bwt_heads, std::vector<ulint> bwt_run_lengths, std::string text) {
    std::cout << "Testing MoveFL" << std::endl;
    MoveFL<> move_fl(bwt_heads, bwt_run_lengths);

    std::cout << "Size: " << move_fl.size() << std::endl;
    std::cout << "Move Runs: " << move_fl.move_runs() << std::endl;
    std::cout << "Permutation Runs: " << move_fl.permutation_runs() << std::endl;
    using Position = typename MoveFL<>::Position;
    Position pos = move_fl.get_position();
    for (size_t i = 0; i < move_fl.size(); ++i) {
        move_fl.next();
        pos = move_fl.get_position();
    }
    assert(pos.interval == 0);
    assert(pos.offset == 0);
    std::cout << "Successfully moved to start" << std::endl;

    std::string recovered_text = "";
    move_fl.first();
    move_fl.FL();
    for (size_t i = 1; i < move_fl.size(); ++i) {
        recovered_text += (char)move_fl.get_character();
        move_fl.FL();
    }
    std::cout << "Recovered text: " << recovered_text << std::endl;
    std::cout << "Original text:  " << text << std::endl;
    assert(recovered_text.compare(text) == 0);
    std::cout << "Successfully recovered text" << std::endl;
}

void test_runperm_lf(std::vector<uchar> bwt_heads, std::vector<ulint> bwt_run_lengths, std::string text) {
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
    
    for (size_t i = 0; i < runperm_lf.size(); ++i) {
        runperm_lf.LF();
        pos = runperm_lf.get_position();
    }
    assert(pos.interval == 0);
    assert(pos.offset == 0);
    std::cout << "Successfully moved to start" << std::endl;

    std::string recovered_text(text.size(), '\0');
    runperm_lf.first();
    for (size_t i = 1; i < runperm_lf.size(); ++i) {
        recovered_text[text.size() - i] = (char)runperm_lf.get_character();
        runperm_lf.LF();
    }
    std::cout << "Recovered text: " << recovered_text << std::endl;
    std::cout << "Original text:  " << text << std::endl;
    assert(recovered_text.compare(text) == 0);
    std::cout << "Successfully recovered text" << std::endl;
}

void test_runperm_fl(std::vector<uchar> bwt_heads, std::vector<ulint> bwt_run_lengths, std::string text) {
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
    
    for (size_t i = 0; i < runperm_fl.size(); ++i) {
        runperm_fl.FL();
        pos = runperm_fl.get_position();
    }
    assert(pos.interval == 0);
    assert(pos.offset == 0);
    std::cout << "Successfully moved to start" << std::endl;

    runperm_fl.first();
    runperm_fl.FL();
    std::string recovered_text = "";
    for (size_t i = 1; i < runperm_fl.size(); ++i) {
        recovered_text += (char)runperm_fl.get_character();
        runperm_fl.FL();
    }
    std::cout << "Recovered text: " << recovered_text << std::endl;
    std::cout << "Original text:  " << text << std::endl;
    assert(recovered_text.compare(text) == 0);
    std::cout << "Successfully recovered text" << std::endl;
}

void test_runperm_phi(std::vector<uchar> bwt_heads, std::vector<ulint> bwt_run_lengths, std::vector<ulint> sa) {
    std::cout << "Testing RunPermPhi" << std::endl;
    
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

    auto [phi_lengths, phi_interval_permutations, domain] = rlbwt_to_phi(bwt_heads, bwt_run_lengths);
    RunPermPhi<RunData> runperm_phi(phi_lengths, phi_interval_permutations, domain, run_data);
    std::cout << "Size: " << runperm_phi.size() << std::endl;
    std::cout << "Move Runs: " << runperm_phi.move_runs() << std::endl;
    std::cout << "Permutation Runs: " << runperm_phi.permutation_runs() << std::endl;

    runperm_phi.first();
    using Position = typename RunPermPhi<RunData>::Position;
    Position pos = runperm_phi.get_position();
    for (size_t i = 0; i < runperm_phi.size(); ++i) {
        runperm_phi.next();
        pos = runperm_phi.get_position();
    }
    assert(pos.interval == 0);
    assert(pos.offset == 0);
    std::cout << "Successfully moved to start" << std::endl;

    runperm_phi.last();
    runperm_phi.Phi();
    std::vector<ulint> sa_recovered(sa.size());
    for (size_t i = 0; i < sa_recovered.size(); ++i) {
        sa_recovered[sa.size() - i - 1] = runperm_phi.SA();
        runperm_phi.Phi();
        pos = runperm_phi.get_position();
    }
    std::cout << "Recovered SA: ";
    print_vector(sa_recovered);
    std::cout << "Original SA:  ";
    print_vector(sa);
    for (size_t i = 0; i < sa_recovered.size(); ++i) {
        assert(sa_recovered[i] == sa[i]);
    }
    std::cout << "Successfully moved to end" << std::endl;
}

void test_move_phi(std::vector<uchar> bwt_heads, std::vector<ulint> bwt_run_lengths, std::vector<ulint> sa) {
    std::cout << "Testing MovePhi" << std::endl;
    auto [phi_lengths, phi_interval_permutations, domain] = rlbwt_to_phi(bwt_heads, bwt_run_lengths);
    MovePhi move_phi(phi_lengths, phi_interval_permutations, domain);
    std::cout << "Size: " << move_phi.size() << std::endl;
    std::cout << "Move Runs: " << move_phi.move_runs() << std::endl;
    std::cout << "Permutation Runs: " << move_phi.permutation_runs() << std::endl;
    using Position = typename MovePhi::Position;
    Position pos = move_phi.get_position();
    for (size_t i = 0; i < move_phi.size(); ++i) {
        move_phi.next();
        pos = move_phi.get_position();
    }
    assert(pos.interval == 0);
    assert(pos.offset == 0);
    std::cout << "Successfully moved to start" << std::endl;

    move_phi.last();
    move_phi.Phi();
    std::vector<ulint> sa_recovered(sa.size());
    for (size_t i = 0; i < move_phi.size(); ++i) {
        sa_recovered[sa.size() - i - 1] = move_phi.SA();
        move_phi.Phi();
        pos = move_phi.get_position();
    }
    std::cout << "Recovered SA: ";
    print_vector(sa_recovered);
    std::cout << "Original SA:  ";
    print_vector(sa);
    for (size_t i = 0; i < sa_recovered.size(); ++i) {
        assert(sa_recovered[i] == sa[i]);
    }
    std::cout << "Successfully moved to end" << std::endl;
}

void test_runperm_invphi(std::vector<uchar> bwt_heads, std::vector<ulint> bwt_run_lengths, std::vector<ulint> sa) {
    std::cout << "Testing RunPermInvPhi" << std::endl;
    
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

    auto [invphi_lengths, invphi_interval_permutations, domain] = rlbwt_to_invphi(bwt_heads, bwt_run_lengths);
    RunPermInvPhi<RunData> runperm_invphi(invphi_lengths, invphi_interval_permutations, domain, run_data);
    std::cout << "Size: " << runperm_invphi.size() << std::endl;
    std::cout << "Move Runs: " << runperm_invphi.move_runs() << std::endl;
    std::cout << "Permutation Runs: " << runperm_invphi.permutation_runs() << std::endl;

    runperm_invphi.first();
    using Position = typename RunPermPhi<RunData>::Position;
    Position pos = runperm_invphi.get_position();
    for (size_t i = 0; i < runperm_invphi.size(); ++i) {
        runperm_invphi.next();
        pos = runperm_invphi.get_position();
    }
    assert(pos.interval == 0);
    assert(pos.offset == 0);
    std::cout << "Successfully moved to start" << std::endl;

    runperm_invphi.last();
    std::vector<ulint> sa_recovered(sa.size());
    for (size_t i = 0; i < sa_recovered.size(); ++i) {
        sa_recovered[i] = runperm_invphi.SA();
        runperm_invphi.InvPhi();
        pos = runperm_invphi.get_position();
    }
    std::cout << "Recovered SA: ";
    print_vector(sa_recovered);
    std::cout << "Original SA:  ";
    print_vector(sa);
    for (size_t i = 0; i < sa_recovered.size(); ++i) {
        assert(sa_recovered[i] == sa[i]);
    }
    std::cout << "Successfully recovered SA" << std::endl;
}

void test_move_invphi(std::vector<uchar> bwt_heads, std::vector<ulint> bwt_run_lengths, std::vector<ulint> sa) {
    std::cout << "Testing MoveInvPhi" << std::endl;
    auto [invphi_lengths, invphi_interval_permutations, domain] = rlbwt_to_invphi(bwt_heads, bwt_run_lengths);
    MoveInvPhi move_invphi(invphi_lengths, invphi_interval_permutations, domain);
    std::cout << "Size: " << move_invphi.size() << std::endl;
    std::cout << "Move Runs: " << move_invphi.move_runs() << std::endl;
    std::cout << "Permutation Runs: " << move_invphi.permutation_runs() << std::endl;
    using Position = typename MoveInvPhi::Position;
    Position pos = move_invphi.get_position();
    for (size_t i = 0; i < move_invphi.size(); ++i) {
        move_invphi.next();
        pos = move_invphi.get_position();
    }
    assert(pos.interval == 0);
    assert(pos.offset == 0);
    std::cout << "Successfully moved to start" << std::endl;

    move_invphi.last();
    std::vector<ulint> sa_recovered(sa.size());
    for (size_t i = 0; i < move_invphi.size(); ++i) {
        sa_recovered[i] = move_invphi.SA();
        move_invphi.InvPhi();
        pos = move_invphi.get_position();
    }
    std::cout << "Recovered SA: ";
    print_vector(sa_recovered);
    std::cout << "Original SA:  ";
    print_vector(sa);
    for (size_t i = 0; i < sa_recovered.size(); ++i) {
        assert(sa_recovered[i] == sa[i]);
    }
    std::cout << "Successfully recovered SA" << std::endl;
}

int main() {
    // TEXT: GATTACATGATTACATAGATTACATT$
    // BWT:  TTTTTCCCGGGAAAT$ATTTTAAAAAA
    // RLBWT: TCGAT$ATA
    std::string text = "GATTACATGATTACATAGATTACATT";
    std::vector<ulint> sa = {26, 12, 4, 21, 16, 14, 6, 23, 9, 1, 18, 13, 5, 22, 8, 0, 17, 25, 11, 3, 20, 15, 7, 24, 10, 2, 19};
    std::vector<uchar> bwt_heads =       {'T','C','G','A','T', 1 ,'A','T','A'};
    std::vector<ulint> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };
    test_move_lf(bwt_heads, bwt_run_lengths, text);
    std::cout << std::endl;
    test_move_fl(bwt_heads, bwt_run_lengths, text);
    std::cout << std::endl;
    test_runperm_lf(bwt_heads, bwt_run_lengths, text);
    std::cout << std::endl;
    test_runperm_fl(bwt_heads, bwt_run_lengths, text);
    std::cout << std::endl;
    test_runperm_phi(bwt_heads, bwt_run_lengths, sa);
    std::cout << std::endl;
    test_move_phi(bwt_heads, bwt_run_lengths, sa);
    std::cout << std::endl;
    test_runperm_invphi(bwt_heads, bwt_run_lengths, sa);
    std::cout << std::endl;
    test_move_invphi(bwt_heads, bwt_run_lengths, sa);
    std::cout << std::endl;
    return 0;
}