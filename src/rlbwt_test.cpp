#include <iostream>
#include <cassert>
#include <sstream>
#include <optional>
#include <vector>
#include <chrono>
#include <functional>

#include "rlbwt.hpp"
#include "runperm.hpp"

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

    std::cout << "Domain: " << move_lf.domain() << std::endl;
    std::cout << "Move Runs: " << move_lf.move_runs() << std::endl;
    std::cout << "Permutation Runs: " << move_lf.permutation_runs() << std::endl;
    using Position = typename MoveLF<>::Position;
    auto pos = move_lf.first();
    for (size_t i = 0; i < move_lf.domain(); ++i) {
        pos = move_lf.next(pos);
    }
    assert(pos.interval == 0);
    assert(pos.offset == 0);
    std::cout << "Successfully moved to start" << std::endl;

    std::string recovered_text(text.size(), '\0');
    pos = move_lf.first();
    for (size_t i = 1; i < move_lf.domain(); ++i) {
        recovered_text[text.size() - i] = (char)move_lf.get_character(pos);
        pos = move_lf.next(pos);
    }
    std::cout << "Recovered text: " << recovered_text << std::endl;
    std::cout << "Original text:  " << text << std::endl;
    assert(recovered_text.compare(text) == 0);
    std::cout << "Successfully recovered text" << std::endl;
}

void test_move_fl(std::vector<uchar> bwt_heads, std::vector<ulint> bwt_run_lengths, std::string text) {
    std::cout << "Testing MoveFL" << std::endl;
    MoveFL<> move_fl(bwt_heads, bwt_run_lengths);

    std::cout << "Domain: " << move_fl.domain() << std::endl;
    std::cout << "Move Runs: " << move_fl.move_runs() << std::endl;
    std::cout << "Permutation Runs: " << move_fl.permutation_runs() << std::endl;
    using Position = typename MoveFL<>::Position;
    Position pos = move_fl.first();
    for (size_t i = 0; i < move_fl.domain(); ++i) {
        pos = move_fl.next(pos);
    }
    assert(pos.interval == 0);
    assert(pos.offset == 0);
    std::cout << "Successfully moved to start" << std::endl;

    std::string recovered_text = "";
    pos = move_fl.first();
    pos = move_fl.FL(pos);
    for (size_t i = 1; i < move_fl.domain(); ++i) {
        recovered_text += (char)move_fl.get_character(pos);
        pos = move_fl.FL(pos);
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
    
    std::cout << "Domain: " << runperm_lf.domain() << std::endl;
    std::cout << "Move Runs: " << runperm_lf.move_runs() << std::endl;
    std::cout << "Permutation Runs: " << runperm_lf.permutation_runs() << std::endl;
    
    runperm_lf.first();
    using Position = typename RunPermLF<RunData>::Position;
    Position pos = runperm_lf.first();
    
    for (size_t i = 0; i < runperm_lf.domain(); ++i) {
        pos = runperm_lf.LF(pos);
    }
    assert(pos.interval == 0);
    assert(pos.offset == 0);
    std::cout << "Successfully moved to start" << std::endl;

    std::string recovered_text(text.size(), '\0');
    pos = runperm_lf.first();
    for (size_t i = 1; i < runperm_lf.domain(); ++i) {
        recovered_text[text.size() - i] = (char)runperm_lf.get_character(pos);
        pos = runperm_lf.LF(pos);
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
    
    std::cout << "Domain: " << runperm_fl.domain() << std::endl;
    std::cout << "Move Runs: " << runperm_fl.move_runs() << std::endl;
    std::cout << "Permutation Runs: " << runperm_fl.permutation_runs() << std::endl;
    
    using Position = typename RunPermFL<RunData>::Position;
    Position pos = runperm_fl.first();
    
    for (size_t i = 0; i < runperm_fl.domain(); ++i) {
        pos = runperm_fl.FL(pos);
    }
    assert(pos.interval == 0);
    assert(pos.offset == 0);
    std::cout << "Successfully moved to start" << std::endl;

    pos = runperm_fl.first();
    pos = runperm_fl.FL(pos);
    std::string recovered_text = "";
    for (size_t i = 1; i < runperm_fl.domain(); ++i) {
        recovered_text += (char)runperm_fl.get_character(pos);
        pos = runperm_fl.FL(pos);
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
    std::cout << "Domain: " << runperm_phi.domain() << std::endl;
    std::cout << "Move Runs: " << runperm_phi.move_runs() << std::endl;
    std::cout << "Permutation Runs: " << runperm_phi.permutation_runs() << std::endl;

    using Position = typename RunPermPhi<RunData>::Position;
    Position pos = runperm_phi.first();
    for (size_t i = 0; i < runperm_phi.domain(); ++i) {
        pos = runperm_phi.next(pos);
    }
    assert(pos.interval == 0);
    assert(pos.offset == 0);
    std::cout << "Successfully moved to start" << std::endl;

    pos = runperm_phi.last();
    pos = runperm_phi.Phi(pos);
    std::vector<ulint> sa_recovered(sa.size());
    for (size_t i = 0; i < sa_recovered.size(); ++i) {
        sa_recovered[sa.size() - i - 1] = runperm_phi.SA(pos);
        pos = runperm_phi.next(pos);
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
    std::cout << "Domain: " << move_phi.domain() << std::endl;
    std::cout << "Move Runs: " << move_phi.move_runs() << std::endl;
    std::cout << "Permutation Runs: " << move_phi.permutation_runs() << std::endl;
    using Position = typename MovePhi::Position;
    Position pos = move_phi.first();
    for (size_t i = 0; i < move_phi.domain(); ++i) {
        pos = move_phi.next(pos);
    }
    assert(pos.interval == 0);
    assert(pos.offset == 0);
    std::cout << "Successfully moved to start" << std::endl;

    pos = move_phi.last();
    pos = move_phi.Phi(pos);
    std::vector<ulint> sa_recovered(sa.size());
    for (size_t i = 0; i < move_phi.domain(); ++i) {
        sa_recovered[sa.size() - i - 1] = move_phi.SA(pos);
        pos = move_phi.Phi(pos);
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
    std::cout << "Domain: " << runperm_invphi.domain() << std::endl;
    std::cout << "Move Runs: " << runperm_invphi.move_runs() << std::endl;
    std::cout << "Permutation Runs: " << runperm_invphi.permutation_runs() << std::endl;

    runperm_invphi.first();
    using Position = typename RunPermPhi<RunData>::Position;
    Position pos = runperm_invphi.first();
    for (size_t i = 0; i < runperm_invphi.domain(); ++i) {
        pos = runperm_invphi.next(pos);
    }
    assert(pos.interval == 0);
    assert(pos.offset == 0);
    std::cout << "Successfully moved to start" << std::endl;

    pos = runperm_invphi.last();
    std::vector<ulint> sa_recovered(sa.size());
    for (size_t i = 0; i < sa_recovered.size(); ++i) {
        sa_recovered[i] = runperm_invphi.SA(pos);
        pos = runperm_invphi.InvPhi(pos);
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
    std::cout << "Domain: " << move_invphi.domain() << std::endl;
    std::cout << "Move Runs: " << move_invphi.move_runs() << std::endl;
    std::cout << "Permutation Runs: " << move_invphi.permutation_runs() << std::endl;
    using Position = typename MoveInvPhi::Position;
    Position pos = move_invphi.first();
    for (size_t i = 0; i < move_invphi.domain(); ++i) {
        pos = move_invphi.next(pos);
    }
    assert(pos.interval == 0);
    assert(pos.offset == 0);
    std::cout << "Successfully moved to start" << std::endl;

    pos = move_invphi.last();
    std::vector<ulint> sa_recovered(sa.size());
    for (size_t i = 0; i < move_invphi.domain(); ++i) {
        sa_recovered[i] = move_invphi.SA(pos);
        pos = move_invphi.InvPhi(pos);
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