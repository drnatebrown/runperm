#include "runperm.hpp"
#include "rlbwt.hpp"

typedef unsigned long int ulint;

/***** HELPER FUNCTIONS *****/
template<typename RunPermType>
void print_run_data(RunPermType& rp) {
    using RunCols = typename RunPermType::RunCols;
    std::cout << "Intervals:" << std::endl;
    auto pos = rp.first();
    for (ulint i = 0; i < rp.move_runs(); ++i) {
        std::cout << "Interval: " << pos.interval << ", Length: " << rp.get_length(pos.interval) << " --> ";
        std::cout << "Run Data: " << rp.template get<RunCols::VAL1>(pos) << ", " << rp.template get<RunCols::VAL2>(pos);
        if constexpr (std::is_same_v<RunPermType, RunPerm<RunCols, true, true>>) {
            std::cout << "Absolute Position: " << pos.idx;
        }
        if constexpr (std::is_same_v<RunPermType, RunPermLF<RunCols>> || std::is_same_v<RunPermType, RunPermFL<RunCols>>) {
            if (rp.get_character(pos.interval) == 0) {
                std::cout << ", Character: $";
            } else {
                std::cout << ", Character: " << rp.get_character(pos.interval);
            }
        }
        std::cout << std::endl;
        pos = rp.down(pos);
    }
    assert(pos == rp.first());

    std::cout << "\nPermutation Order:" << std::endl;
    // Permutation Order:
    pos = rp.first();
    for (size_t i = 0; i < rp.domain(); ++i) {
        std::cout << "Position: " << pos.interval << ", " << pos.offset << " --> ";
        std::cout << "Run Data: " << rp.template get<RunCols::VAL1>(pos) << ", " << rp.template get<RunCols::VAL2>(pos);
        if constexpr (std::is_same_v<RunPermType, RunPerm<RunCols, true, true>>) {
            std::cout << "Absolute Position: " << pos.idx;
        }
        if constexpr (std::is_same_v<RunPermType, RunPermLF<RunCols>> || std::is_same_v<RunPermType, RunPermFL<RunCols>>) {
            if (rp.get_character(pos.interval) == 0) {
                std::cout << ", Character: $";
            } else {
                std::cout << ", Character: " << rp.get_character(pos.interval);
            }
        }
        std::cout << std::endl;
        pos = rp.next(pos);
    }
    assert(pos == rp.first());
}

template<typename MovePermType>
void print_move_permutation(MovePermType& mp) {
    std::cout << "Intervals:" << std::endl;
    auto pos = mp.first();
    for (ulint i = 0; i < mp.move_runs(); ++i) {
        std::cout << "Interval: " << pos.interval << ", Length: " << mp.get_length(pos.interval);
        if constexpr (std::is_same_v<MovePermType, MovePermAbsolute>) {
            std::cout << ", Absolute Position: " << pos.idx;
        }
        if constexpr (std::is_same_v<MovePermType, MoveLF<>> || std::is_same_v<MovePermType, MoveFL<>>) {
            if (mp.get_character(pos.interval) == 0) {
                std::cout << ", Character: $";
            } else {
                std::cout << ", Character: " << mp.get_character(pos.interval);
            }
        }
        std::cout << std::endl;
        pos = mp.down(pos);
    }
    assert(pos == mp.first());

    std::cout << "\nPermutation Order:" << std::endl;
    pos = mp.first();
    for (size_t i = 0; i < mp.domain(); ++i) {
        std::cout << "Position: " << pos.interval << ", " << pos.offset;
        if constexpr (std::is_same_v<MovePermType, MovePermAbsolute>) {
            std::cout << ", Absolute Position: " << pos.idx;
        }
        if constexpr (std::is_same_v<MovePermType, MoveLF<>> || std::is_same_v<MovePermType, MoveFL<>>) {
            if (mp.get_character(pos.interval) == 0) {
                std::cout << ", Character: $";
            } else {
                std::cout << ", Character: " << mp.get_character(pos.interval);
            }
        }
        std::cout << std::endl;
        pos = mp.next(pos);
    }
    assert(pos == mp.first());
}

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

/************** EXAMPLES **************/
std::vector<std::string> example_names = {
    "RunPerm",
    "MovePerm",
    "MovePerm with splitting",
    "LF from RLBWT",
    "Text from RLBWT using FL",
    "SA from RLBWT using Inverse Phi"
};

// Create a run-length permutation
void example1() {
    std::cout << "Example 1: " << example_names[0] << std::endl;

    std::vector<ulint> lengths =     { 2, 3, 1,  2, 2,  1, 1,  1, 3};        // Length of runs
    std::vector<ulint> permutation = { 1, 9, 3, 12, 4, 14, 0, 15, 6};       // Permutation of runs
    ulint domain = 16; // Total domain size

    // Some example data columns to store alongside these runs.
    // Must always include COUNT as the last entry to signal number of fields.
    enum class RunCols {
        VAL1,
        VAL2,
        COUNT
    };

    using RunColsTuple = DataTuple<RunCols>;
    std::vector<RunColsTuple> run_data(lengths.size()); // Some data with tuples per run
    // Fill with dummy data
    for (size_t i = 0; i < lengths.size(); ++i) {
        run_data[i] = {i, i + 1};
    }

    // Basic construction
    RunPerm<RunCols> rp(lengths, permutation, domain, run_data);
    print_run_data(rp);
}

// Create a move permutation
void example2() {
    std::cout << "Example 2: " << example_names[1] << std::endl;
    // Can create from an exact permutation
    std::vector<ulint> permutation = {1, 2, 9, 10, 11, 3, 12, 13, 4, 5, 14, 0, 15, 6, 7, 8};
    // Same as MovePerm<> --> Know interval/offset at any point during navigation
    std::cout << "Move Permutation (Relative):" << std::endl;
    MovePermRelative mp_relative(permutation);
    print_move_permutation(mp_relative);

    std::cout << "\n\nMove Permutation (Absolute):" << std::endl;
    // Or, can create from lengths and interval permutation
    auto [lengths, interval_permutation] = get_permutation_intervals(permutation);
    ulint domain = permutation.size();
    // MovePermAbsolute also stores the absolute position in the permutation
    MovePermAbsolute mp_absolute(lengths, interval_permutation, domain);
    print_move_permutation(mp_absolute);
}

// Create a move permutation, with splitting
void example3() {
    std::cout << "Example 3: " << example_names[2] << std::endl;
    std::vector<ulint> lengths = {2, 1, 8};
    std::vector<ulint> interval_permutation = {9, 0, 1};
    ulint domain = 11;

    // Original move permutation
    std::cout << "Original Move Permutation (Relative):" << std::endl;
    MovePermRelative mp_relative(lengths, interval_permutation, domain);
    print_move_permutation(mp_relative);

    std::cout << "\n\nSplitting..." << std::endl;
    std::cout << "Runs, r: " << lengths.size() << std::endl;
    std::cout << "Domain, n: " << domain << std::endl;
    std::cout << "n/r: " << static_cast<double>(domain) / static_cast<double>(lengths.size()) << std::endl;

    SplitParams split_params;
    split_params.length_capping_factor = 1;

    std::cout << "Length Capping Factor: " << split_params.length_capping_factor.value() << std::endl;
    std::cout << "Capped Length (n/r * length_capping_factor): " << static_cast<ulint>(std::ceil(static_cast<double>(domain) / static_cast<double>(lengths.size()) * split_params.length_capping_factor.value())) << std::endl;

    std::cout << "\nMove Permutation (Relative):" << std::endl;
    MovePermRelative mp_relative_split(lengths, interval_permutation, domain, split_params);
    print_move_permutation(mp_relative_split);

    std::cout << "\n\nMove Permutation (Absolute):" << std::endl;
    MovePermAbsolute mp_absolute_split(lengths, interval_permutation, domain, split_params);
    print_move_permutation(mp_absolute_split);
}

// Create objects for LF from RLBWT
void example4() {
    std::cout << "Example 4: " << example_names[3] << std::endl;
    // TEXT: GATTACATGATTACATAGATTACATT$
    // BWT:  TTTTTCCCGGGAAAT$ATTTTAAAAAA
    // RLBWT as run heads (characters) and run lengths, with 0 as terminator
    std::vector<uchar> bwt_heads =       {'T','C','G','A','T', 0 ,'A','T','A'};
    std::vector<ulint> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };

    // LF using Move-only interface
    std::cout << "LF using Move-only interface" << std::endl;
    MoveLF<> move_lf(bwt_heads, bwt_run_lengths);
    print_move_permutation(move_lf);

    // LF with RunPerm + run data columns
    std::cout << "\n\nLF with RunPerm + run data columns" << std::endl;
    enum class RunCols { VAL1, VAL2, COUNT };
    using LFRunData = DataTuple<RunCols>;
    std::vector<LFRunData> lf_run_data(bwt_heads.size()); // Insert with some type of data
    for (size_t i = 0; i < bwt_heads.size(); ++i) {
        lf_run_data[i] = {i, i + 1};
    }
    RunPermLF<RunCols> runperm_lf(bwt_heads, bwt_run_lengths, lf_run_data);
    print_run_data(runperm_lf);
}

// Extract text from RLBWT using FL
void example5() {
    std::cout << "Example 5: " << example_names[4] << std::endl;
    // TEXT: GATTACATGATTACATAGATTACATT$
    // BWT:  TTTTTCCCGGGAAAT$ATTTTAAAAAA
    // RLBWT as run heads (characters) and run lengths, with 0 as terminator
    std::vector<uchar> bwt_heads =       {'T','C','G','A','T', 0 ,'A','T','A'};
    std::vector<ulint> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };
    MoveFL<> move_fl(bwt_heads, bwt_run_lengths);

    std::cout << "Input RLBWT: (T, 5), (C, 3), (G, 3), (A, 3), (T, 1), ($, 1), (A, 1), (T, 4), (A, 6)" << std::endl;
    
    auto pos = move_fl.first();
    pos = move_fl.FL(pos);
    std::string recovered_text = "";
    for (size_t i = 0; i < move_fl.domain(); ++i) {
        recovered_text += (move_fl.get_character(pos) == 0) ? '$' : (char)move_fl.get_character(pos);
        pos = move_fl.FL(pos);
    }
    std::cout << "Recovered text: " << recovered_text << std::endl;
}

// Extract SA from RLBWT using Inverse Phi
void example6() {
    std::vector<uchar> bwt_heads =       {'T','C','G','A','T', 0 ,'A','T','A'};
    std::vector<ulint> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };

    std::cout << "Input RLBWT: (T, 5), (C, 3), (G, 3), (A, 3), (T, 1), ($, 1), (A, 1), (T, 4), (A, 6)" << std::endl;

    auto [phi_lengths, phi_interval_permutations, domain] = rlbwt_to_phi(bwt_heads, bwt_run_lengths);
    MoveInvPhi move_invphi(phi_lengths, phi_interval_permutations, domain);
    auto pos = move_invphi.last();
    std::vector<ulint> sa_recovered(domain);
    for (size_t i = 0; i < domain; ++i) {
        sa_recovered[i] = move_invphi.SA(pos);
        pos = move_invphi.InvPhi(pos);
    }
    std::cout << "Recovered SA: ";
    print_vector(sa_recovered);
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <example_name>" << std::endl;
        std::cerr << "Examples: " << std::endl;
        size_t i = 1;
        for (const auto& name : example_names) {
            std::cerr << i << ": " << name << std::endl;
            ++i;
        }
        return 1;
    }
    size_t example_number = std::stoi(argv[1]);
    switch (example_number) {
        case 1:
            example1();
            break;
        case 2:
            example2();
            break;
        case 3:
            example3();
            break;
        case 4:
            example4();
            break;
        case 5:
            example5();
            break;
        case 6:
            example6();
            break;
        default:
            std::cerr << "Invalid example number" << std::endl;
            return 1;
    }
    return 0;
}