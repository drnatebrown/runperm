#include "orbit/runperm.hpp"
#include "orbit/rlbwt.hpp"

typedef unsigned long int ulint;
typedef unsigned char uchar;

/***** HELPER FUNCTIONS *****/
template<typename run_perm_type>
void print_run_data(run_perm_type& rp) {
    using run_cols = typename run_perm_type::data_columns;
    std::cout << "Intervals:" << std::endl;
    auto pos = rp.first();
    for (ulint i = 0; i < rp.intervals(); ++i) {
        std::cout << "Interval: " << pos.interval << ", Length: " << rp.get_length(pos.interval) << " --> ";
        std::cout << "Run Data: " << rp.template get<run_cols::VAL1>(pos) << ", " << rp.template get<run_cols::VAL2>(pos);
        if constexpr (std::is_same_v<run_perm_type, orbit::runperm<run_cols, true, true>>) {
            std::cout << "Absolute Position: " << pos.idx;
        }
        if constexpr (std::is_same_v<run_perm_type, orbit::rlbwt::runperm_lf<run_cols>> || std::is_same_v<run_perm_type, orbit::rlbwt::runperm_fl<run_cols>>) {
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
        std::cout << "Run Data: " << rp.template get<run_cols::VAL1>(pos) << ", " << rp.template get<run_cols::VAL2>(pos);
        if constexpr (std::is_same_v<run_perm_type, orbit::runperm<run_cols, true, true>>) {
            std::cout << "Absolute Position: " << pos.idx;
        }
        if constexpr (std::is_same_v<run_perm_type, orbit::rlbwt::runperm_lf<run_cols>> || std::is_same_v<run_perm_type, orbit::rlbwt::runperm_fl<run_cols>>) {
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

template<typename move_perm_type>
void print_move_permutation(move_perm_type& mp) {
    std::cout << "Intervals:" << std::endl;
    auto pos = mp.first();
    for (ulint i = 0; i < mp.intervals(); ++i) {
        std::cout << "Interval: " << pos.interval << ", Length: " << mp.get_length(pos.interval);
        if constexpr (std::is_same_v<move_perm_type, orbit::moveperm_absolute>) {
            std::cout << ", Absolute Position: " << pos.idx;
        }
        if constexpr (std::is_same_v<move_perm_type, orbit::rlbwt::move_lf<>> || std::is_same_v<move_perm_type, orbit::rlbwt::move_fl<>>) {
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
        if constexpr (std::is_same_v<move_perm_type, orbit::moveperm_absolute>) {
            std::cout << ", Absolute Position: " << pos.idx;
        }
        if constexpr (std::is_same_v<move_perm_type, orbit::rlbwt::move_lf<>> || std::is_same_v<move_perm_type, orbit::rlbwt::move_fl<>>) {
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
    "Permutation",
    "Move Permutation",
    "Move Permutation with splitting",
    "LF from RLBWT",
    "Text from RLBWT using FL",
    "SA from RLBWT using Inverse Phi"
};

// Create a run-length permutation
void example1() {
    std::cout << "Example 1: " << example_names[0] << std::endl;

    std::vector<ulint> lengths =     { 2, 3, 1,  2, 2,  1, 1,  1, 3};        // Length of runs
    std::vector<ulint> permutation = { 1, 9, 3, 12, 4, 14, 0, 15, 6};       // Permutation of runs

    // Some example data columns to store alongside these runs.
    DEFINE_COLUMNS(run_cols, VAL1, VAL2);
    // The DEFINE_COLUMNS(enum_name, ...) macro above is equivalent to:
    // enum class run_cols {
    //     VAL1,
    //     VAL2,
    //     COUNT
    // };
    // !!! The COUNT enumerator is automatically added by the DEFINE_COLUMNS macro,
    // !!! but must be included manually in the enum definition.

    // Defines std::array<ulint, static_cast<size_t>(run_cols::COUNT)>
    using run_cols_tuple = orbit::columns_tuple<run_cols>;
    std::vector<run_cols_tuple> run_data(lengths.size()); // Some data with tuples per run
    // Fill with dummy data
    for (size_t i = 0; i < lengths.size(); ++i) {
        run_data[i] = {i, i + 1};
    }

    // Basic construction
    orbit::runperm<run_cols> rp(lengths, permutation, run_data);
    print_run_data(rp);
}

// Create a move permutation
void example2() {
    std::cout << "Example 2: " << example_names[1] << std::endl;
    // Can create from an exact permutation
    std::vector<ulint> permutation = {1, 2, 9, 10, 11, 3, 12, 13, 4, 5, 14, 0, 15, 6, 7, 8};
    // Same as MovePerm<> --> Know interval/offset at any point during navigation
    std::cout << "Move Permutation (Relative):" << std::endl;
    orbit::moveperm_relative mp_relative(permutation);
    print_move_permutation(mp_relative);

    std::cout << "\n\nMove Permutation (Absolute):" << std::endl;
    // Or, can create from lengths and interval permutation
    auto [lengths, images] = orbit::get_permutation_intervals(permutation);
    // moveperm_absolute also stores the absolute position in the permutation
    orbit::moveperm_absolute mp_absolute(lengths, images);
    print_move_permutation(mp_absolute);
}

// Create a move permutation, with splitting
void example3() {
    std::cout << "Example 3: " << example_names[2] << std::endl;
    std::vector<ulint> lengths = {2, 1, 8};
    std::vector<ulint> images = {9, 0, 1};

    // Original move permutation
    std::cout << "Original Move Permutation (Relative):" << std::endl;
    orbit::moveperm_relative mp_relative(lengths, images);
    print_move_permutation(mp_relative);

    std::cout << "\n\nSplitting..." << std::endl;
    std::cout << "Runs, r: " << lengths.size() << std::endl;
    std::cout << "Domain, n: " << mp_relative.domain() << std::endl;
    std::cout << "n/r: " << static_cast<double>(mp_relative.domain()) / static_cast<double>(lengths.size()) << std::endl;

    orbit::split_params sp;
    sp.length_capping = 1;

    std::cout << "Length Capping Factor: " << sp.length_capping.value() << std::endl;
    std::cout << "Capped Length (n/r * length_capping_factor): " << static_cast<ulint>(std::ceil(static_cast<double>(mp_relative.domain()) / static_cast<double>(lengths.size()) * sp.length_capping.value())) << std::endl;
    std::cout << "Round up to use all log2(n/r * length_capping_factor) bits: " << orbit::max_val(orbit::bit_width(static_cast<ulint>(std::ceil(static_cast<double>(mp_relative.domain()) / static_cast<double>(lengths.size()) * sp.length_capping.value())))) << std::endl;

    std::cout << "\nMove Permutation (Relative):" << std::endl;
    orbit::moveperm_relative mp_relative_split(lengths, images, sp);
    print_move_permutation(mp_relative_split);

    std::cout << "\n\nMove Permutation (Absolute):" << std::endl;
    orbit::moveperm_absolute mp_absolute_split(lengths, images, sp);
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
    orbit::rlbwt::move_lf<> move_lf(bwt_heads, bwt_run_lengths);
    print_move_permutation(move_lf);

    // LF with user data columns
    std::cout << "\n\nLF with with user data columns columns" << std::endl;
    // Alternative to example 1, without using the macro
    enum class run_cols {
        VAL1,
        VAL2,
        COUNT
    };
    using run_cols_tuple = orbit::columns_tuple<run_cols>;
    std::vector<run_cols_tuple> lf_run_data(bwt_heads.size()); // Insert with some type of data
    for (size_t i = 0; i < bwt_heads.size(); ++i) {
        lf_run_data[i] = {i, i + 1};
    }
    orbit::rlbwt::runperm_lf<run_cols> runperm_lf(bwt_heads, bwt_run_lengths, lf_run_data);
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
    orbit::rlbwt::move_fl<> move_fl(bwt_heads, bwt_run_lengths);

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
    std::cout << "Example 6: " << example_names[5] << std::endl;
    std::vector<uchar> bwt_heads =       {'T','C','G','A','T', 0 ,'A','T','A'};
    std::vector<ulint> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };
    size_t domain = std::accumulate(bwt_run_lengths.begin(), bwt_run_lengths.end(), 0);

    std::cout << "Input RLBWT: (T, 5), (C, 3), (G, 3), (A, 3), (T, 1), ($, 1), (A, 1), (T, 4), (A, 6)" << std::endl;

    auto permutation = orbit::rlbwt::rlbwt_to_phi(bwt_heads, bwt_run_lengths);
    orbit::rlbwt::move_phi_inv<> move_phi_inv(permutation);
    auto pos = move_phi_inv.last();
    std::vector<ulint> sa_recovered(domain);
    for (size_t i = 0; i < domain; ++i) {
        sa_recovered[i] = move_phi_inv.SA(pos);
        pos = move_phi_inv.phi_inv(pos);
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