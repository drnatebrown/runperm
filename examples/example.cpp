#include "example.hpp"
#include "orbit/permutation.hpp"
#include "orbit/rlbwt.hpp"

typedef unsigned long int ulint;
typedef unsigned char uchar;

/************** EXAMPLES **************/
std::vector<std::string> example_names = {
    "Permutation",
    "Permutation with Data Columns",
    "Permutation with Splitting",
    "LF from RLBWT",
    "Text from RLBWT using FL",
    "SA from RLBWT using Inverse Phi"
};

// Create a move permutation
void example1() {
    std::cout << "Example 1: " << example_names[1] << std::endl;
    // Can create from an exact permutation
    std::vector<ulint> permutation = {1, 2, 9, 10, 11, 3, 12, 13, 4, 5, 14, 0, 15, 6, 7, 8};

    // Relative positions mean we know the interval/offset at any point during navigation
    std::cout << "Move Permutation (Relative):" << std::endl;
    orbit::permutation<> mp_relative(permutation);
    print_move_permutation(mp_relative);

    // Absolute positions mean we know the absolute position in the permutation domain at any point during navigation
    std::cout << "\n\nMove Permutation (Absolute):" << std::endl;
    // Or, can create from lengths and interval permutation
    std::vector<ulint> lengths =     { 2, 3, 1,  2, 2,  1, 1,  1, 3};  // Length of contigiously permuted runs
    std::vector<ulint> images  =    { 1, 9, 3, 12, 4, 14, 0, 15, 6};   // Permutation of first values of contigiously permuted runs
    orbit::permutation_absolute<> mp_absolute(lengths, images);
    print_move_permutation(mp_absolute);
}

// Create a run-length permutation
void example2() {
    std::cout << "Example 2: " << example_names[1] << std::endl;

    std::vector<ulint> lengths =     { 2, 3, 1,  2, 2,  1, 1,  1, 3};  // Length of contigiously permuted runs
    std::vector<ulint> images  =    { 1, 9, 3, 12, 4, 14, 0, 15, 6};   // Permutation of first values of contigiously permuted runs

    // Some example data columns to store alongside these runs.
    DEFINE_ORBIT_COLUMNS(run_cols, VAL1, VAL2);
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
    orbit::permutation<run_cols> rp(lengths, images, run_data);
    print_run_data(rp);
}

// Create a permutation, with splitting
void example3() {
    std::cout << "Example 3: " << example_names[2] << std::endl;
    std::vector<ulint> lengths = {2, 1, 8};
    std::vector<ulint> images =  {9, 0, 1};

    // Original move permutation
    std::cout << "Original Permutation (Relative):" << std::endl;
    orbit::permutation<> mp_relative(lengths, images);
    print_move_permutation(mp_relative);

    std::cout << "\n\nSplitting..." << std::endl;
    std::cout << "Runs, r: " << lengths.size() << std::endl;
    std::cout << "Domain, n: " << mp_relative.domain() << std::endl;
    std::cout << "n/r: " << static_cast<double>(mp_relative.domain()) / static_cast<double>(lengths.size()) << std::endl;

    orbit::split_params sp;
    sp.length_capping = 1;

    std::cout << "Length Capping Factor: " << sp.length_capping.value() << std::endl;
    std::cout << "Capped Length (n/r * length_capping_factor): " << calculate_capped_length(mp_relative.domain(), lengths.size(), sp.length_capping.value()) << std::endl;
    std::cout << "Round up to use all log2(n/r * length_capping_factor) bits: " << calculate_max_val(mp_relative.domain(), lengths.size(), sp.length_capping.value()) << std::endl;

    std::cout << "\nSplit Permutation (Relative):" << std::endl;
    orbit::permutation<> mp_relative_split(lengths, images, sp);
    print_move_permutation(mp_relative_split);

    std::cout << "\n\nSplit Permutation (Absolute):" << std::endl;
    orbit::permutation_absolute<> mp_absolute_split(lengths, images, sp);
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
    std::cout << "LF using Orbit Permutation" << std::endl;
    orbit::rlbwt::lf_permutation<> move_lf(bwt_heads, bwt_run_lengths);
    print_move_permutation(move_lf);

    // LF with user data columns
    std::cout << "\n\nLF with user data columns" << std::endl;
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
    orbit::rlbwt::lf_permutation<run_cols> lf(bwt_heads, bwt_run_lengths, lf_run_data);
    print_run_data(lf);
}

// Extract text from RLBWT using FL
void example5() {
    std::cout << "Example 5: " << example_names[4] << std::endl;
    // TEXT: GATTACATGATTACATAGATTACATT$
    // BWT:  TTTTTCCCGGGAAAT$ATTTTAAAAAA
    // RLBWT as run heads (characters) and run lengths, with 0 as terminator
    std::vector<uchar> bwt_heads =       {'T','C','G','A','T', 0 ,'A','T','A'};
    std::vector<ulint> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };
    orbit::rlbwt::fl_permutation<> move_fl(bwt_heads, bwt_run_lengths);

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

    std::cout << "Input RLBWT: (T, 5), (C, 3), (G, 3), (A, 3), (T, 1), ($, 1), (A, 1), (T, 4), (A, 6)" << std::endl;

    orbit::rlbwt::phi_inv_permutation<> move_phi_inv(bwt_heads, bwt_run_lengths);
    auto pos = move_phi_inv.last();
    std::vector<ulint> sa_recovered(move_phi_inv.domain());
    for (size_t i = 0; i < move_phi_inv.domain(); ++i) {
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