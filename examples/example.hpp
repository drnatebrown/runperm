#ifndef EXAMPLES_HPP
#define EXAMPLES_HPP

#include "orbit/permutation.hpp"
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
        if constexpr (std::is_same_v<run_perm_type, orbit::permutation_absolute<run_cols>>) {
            std::cout << "Absolute Position: " << pos.idx;
        }
        if constexpr (std::is_same_v<run_perm_type, orbit::rlbwt::lf_permutation<run_cols>>) {
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
        if constexpr (std::is_same_v<run_perm_type, orbit::permutation_absolute<run_cols>>) {
            std::cout << "Absolute Position: " << pos.idx;
        }
        if constexpr (std::is_same_v<run_perm_type, orbit::rlbwt::lf_permutation<run_cols>>) {
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
        if constexpr (std::is_same_v<move_perm_type, orbit::permutation_absolute<>>) {
            std::cout << ", Absolute Position: " << pos.idx;
        }
        if constexpr (std::is_same_v<move_perm_type, orbit::rlbwt::lf_permutation<>>) {
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
        if constexpr (std::is_same_v<move_perm_type, orbit::permutation_absolute<>>) {
            std::cout << ", Absolute Position: " << pos.idx;
        }
        if constexpr (std::is_same_v<move_perm_type, orbit::rlbwt::lf_permutation<>>) {
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

ulint calculate_capped_length(ulint domain, ulint lengths, ulint length_capping) {
    return static_cast<ulint>(std::ceil(static_cast<double>(domain) / static_cast<double>(lengths) * length_capping));
}

ulint calculate_max_val(ulint domain, ulint lengths, ulint length_capping) {
    return orbit::max_val(orbit::bit_width(calculate_capped_length(domain, lengths, length_capping)));
}

#endif // EXAMPLES_HPP