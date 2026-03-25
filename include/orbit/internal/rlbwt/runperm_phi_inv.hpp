#ifndef _RUNPERM_phi_inv_HPP
#define _RUNPERM_phi_inv_HPP

#include "orbit/common.hpp"
#include "orbit/internal/runperm/runperm_impl.hpp"
#include "orbit/internal/rlbwt/runperm_lf.hpp"

namespace orbit::rlbwt {

// Always use absolute positions for inv_phi
template<typename data_columns_t,
         bool integrated_move_structure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool exponential_search = DEFAULT_EXPONENTIAL_SEARCH,
         template<typename> class table_t = move_vector>
class runperm_phi_inv : public runperm_impl<data_columns_t, integrated_move_structure, true, exponential_search, move_columns, move_structure, table_t> {
    using base = runperm_impl<data_columns_t, integrated_move_structure, true, exponential_search, move_columns, move_structure, table_t>;
public:
    using base::base;
    using base::operator=;
    using position = typename base::position;

    position phi_inv(position pos) {
        return base::next(pos);
    }

    position phi_inv(position pos, ulint steps) {
        return base::next(pos, steps);
    }

    ulint SA(position pos) {
        return pos.idx;
    }
};

// Always use absolute positions for inv_phi
template<bool exponential_search = DEFAULT_EXPONENTIAL_SEARCH>
class move_phi_inv : public moveperm_impl<true, exponential_search, move_columns, move_structure, move_vector> {
    using base = moveperm_impl<true, exponential_search, move_columns, move_structure, move_vector>;
public:
    using base::base;
    using base::operator=;
    using position = typename base::position;


    move_phi_inv(const std::vector<uchar>& rlbwt_heads,
        const std::vector<ulint>& rlbwt_run_lengths,
        const split_params& sp = split_params())
    : base([&] {
        size_t domain = 0;
        ulint max_length = 0;
        auto [phi_inv_lengths, phi_inv_img_rank_inv] =
            rlbwt_to_phi_inv_img_rank_inv<>(rlbwt_heads, rlbwt_run_lengths,
                                            &domain, &max_length);
        return permutation_impl<>::from_lengths_and_img_rank_inv(
            phi_inv_lengths, phi_inv_img_rank_inv, domain, max_length, sp);
    }()) {}

    position phi_inv(position pos) {
        return base::next(pos);
    }

    position phi_inv(position pos, ulint steps) {
        return base::next(pos, steps);
    }

    ulint SA(position pos) {
        return pos.idx;
    }
};

} // namespace orbit::rlbwt

#endif /* end of include guard: _RUNPERM_phi_inv_HPP */