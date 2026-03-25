#ifndef _RUNPERM_PHI_HPP
#define _RUNPERM_PHI_HPP

#include "orbit/common.hpp"
#include "orbit/internal/runperm/runperm_impl.hpp"
#include "orbit/internal/rlbwt/phi_helpers.hpp"

namespace orbit::rlbwt {

// Always use absolute positions for Phi
template<typename data_columns_t,
         bool integrated_move_structure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool exponential_search = DEFAULT_EXPONENTIAL_SEARCH,
         template<typename> class table_t = move_vector>
class runperm_phi : public runperm_impl<data_columns_t, integrated_move_structure, true, exponential_search, move_columns, move_structure, table_t> {
    using base = runperm_impl<data_columns_t, integrated_move_structure, true, exponential_search, move_columns, move_structure, table_t>;
public:
    using base::base;
    using base::operator=;
    using position = typename base::position;

    position phi(position pos) {
        return base::next(pos);
    }

    position phi(position pos, ulint steps) {
        return base::next(pos, steps);
    }

    ulint SA(position pos) {
        return pos.idx;
    }
};

// Always use absolute positions for Phi
template<bool exponential_search = DEFAULT_EXPONENTIAL_SEARCH>
class move_phi : public moveperm_impl<true, exponential_search, move_columns, move_structure, move_vector> {
    using base = moveperm_impl<true, exponential_search, move_columns, move_structure, move_vector>;
public:
    using base::base;
    using base::operator=;
    using position = typename base::position;

    move_phi(const std::vector<uchar>& rlbwt_heads,
        const std::vector<ulint>& rlbwt_run_lengths,
        const split_params& split_params = split_params())
    : base([&] {
        size_t domain = 0;
        ulint max_length = 0;
        auto [phi_lengths, phi_img_rank_inv] =
            rlbwt_to_phi_img_rank_inv<>(rlbwt_heads, rlbwt_run_lengths,
                                        &domain, &max_length);
        return permutation_impl<>::from_lengths_and_img_rank_inv(
            phi_lengths, phi_img_rank_inv, domain, max_length, split_params);
      }()) {}

    position phi(position pos) {
        return base::next(pos);
    }

    position phi(position pos, ulint steps) {
        return base::next(pos, steps);
    }

    ulint SA(position pos) {
        return pos.idx;
    }
};

} // namespace orbit::rlbwt

#endif /* end of include guard: _RUNPERM_PHI_HPP */