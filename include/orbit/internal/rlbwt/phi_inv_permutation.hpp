#ifndef _PERMUTATION_PHI_INV_HPP
#define _PERMUTATION_PHI_INV_HPP

#include "orbit/common.hpp"
#include "orbit/internal/perm/permutation_impl.hpp"
#include "orbit/internal/rlbwt/phi_helpers.hpp"

namespace orbit::rlbwt {

// Always use absolute positions for phi_inv
template<typename data_columns_t = empty_data_columns,
         bool integrated_move_structure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool exponential_search = DEFAULT_EXPONENTIAL_SEARCH,
         template<typename> class table_t = move_vector>
class phi_inv_permutation_impl : public permutation_impl<data_columns_t, integrated_move_structure, true, exponential_search, move_columns, move_structure, table_t> {
    using base = permutation_impl<data_columns_t, integrated_move_structure, true, exponential_search, move_columns, move_structure, table_t>;
public:
    using base::base;
    using base::operator=;
    using position = typename base::position;

    template<class dc = data_columns_t,
             std::enable_if_t<std::is_same_v<dc, empty_data_columns>, int> = 0>
    phi_inv_permutation_impl(const std::vector<uchar>& rlbwt_heads,
        const std::vector<ulint>& rlbwt_run_lengths,
        const split_params& sp = split_params())
    : base([&] {
        size_t domain = 0;
        ulint max_length = 0;
        auto [phi_inv_lengths, phi_inv_img_rank_inv] =
            rlbwt_to_phi_inv_img_rank_inv<>(rlbwt_heads, rlbwt_run_lengths,
                                            &domain, &max_length);
        return interval_encoding_impl<>::from_lengths_and_img_rank_inv(
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

// Always use absolute positions for phi_inv
template<bool exponential_search = DEFAULT_EXPONENTIAL_SEARCH,
         template<typename> class table_t = move_vector>
using phi_inv_move_impl = phi_inv_permutation_impl<empty_data_columns, false, exponential_search, table_t>;

} // namespace orbit::rlbwt

#endif /* end of include guard: _PERMUTATION_PHI_INV_HPP */