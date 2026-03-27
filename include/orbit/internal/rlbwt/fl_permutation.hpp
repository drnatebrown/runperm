#ifndef _RLBWT_FL_HPP
#define _RLBWT_FL_HPP

#include "orbit/common.hpp"
#include "orbit/internal/ds/alphabet.hpp"
#include "orbit/internal/rlbwt/specializations/rlbwt_permutation.hpp"
#include "orbit/internal/rlbwt/specializations/rlbwt_interval_encoding.hpp"
#include "orbit/internal/rlbwt/rlbwt_helpers.hpp"

namespace orbit::rlbwt {

template<typename data_columns_t = empty_data_columns,
         bool integrated_move_structure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         bool exponential_search = DEFAULT_EXPONENTIAL_SEARCH,
         typename alphabet_t = nucleotide,
         template<typename> class table_t = move_vector>
class fl_permutation_impl : public rlbwt_permutation<fl_permutation_impl<data_columns_t, integrated_move_structure, store_absolute_positions, exponential_search, alphabet_t, table_t>,
                         data_columns_t, integrated_move_structure, store_absolute_positions, exponential_search, alphabet_t, table_t> {
    using base = rlbwt_permutation<fl_permutation_impl<data_columns_t, integrated_move_structure, store_absolute_positions, exponential_search, alphabet_t, table_t>,
                         data_columns_t, integrated_move_structure, store_absolute_positions, exponential_search, alphabet_t, table_t>;
    using rlbwt_interval_encoding = typename base::rlbwt_interval_encoding;
public:
    using data_columns = typename base::data_columns;
    using base::base;
    using base::operator=;
    using position = typename base::position;

    rlbwt_interval_encoding find_interval_encoding(
        const std::vector<uchar>& rlbwt_heads,
        const std::vector<ulint>& rlbwt_run_lengths,
        const split_params& sp
    ) {
        return rlbwt_interval_encoding::fl_interval_encoding(rlbwt_heads, rlbwt_run_lengths, sp);
    }

    position FL(position pos) {
        return base::next(pos);
    }

    position FL(position pos, ulint steps) {
        return base::next(pos, steps);
    }
};

template<bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         bool exponential_search = DEFAULT_EXPONENTIAL_SEARCH,
         typename alphabet_t = nucleotide,
         template<typename> class table_t = move_vector>
using fl_move_impl = fl_permutation_impl<empty_data_columns, false, store_absolute_positions, exponential_search, alphabet_t, table_t>;

} // namespace orbit::rlbwt

#endif /* end of include guard: _RLBWT_FL_HPP */