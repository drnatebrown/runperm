#ifndef _RLBWT_LF_HPP
#define _RLBWT_LF_HPP

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
         typename base_columns_t = rlbwt_columns,
         template<typename, template<typename> class> class move_structure_t = rlbwt_move_structure,
         template<typename> class table_t = move_vector>
class lf_permutation_impl : public rlbwt_permutation<lf_permutation_impl<data_columns_t, integrated_move_structure, store_absolute_positions, exponential_search, alphabet_t, base_columns_t, move_structure_t, table_t>,
                         data_columns_t, integrated_move_structure, store_absolute_positions, exponential_search, alphabet_t, base_columns_t, move_structure_t, table_t> {
    using base = rlbwt_permutation<lf_permutation_impl, data_columns_t, integrated_move_structure, store_absolute_positions, exponential_search, alphabet_t, base_columns_t, move_structure_t, table_t>;
    using rlbwt_interval_encoding_t = typename base::rlbwt_interval_encoding_t;
public:
    using data_columns = typename base::data_columns;
    using base::base;
    using base::operator=;
    using position = typename base::position;

    template<typename container1_t, typename container2_t>
    rlbwt_interval_encoding_t find_interval_encoding(
        const container1_t& rlbwt_heads,
        const container2_t& rlbwt_run_lengths,
        const split_params& sp
    ) {
        return rlbwt_interval_encoding_t::lf_interval_encoding(rlbwt_heads, rlbwt_run_lengths, sp);
    }

    position LF(position pos) {
        return base::next(pos);
    }

    position LF(position pos, ulint steps) {
        return base::next(pos, steps);
    }
};

template<bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         bool exponential_search = DEFAULT_EXPONENTIAL_SEARCH,
         typename alphabet_t = nucleotide,
         typename base_columns_t = rlbwt_columns,
         template<typename, template<typename> class> class move_structure_t = rlbwt_move_structure,
         template<typename> class table_t = move_vector>
using lf_move_impl = lf_permutation_impl<empty_data_columns, false, store_absolute_positions, exponential_search, alphabet_t, base_columns_t, move_structure_t, table_t>;

template<typename alphabet_t=nucleotide>
using lf_move_impl_default = lf_move_impl<DEFAULT_STORE_ABSOLUTE_POSITIONS, DEFAULT_EXPONENTIAL_SEARCH, alphabet_t, rlbwt_columns, rlbwt_move_structure, move_vector>;

} // namespace orbit::rlbwt

#endif /* end of include guard: _RLBWT_LF_HPP */