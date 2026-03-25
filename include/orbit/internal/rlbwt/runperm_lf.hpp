#ifndef _RLBWT_LF_HPP
#define _RLBWT_LF_HPP

#include "orbit/common.hpp"
#include "orbit/internal/rlbwt/specializations/runperm_rlbwt.hpp"
#include "orbit/internal/ds/alphabet.hpp"
#include "orbit/internal/rlbwt/specializations/rlbwt_interval_encoding.hpp"
#include "orbit/internal/rlbwt/rlbwt_helpers.hpp"

namespace orbit::rlbwt {
// TODO use int_vector here

template<typename data_columns_t,
         bool integrated_move_structure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         bool exponential_search = DEFAULT_EXPONENTIAL_SEARCH,
         typename alphabet_t = nucleotide,
         template<typename> class table_t = move_vector>
class runperm_lf_impl : public runperm_rlbwt<runperm_lf_impl<data_columns_t, integrated_move_structure, store_absolute_positions, exponential_search, alphabet_t, table_t>,
                         data_columns_t, integrated_move_structure, store_absolute_positions, exponential_search, alphabet_t, table_t> {
    using base = runperm_rlbwt<runperm_lf_impl, data_columns_t, integrated_move_structure, store_absolute_positions, exponential_search, alphabet_t, table_t>;
    using rlbwt_interval_encoding = typename base::rlbwt_interval_encoding;
public:
    using data_columns = typename base::data_columns;
    using base::base;
    using base::operator=;
    using position = typename base::position;

    // TODO use int_vector and container templates here
    rlbwt_interval_encoding find_interval_encoding(
        const std::vector<uchar>& rlbwt_heads,
        const std::vector<ulint>& rlbwt_run_lengths,
        const split_params& sp
    ) {
        return rlbwt_interval_encoding::lf_interval_encoding(rlbwt_heads, rlbwt_run_lengths, sp);
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
         template<typename> class table_t = move_vector>
class move_lf_impl : public moveperm_rlbwt<runperm_lf_impl<empty_data_columns, false, store_absolute_positions, exponential_search, alphabet_t, table_t>, 
                      store_absolute_positions, exponential_search, alphabet_t, table_t> {
    using base = moveperm_rlbwt<runperm_lf_impl<empty_data_columns, false, store_absolute_positions, exponential_search, alphabet_t, table_t>,
                 store_absolute_positions, exponential_search, alphabet_t, table_t>;
public:
    using base::base;
    using base::operator=;
    using position = typename base::position;

    position LF(position pos) {
        return base::next(pos);
    }

    position LF(position pos, ulint steps) {
        return base::next(pos, steps);
    }
};

template<typename alphabet_t=nucleotide>
using move_lf_impl_default = move_lf_impl<DEFAULT_STORE_ABSOLUTE_POSITIONS, DEFAULT_EXPONENTIAL_SEARCH, alphabet_t, move_vector>;

} // namespace orbit::rlbwt

#endif /* end of include guard: _RLBWT_LF_HPP */