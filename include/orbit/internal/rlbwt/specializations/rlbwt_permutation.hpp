#ifndef _RLBWT_PERMUTATION_HPP
#define _RLBWT_PERMUTATION_HPP

#include "orbit/common.hpp"
#include "orbit/internal/perm/permutation_impl.hpp"
#include "orbit/internal/ds/alphabet.hpp"
#include "orbit/internal/rlbwt/specializations/rlbwt_columns.hpp"
#include "orbit/internal/rlbwt/specializations/rlbwt_structure.hpp"
#include "orbit/internal/rlbwt/specializations/rlbwt_row.hpp"
#include "orbit/internal/rlbwt/specializations/rlbwt_interval_encoding.hpp"

namespace orbit::rlbwt {

template<typename derived,
         typename data_columns_t = empty_data_columns,
         bool integrated_move_structure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         bool exponential_search = DEFAULT_EXPONENTIAL_SEARCH,
         typename alphabet_t = nucleotide,
         typename base_columns_t = rlbwt_columns,
         template<typename, template<typename> class> class move_structure_t = rlbwt_move_structure,
         template<typename> class table_t = move_vector>
class rlbwt_permutation : public permutation_impl<data_columns_t, integrated_move_structure, store_absolute_positions, exponential_search, base_columns_t, move_structure_t, table_t> {
    using base = permutation_impl<data_columns_t, integrated_move_structure, store_absolute_positions, exponential_search, base_columns_t, move_structure_t, table_t>;
protected:
    using base_columns = typename base::base_columns;
    using move_structure_perm = typename base::move_structure_perm;
    // TODO: Add invertible template parameter
    public:
    // Sets Columns, ColsTraits, and NumCols
    MOVE_CLASS_TRAITS(base_columns)
    using data_columns = typename base::data_columns;
    using data_tuple = typename base::data_tuple;
    using position = typename base::position;
    using rlbwt_interval_encoding_t = rlbwt_interval_encoding_impl<cols_traits::INVERTIBLE, int_vector_aligned, alphabet_t>;

    rlbwt_permutation() = default;

    // Gated constructors for empty data columns
    template<typename rlbwt_interval_encoding_impl_t,
             typename dc = data_columns_t,
             std::enable_if_t<std::is_same_v<dc, empty_data_columns>, int> = 0>
    rlbwt_permutation(const rlbwt_interval_encoding_impl_t &enc) 
    : rlbwt_permutation(enc, std::vector<data_tuple>(enc.intervals())) {}

    template<typename dc = data_columns_t,
             std::enable_if_t<std::is_same_v<dc, empty_data_columns>, int> = 0>
    rlbwt_permutation(const std::vector<uchar>& bwt, const split_params& sp = {})
    : rlbwt_permutation([&]{
        auto [heads, lens] = bwt_to_rlbwt(bwt);
        return rlbwt_interval_encoding_t(heads, lens, sp);
    }()) {}

    template<typename container1_t, typename container2_t,
             typename dc = data_columns_t,
             std::enable_if_t<std::is_same_v<dc, empty_data_columns>, int> = 0>
    rlbwt_permutation(const container1_t &rlbwt_heads, const container2_t &rlbwt_run_lengths, 
                  const split_params& sp = split_params()) 
    : rlbwt_permutation(rlbwt_heads, rlbwt_run_lengths, sp, std::vector<data_tuple>(rlbwt_heads.size())) {}

    // Regular non-gated constructors for user defined data columns
    template<typename rlbwt_interval_encoding_impl_t>
    rlbwt_permutation(const rlbwt_interval_encoding_impl_t& enc, const std::vector<data_tuple> &run_data) {
        static_assert(std::is_same_v<alphabet_t, typename rlbwt_interval_encoding_impl_t::alphabet_tag>, "alphabet_t must be the same as the alphabet type used to create the interval encoding");
        alphabet_ = enc.get_alphabet();
        base::build_from_interval_encoding(enc, run_data);
    }

    template<typename container1_t, typename container2_t>
    rlbwt_permutation(const container1_t &rlbwt_heads, const container2_t &rlbwt_run_lengths, const std::vector<data_tuple> &run_data)
        : rlbwt_permutation(rlbwt_heads, rlbwt_run_lengths, NO_SPLITTING, run_data) {}

    template<typename container1_t, typename container2_t>
    rlbwt_permutation(const container1_t &rlbwt_heads, const container2_t &rlbwt_run_lengths, const split_params& sp, const std::vector<data_tuple> &run_data)
        : rlbwt_permutation(rlbwt_heads, rlbwt_run_lengths, sp,
            [&run_data](ulint orig_interval, ulint orig_interval_length, ulint new_offset_from_orig_start, ulint new_length) {
                return run_data[orig_interval];
            }
        ){}

    template<typename container1_t, typename container2_t>
    rlbwt_permutation(const container1_t &rlbwt_heads, const container2_t &rlbwt_run_lengths, const split_params& sp, const std::vector<data_tuple> &run_data, std::function<data_tuple(ulint, ulint, ulint, ulint)> get_run_cols_data) {
        build_from_rlbwt(rlbwt_heads, rlbwt_run_lengths, sp, &run_data, std::move(get_run_cols_data));
    }

    template<typename container1_t, typename container2_t>
    rlbwt_permutation(const container1_t &rlbwt_heads, const container2_t &rlbwt_run_lengths, const split_params& sp, std::function<data_tuple(ulint, ulint, ulint, ulint)> get_run_cols_data) {
        build_from_rlbwt(rlbwt_heads, rlbwt_run_lengths, sp, nullptr, std::move(get_run_cols_data));
    }

    static rlbwt_permutation from_structure(packed_vector<base_columns> &&structure, const size_t domain, const size_t runs) {
        return rlbwt_permutation(std::move(structure), domain, runs);
    }

    static rlbwt_permutation from_structure(packed_vector<base_columns> &&structure, std::vector<data_tuple> &run_data, const size_t domain, const size_t runs) {
        return rlbwt_permutation(std::move(structure), run_data, domain, runs);
    }

    static rlbwt_permutation from_move_structure(move_structure_perm &&ms) {
        return rlbwt_permutation(std::move(ms));
    }

    static rlbwt_permutation from_move_structure(move_structure_perm &&ms, std::vector<data_tuple> &run_data) {
        return rlbwt_permutation(std::move(ms), run_data);
    }

    uchar get_character(ulint interval) {
        return alphabet_.unmap_char(base::template get_base_column<base_columns::CHARACTER>(interval));
    }
    uchar get_character(position pos) {
        return alphabet_.unmap_char(base::template get_base_column<base_columns::CHARACTER>(pos.interval));
    }

    // Returns row/offset of largest idx before or at position run with matching character
    std::optional<position> pred_char(position position, uchar c) {
        if (get_character(position) == c) {
            return position;
        }

        while (get_character(position) != c)
        {
            if (position.interval == 0) return std::nullopt;
            --position.interval;
        }
        position.offset = base::move_structure.get_length(position.interval) - 1;
        if constexpr (store_absolute_positions) {
            position.idx = base::move_structure.get_start(position.interval) + position.offset;
        }
        return position;
    }

    // Returns row/offset of smallest idx after or at position run with matching character
    std::optional<position> succ_char(position position, uchar c) {
        if (get_character(position) == c) {
            return position;
        }

        while (get_character(position) != c)
        {
            if (position.interval == base::move_structure.runs() - 1) return std::nullopt;
            ++position.interval;
        }
        position.offset = 0;
        if constexpr (store_absolute_positions) {
            position.idx = base::move_structure.get_start(position.interval);
        }
        return position;
    }

    std::vector<uchar> get_alphabet() const {
        std::vector<uchar> sigma;
        sigma.reserve(static_cast<size_t>(alphabet_.size()));
        for (size_t i = 0; i < static_cast<size_t>(alphabet_.size()); ++i)
            sigma.push_back(alphabet_.unmap_char(static_cast<uchar>(i)));
        return sigma;
    }

    size_t serialize(std::ostream& out) {
        return base::serialize(out) + alphabet_.serialize(out);
    }

    void load(std::istream& in) {
        base::load(in);
        alphabet_.load(in);
    }

protected:
    alphabet_t alphabet_;

    template<typename container1_t, typename container2_t>
    void build_from_rlbwt(const container1_t &rlbwt_heads, const container2_t &rlbwt_run_lengths, const split_params& sp, const std::vector<data_tuple>* run_data, std::function<data_tuple(ulint, ulint, ulint, ulint)> get_run_cols_data)
    {
        assert(rlbwt_heads.size() == rlbwt_run_lengths.size());
        auto enc = find_interval_encoding(rlbwt_heads, rlbwt_run_lengths, sp);
        alphabet_ = enc.get_alphabet();
        base::build_from_interval_encoding_callback(enc, rlbwt_run_lengths, run_data, std::move(get_run_cols_data));
    }

    // Special logic for LF or FL
    template<typename container1_t, typename container2_t>
    rlbwt_interval_encoding_t find_interval_encoding(
        const container1_t& rlbwt_heads,
        const container2_t& rlbwt_run_lengths,
        const split_params& sp
    ) {
        return static_cast<derived*>(this)->find_interval_encoding(rlbwt_heads, rlbwt_run_lengths, sp);
    }
};

template<typename derived,
         bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         bool exponential_search = DEFAULT_EXPONENTIAL_SEARCH,
         typename alphabet_t = nucleotide,
         typename base_columns_t = rlbwt_columns,
         template<typename, template<typename> class> class move_structure_t = rlbwt_move_structure,
         template<typename> class table_t = move_vector>
using rlbwt_move = rlbwt_permutation<derived, empty_data_columns, false, store_absolute_positions, exponential_search, alphabet_t, base_columns_t, move_structure_t, table_t>;

} // namespace orbit::rlbwt

#endif /* end of include guard: _RLBWT_PERMUTATION_HPP */