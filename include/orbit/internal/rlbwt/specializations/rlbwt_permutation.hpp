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
         template<typename> class table_t = move_vector>
class rlbwt_permutation : public permutation_impl<data_columns_t, integrated_move_structure, store_absolute_positions, exponential_search, rlbwt_columns, rlbwt_move_structure, table_t> {
    using base = permutation_impl<data_columns_t, integrated_move_structure, store_absolute_positions, exponential_search, rlbwt_columns, rlbwt_move_structure, table_t>;
protected:
    using base_columns = typename base::base_columns;
    using move_structure_perm = typename base::move_structure_perm;
    using rlbwt_interval_encoding = rlbwt_interval_encoding_impl<int_vector_aligned, alphabet_t>;
public:
    // Sets Columns, ColsTraits, and NumCols
    MOVE_CLASS_TRAITS(data_columns_t)
    using data_columns = typename base::data_columns;
    using data_tuple = typename base::data_tuple;
    using position = typename base::position;

    rlbwt_permutation() = default;

    // Gated constructors for empty data columns
    template<typename rlbwt_interval_encoding_t,
             typename dc = data_columns_t,
             std::enable_if_t<std::is_same_v<dc, empty_data_columns>, int> = 0>
    rlbwt_permutation(const rlbwt_interval_encoding_t &enc) 
    : rlbwt_permutation(enc, std::vector<data_tuple>(enc.intervals())) {}

    template<typename dc = data_columns_t,
             std::enable_if_t<std::is_same_v<dc, empty_data_columns>, int> = 0>
    rlbwt_permutation(const std::vector<uchar>& bwt, const split_params& sp = {})
    : rlbwt_permutation([&]{
        auto [heads, lens] = bwt_to_rlbwt(bwt);
        return rlbwt_interval_encoding(heads, lens, sp);
    }()) {}

    template<typename container1_t, typename container2_t,
             typename dc = data_columns_t,
             std::enable_if_t<std::is_same_v<dc, empty_data_columns>, int> = 0>
    rlbwt_permutation(const container1_t &rlbwt_heads, const container2_t &rlbwt_run_lengths, 
                  const split_params& sp = split_params()) 
    : rlbwt_permutation(rlbwt_heads, rlbwt_run_lengths, sp, std::vector<data_tuple>(rlbwt_heads.size())) {}

    // Regular non-gated constructors for user defined data columns
    template<typename rlbwt_interval_encoding_t>
    rlbwt_permutation(const rlbwt_interval_encoding_t& enc, const std::vector<data_tuple> &run_data) {
        static_assert(std::is_same_v<alphabet_t, typename rlbwt_interval_encoding_t::alphabet_tag>, "alphabet_t must be the same as the alphabet type used to create the interval encoding");

        base::split_params_ = enc.get_split_params();
        alphabet_ = enc.get_alphabet();
        packed_vector<base_columns> base_structure = base::move_structure_base::find_structure(enc);
        if (run_data.size() == enc.intervals()) {
            base::populate_structure(std::move(base_structure), run_data, enc.domain(), enc.runs());
        }
        else if (run_data.size() == enc.runs()) {
            throw std::invalid_argument("Run data size is same as number of runs, not intervals after splitting; avoid splitting, manually split run data, or use permutation copy split.");
        } else {
            throw std::invalid_argument("Run data size must be the same as the number of intervals (user defined splits) or number of runs (no splitting).");
        }
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

    // TODO reduce code duplication with the constructors next
    template<typename container1_t, typename container2_t>
    rlbwt_permutation(const container1_t &rlbwt_heads, const container2_t &rlbwt_run_lengths, const split_params& sp, const std::vector<data_tuple> &run_data, std::function<data_tuple(ulint, ulint, ulint, ulint)> get_run_cols_data) {
        assert(rlbwt_heads.size() == rlbwt_run_lengths.size());

        rlbwt_interval_encoding enc = find_interval_encoding(rlbwt_heads, rlbwt_run_lengths, sp);
        alphabet_ = enc.get_alphabet();
        base::split_params_ = sp;
        packed_vector<base_columns> base_structure = base::move_structure_base::find_structure(enc);

        if (sp == NO_SPLITTING) {
            base::populate_structure(std::move(base_structure), run_data, enc.domain(), enc.runs());
        } else {
            std::vector<data_tuple> final_run_data = base::extend_run_data(rlbwt_run_lengths, base_structure, enc.domain(), get_run_cols_data);
            base::populate_structure(std::move(base_structure), final_run_data, enc.domain(), enc.runs());
        }
    }

    template<typename container1_t, typename container2_t>
    rlbwt_permutation(const container1_t &rlbwt_heads, const container2_t &rlbwt_run_lengths, const split_params& sp, std::function<data_tuple(ulint, ulint, ulint, ulint)> get_run_cols_data) {
        assert(rlbwt_heads.size() == rlbwt_run_lengths.size());

        rlbwt_interval_encoding enc = find_interval_encoding(rlbwt_heads, rlbwt_run_lengths, sp);
        alphabet_ = enc.get_alphabet();
        base::split_params_ = sp;
        packed_vector<base_columns> base_structure = base::move_structure_base::find_structure(enc);

        /* extend_run_data is required when find_structure applies splitting:
           base_structure may have more rows than run_data; we copy run_data[orig_interval] for each split row */
        std::vector<data_tuple> final_run_data = base::extend_run_data(rlbwt_run_lengths, base_structure, enc.domain(), get_run_cols_data);
        base::populate_structure(std::move(base_structure), final_run_data, enc.domain(), enc.runs());
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

    // Special logic for LF or FL
    rlbwt_interval_encoding find_interval_encoding(
        const std::vector<uchar>& rlbwt_heads,
        const std::vector<ulint>& rlbwt_run_lengths,
        const split_params& sp
    ) {
        return static_cast<derived*>(this)->find_interval_encoding(rlbwt_heads, rlbwt_run_lengths, sp);
    }
};

template<typename derived,
         bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         bool exponential_search = DEFAULT_EXPONENTIAL_SEARCH,
         typename alphabet_t = nucleotide,
         template<typename> class table_t = move_vector>
using rlbwt_move = rlbwt_permutation<derived, empty_data_columns, false, store_absolute_positions, exponential_search, alphabet_t, table_t>;

} // namespace orbit::rlbwt

#endif /* end of include guard: _RLBWT_PERMUTATION_HPP */