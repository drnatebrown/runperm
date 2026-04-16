#ifndef _RLBWT_STRUCTURE_HPP
#define _RLBWT_STRUCTURE_HPP

#include "orbit/common.hpp"
#include "orbit/internal/move/move_structure_impl.hpp"
#include "orbit/internal/rlbwt/specializations/rlbwt_columns.hpp"
#include "orbit/internal/rlbwt/specializations/rlbwt_interval_encoding.hpp"

namespace orbit::rlbwt {

// Shared implementation for both RLBWT column types
template<typename rlbwt_columns_t = rlbwt_columns, template<typename> class table_t = move_vector>
class rlbwt_move_structure : public move_structure<rlbwt_columns_t, table_t> {
    using base = move_structure<rlbwt_columns_t, table_t>;  
public:
    // Sets Columns, ColsTraits, and NumCols
    MOVE_CLASS_TRAITS(typename base::columns)
    using position = typename cols_traits::position;

    rlbwt_move_structure() = default;

    rlbwt_move_structure(const std::vector<uchar>& head_chars, const std::vector<ulint>& lengths, const std::vector<ulint>& images, const ulint domain, const uchar sigma, const split_params& sp = split_params())
    : base(find_structure(head_chars, lengths, images, domain, sigma, sp), domain, lengths.size()) {}

    template<typename rlbwt_permutation_t>
    rlbwt_move_structure(const rlbwt_permutation_t& permutation) {
        this->n = permutation.domain();
        this->r = permutation.runs();
        this->table = find_structure(permutation);
    }

    rlbwt_move_structure(packed_vector<columns> &&structure, const size_t domain, ulint runs) : base(std::move(structure), domain, runs) {}

    static packed_vector<columns> find_structure(const std::vector<ulint>& lengths, const std::vector<ulint>& images, const ulint domain, const uchar char_width, const split_params& sp = split_params()) = delete;

    // When the permutation is not already computed
    static packed_vector<columns> find_structure(const std::vector<uchar>& head_chars, const std::vector<ulint>& lengths, const std::vector<ulint>& images, const uchar sigma, const split_params& sp = split_params()) {
        assert(head_chars.size() == lengths.size());

        interval_encoding_impl<> enc(lengths, images, sp);
        assert(enc.runs() == lengths.size());

        packed_vector<columns> structure(enc.intervals(), get_move_widths(enc.domain(), enc.intervals(), enc.max_length(), sigma));
        base::populate_structure(structure, enc);
        set_characters(structure, head_chars, lengths, enc.domain());

        return structure; 
    }

    template<typename interval_encoding_t>
    static packed_vector<columns> find_structure(const std::vector<uchar>& rlbwt_chars, const interval_encoding_t& enc, const uchar sigma) {
        assert(rlbwt_chars.size() == enc.intervals());

        // Also initialize with the character width
        packed_vector<columns> structure(enc.intervals(), get_move_widths(enc.domain(), enc.intervals(), enc.max_length(), sigma));
        base::populate_structure(structure, enc);
        // Set the character field
        set_characters(structure, rlbwt_chars);

        return structure; 
    }
    
    template<typename rlbwt_permutation_t>
    static packed_vector<columns> find_structure(const rlbwt_permutation_t& permutation) {
        assert(permutation.get_heads().size() == permutation.intervals());

        // Also initialize with the character width
        packed_vector<columns> structure(permutation.intervals(), get_move_widths(permutation.domain(), permutation.intervals(), permutation.max_length(), permutation.sigma()));
        base::populate_structure(structure, permutation);
        // Set the character field
        set_characters(structure, permutation.get_heads());

        return structure; 
    }

    uchar get_character(size_t i) const {
        return base::table.template get<cols_traits::CHARACTER>(i);
    }
    uchar get_character(position pos) const {
        return get_character(pos.interval);
    }

private:
    static std::array<uchar, num_cols> get_move_widths(const ulint domain, const ulint runs, const ulint max_length, const uchar sigma) {
        auto widths = base::get_move_widths(domain, runs, max_length);
        widths[static_cast<size_t>(cols_traits::CHARACTER)] = bit_width(sigma - 1);
        return widths;
    }

    // Copy characters over while accounting for splitting
    template<typename collection_t>
    static void set_characters(packed_vector<columns>& structure, const collection_t& rlbwt_chars) {
        assert(rlbwt_chars.size() == structure.size());

        for (size_t i = 0; i < rlbwt_chars.size(); ++i) {
            structure.template set<to_cols(cols_traits::CHARACTER)>(i, static_cast<ulint>(rlbwt_chars[i]));
        }
    }

    static void set_characters(packed_vector<columns>& structure, const std::vector<uchar>& rlbwt_chars, const std::vector<ulint>& lengths, const ulint domain) {
        assert(rlbwt_chars.size() == lengths.size());
        
        auto get_structure_length = [&structure, &domain](size_t row) {
            if constexpr (!cols_traits::RELATIVE) {
                return (row == structure.size() - 1) 
                    ? domain - structure.template get<cols_traits::PRIMARY>(row)
                    : structure.template get<cols_traits::PRIMARY>(row + 1) - structure.template get<cols_traits::PRIMARY>(row);
            } else {
                return structure.template get<cols_traits::PRIMARY>(row);
            }
        };

        size_t rlbwt_run_start = 0;
        size_t structure_run_start = 0;
        size_t structure_run_idx = 0;
        for (size_t i = 0; i < rlbwt_chars.size(); ++i) {
            while (structure_run_idx < structure.size() && structure_run_start < rlbwt_run_start + lengths[i]) {
                structure.template set<to_cols(cols_traits::CHARACTER)>(structure_run_idx, rlbwt_chars[i]);
                structure_run_start += get_structure_length(structure_run_idx++);
            }
            rlbwt_run_start += lengths[i];
        }
    }
};

} // namespace orbit::rlbwt

#endif /* end of include guard: _RLBWT_STRUCTURE_HPP */