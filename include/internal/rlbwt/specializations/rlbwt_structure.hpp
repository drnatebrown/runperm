#ifndef _RLBWT_STRUCTURE_HPP
#define _RLBWT_STRUCTURE_HPP

#include "internal/common.hpp"
#include "internal/move/move_structure.hpp"
#include "internal/rlbwt/specializations/rlbwt_columns.hpp"

// Shared implementation for both RLBWT column types
template<typename RLBWTColsType = RLBWTCols, template<typename> class TableType = MoveVector>
class RLBWTMoveStructure : public MoveStructure<RLBWTColsType, TableType> {
    using Base = MoveStructure<RLBWTColsType, TableType>;  
public:
    // Sets Columns, ColsTraits, and NumCols
    MOVE_CLASS_TRAITS(typename Base::Columns)
    using Position = typename ColsTraits::Position;

    RLBWTMoveStructure() = default;

    RLBWTMoveStructure(const std::vector<uchar>& bwt_chars, const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const ulint domain, const uchar char_width, const SplitParams& split_params = SplitParams())
    : Base::table(find_structure(bwt_chars, lengths, interval_permutation, domain, char_width, split_params)), Base::n(domain), Base::r(Base::table.size()) {}

    RLBWTMoveStructure(PackedVector<Columns> &&structure, const ulint domain) : Base(std::move(structure), domain) {}

    static PackedVector<Columns> find_structure(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const ulint domain, const uchar char_width, const SplitParams& split_params = SplitParams()) = delete;

    // Duplicates from MoveStructure, but adds bwt_chars to widths
    // Pass with the width of the character field as well as the characters themselves
    // Comments for MoveStructure duplicates are not repeated here
    static PackedVector<Columns> find_structure(const std::vector<uchar>& rlbwt_chars, const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const ulint domain, const uchar char_width, const SplitParams& split_params = SplitParams()) {
        /* === Same as MoveStructure === */
        assert(lengths.size() == interval_permutation.size());
        const std::vector<ulint>* final_lengths = &lengths;
        const std::vector<ulint>* final_interval_permutation = &interval_permutation;
        ulint max_length = 0;
        SplitResult split_result;
        Base::apply_splitting(final_lengths, final_interval_permutation, max_length, split_params, split_result);
        /* === End of MoveStructure === */

        // Also initialize with the character width
        PackedVector<Columns> structure(final_lengths->size(), get_move_widths(domain, final_lengths->size(), max_length, char_width));
        
        /* === Same as MoveStructure === */
        auto sorted_indices = Base::compute_sorted_indices(*final_interval_permutation);
        Base::populate_structure(structure, *final_lengths, *final_interval_permutation, sorted_indices);
        /* === End of MoveStructure === */

        // Set the character field
        set_characters(structure, rlbwt_chars, lengths, domain);

        return structure; 
    }

    uchar get_character(size_t i) const {
        return Base::table.template get<ColsTraits::CHARACTER>(i);
    }
    uchar get_character(Position pos) const {
        return Base::get_character(pos.interval);
    }

private:
    static std::array<uchar, NumCols> get_move_widths(const ulint domain, const ulint runs, const ulint max_length, const uchar char_width) {
        auto widths = Base::get_move_widths(domain, runs, max_length);
        widths[static_cast<size_t>(ColsTraits::CHARACTER)] = char_width;
        return widths;
    }

    // Copy characters over while accounting for splitting
    static void set_characters(PackedVector<Columns>& structure, const std::vector<uchar>& rlbwt_chars, const std::vector<ulint>& lengths, const ulint domain) {
        assert(rlbwt_chars.size() == lengths.size());
        
        auto get_structure_length = [&](size_t idx) {
            if constexpr (ColsTraits::RELATIVE) {
                return (idx == structure.size() - 1) 
                    ? domain - structure.template get<ColsTraits::PRIMARY>(idx)
                    : structure.template get<ColsTraits::PRIMARY>(idx + 1) - structure.template get<ColsTraits::PRIMARY>(idx);
            } else {
                return structure.template get<ColsTraits::PRIMARY>(idx);
            }
        };

        size_t rlbwt_run_start = 0;

        size_t structure_run_start = 0;
        size_t structure_run_idx = 0;
        for (size_t i = 0; i < rlbwt_chars.size(); ++i) {
            while (structure_run_idx < structure.size() && structure_run_start < rlbwt_run_start + lengths[i]) {
                structure.template set<ColsTraits::CHARACTER>(structure_run_idx, rlbwt_chars[i]);
                ++structure_run_idx;
                structure_run_start += get_structure_length(structure_run_idx);
            }
            rlbwt_run_start += lengths[i];
        }
    }
};

#endif /* end of include guard: _RLBWT_STRUCTURE_HPP */