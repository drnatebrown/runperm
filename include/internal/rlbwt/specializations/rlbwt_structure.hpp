#ifndef _RLBWT_STRUCTURE_HPP
#define _RLBWT_STRUCTURE_HPP

#include "common.hpp"
#include "internal/move/move_structure.hpp"
#include "internal/rlbwt/specializations/rlbwt_columns.hpp"
#include "internal/rlbwt/specializations/rlbwt_permutation.hpp"

// Shared implementation for both RLBWT column types
template<typename RLBWTColsType = RLBWTCols, template<typename> class TableType = MoveVector>
class RLBWTMoveStructure : public MoveStructure<RLBWTColsType, TableType> {
    using Base = MoveStructure<RLBWTColsType, TableType>;  
public:
    // Sets Columns, ColsTraits, and NumCols
    MOVE_CLASS_TRAITS(typename Base::Columns)
    using Position = typename ColsTraits::Position;

    RLBWTMoveStructure() = default;

    RLBWTMoveStructure(const std::vector<uchar>& head_chars, const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const ulint domain, const uchar sigma, const SplitParams& split_params = SplitParams())
    : Base(find_structure(head_chars, lengths, interval_permutation, domain, sigma, split_params), domain, lengths.size()) {}

    template<typename RLBWTPermutationType>
    RLBWTMoveStructure(const RLBWTPermutationType& permutation) {
        this->n = permutation.domain();
        this->r = permutation.runs();
        this->table = find_structure(permutation);
    }

    RLBWTMoveStructure(PackedVector<Columns> &&structure, const size_t domain, ulint runs) : Base(std::move(structure), domain, runs) {}

    static PackedVector<Columns> find_structure(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const ulint domain, const uchar char_width, const SplitParams& split_params = SplitParams()) = delete;

    // When the permutation is not already computed
    static PackedVector<Columns> find_structure(const std::vector<uchar>& head_chars, const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const ulint domain, const uchar sigma, const SplitParams& split_params = SplitParams()) {
        assert(head_chars.size() == lengths.size());

        PermutationImpl<> permutation(lengths, interval_permutation, split_params);
        assert(permutation.runs() == lengths.size());
        assert(permutation.domain() == domain);

        PackedVector<Columns> structure(permutation.intervals(), get_move_widths(domain, permutation.intervals(), permutation.max_length(), sigma));
        Base::populate_structure(structure, permutation);
        set_characters(structure, head_chars, lengths, domain);

        return structure; 
    }

    template<typename PermutationType>
    static PackedVector<Columns> find_structure(const std::vector<uchar>& rlbwt_chars, const PermutationType& permutation, const uchar sigma) {
        assert(rlbwt_chars.size() == permutation.intervals());

        // Also initialize with the character width
        PackedVector<Columns> structure(permutation.intervals(), get_move_widths(permutation.domain(), permutation.intervals(), permutation.max_length(), sigma));
        Base::populate_structure(structure, permutation);
        // Set the character field
        set_characters(structure, rlbwt_chars);

        return structure; 
    }
    
    template<typename RLBWTPermutationType>
    static PackedVector<Columns> find_structure(const RLBWTPermutationType& permutation) {
        assert(permutation.get_heads().size() == permutation.intervals());

        // Also initialize with the character width
        PackedVector<Columns> structure(permutation.intervals(), get_move_widths(permutation.domain(), permutation.intervals(), permutation.max_length(), permutation.sigma()));
        Base::populate_structure(structure, permutation);
        // Set the character field
        set_characters(structure, permutation.get_heads());

        return structure; 
    }

    uchar get_character(size_t i) const {
        return Base::table.template get<ColsTraits::CHARACTER>(i);
    }
    uchar get_character(Position pos) const {
        return get_character(pos.interval);
    }

private:
    static std::array<uchar, NumCols> get_move_widths(const ulint domain, const ulint runs, const ulint max_length, const uchar sigma) {
        auto widths = Base::get_move_widths(domain, runs, max_length);
        widths[static_cast<size_t>(ColsTraits::CHARACTER)] = bit_width(sigma - 1);
        return widths;
    }

    // Copy characters over while accounting for splitting
    template<class Collection>
    static void set_characters(PackedVector<Columns>& structure, const Collection& rlbwt_chars) {
        assert(rlbwt_chars.size() == structure.size());

        for (size_t i = 0; i < rlbwt_chars.size(); ++i) {
            structure.template set<to_cols(ColsTraits::CHARACTER)>(i, static_cast<ulint>(rlbwt_chars[i]));
        }
    }

    static void set_characters(PackedVector<Columns>& structure, const std::vector<uchar>& rlbwt_chars, const std::vector<ulint>& lengths, const ulint domain) {
        assert(rlbwt_chars.size() == lengths.size());
        
        auto get_structure_length = [&structure, &domain](size_t row) {
            if constexpr (!ColsTraits::RELATIVE) {
                return (row == structure.size() - 1) 
                    ? domain - structure.template get<ColsTraits::PRIMARY>(row)
                    : structure.template get<ColsTraits::PRIMARY>(row + 1) - structure.template get<ColsTraits::PRIMARY>(row);
            } else {
                return structure.template get<ColsTraits::PRIMARY>(row);
            }
        };

        size_t rlbwt_run_start = 0;
        size_t structure_run_start = 0;
        size_t structure_run_idx = 0;
        for (size_t i = 0; i < rlbwt_chars.size(); ++i) {
            while (structure_run_idx < structure.size() && structure_run_start < rlbwt_run_start + lengths[i]) {
                structure.template set<to_cols(ColsTraits::CHARACTER)>(structure_run_idx, rlbwt_chars[i]);
                structure_run_start += get_structure_length(structure_run_idx++);
            }
            rlbwt_run_start += lengths[i];
        }
    }
};

#endif /* end of include guard: _RLBWT_STRUCTURE_HPP */