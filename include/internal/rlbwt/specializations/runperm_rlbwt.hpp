#ifndef _INTERNAL_RUNPERM_RLBWT_HPP
#define _INTERNAL_RUNPERM_RLBWT_HPP

#include "internal/common.hpp"
#include "internal/runperm/runperm.hpp"
#include "internal/ds/alphabet.hpp"
#include "internal/rlbwt/specializations/rlbwt_columns.hpp"
#include "internal/rlbwt/specializations/rlbwt_structure.hpp"
#include "internal/rlbwt/specializations/rlbwt_row.hpp"

template<typename Derived,
         typename RunColsType,
         bool IntegratedMoveStructure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         typename AlphabetType = Nucleotide,
         template<typename> class TableType = MoveVector>
class RunPermRLBWT : public RunPermImpl<RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, RLBWTCols, RLBWTMoveStructure, TableType> {
    using Base = RunPermImpl<RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, RLBWTCols, RLBWTMoveStructure, TableType>;
protected:
    using BaseColumns = typename Base::BaseColumns;
    using MoveStructurePerm = typename Base::MoveStructurePerm;

public:
    // Sets Columns, ColsTraits, and NumCols
    MOVE_CLASS_TRAITS(RunColsType)
    using RunData = typename Base::RunData;
    using Position = typename Base::Position;

    RunPermRLBWT() = default;

    RunPermRLBWT(const std::vector<uchar> &rlbwt_heads, const std::vector<ulint> &rlbwt_run_lengths, const std::vector<RunData> &run_data) {
        assert(rlbwt_heads.size() == rlbwt_run_lengths.size());
        assert(run_data.size() == rlbwt_heads.size());
        Base::orig_intervals = rlbwt_heads.size();
        Base::position = Position();

        alphabet = AlphabetType();
        ulint num_chars;
        PackedVector<BaseColumns> base_structure;
        find_permutation_and_alphabet(rlbwt_heads, rlbwt_run_lengths, alphabet, num_chars, base_structure);

        Base::populate_structure(std::move(base_structure), run_data, num_chars);
    }

    RunPermRLBWT(const std::vector<uchar> &rlbwt_heads, const std::vector<ulint> &rlbwt_run_lengths, const SplitParams &split_params, const std::vector<RunData> &run_data)
        : RunPermRLBWT(rlbwt_heads, rlbwt_run_lengths, split_params,
            [&run_data](ulint orig_interval, ulint orig_interval_length, ulint new_offset_from_orig_start, ulint new_length) {
                return run_data[orig_interval];
            }
        ){}

    RunPermRLBWT(const std::vector<uchar> &rlbwt_heads, const std::vector<ulint> &rlbwt_run_lengths, const SplitParams &split_params, std::function<RunData(ulint, ulint, ulint, ulint)> get_run_cols_data) {
        assert(rlbwt_heads.size() == rlbwt_run_lengths.size());
        Base::orig_intervals = rlbwt_heads.size();
        Base::position = Position();

        alphabet = AlphabetType();
        ulint num_chars;
        PackedVector<BaseColumns> base_structure;
        find_permutation_and_alphabet(rlbwt_heads, rlbwt_run_lengths, alphabet, num_chars, base_structure);

        std::vector<RunData> final_run_data = Base::extend_run_data(rlbwt_run_lengths, num_chars, base_structure, get_run_cols_data);
        Base::populate_structure(std::move(base_structure), final_run_data, num_chars);
    }

    RunPermRLBWT(PackedVector<Columns> &&structure, const ulint domain) : Base::move_structure(std::move(structure), domain), Base::position(Position()), Base::orig_intervals(structure.size()) {
        static_assert(IntegratedMoveStructure, "Cannot construct RunPermRLBWT with pre-computed permutation structure if not integrating user data with move structure");
    }

    RunPermRLBWT(PackedVector<BaseColumns> &&structure, std::vector<RunData> &run_data, const ulint domain) : Base::move_structure(std::move(structure), domain), Base::position(Position()), Base::orig_intervals(structure.size()) {
        static_assert(!IntegratedMoveStructure, "Cannot construct RunPermRLBWT with pre-computed permutation structure if integrating user data with move structure");
        populate_run_data(std::move(structure), run_data, domain);
    }

    RunPermRLBWT(MoveStructurePerm &&ms, const ulint domain) : Base::move_structure(std::move(ms)), Base::position(Position()), Base::orig_intervals(Base::move_structure.size()) {
        static_assert(IntegratedMoveStructure, "Cannot construct RunPermRLBWT with pre-computed move structure if not integrating user data with move structure");
    }

    RunPermRLBWT(MoveStructurePerm &&ms, std::vector<RunData> &run_data, const ulint domain) : Base::move_structure(std::move(ms)), Base::position(Position()), Base::orig_intervals(Base::move_structure.size()) {
        static_assert(!IntegratedMoveStructure, "Cannot construct RunPermRLBWT with pre-computed move structure if integrating user data with move structure");
        auto run_cols_widths = get_run_cols_widths(run_data);
        fill_seperated_data(run_data, run_cols_widths);
    }

    uchar get_character() {
        return alphabet.unmap_char(Base::template get_base_column<BaseColumns::CHARACTER>());
    }

    size_t serialize(std::ostream& out) {
        return Base::serialize(out) + alphabet.serialize(out);
    }

    void load(std::istream& in) {
        Base::load(in);
        alphabet.load(in);
    }

protected:
    AlphabetType alphabet;

    // Special logic for LF or FL
    void find_permutation_and_alphabet(
        const std::vector<uchar>& rlbwt_heads,
        const std::vector<ulint>& rlbwt_run_lengths,
        AlphabetType& alphabet,
        ulint& num_chars,
        PackedVector<BaseColumns>& base_structure
    ) {
        return static_cast<Derived*>(this)->find_permutation_and_alphabet(
            rlbwt_heads, rlbwt_run_lengths, alphabet, num_chars, base_structure);
    }
};

// A wrapper around RunPermRLBWT without any run data, essentially just a MoveStructure for RLBWT
template<typename Derived,
         bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         typename AlphabetType = Nucleotide,
         template<typename> class TableType = MoveVector>
class MovePermRLBWT {
private:
    using RunPermRLBWTType = RunPermRLBWT<Derived, EmptyRunCols, false, StoreAbsolutePositions, AlphabetType, TableType>;
    RunPermRLBWTType run_perm_rlbwt;
    
public:
    using Position = typename RunPermRLBWTType::Position;
    
    MovePermRLBWT() = default;

    MovePermRLBWT(const std::vector<uchar> &bwt, SplitParams split_params = SplitParams()) {
        auto [rlbwt_heads, rlbwt_run_lengths] = bwt_to_rlbwt(bwt);
        std::vector<std::array<ulint, 0>> empty_run_data(rlbwt_heads.size());
        run_perm_rlbwt = RunPermRLBWTType(rlbwt_heads, rlbwt_run_lengths, split_params, empty_run_data);
    }
    
    // Constructor from RLBWT data
    MovePermRLBWT(const std::vector<uchar> &rlbwt_heads, 
                  const std::vector<ulint> &rlbwt_run_lengths, 
                  SplitParams split_params = SplitParams()) {
        std::vector<std::array<ulint, 0>> empty_run_data(rlbwt_heads.size());
        run_perm_rlbwt = RunPermRLBWTType(rlbwt_heads, rlbwt_run_lengths, split_params, empty_run_data);
    }
    
    void first() { run_perm_rlbwt.first(); }
    void last() { run_perm_rlbwt.last(); }
    void move(Position pos) const { run_perm_rlbwt.move(pos); }
    Position get_position() const { return run_perm_rlbwt.get_position(); }
    void next() { run_perm_rlbwt.next(); }
    
    ulint size() const { return run_perm_rlbwt.size(); }
    ulint move_runs() const { return run_perm_rlbwt.move_runs(); }
    ulint permutation_runs() const { return run_perm_rlbwt.permutation_runs(); }
    
    // RLBWT-specific method
    uchar get_character() { return run_perm_rlbwt.get_character(); }
    
    size_t serialize(std::ostream& out) { return run_perm_rlbwt.serialize(out); }
    void load(std::istream& in) { run_perm_rlbwt.load(in); }
};

#endif /* end of include guard: _INTERNAL_RUNPERM_RLBWT_HPP */