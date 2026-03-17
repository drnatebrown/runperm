#ifndef _INTERNAL_RUNPERM_RLBWT_HPP
#define _INTERNAL_RUNPERM_RLBWT_HPP

#include "common.hpp"
#include "internal/runperm/runperm.hpp"
#include "internal/ds/alphabet.hpp"
#include "internal/rlbwt/specializations/rlbwt_columns.hpp"
#include "internal/rlbwt/specializations/rlbwt_structure.hpp"
#include "internal/rlbwt/specializations/rlbwt_row.hpp"
#include "internal/rlbwt/specializations/rlbwt_permutation.hpp"

template<typename Derived,
         typename RunColsType,
         bool IntegratedMoveStructure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         bool ExponentialSearch = DEFAULT_EXPONENTIAL_SEARCH,
         typename AlphabetType = Nucleotide,
         template<typename> class TableType = MoveVector>
class RunPermRLBWT : public RunPermImpl<RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, ExponentialSearch, RLBWTCols, RLBWTMoveStructure, TableType> {
    using Base = RunPermImpl<RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, ExponentialSearch, RLBWTCols, RLBWTMoveStructure, TableType>;
protected:
    using BaseColumns = typename Base::BaseColumns;
    using MoveStructurePerm = typename Base::MoveStructurePerm;
    using RLBWTPermutation = RLBWTPermutationImpl<IntVectorAligned, AlphabetType>;
public:
    // Sets Columns, ColsTraits, and NumCols
    MOVE_CLASS_TRAITS(RunColsType)
    using RunCols = typename Base::RunCols;
    using RunData = typename Base::RunData;
    using Position = typename Base::Position;

    RunPermRLBWT() = default;

    template<typename RLBWTPermutationType>
    RunPermRLBWT(const RLBWTPermutationType& permutation, const std::vector<RunData> &run_data) {
        static_assert(std::is_same_v<AlphabetType, typename RLBWTPermutationType::AlphabetTag>, "AlphabetType must be the same as the alphabet type used to create the permutation");

        Base::split_params_ = permutation.get_split_params();
        alphabet = permutation.get_alphabet();
        PackedVector<BaseColumns> base_structure = Base::MoveStructureBase::find_structure(permutation);
        if (run_data.size() == permutation.intervals()) {
            Base::populate_structure(std::move(base_structure), run_data, permutation.domain(), permutation.runs());
        }
        else if (run_data.size() == permutation.runs()) {
            throw std::invalid_argument("Run data size is same as number of runs, not intervals after splitting; avoid splitting, manually split run data, or use permutation copy split.");
        } else {
            throw std::invalid_argument("Run data size must be the same as the number of intervals (user defined splits) or number of runs (no splitting).");
        }
    }

    template<class Container1, class Container2>
    RunPermRLBWT(const Container1 &rlbwt_heads, const Container2 &rlbwt_run_lengths, const std::vector<RunData> &run_data)
        : RunPermRLBWT(rlbwt_heads, rlbwt_run_lengths, NO_SPLITTING, run_data) {}

    template<class Container1, class Container2>
    RunPermRLBWT(const Container1 &rlbwt_heads, const Container2 &rlbwt_run_lengths, const SplitParams &split_params, const std::vector<RunData> &run_data)
        : RunPermRLBWT(rlbwt_heads, rlbwt_run_lengths, split_params,
            [&run_data](ulint orig_interval, ulint orig_interval_length, ulint new_offset_from_orig_start, ulint new_length) {
                return run_data[orig_interval];
            }
        ){}

    // TODO reduce code duplication with the constructors next
    template<class Container1, class Container2>
    RunPermRLBWT(const Container1 &rlbwt_heads, const Container2 &rlbwt_run_lengths, const SplitParams &split_params, const std::vector<RunData> &run_data, std::function<RunData(ulint, ulint, ulint, ulint)> get_run_cols_data) {
        assert(rlbwt_heads.size() == rlbwt_run_lengths.size());

        RLBWTPermutation permutation = find_permutation(rlbwt_heads, rlbwt_run_lengths, split_params);
        alphabet = permutation.get_alphabet();
        Base::split_params_ = split_params;
        PackedVector<BaseColumns> base_structure = Base::MoveStructureBase::find_structure(permutation);

        if (split_params == NO_SPLITTING) {
            Base::populate_structure(std::move(base_structure), run_data, permutation.domain(), permutation.runs());
        } else {
            std::vector<RunData> final_run_data = Base::extend_run_data(rlbwt_run_lengths, base_structure, permutation.domain(), get_run_cols_data);
            Base::populate_structure(std::move(base_structure), final_run_data, permutation.domain(), permutation.runs());
        }
    }

    template<class Container1, class Container2>
    RunPermRLBWT(const Container1 &rlbwt_heads, const Container2 &rlbwt_run_lengths, const SplitParams &split_params, std::function<RunData(ulint, ulint, ulint, ulint)> get_run_cols_data) {
        assert(rlbwt_heads.size() == rlbwt_run_lengths.size());

        RLBWTPermutation permutation = find_permutation(rlbwt_heads, rlbwt_run_lengths, split_params);
        alphabet = permutation.get_alphabet();
        Base::split_params_ = split_params;
        PackedVector<BaseColumns> base_structure = Base::MoveStructureBase::find_structure(permutation);

        /* extend_run_data is required when find_structure applies splitting:
           base_structure may have more rows than run_data; we copy run_data[orig_interval] for each split row */
        std::vector<RunData> final_run_data = Base::extend_run_data(rlbwt_run_lengths, base_structure, permutation.domain(), get_run_cols_data);
        Base::populate_structure(std::move(base_structure), final_run_data, permutation.domain(), permutation.runs());
    }

    static RunPermRLBWT from_structure(PackedVector<BaseColumns> &&structure, const size_t domain, const size_t runs) {
        return RunPermRLBWT(std::move(structure), domain, runs);
    }

    static RunPermRLBWT from_structure(PackedVector<BaseColumns> &&structure, std::vector<RunData> &run_data, const size_t domain, const size_t runs) {
        return RunPermRLBWT(std::move(structure), run_data, domain, runs);
    }

    static RunPermRLBWT from_move_structure(MoveStructurePerm &&ms) {
        return RunPermRLBWT(std::move(ms));
    }

    static RunPermRLBWT from_move_structure(MoveStructurePerm &&ms, std::vector<RunData> &run_data) {
        return RunPermRLBWT(std::move(ms), run_data);
    }

    uchar get_character(ulint interval) {
        return alphabet.unmap_char(Base::template get_base_column<BaseColumns::CHARACTER>(interval));
    }
    uchar get_character(Position pos) {
        return alphabet.unmap_char(Base::template get_base_column<BaseColumns::CHARACTER>(pos.interval));
    }

    std::vector<uchar> get_alphabet() const {
        std::vector<uchar> sigma;
        sigma.reserve(static_cast<size_t>(alphabet.size()));
        for (size_t i = 0; i < static_cast<size_t>(alphabet.size()); ++i)
            sigma.push_back(alphabet.unmap_char(static_cast<uchar>(i)));
        return sigma;
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
    RLBWTPermutation find_permutation(
        const std::vector<uchar>& rlbwt_heads,
        const std::vector<ulint>& rlbwt_run_lengths,
        const SplitParams& split_params
    ) {
        return static_cast<Derived*>(this)->find_permutation(rlbwt_heads, rlbwt_run_lengths, split_params);
    }
};

// A wrapper around RunPermRLBWT without any run data, essentially just a MoveStructure for RLBWT
template<typename Derived,
         bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         bool ExponentialSearch = DEFAULT_EXPONENTIAL_SEARCH,
         typename AlphabetType = Nucleotide,
         template<typename> class TableType = MoveVector>
class MovePermRLBWT {
private:
    using RunPermRLBWTType = RunPermRLBWT<Derived, EmptyRunCols, false, StoreAbsolutePositions, ExponentialSearch, AlphabetType, TableType>;
    RunPermRLBWTType run_perm_rlbwt;
    
public:
    using Position = typename RunPermRLBWTType::Position;
    
    MovePermRLBWT() = default;

    template<typename RLBWTPermutationType>
    MovePermRLBWT(const RLBWTPermutationType &permutation) {
        std::vector<std::array<ulint, 0>> empty_run_data(permutation.intervals());
        run_perm_rlbwt = RunPermRLBWTType(permutation, empty_run_data);
    }

    MovePermRLBWT(const std::vector<uchar> &bwt, SplitParams split_params = SplitParams()) {
        auto [rlbwt_heads, rlbwt_run_lengths] = bwt_to_rlbwt(bwt);
        std::vector<std::array<ulint, 0>> empty_run_data(rlbwt_heads.size());
        run_perm_rlbwt = RunPermRLBWTType(rlbwt_heads, rlbwt_run_lengths, split_params, empty_run_data);
    }
    
    // Constructor from RLBWT data
    template<class Container1, class Container2>
    MovePermRLBWT(const Container1 &rlbwt_heads, const Container2 &rlbwt_run_lengths, 
                  SplitParams split_params = SplitParams()) {
        std::vector<std::array<ulint, 0>> empty_run_data(rlbwt_heads.size());
        run_perm_rlbwt = RunPermRLBWTType(rlbwt_heads, rlbwt_run_lengths, split_params, empty_run_data);
    }
    
    Position first() { return run_perm_rlbwt.first(); }
    void last() { run_perm_rlbwt.last(); }
    ulint get_length(ulint interval) const { return run_perm_rlbwt.get_length(interval); }
    ulint get_length(Position pos) const { return run_perm_rlbwt.get_length(pos); }
    Position next(Position pos) { return run_perm_rlbwt.next(pos); }
    Position next(Position pos, ulint steps) { return run_perm_rlbwt.next(pos, steps); }
    Position up(Position pos) { return run_perm_rlbwt.up(pos); }
    Position down(Position pos) { return run_perm_rlbwt.down(pos); }
    
    ulint domain() const { return run_perm_rlbwt.domain(); }
    ulint runs() const { return run_perm_rlbwt.runs(); }
    ulint intervals() const { return run_perm_rlbwt.intervals(); }
    
    // RLBWT-specific method
    uchar get_character(ulint interval) { return run_perm_rlbwt.get_character(interval); }
    uchar get_character(Position pos) { return run_perm_rlbwt.get_character(pos); }
    
    size_t serialize(std::ostream& out) { return run_perm_rlbwt.serialize(out); }
    void load(std::istream& in) { run_perm_rlbwt.load(in); }
};

#endif /* end of include guard: _INTERNAL_RUNPERM_RLBWT_HPP */