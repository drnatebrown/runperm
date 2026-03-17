#ifndef _RLBWT_LF_HPP
#define _RLBWT_LF_HPP

#include "common.hpp"
#include "internal/rlbwt/specializations/runperm_rlbwt.hpp"
#include "internal/ds/alphabet.hpp"
#include "internal/rlbwt/specializations/rlbwt_permutation.hpp"
#include "internal/rlbwt/rlbwt_helpers.hpp"

// TODO use IntVector here

template<typename RunColsType,
         bool IntegratedMoveStructure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         bool ExponentialSearch = DEFAULT_EXPONENTIAL_SEARCH,
         typename AlphabetType = Nucleotide,
         template<typename> class TableType = MoveVector>
class RunPermLFImpl : public RunPermRLBWT<RunPermLFImpl<RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, ExponentialSearch, AlphabetType, TableType>,
                         RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, ExponentialSearch, AlphabetType, TableType> {
    using Base = RunPermRLBWT<RunPermLFImpl, RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, ExponentialSearch, AlphabetType, TableType>;
    using BaseColumns = typename Base::BaseColumns;
    using RLBWTPermutation = typename Base::RLBWTPermutation;
public:
    using Base::Base;
    using Base::operator=;
    using Position = typename Base::Position;

    // TODO use IntVector and container templates here
    RLBWTPermutation find_permutation(
        const std::vector<uchar>& rlbwt_heads,
        const std::vector<ulint>& rlbwt_run_lengths,
        const SplitParams& split_params
    ) {
        return RLBWTPermutation::lf_permutation(rlbwt_heads, rlbwt_run_lengths, split_params);
    }

    Position LF(Position pos) {
        return Base::next(pos);
    }

    Position LF(Position pos, ulint steps) {
        return Base::next(pos, steps);
    }
};

template<bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         bool ExponentialSearch = DEFAULT_EXPONENTIAL_SEARCH,
         typename AlphabetType = Nucleotide,
         template<typename> class TableType = MoveVector>
class MoveLFImpl : public MovePermRLBWT<RunPermLFImpl<EmptyRunCols, false, StoreAbsolutePositions, ExponentialSearch, AlphabetType, TableType>, 
                      StoreAbsolutePositions, ExponentialSearch, AlphabetType, TableType> {
    using Base = MovePermRLBWT<RunPermLFImpl<EmptyRunCols, false, StoreAbsolutePositions, ExponentialSearch, AlphabetType, TableType>,
                 StoreAbsolutePositions, ExponentialSearch, AlphabetType, TableType>;
public:
    using Base::Base;
    using Base::operator=;
    using Position = typename Base::Position;

    Position LF(Position pos) {
        return Base::next(pos);
    }

    Position LF(Position pos, ulint steps) {
        return Base::next(pos, steps);
    }
};

template<typename Alphabet=Nucleotide>
using MoveLFImplDefault = MoveLFImpl<DEFAULT_STORE_ABSOLUTE_POSITIONS, DEFAULT_EXPONENTIAL_SEARCH, Alphabet, MoveVector>;

#endif /* end of include guard: _RLBWT_LF_HPP */