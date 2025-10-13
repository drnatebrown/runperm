// Convenience Header for RunPermLF, RunPermFL, MoveLF, MoveFL
#ifndef _RLBWT_HPP
#define _RLBWT_HPP

#include "internal/common.hpp"
// #include "internal/rlbwt/runperm_fl.hpp"
#include "internal/rlbwt/runperm_lf.hpp"

// === RunPermLF ===
template<typename RunColsType,
         bool IntegratedMoveStructure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         typename AlphabetType = Nucleotide>
class RunPermLF : public RunPermLFImpl<RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, AlphabetType, MoveVector> {
    using Base = RunPermLFImpl<RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, AlphabetType, MoveVector>;
public:
    using Base::Base;
    using Base::operator=;
};

template<typename RunColsType, typename AlphabetType = Nucleotide>
using RunPermLFIntegrated = RunPermLF<RunColsType, true, DEFAULT_STORE_ABSOLUTE_POSITIONS, AlphabetType>;
template<typename RunColsType, typename AlphabetType = Nucleotide>
using RunPermLFIntegratedAbsolute = RunPermLF<RunColsType, true, true, AlphabetType>;
template<typename RunColsType, typename AlphabetType = Nucleotide>
using RunPermLFSeperated = RunPermLF<RunColsType, false, DEFAULT_STORE_ABSOLUTE_POSITIONS, AlphabetType>;
template<typename RunColsType, typename AlphabetType = Nucleotide>
using RunPermLFSeperatedAbsolute = RunPermLF<RunColsType, false, true, AlphabetType>;

// === MoveLF ===

template<bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         typename AlphabetType = Nucleotide>
class MoveLF : public MoveLFImpl<StoreAbsolutePositions, AlphabetType, MoveVector> {
    using Base = MoveLFImpl<StoreAbsolutePositions, AlphabetType, MoveVector>;
public:
    using Base::Base;
    using Base::operator=;
};

#endif /* end of include guard: _PUBLIC_RUNPERM_RLBWT_HPP */