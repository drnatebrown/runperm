// Convenience Header for RunPermLF, RunPermFL, MoveLF, MoveFL
#ifndef _RLBWT_HPP
#define _RLBWT_HPP

#include "internal/common.hpp"
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
using RunPermLFSeperated = RunPermLF<RunColsType, false, false, AlphabetType>; // Default
template<typename RunColsType, typename AlphabetType = Nucleotide>
using RunPermLFIntegrated = RunPermLF<RunColsType, true, false, AlphabetType>;
template<typename RunColsType, typename AlphabetType = Nucleotide>
using RunPermLFSeperatedAbsolute = RunPermLF<RunColsType, false, true, AlphabetType>;
template<typename RunColsType, typename AlphabetType = Nucleotide>
using RunPermLFIntegratedAbsolute = RunPermLF<RunColsType, true, true, AlphabetType>;

// === MoveLF ===
template<bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         typename AlphabetType = Nucleotide>
class MoveLF : public MoveLFImpl<StoreAbsolutePositions, AlphabetType, MoveVector> {
    using Base = MoveLFImpl<StoreAbsolutePositions, AlphabetType, MoveVector>;
public:
    using Base::Base;
    using Base::operator=;
};

#include "internal/rlbwt/runperm_fl.hpp"

// === RunPermFL ===
template<typename RunColsType,
         bool IntegratedMoveStructure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         typename AlphabetType = Nucleotide>
class RunPermFL : public RunPermFLImpl<RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, AlphabetType, MoveVector> {
    using Base = RunPermFLImpl<RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, AlphabetType, MoveVector>;
public:
    using Base::Base;
    using Base::operator=;
};

template<typename RunColsType, typename AlphabetType = Nucleotide>
using RunPermFLSeperated = RunPermFL<RunColsType, false, false, AlphabetType>; // Default
template<typename RunColsType, typename AlphabetType = Nucleotide>
using RunPermFLIntegrated = RunPermFL<RunColsType, true, false, AlphabetType>;
template<typename RunColsType, typename AlphabetType = Nucleotide>
using RunPermFLSeperatedAbsolute = RunPermFL<RunColsType, false, true, AlphabetType>;
template<typename RunColsType, typename AlphabetType = Nucleotide>
using RunPermFLIntegratedAbsolute = RunPermFL<RunColsType, true, true, AlphabetType>;

// === MoveFL ===
template<bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         typename AlphabetType = Nucleotide>
class MoveFL : public MoveFLImpl<StoreAbsolutePositions, AlphabetType, MoveVector> {
    using Base = MoveFLImpl<StoreAbsolutePositions, AlphabetType, MoveVector>;
public:
    using Base::Base;
    using Base::operator=;
};

#endif /* end of include guard: _PUBLIC_RUNPERM_RLBWT_HPP */