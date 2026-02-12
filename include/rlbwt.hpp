// Convenience Header for RunPermLF, RunPermFL, MoveLF, MoveFL
#ifndef _RLBWT_HPP
#define _RLBWT_HPP

#include "internal/common.hpp"
#include "internal/rlbwt/runperm_lf.hpp"
#include "internal/rlbwt/runperm_fl.hpp"
#include "internal/rlbwt/runperm_phi.hpp"
#include "internal/rlbwt/runperm_invphi.hpp"

// === RunPermLF ===
template<typename RunColsType,
         bool IntegratedMoveStructure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         typename AlphabetType = Nucleotide>
class RunPermLF : public RunPermLFImpl<RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, DEFAULT_EXPONENTIAL_SEARCH, AlphabetType, MoveVector> {
    using Base = RunPermLFImpl<RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, DEFAULT_EXPONENTIAL_SEARCH, AlphabetType, MoveVector>;
public:
    using Base::Base;
    using Base::operator=;
};

// Takes std::vector<uchar> bwt_and std::vector<ulint> as bwt_heads and bwt_run_lengths as input in place of lengths and interval permutations
// Otherwise, same as RunPerm
// Also implements LF, LF(steps), get_character(), get_character(row)

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
class MoveLF : public MoveLFImpl<StoreAbsolutePositions, DEFAULT_EXPONENTIAL_SEARCH, AlphabetType, MoveVector> {
    using Base = MoveLFImpl<StoreAbsolutePositions, DEFAULT_EXPONENTIAL_SEARCH, AlphabetType, MoveVector>;
public:
    using Base::Base;
    using Base::operator=;
};

// See above

// === RunPermFL ===
template<typename RunColsType,
         bool IntegratedMoveStructure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         typename AlphabetType = Nucleotide>
class RunPermFL : public RunPermFLImpl<RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, DEFAULT_EXPONENTIAL_SEARCH, AlphabetType, MoveVector> {
    using Base = RunPermFLImpl<RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, DEFAULT_EXPONENTIAL_SEARCH, AlphabetType, MoveVector>;
public:
    using Base::Base;
    using Base::operator=;
};

// Takes std::vector<uchar> bwt_and std::vector<ulint> as bwt_heads and bwt_run_lengths as input in place of lengths and interval permutations
// Otherwise, same as RunPerm
// Also implements FL(pos), FL(pos, steps), get_character(pos), get_character(interval)

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
class MoveFL : public MoveFLImpl<StoreAbsolutePositions, DEFAULT_EXPONENTIAL_SEARCH, AlphabetType, MoveVector> {
    using Base = MoveFLImpl<StoreAbsolutePositions, DEFAULT_EXPONENTIAL_SEARCH, AlphabetType, MoveVector>;
public:
    using Base::Base;
    using Base::operator=;
};

// See above

// === RunPermPhi ===
// Need to call rlbwt_to_phi(rlbwt_heads, rlbwt_run_lengths) to get lengths and interval permutations
// Otherwise, same as RunPerm
// RunPhiPerm is just a named specilization which sets StoreAbsolutePositions to true
// Also implements Phi(pos), Phi(pos, steps), SA(pos)

// === MovePhi ===
// See above

// === RunPermInvPhi ===
// Need to call rlbwt_to_invphi(rlbwt_heads, rlbwt_run_lengths) to get lengths and interval permutations
// Otherwise, same as RunPerm -->
// RunInvPhiPerm is just a named specilization which sets StoreAbsolutePositions to true
// Also implements InvPhi(pos), InvPhi(pos, steps), SA(pos)

// === MoveInvPhi ===
// See above

#endif /* end of include guard: _PUBLIC_RUNPERM_RLBWT_HPP */