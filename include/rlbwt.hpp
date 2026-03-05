// Convenience header for RLBWT-based permutations.
//
// This header exposes high-level wrappers for working with run-length BWT
// permutations:
//   - LF / FL navigation over an RLBWT (`MoveLF`, `MoveFL`, `RunPermLF`, `RunPermFL`)
//   - Phi / InvPhi permutations for suffix-array based locate queries
//
// Usage overview:
//   - You start from a run-length BWT, given as:
//       std::vector<uchar> bwt_heads;       // run heads (characters)
//       std::vector<ulint> bwt_run_lengths; // run lengths
//   - For "move-only" access (no user run data), use:
//       MoveLF<> lf(bwt_heads, bwt_run_lengths);
//       MoveFL<> fl(bwt_heads, bwt_run_lengths);
//   - For RLBWT permutations that also carry user-defined run data, use:
//       enum class RunCols { FIELD1, FIELD2, COUNT };
//       where you must include COUNT as the last entry to signal number of fields.
//       using RunData = DataTuple<RunCols>;
//       std::vector<RunData> run_data(bwt_heads.size());
//       RunPermLF<RunCols> lf_rp(bwt_heads, bwt_run_lengths, run_data);
//       RunPermFL<RunCols> fl_rp(bwt_heads, bwt_run_lengths, run_data);
//
//   - Phi / InvPhi structures are built from an RLBWT via helpers
//     declared in the internal headers (see README for examples):
//       auto [phi_lengths, phi_interval_perm, n] =
//           rlbwt_to_phi(bwt_heads, bwt_run_lengths);
//       auto [inv_lengths, inv_interval_perm, n2] =
//           rlbwt_to_invphi(bwt_heads, bwt_run_lengths);
//
//     and then wrapped by:
//       MovePhi phi(phi_lengths, phi_interval_perm, n);
//       MoveInvPhi invphi(inv_lengths, inv_interval_perm, n2);
//       RunPermPhi<RunCols> rp_phi(phi_lengths, phi_interval_perm, n, run_data);
//       RunPermInvPhi<RunCols> rp_inv(inv_lengths, inv_interval_perm, n2, run_data);
//
// Aliases:
//   - *Separated* vs *Integrated* controls whether run data is stored in a
//     separate packed vector (Separated, default) or integrated into the
//     move structure rows (Integrated, better locality but larger rows).
//   - *Absolute* variants additionally store absolute positions for each
//     entry, which enables direct index lookups at the cost of extra space.
//
// See the README and tests under `tests/unit/rlbwt/` and
// `tests/integration/rlbwt_test.cpp` for concrete usage.
// Expanded documentation below:
#ifndef _RLBWT_HPP
#define _RLBWT_HPP

#include "internal/common.hpp"
#include "internal/rlbwt/runperm_lf.hpp"
#include "internal/rlbwt/runperm_fl.hpp"
#include "internal/rlbwt/runperm_phi.hpp"
#include "internal/rlbwt/runperm_invphi.hpp"

// === RunPermLF ===
template<typename DataColumns,
         bool IntegratedMoveStructure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         typename AlphabetType = Nucleotide>
class RunPermLF : public RunPermLFImpl<DataColumns, IntegratedMoveStructure, StoreAbsolutePositions, DEFAULT_EXPONENTIAL_SEARCH, AlphabetType, MoveVector> {
    using Base = RunPermLFImpl<DataColumns, IntegratedMoveStructure, StoreAbsolutePositions, DEFAULT_EXPONENTIAL_SEARCH, AlphabetType, MoveVector>;
public:
    using RunDataColumns = typename Base::RunCols;
    using RunDataTuple = typename Base::RunData;
    
    using Base::Base;
    using Base::operator=;

    using Position = typename Base::Position;
};

// Takes std::vector<uchar> bwt_and std::vector<ulint> as bwt_heads and bwt_run_lengths as input in place of lengths and interval permutations
// Otherwise, same as RunPerm
// Also implements LF, LF(steps), get_character(), get_character(row)

template<typename DataColumns, typename AlphabetType = Nucleotide>
using RunPermLFSeparated = RunPermLF<DataColumns, false, false, AlphabetType>; // Default
template<typename DataColumns, typename AlphabetType = Nucleotide>
using RunPermLFIntegrated = RunPermLF<DataColumns, true, false, AlphabetType>;
template<typename DataColumns, typename AlphabetType = Nucleotide>
using RunPermLFSeparatedAbsolute = RunPermLF<DataColumns, false, true, AlphabetType>;
template<typename DataColumns, typename AlphabetType = Nucleotide>
using RunPermLFIntegratedAbsolute = RunPermLF<DataColumns, true, true, AlphabetType>;

// === MoveLF ===
template<bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         typename AlphabetType = Nucleotide>
class MoveLF : public MoveLFImpl<StoreAbsolutePositions, DEFAULT_EXPONENTIAL_SEARCH, AlphabetType, MoveVector> {
    using Base = MoveLFImpl<StoreAbsolutePositions, DEFAULT_EXPONENTIAL_SEARCH, AlphabetType, MoveVector>;
public:
    using Base::Base;
    using Base::operator=;

    using Position = typename Base::Position;
};

// See above

// === RunPermFL ===
template<typename DataColumns,
         bool IntegratedMoveStructure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         typename AlphabetType = Nucleotide>
class RunPermFL : public RunPermFLImpl<DataColumns, IntegratedMoveStructure, StoreAbsolutePositions, DEFAULT_EXPONENTIAL_SEARCH, AlphabetType, MoveVector> {
    using Base = RunPermFLImpl<DataColumns, IntegratedMoveStructure, StoreAbsolutePositions, DEFAULT_EXPONENTIAL_SEARCH, AlphabetType, MoveVector>;
public:
    using RunDataColumns = typename Base::RunCols;
    using RunDataTuple = typename Base::RunData;

    using Base::Base;
    using Base::operator=;

    using Position = typename Base::Position;
};

// Takes std::vector<uchar> bwt_and std::vector<ulint> as bwt_heads and bwt_run_lengths as input in place of lengths and interval permutations
// Otherwise, same as RunPerm
// Also implements FL(pos), FL(pos, steps), get_character(pos), get_character(interval)

template<typename DataColumns, typename AlphabetType = Nucleotide>
using RunPermFLSeparated = RunPermFL<DataColumns, false, false, AlphabetType>; // Default
template<typename DataColumns, typename AlphabetType = Nucleotide>
using RunPermFLIntegrated = RunPermFL<DataColumns, true, false, AlphabetType>;
template<typename DataColumns, typename AlphabetType = Nucleotide>
using RunPermFLSeparatedAbsolute = RunPermFL<DataColumns, false, true, AlphabetType>;
template<typename DataColumns, typename AlphabetType = Nucleotide>
using RunPermFLIntegratedAbsolute = RunPermFL<DataColumns, true, true, AlphabetType>;

// === MoveFL ===
template<bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         typename AlphabetType = Nucleotide>
class MoveFL : public MoveFLImpl<StoreAbsolutePositions, DEFAULT_EXPONENTIAL_SEARCH, AlphabetType, MoveVector> {
    using Base = MoveFLImpl<StoreAbsolutePositions, DEFAULT_EXPONENTIAL_SEARCH, AlphabetType, MoveVector>;
public:
    using Base::Base;
    using Base::operator=;

    using Position = typename Base::Position;
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