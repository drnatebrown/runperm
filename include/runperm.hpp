// Convenience Header for RunPerm and MovePerm
#ifndef _PUBLIC_RUNPERM_HPP
#define _PUBLIC_RUNPERM_HPP

#include "internal/common.hpp"
#include "internal/runperm/runperm.hpp"



// =============================== RunPerm ===============================

// Actual implementation, see documentation below
// Advanced users can use the full implementation in include/internal/runperm/runperm.hpp
template<typename RunColsType,
         bool IntegratedMoveStructure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS>
class RunPerm : public RunPermImpl<RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, MoveCols, MoveStructure, MoveVector> {
    using Base = RunPermImpl<RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, MoveCols, MoveStructure, MoveVector>;
public:
    using Base::RunData;
    using Base::MoveStructurePerm;
    using Base::Position;

    using Base::Base;
    using Base::operator=;
};

/* === Basic types === */
template<typename RunColsType>
using RunPermSeperated = RunPerm<RunColsType, false, false>; // Same as RunPerm<RunColsType>, the default
template<typename RunColsType>
using RunPermIntegrated = RunPerm<RunColsType, true, false>;
template<typename RunColsType>
using RunPermSeperatedAbsolute = RunPerm<RunColsType, false, true>;
template<typename RunColsType>
using RunPermIntegratedAbsolute = RunPerm<RunColsType, true, true>;

/* === Simplified interface for basic users, see full RunPermImpl in include/internal/runperm/runperm.hpp for more template parameters ===
*   template<typename RunColsType,
*         bool IntegratedMoveStructure = DEFAULT_INTEGRATED_MOVE_STRUCTURE, // Whether to integrate the run data alongside the move structure for cache locality, default is false
*         bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS> // Whether to store absolute positions instead of interval/offset to support idx lookups, default is false
*  class RunPerm : public RunPermHelper<RunColsType, IntegratedMoveStructure, StoreAbsolutePositions> {
*   public:
*    // === Implemented types ===
*    using RunData = typename Base::RunData; // std::array<ulint, RunColsType::COUNT>
*    using Position = typename MoveStructure::Position; // see set position method below
*    using MoveStructurePerm = typename Base::MoveStructurePerm; // MoveStructurePerm as determined by the template parameters
*
*    // === Basic Constructors ===
*    // lengths -> length of each interval which permutes contiguously
*    // interval_permutation -> permutation position of the first position of each interval
*    // domain -> domain of the permutation, i.e. a permutatation over 1..n has domain n
*    // run_data -> run data for each interval, the size of this vector should be the same as the number of intervals
*    // split_params -> parameters for splitting the move structure, (max_allowed_length, balancing_factor)
*    RunPerm(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const ulint domain, const std::vector<RunData> &run_data);
*    RunPerm(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const ulint domain, const SplitParams &split_params, const std::vector<RunData> &run_data);
*
*    // === Advanced Constructors ===
*    // get_run_cols_data -> function to get run data for each interval, (orig_interval, orig_interval_length, new_offset_from_orig_start, new_length) -> run data
*    // structure -> pre-computed move structure stored in PackedVector
*    // ms -> pre-computed move structure stored in MoveStructure
*    RunPerm(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const ulint domain, const SplitParams &split_params, std::function<RunData(ulint, ulint, ulint, ulint)> get_run_cols_data);
*    RunPerm(PackedVector<MoveCols> &&structure, const ulint domain);
*    RunPerm(PackedVector<MoveCols> &&structure, std::vector<RunData> &run_data, const ulint domain);
*    RunPerm(MoveStructure &&ms, const ulint domain);
*    RunPerm(MoveStructure &&ms, std::vector<RunData> &run_data, const ulint domain);
*
*    // === Navigation methods ===
*   // get_run_cols_data -> function to get run data for each interval, (orig_interval, orig_interval_length, new_offset_from_orig_start, new_length) -> run data
*    // structure -> pre-computed move structure stored in PackedVector
*    // ms -> pre-computed move structure stored in MoveStructure
*    RunPerm(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const ulint domain, const SplitParams &split_params, std::function<RunData(ulint, ulint, ulint, ulint)> get_run_cols_data);
*    RunPerm(PackedVector<MoveCols> &&structure, const ulint domain);
*    RunPerm(PackedVector<MoveCols> &&structure, std::vector<RunData> &run_data, const ulint domain);
*    RunPerm(MoveStructure &&ms, const ulint domain);
*    RunPerm(MoveStructure &&ms, std::vector<RunData> &run_data, const ulint domain);
*
*    // === Navigation methods ===
*    void next(); // Apply permutation
*    void next(ulint steps); // Apply permutation multiple times
*    bool up(); // Move up one interval
*    bool down(); // Move down one interval
*    void first(); // Move to first interval
*    void last(); // Move to last interval
*    
*    // === Position methods ===
*    Position get_position() const; // Get current position
*    void set_position(Position pos); // Set current position, (interval, offset) for relative positions, (interval, offset, idx) for absolute positions
*
*    // === Query methods ===
*    ulint size() const; // Get size of permutation
*    ulint move_runs() const; // Get number of runs/intervals in move structure
*    ulint permutation_runs() const; // Get number of runs/intervals in original permutation
*
*    // === Run data access ===
*    template<RunColsType Col>
*    ulint get() const; // Get value of run data column
*    template<RunColsType Col>
*    ulint get(size_t i) const; // Get value of run data column for interval i
*    ulint get_length() const; // Get length of current interval
*    ulint get_length(size_t i) const; // Get length of interval i
*
*    // === Search methods ===
*    template<RunColsType Col>
*    std::optional<Position> pred(ulint val); // Get position of largest idx before or at position run with matching run data value
*    template<RunColsType Col>
*    std::optional<Position> succ(ulint val); // Get position of smallest idx after or at position run with matching run data value
*
*    // === Serialization ===
*    size_t serialize(std::ostream& os); // Serialize data structure to ostream
*    void load(std::istream& is); // Load data structure from istream
*};
*/



// =============================== MovePerm ===============================

// Actual implementation, see documentation below
// Advanced users can use the full implementation in include/internal/runperm/runperm.hpp
template<bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS>
class MovePerm : public MovePermImpl<StoreAbsolutePositions, MoveCols, MoveStructure, MoveVector> {
    using Base = MovePermImpl<StoreAbsolutePositions, MoveCols, MoveStructure, MoveVector>;
public:
    using Base::Base;
    using Base::operator=;
};

// Basic types
using MovePermAbsolute = MovePerm<true>;
using MovePermRelative = MovePerm<false>; // Same as MovePerm<>, the default

/* === Simplified interface for basic users, see full MovePermImpl in include/internal/runperm/runperm.hpp for more template parameters ===
/* NOTE: All methods here share the documentation above for RunPerm
* template<bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS>
* class MovePerm {
*
* public:
*    // === Implemented types ===
*    using Position = typename RunPermType::Position;
*    
*    // === Constructors ===
*    MovePermImpl(std::vector<ulint>& permutation, SplitParams split_params = SplitParams()); // Constructor from permutation vector
*    MovePermImpl(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const ulint domain, SplitParams split_params = SplitParams()); // See RunPerm
*    
*    // === Navigation methods ===
*    void next(); // Apply permutation
*    void next(ulint steps); // Apply permutation multiple times
*    bool up(); // Move up one interval
*    bool down(); // Move down one interval
*    void first(); // Move to first interval
*    void last(); // Move to last interval
*
*    // === Query methods ===
*    ulint size() const; // Get size of permutation
*    ulint move_runs() const; // Get number of runs/intervals in move structure
*    ulint permutation_runs() const; // Get number of runs/intervals in original permutation
*
*    // === Run data access ===
*    ulint get_length() const; // Get length of current interval
*    ulint get_length(size_t i) const; // Get length of interval i
*
*    // === Serialization ===
*    size_t serialize(std::ostream& os); // Serialize data structure to ostream
*    void load(std::istream& is); // Load data structure from istream
* };
*/

#endif /* end of include guard: _PUBLIC_RUNPERM_HPP */