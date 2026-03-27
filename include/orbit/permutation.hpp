// Convenience header for run-length permutations (permutation / move_permutation).
//
// This header exposes high-level wrappers around the internal move structure
// for generic "runny" permutations:
//   - permutation: run-length encoded permutations with user-defined run data.
//   - move_permutation: run-length encoded permutations without extra data.
//
// Typical usage:
//   - Start from a permutation over [0, n) that is "runny", or from its
//     run decomposition:
//       std::vector<ulint> lengths;          // lengths of contiguous intervals
//       std::vector<ulint> interval_perm;    // perm position of each interval head
//   - Or, from its interval encoding given by the permutation class:
//       orbit::interval_encoding enc;
//   - Define an enum for run data columns and a columns_tuple alias:
//       enum class data_columns { VAL1, VAL2, COUNT };
//       // Must include COUNT as the last entry to signal number of fields.
//       // alternatively, use the DEFINE_COLUMNS macro:
//       DEFINE_COLUMNS(data_columns, VAL1, VAL2)
//       using data_tuple = orbit::columns_tuple<data_columns>;
//       std::vector<data_tuple> data(lengths.size()); // Some data with tuples per run
//   - Construct a permutation instance:
//       permutation<data_columns> p(lengths, interval_perm, domain, data);
//
// Aliases:
//   - permutation_separated / permutation_integrated select whether run data is stored
//     in a separate packed vector (Separated, default) or integrated into the
//     move structure rows (Integrated, better locality but larger rows).
//   - permutation_separated_absolute / permutation_integrated_absolute additionally store
//     absolute positions, enabling direct index lookups at extra space cost.
//   - move_permutation_relative / move_permutation_absolute provide analogous choices for the
//     move-only interface without run data.
//   - runperm / runperm_separated / runperm_integrated / runperm_separated_absolute / runperm_integrated_absolute
//     are alpha release aliases for permutation / permutation_separated / permutation_integrated / permutation_separated_absolute / permutation_integrated_absolute
//   - moveperm / moveperm_absolute / moveperm_relative
//     are alpha release aliases for move_permutation / move_permutation_absolute / move_permutation_relative
//
// See the README and tests under `tests/unit/perm/` and
// `tests/integration/permutation_test.cpp` for concrete usage.
// Expanded documentation below:
#ifndef _PUBLIC_PERMUTATION_HPP
#define _PUBLIC_PERMUTATION_HPP

#include "orbit/common.hpp"
#include "orbit/internal/move/interval_encoding_impl.hpp"
#include "orbit/internal/perm/permutation_impl.hpp"

namespace orbit {

// =============================== permutation ===============================

// Actual implementation, see documentation below
// Advanced users can use the full implementation in include/internal/perm/permutation_impl.hpp
template<typename data_columns_t = empty_data_columns,
         bool integrated_move_structure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS>
class permutation : public permutation_impl<data_columns_t, integrated_move_structure, store_absolute_positions, DEFAULT_EXPONENTIAL_SEARCH, move_columns, move_structure, move_vector> {
    using base = permutation_impl<data_columns_t, integrated_move_structure, store_absolute_positions, DEFAULT_EXPONENTIAL_SEARCH, move_columns, move_structure, move_vector>;
public:
    using typename base::data_columns;
    using typename base::data_tuple;
    using typename base::position;

    using base::base;
    using base::operator=;
};

/* === Basic types === */
template<typename data_columns_t = empty_data_columns>
using permutation_separated = permutation<data_columns_t, false, false>; // Same as permutation<data_columns_t>, the default
template<typename data_columns_t = empty_data_columns>
using permutation_integrated = permutation<data_columns_t, true, false>;
template<typename data_columns_t = empty_data_columns>
using permutation_absolute = permutation<data_columns_t, false, true>;
template<typename data_columns_t = empty_data_columns>
using permutation_separated_absolute = permutation<data_columns_t, false, true>;
template<typename data_columns_t = empty_data_columns>
using permutation_integrated_absolute = permutation<data_columns_t, true, true>;

// Alpha release aliases
template<typename data_columns_t = empty_data_columns,
         bool integrated_move_structure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS>
using runperm = permutation<data_columns_t, integrated_move_structure, store_absolute_positions>;
template<typename data_columns_t = empty_data_columns>
using runperm_separated = permutation_separated<data_columns_t>;
template<typename data_columns_t = empty_data_columns>
using runperm_integrated = permutation_integrated<data_columns_t>;
template<typename data_columns_t = empty_data_columns>
using runperm_absolute = permutation_absolute<data_columns_t>;
template<typename data_columns_t = empty_data_columns>
using runperm_separated_absolute = permutation_separated_absolute<data_columns_t>;
template<typename data_columns_t = empty_data_columns>
using runperm_integrated_absolute = permutation_integrated_absolute<data_columns_t>;

/* === Simplified interface for basic users, see full permutation_impl in include/internal/perm/permutation_impl.hpp for more template parameters ===
*   template<typename data_columns_t,
*         bool integrated_move_structure = DEFAULT_INTEGRATED_MOVE_STRUCTURE, // Whether to integrate the run data alongside the move structure for cache locality, default is false
*         bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS> // Whether to store absolute positions instead of interval/offset to support idx lookups, default is false
*  class permutation : public permutation_impl<data_columns_t, integrated_move_structure, store_absolute_positions> {
*   public:
*    // === Implemented types ===
*    using data_tuple = typename base::data_tuple; // std::array<ulint, data_columns::COUNT>
*    using position = typename move_structure::position; // see set position method below
*    using move_structure_perm = typename base::move_structure_perm; // move_structure_perm as determined by the template parameters
*
*    // === Basic Constructors ===
*    // lengths -> length of each interval which permutes contiguously
*    // images -> permutation position of the first position of each interval
*    // domain -> domain of the permutation, i.e. a permutatation over 1..n has domain n
*    // run_data -> run data for each interval, the size of this vector should be the same as the number of intervals
*    // split_params -> parameters for splitting the move structure, (max_allowed_length, balancing_factor)
     // !!! NOTE: If using splitting, user data will be copied by default. See Splitting Constructor if this is not desired.
     // !!!       For this reason, splitting is turned off by default when using the simplified interface.
*    permutation(const std::vector<ulint>& lengths, const std::vector<ulint>& images, const ulint domain, const std::vector<data_tuple> &run_data, split_params split_params = NO_SPLITTING);
*    permutation(const std::vector<ulint>& lengths, const std::vector<ulint>& images, const ulint domain, const split_params &split_params, const std::vector<data_tuple> &run_data);
*
*    // === Splitting Constructor ===
*    // interval_encoding -> object built by passing lengths, images to interval_encoding constructor
*    // run_data -> run data for each interval, the size of this vector should be the same as the number of intervals after splitting by modifiying based on the permutation object
*    permutation(const interval_encoding &enc, const std::vector<data_tuple> &run_data);
*
*    // === Advanced Constructors ===
*    // get_run_cols_data -> function to get run data for each interval, (orig_interval, orig_interval_length, new_offset_from_orig_start, new_length) -> run data
*    // structure -> pre-computed move structure stored in packed_vector
*    // ms -> pre-computed move structure stored in move_structure
*    permutation(const std::vector<ulint>& lengths, const std::vector<ulint>& images, const ulint domain, const SplitParams &split_params, std::function<RunData(ulint, ulint, ulint, ulint)> get_run_cols_data);
*    permutation from_structure(packed_vector<columns> &&structure, const size_t domain, const size_t runs);
*    permutation from_structure(packed_vector<base_columns> &&structure, std::vector<data_tuple> &run_data, const size_t domain, const size_t runs);
*    permutation from_move_structure(move_structure_perm &&ms);
*    permutation from_move_structure(move_structure_perm &&ms, std::vector<data_tuple> &run_data);
*
*    // === Navigation methods ===
*    position next(position pos); // Apply permutation
*    void next(position pos, ulint steps); // Apply permutation multiple times
*    position up(position pos); // Move up one interval, circularly if at top
*    position down(position pos); // Move down one interval, circularly if at bottom
*    position first(); // Move to first interval
*    position last(); // Move to last interval
*
*    // === Query methods ===
*    ulint domain() const; // Get size of permutation
*    ulint intervals() const; // Get number of intervals in move structure (could be greater than runs due to splitting)
*    ulint runs() const; // Get number of runs in original permutation
*
*    // === Run data access ===
*    template<data_columns_t col>
*    ulint get(position pos) const; // Get value of run data column
*    template<data_columns_t col>
*    ulint get(size_t i) const; // Get value of run data column for interval i
*    ulint get_length(position pos) const; // Get length of interval containing position
*    ulint get_length(size_t i) const; // Get length of interval i
*
*    // === Search methods ===
*    template<data_columns_t col>
*    std::optional<position> pred(position pos, ulint val); // Get position of largest idx before or at position run with matching run data value
*    template<data_columns_t col>
*    std::optional<position> succ(position pos, ulint val); // Get position of smallest idx after or at position run with matching run data value
*
*    // === Serialization ===
*    size_t serialize(std::ostream& os); // Serialize data structure to ostream
*    void load(std::istream& is); // Load data structure from istream
*};
*/


// =============================== move_permutation ===============================

// Actual implementation, see documentation above
// Helper alias to use empty data columns
template<bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS>
using move_permutation = move_permutation_impl<store_absolute_positions>;

// Basic types
using move_permutation_absolute = move_permutation<true>;
using move_permutation_relative = move_permutation<false>; // Same as move_permutation<>, the default

// Alpha release aliases
template<bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS>
using moveperm = move_permutation<store_absolute_positions>;
using moveperm_absolute = move_permutation_absolute;
using moveperm_relative = move_permutation_relative;

} // namespace orbit

#endif /* end of include guard: _PUBLIC_PERMUTATION_HPP */