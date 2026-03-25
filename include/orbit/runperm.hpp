// Convenience header for run-length permutations (runperm / moveperm).
//
// This header exposes high-level wrappers around the internal move structure
// for generic "runny" permutations:
//   - runperm: run-length encoded permutations with user-defined run data.
//   - moveperm: run-length encoded permutations without extra data.
//
// Typical usage:
//   - Start from a permutation over [0, n) that is "runny", or from its
//     run decomposition:
//       std::vector<ulint> lengths;          // lengths of contiguous intervals
//       std::vector<ulint> interval_perm;    // perm position of each interval head
//   - Or, from its interval encoding given by the permutation class:
//       orbit::permutation perm;
//   - Define an enum for run data columns and a columns_tuple alias:
//       enum class data_columns { VAL1, VAL2, COUNT };
//       // Must include COUNT as the last entry to signal number of fields.
//       // alternatively, use the DEFINE_COLUMNS macro:
//       DEFINE_COLUMNS(data_columns, VAL1, VAL2)
//       using data_tuple = orbit::columns_tuple<data_columns>;
//       std::vector<data_tuple> data(lengths.size()); // Some data with tuples per run
//   - Construct a runperm instance:
//       runperm<data_columns> rp(lengths, interval_perm, domain, data);
//
// Aliases:
//   - runperm_separated / runperm_integrated select whether run data is stored
//     in a separate packed vector (Separated, default) or integrated into the
//     move structure rows (Integrated, better locality but larger rows).
//   - runperm_separated_absolute / runperm_integrated_absolute additionally store
//     absolute positions, enabling direct index lookups at extra space cost.
//   - moveperm_relative / moveperm_absolute provide analogous choices for the
//     move-only interface without run data.
//
// See the README and tests under `tests/unit/runperm/` and
// `tests/integration/runperm_test.cpp` for concrete usage.
// Expanded documentation below:
#ifndef _PUBLIC_RUNPERM_HPP
#define _PUBLIC_RUNPERM_HPP

#include "orbit/common.hpp"
#include "orbit/internal/move/permutation_impl.hpp"
#include "orbit/internal/runperm/runperm_impl.hpp"

namespace orbit {


// =============================== runperm ===============================

// Actual implementation, see documentation below
// Advanced users can use the full implementation in include/internal/runperm/runperm.hpp
template<typename data_columns_t,
         bool integrated_move_structure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS>
class runperm : public runperm_impl<data_columns_t, integrated_move_structure, store_absolute_positions, DEFAULT_EXPONENTIAL_SEARCH, move_columns, move_structure, move_vector> {
    using base = runperm_impl<data_columns_t, integrated_move_structure, store_absolute_positions, DEFAULT_EXPONENTIAL_SEARCH, move_columns, move_structure, move_vector>;
public:
    using data_columns = typename base::data_columns;
    using data_tuple = typename base::data_tuple;
    using typename base::position;
    using typename base::move_structure_perm;

    using base::base;
    using base::operator=;
};

/* === Basic types === */
template<typename data_columns_t>
using runperm_separated = runperm<data_columns_t, false, false>; // Same as runperm<data_columns_t>, the default
template<typename data_columns_t>
using runperm_integrated = runperm<data_columns_t, true, false>;
template<typename data_columns_t>
using runperm_separated_absolute = runperm<data_columns_t, false, true>;
template<typename data_columns_t>
using runperm_integrated_absolute = runperm<data_columns_t, true, true>;

/* === Simplified interface for basic users, see full RunPermImpl in include/internal/runperm/runperm.hpp for more template parameters ===
*   template<typename data_columns_t,
*         bool integrated_move_structure = DEFAULT_INTEGRATED_MOVE_STRUCTURE, // Whether to integrate the run data alongside the move structure for cache locality, default is false
*         bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS> // Whether to store absolute positions instead of interval/offset to support idx lookups, default is false
*  class runperm : public runperm_impl<data_columns_t, integrated_move_structure, store_absolute_positions> {
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
*    runperm(const std::vector<ulint>& lengths, const std::vector<ulint>& images, const ulint domain, const std::vector<data_tuple> &run_data, split_params split_params = NO_SPLITTING);
*    runperm(const std::vector<ulint>& lengths, const std::vector<ulint>& images, const ulint domain, const split_params &split_params, const std::vector<data_tuple> &run_data);
*
*    // === Splitting Constructor ===
*    // permutation -> object built by passing lengths, images to permutation constructor
*    // run_data -> run data for each interval, the size of this vector should be the same as the number of intervals after splitting by modifiying based on the permutation object
*    runperm(const permutation &perm, const std::vector<data_tuple> &run_data);
*
*    // === Advanced Constructors ===
*    // get_run_cols_data -> function to get run data for each interval, (orig_interval, orig_interval_length, new_offset_from_orig_start, new_length) -> run data
*    // structure -> pre-computed move structure stored in packed_vector
*    // ms -> pre-computed move structure stored in move_structure
*    RunPerm(const std::vector<ulint>& lengths, const std::vector<ulint>& images, const ulint domain, const SplitParams &split_params, std::function<RunData(ulint, ulint, ulint, ulint)> get_run_cols_data);
*    runperm from_structure(packed_vector<columns> &&structure, const size_t domain, const size_t runs);
*    runperm from_structure(packed_vector<base_columns> &&structure, std::vector<data_tuple> &run_data, const size_t domain, const size_t runs);
*    runperm from_move_structure(move_structure_perm &&ms);
*    runperm from_move_structure(move_structure_perm &&ms, std::vector<data_tuple> &run_data);
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



// =============================== MovePerm ===============================

// Actual implementation, see documentation below
// Advanced users can use the full implementation in include/internal/runperm/runperm.hpp
template<bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS>
class moveperm : public moveperm_impl<store_absolute_positions, DEFAULT_EXPONENTIAL_SEARCH, move_columns, move_structure, move_vector> {
    using base = moveperm_impl<store_absolute_positions, DEFAULT_EXPONENTIAL_SEARCH, move_columns, move_structure, move_vector>;
public:
    using typename base::position;

    using base::base;
    using base::operator=;
};

// Basic types
using moveperm_absolute = moveperm<true>;
using moveperm_relative = moveperm<false>; // Same as moveperm<>, the default

/* === Simplified interface for basic users, see full MovePermImpl in include/internal/runperm/runperm.hpp for more template parameters ===
* NOTE: All methods here share the documentation above for RunPerm
* template<bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS> // Whether to store absolute positions instead of interval/offset to support idx lookups, default is false
* class MovePerm {
*
* public:
*    // === Implemented types ===
*    using position = typename RunPermType::position;
*    
*    // === Constructors ===
*    MovePermImpl(std::vector<ulint>& permutation, split_params split_params = split_params()); // Constructor from permutation vector
*    MovePermImpl(const std::vector<ulint>& lengths, const std::vector<ulint>& images, split_params split_params = split_params()); // Splitting on by default
*    MovePermImpl(const permutation &perm); // Constructor from permutation object
*    
*    // === Navigation methods ===
*    position next(position pos); // Apply permutation
*    position next(position pos, ulint steps); // Return position of applying permutation multiple times
*    position up(position pos); // Return position of up one interval, circularly if at top
*    position down(position pos); // Return position of down one interval, circularly if at bottom
*    position first(); // Return first interval
*    position last(); // Return last interval
*
*    // === Query methods ===
*    ulint domain() const; // Get size of permutation
*    ulint intervals() const; // Get number of intervals in move structure (could be greater than runs due to splitting)
*    ulint runs() const; // Get number of runs in original permutation
*
*    // === Run data access ===
*    ulint get_length(position pos) const; // Get length of interval containing position
*    ulint get_length(size_t i) const; // Get length of interval i
*
*    // === Serialization ===
*    size_t serialize(std::ostream& os); // Serialize data structure to ostream
*    void load(std::istream& is); // Load data structure from istream
* };
*/

} // namespace orbit

#endif /* end of include guard: _PUBLIC_RUNPERM_HPP */