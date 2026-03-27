// Convenience header for RLBWT-based permutations.
//
// This header exposes high-level wrappers for working with run-length BWT
// permutations:
//   - LF / FL navigation over an RLBWT (`lf_permutation`, `fl_permutation`, `lf_move`, `fl_move`)
//   - phi / phi_inv permutations for suffix-array based locate queries
//   (phi_permutation, phi_inv_permutation, phi_move, phi_inv_move)
//
// Usage overview (RLBWT):
//   - You start from a run-length BWT, given as:
//       std::vector<uchar> bwt_heads;       // run heads (characters)
//       std::vector<ulint> bwt_run_lengths; // run lengths
//   - For "move-only" access (no user run data), use:
//       lf_move<> lf(bwt_heads, bwt_run_lengths);
//       fl_move<> fl(bwt_heads, bwt_run_lengths);
//   - For RLBWT permutations that also carry user-defined run data, use:
//       enum class data_columns { FIELD1, FIELD2, COUNT };
//       where you must include COUNT as the last entry to signal number of fields.
//       using data_tuple = DataTuple<data_columns>;
//       std::vector<data_tuple> run_data(bwt_heads.size());
//       lf_permutation<data_columns> lf_rp(bwt_heads, bwt_run_lengths, run_data);
//       fl_permutation<data_columns> fl_rp(bwt_heads, bwt_run_lengths, run_data);
//
//   - phi / phi_inv structures are built from an RLBWT via helpers
//     declared in the internal headers (see README for examples):
//       auto encoding = rlbwt_to_phi(bwt_heads, bwt_run_lengths);
//       auto encoding_inv = rlbwt_to_phi_inv(bwt_heads, bwt_run_lengths);
//
//     and then wrapped by:
//       phi_move phi(encoding);
//       phi_inv_move phi_inv(encoding_inv);
//       phi_permutation<data_columns> phi(encoding, run_data);
//       phi_inv_permutation<data_columns> phi_inv(encoding_inv, run_data);
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

#include "orbit/interval_encoding.hpp"
#include "orbit/common.hpp"
#include "orbit/internal/rlbwt/lf_permutation.hpp"
#include "orbit/internal/rlbwt/fl_permutation.hpp"
#include "orbit/internal/rlbwt/phi_permutation.hpp"
#include "orbit/internal/rlbwt/phi_inv_permutation.hpp"
#include "orbit/internal/rlbwt/rlbwt_helpers.hpp"
#include "orbit/internal/rlbwt/phi_helpers.hpp"

namespace orbit::rlbwt {

template<typename alphabet_t = nucleotide>
using rlbwt_interval_encoding = rlbwt_interval_encoding_impl<int_vector_aligned, alphabet_t>;

interval_encoding_impl<> rlbwt_to_phi(const std::vector<uchar>& bwt_heads, const std::vector<ulint>& bwt_run_lengths, const split_params& sp = split_params()) {
    size_t domain;
    ulint max_length;
    auto [phi_lengths, phi_img_rank_inv] = rlbwt_to_phi_img_rank_inv(bwt_heads, bwt_run_lengths, &domain, &max_length);
    return interval_encoding_impl<>::from_lengths_and_img_rank_inv(phi_lengths, phi_img_rank_inv, domain, max_length, sp);
}

interval_encoding_impl<> rlbwt_to_phi_inv(const std::vector<uchar>& bwt_heads, const std::vector<ulint>& bwt_run_lengths, const split_params& sp = split_params()) {
    size_t domain;
    ulint max_length;
    auto [phi_inv_lengths, phi_inv_img_rank_inv] = rlbwt_to_phi_inv_img_rank_inv(bwt_heads, bwt_run_lengths, &domain, &max_length);
    return interval_encoding_impl<>::from_lengths_and_img_rank_inv(phi_inv_lengths, phi_inv_img_rank_inv, domain, max_length, sp);
}

// =============================== RLBWT Based Permutations for LF/FL ===============================
// === lf_permutation ===
template<typename data_columns_t = empty_data_columns,
         bool integrated_move_structure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         typename alphabet_t = nucleotide>
class lf_permutation : public lf_permutation_impl<data_columns_t, integrated_move_structure, store_absolute_positions, DEFAULT_EXPONENTIAL_SEARCH, alphabet_t, move_vector> {
    using base = lf_permutation_impl<data_columns_t, integrated_move_structure, store_absolute_positions, DEFAULT_EXPONENTIAL_SEARCH, alphabet_t, move_vector>;
public:
    using typename base::data_columns;
    using typename base::data_tuple;
    using typename base::position;
    
    using base::base;
    using base::operator=;

};

// Takes std::vector<uchar> bwt_and std::vector<ulint> as bwt_heads and bwt_run_lengths as input in place of lengths and images
// Otherwise, same as permutation. Can also pass rlbwt_interval_encoding as input instead of bwt_heads and bwt_run_lengths.
// Also implements LF, LF(steps), get_character(), get_character(row)

template<typename data_columns_t = empty_data_columns, typename alphabet_t = nucleotide>
using lf_permutation_absolute = lf_permutation<data_columns_t, false, true, alphabet_t>;
template<typename data_columns_t = empty_data_columns, typename alphabet_t = nucleotide>
using lf_permutation_separated = lf_permutation<data_columns_t, false, false, alphabet_t>; // Default
template<typename data_columns_t = empty_data_columns, typename alphabet_t = nucleotide>
using lf_permutation_integrated = lf_permutation<data_columns_t, true, false, alphabet_t>;
template<typename data_columns_t = empty_data_columns, typename alphabet_t = nucleotide>
using lf_permutation_separated_absolute = lf_permutation<data_columns_t, false, true, alphabet_t>;
template<typename data_columns_t = empty_data_columns, typename alphabet_t = nucleotide>
using lf_permutation_integrated_absolute = lf_permutation<data_columns_t, true, true, alphabet_t>;

// Alpha release aliases
template<typename data_columns_t = empty_data_columns,
         bool integrated_move_structure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         typename alphabet_t = nucleotide>
using runperm_lf = lf_permutation<data_columns_t, integrated_move_structure, store_absolute_positions, alphabet_t>;
template<typename data_columns_t = empty_data_columns, typename alphabet_t = nucleotide>
using runperm_lf_absolute = lf_permutation<data_columns_t, false, true, alphabet_t>;
template<typename data_columns_t = empty_data_columns, typename alphabet_t = nucleotide>
using runperm_lf_separated = lf_permutation<data_columns_t, false, false, alphabet_t>;
template<typename data_columns_t = empty_data_columns, typename alphabet_t = nucleotide>
using runperm_lf_integrated = lf_permutation<data_columns_t, true, false, alphabet_t>;
template<typename data_columns_t = empty_data_columns, typename alphabet_t = nucleotide>
using runperm_lf_separated_absolute = lf_permutation<data_columns_t, false, true, alphabet_t>;
template<typename data_columns_t = empty_data_columns, typename alphabet_t = nucleotide>
using runperm_lf_integrated_absolute = lf_permutation<data_columns_t, true, true, alphabet_t>;

// === lf_move ===
template<bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         typename alphabet_t = nucleotide>
using lf_move = lf_move_impl<store_absolute_positions, DEFAULT_EXPONENTIAL_SEARCH, alphabet_t, move_vector>;
using lf_move_absolute = lf_move<true>;
using lf_move_relative = lf_move<false>;

// Alpha release aliases
template<bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         typename alphabet_t = nucleotide>
using move_lf = lf_move<store_absolute_positions, alphabet_t>;
using move_lf_absolute = lf_move_absolute;
using move_lf_relative = lf_move_relative;

// === fl_permutation ===
template<typename data_columns_t = empty_data_columns,
         bool integrated_move_structure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         typename alphabet_t = nucleotide>
class fl_permutation : public fl_permutation_impl<data_columns_t, integrated_move_structure, store_absolute_positions, DEFAULT_EXPONENTIAL_SEARCH, alphabet_t, move_vector> {
    using base = fl_permutation_impl<data_columns_t, integrated_move_structure, store_absolute_positions, DEFAULT_EXPONENTIAL_SEARCH, alphabet_t, move_vector>;
public:
    using typename base::data_columns;
    using typename base::data_tuple;
    using typename base::position;

    using base::base;
    using base::operator=;
};

// Takes std::vector<uchar> bwt_and std::vector<ulint> as bwt_heads and bwt_run_lengths as input in place of lengths and interval permutations
// Otherwise, same as permutation. Can also pass rlbwt_interval_encoding as input instead of bwt_heads and bwt_run_lengths.
// Also implements FL(pos), FL(pos, steps), get_character(pos), get_character(interval)

template<typename data_columns_t = empty_data_columns, typename alphabet_t = nucleotide>
using fl_permutation_absolute = fl_permutation<data_columns_t, false, true, alphabet_t>;
template<typename data_columns_t = empty_data_columns, typename alphabet_t = nucleotide>
using fl_permutation_separated = fl_permutation<data_columns_t, false, false, alphabet_t>; // Default
template<typename data_columns_t, typename alphabet_t = nucleotide>
using fl_permutation_integrated = fl_permutation<data_columns_t, true, false, alphabet_t>;
template<typename data_columns_t, typename alphabet_t = nucleotide>
using fl_permutation_separated_absolute = fl_permutation<data_columns_t, false, true, alphabet_t>;
template<typename data_columns_t, typename alphabet_t = nucleotide>
using fl_permutation_integrated_absolute = fl_permutation<data_columns_t, true, true, alphabet_t>;

// Alpha release aliases
template<typename data_columns_t = empty_data_columns,
         bool integrated_move_structure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         typename alphabet_t = nucleotide>
using runperm_fl = fl_permutation<data_columns_t, integrated_move_structure, store_absolute_positions, alphabet_t>;
template<typename data_columns_t = empty_data_columns, typename alphabet_t = nucleotide>
using runperm_fl_absolute = fl_permutation<data_columns_t, false, true, alphabet_t>;
template<typename data_columns_t = empty_data_columns, typename alphabet_t = nucleotide>
using runperm_fl_separated = fl_permutation<data_columns_t, false, false, alphabet_t>;
template<typename data_columns_t, typename alphabet_t = nucleotide>
using runperm_fl_integrated = fl_permutation<data_columns_t, true, false, alphabet_t>;
template<typename data_columns_t, typename alphabet_t = nucleotide>
using runperm_fl_separated_absolute = fl_permutation<data_columns_t, false, true, alphabet_t>;
template<typename data_columns_t, typename alphabet_t = nucleotide>
using runperm_fl_integrated_absolute = fl_permutation<data_columns_t, true, true, alphabet_t>;

// === move_fl ===
template<bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         typename alphabet_t = nucleotide>
using fl_move = fl_move_impl<store_absolute_positions, DEFAULT_EXPONENTIAL_SEARCH, alphabet_t, move_vector>;
using fl_move_absolute = fl_move<true>;
using fl_move_relative = fl_move<false>;

// Alpha release aliases
template<bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         typename alphabet_t = nucleotide>
using move_fl = fl_move<store_absolute_positions, alphabet_t>;
using move_fl_absolute = fl_move_absolute;
using move_fl_relative = fl_move_relative;

// =============================== RLBWT Based Permutations for PHI/PHI_INV ===============================
// === phi_permutation ===
// Need to call rlbwt_to_phi(rlbwt_heads, rlbwt_run_lengths) to get interval encoding
// Otherwise, same as permutation.
// phi_permutation is just a named specialization which sets store_absolute_positions to true
// and lets set exponential_search to true or false based on the template parameter.
// Also implements phi(pos), phi(pos, steps), SA(pos)

template<typename data_columns_t = empty_data_columns,
         bool integrated_move_structure = DEFAULT_INTEGRATED_MOVE_STRUCTURE>
using phi_permutation = phi_permutation_impl<data_columns_t, integrated_move_structure, DEFAULT_EXPONENTIAL_SEARCH, move_vector>;
template<typename data_columns_t = empty_data_columns>
using phi_permutation_integrated = phi_permutation<data_columns_t, true>;

// Alpha release aliases
template<typename data_columns_t = empty_data_columns,
         bool integrated_move_structure = DEFAULT_INTEGRATED_MOVE_STRUCTURE>
using runperm_phi = phi_permutation<data_columns_t, integrated_move_structure>;
template<typename data_columns_t = empty_data_columns>
using runperm_phi_integrated = phi_permutation<data_columns_t, true>;

// === phi_move ===
// See above, but no user data.
using phi_move = phi_move_impl<DEFAULT_EXPONENTIAL_SEARCH, move_vector>;

// Alpha release aliases
using move_phi = phi_move;

// === phi_inv_permutation ===
// Need to call rlbwt_to_phi_inv(rlbwt_heads, rlbwt_run_lengths) to get interval encoding
// Otherwise, same as permutation.
// phi_inv_permutation is just a named specialization which sets store_absolute_positions to true
// and lets set exponential_search to true or false based on the template parameter.
// Also implements phi_inv(pos), phi_inv(pos, steps), SA(pos)

template<typename data_columns_t = empty_data_columns,
         bool integrated_move_structure = DEFAULT_INTEGRATED_MOVE_STRUCTURE>
using phi_inv_permutation = phi_inv_permutation_impl<data_columns_t, integrated_move_structure, DEFAULT_EXPONENTIAL_SEARCH, move_vector>;
template<typename data_columns_t = empty_data_columns>
using phi_inv_permutation_integrated = phi_inv_permutation<data_columns_t, true>;

// Alpha release aliases
template<typename data_columns_t = empty_data_columns,
         bool integrated_move_structure = DEFAULT_INTEGRATED_MOVE_STRUCTURE>
using runperm_phi_inv = phi_inv_permutation<data_columns_t, integrated_move_structure>;
template<typename data_columns_t = empty_data_columns>
using runperm_phi_inv_integrated = phi_inv_permutation<data_columns_t, true>;

// === phi_inv_move ===
// See above, but no user data.
using phi_inv_move = phi_inv_move_impl<DEFAULT_EXPONENTIAL_SEARCH, move_vector>;

// Alpha release aliases
using move_phi_inv = phi_inv_move;

} // end namespace orbit::rlbwt

#endif /* end of include guard: _RLBWT_HPP */