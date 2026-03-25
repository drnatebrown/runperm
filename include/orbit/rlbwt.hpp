// Convenience header for RLBWT-based permutations.
//
// This header exposes high-level wrappers for working with run-length BWT
// permutations:
//   - LF / FL navigation over an RLBWT (`move_lf`, `move_fl`, `runperm_lf`, `runperm_fl`)
//   - Phi / phi_inv permutations for suffix-array based locate queries
//
// Usage overview:
//   - You start from a run-length BWT, given as:
//       std::vector<uchar> bwt_heads;       // run heads (characters)
//       std::vector<ulint> bwt_run_lengths; // run lengths
//   - For "move-only" access (no user run data), use:
//       move_lf<> lf(bwt_heads, bwt_run_lengths);
//       move_fl<> fl(bwt_heads, bwt_run_lengths);
//   - For RLBWT permutations that also carry user-defined run data, use:
//       enum class data_columns { FIELD1, FIELD2, COUNT };
//       where you must include COUNT as the last entry to signal number of fields.
//       using data_tuple = DataTuple<data_columns>;
//       std::vector<data_tuple> run_data(bwt_heads.size());
//       runperm_lf<data_columns> lf_rp(bwt_heads, bwt_run_lengths, run_data);
//       runperm_fl<data_columns> fl_rp(bwt_heads, bwt_run_lengths, run_data);
//
//   - Phi / phi_inv structures are built from an RLBWT via helpers
//     declared in the internal headers (see README for examples):
//       auto [phi_lengths, phi_interval_perm, n] =
//           rlbwt_to_phi(bwt_heads, bwt_run_lengths);
//       auto [inv_lengths, inv_interval_perm, n2] =
//           rlbwt_to_phi_inv(bwt_heads, bwt_run_lengths);
//
//     and then wrapped by:
//       MovePhi phi(phi_lengths, phi_interval_perm, n);
//       Movephi_inv phi_inv(inv_lengths, inv_interval_perm, n2);
//       RunPermPhi<data_columns> rp_phi(phi_lengths, phi_interval_perm, n, run_data);
//       RunPermphi_inv<data_columns> rp_inv(inv_lengths, inv_interval_perm, n2, run_data);
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
#include "orbit/internal/rlbwt/runperm_lf.hpp"
#include "orbit/internal/rlbwt/runperm_fl.hpp"
#include "orbit/internal/rlbwt/runperm_phi.hpp"
#include "orbit/internal/rlbwt/runperm_phi_inv.hpp"
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

// === runperm_lf ===
template<typename data_columns_t,
         bool integrated_move_structure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         typename alphabet_t = nucleotide>
class runperm_lf : public runperm_lf_impl<data_columns_t, integrated_move_structure, store_absolute_positions, DEFAULT_EXPONENTIAL_SEARCH, alphabet_t, move_vector> {
    using base = runperm_lf_impl<data_columns_t, integrated_move_structure, store_absolute_positions, DEFAULT_EXPONENTIAL_SEARCH, alphabet_t, move_vector>;
public:
    using data_columns = typename base::data_columns;
    using data_tuple = typename base::data_tuple;
    
    using base::base;
    using base::operator=;

    using position = typename base::position;
};

// Takes std::vector<uchar> bwt_and std::vector<ulint> as bwt_heads and bwt_run_lengths as input in place of lengths and interval permutations
// Otherwise, same as runperm
// Also implements lf, lf(steps), get_character(), get_character(row)

template<typename data_columns_t, typename alphabet_t = nucleotide>
using runperm_lf_separated = runperm_lf<data_columns_t, false, false, alphabet_t>; // Default
template<typename data_columns_t, typename alphabet_t = nucleotide>
using runperm_lf_integrated = runperm_lf<data_columns_t, true, false, alphabet_t>;
template<typename data_columns_t, typename alphabet_t = nucleotide>
using runperm_lf_separated_absolute = runperm_lf<data_columns_t, false, true, alphabet_t>;
template<typename data_columns_t, typename alphabet_t = nucleotide>
using runperm_lf_integrated_absolute = runperm_lf<data_columns_t, true, true, alphabet_t>;

// === move_lf ===
template<bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         typename alphabet_t = nucleotide>
class move_lf : public move_lf_impl<store_absolute_positions, DEFAULT_EXPONENTIAL_SEARCH, alphabet_t, move_vector> {
    using base = move_lf_impl<store_absolute_positions, DEFAULT_EXPONENTIAL_SEARCH, alphabet_t, move_vector>;
public:
    using base::base;
    using base::operator=;

    using position = typename base::position;
};

// See above

// === runperm_fl ===
template<typename data_columns_t,
         bool integrated_move_structure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         typename alphabet_t = nucleotide>
class runperm_fl : public runperm_fl_impl<data_columns_t, integrated_move_structure, store_absolute_positions, DEFAULT_EXPONENTIAL_SEARCH, alphabet_t, move_vector> {
    using base = runperm_fl_impl<data_columns_t, integrated_move_structure, store_absolute_positions, DEFAULT_EXPONENTIAL_SEARCH, alphabet_t, move_vector>;
public:
    using data_columns = typename base::data_columns;
    using data_tuple = typename base::data_tuple;

    using base::base;
    using base::operator=;

    using position = typename base::position;
};

// Takes std::vector<uchar> bwt_and std::vector<ulint> as bwt_heads and bwt_run_lengths as input in place of lengths and interval permutations
// Otherwise, same as RunPerm
// Also implements FL(pos), FL(pos, steps), get_character(pos), get_character(interval)

template<typename data_columns_t, typename alphabet_t = nucleotide>
using runperm_fl_separated = runperm_fl<data_columns_t, false, false, alphabet_t>; // Default
template<typename data_columns_t, typename alphabet_t = nucleotide>
using runperm_fl_integrated = runperm_fl<data_columns_t, true, false, alphabet_t>;
template<typename data_columns_t, typename alphabet_t = nucleotide>
using runperm_fl_separated_absolute = runperm_fl<data_columns_t, false, true, alphabet_t>;
template<typename data_columns_t, typename alphabet_t = nucleotide>
using runperm_fl_integrated_absolute = runperm_fl<data_columns_t, true, true, alphabet_t>;

// === move_fl ===
template<bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         typename alphabet_t = nucleotide>
class move_fl : public move_fl_impl<store_absolute_positions, DEFAULT_EXPONENTIAL_SEARCH, alphabet_t, move_vector> {
    using base = move_fl_impl<store_absolute_positions, DEFAULT_EXPONENTIAL_SEARCH, alphabet_t, move_vector>;
public:
    using base::base;
    using base::operator=;

    using position = typename base::position;
};

// See above

// === runperm_phi ===
// Need to call rlbwt_to_phi(rlbwt_heads, rlbwt_run_lengths) to get interval encoding
// Otherwise, same as runperm
// runperm_phi is just a named specilization which sets store_absolute_positions to true
// also lets set exponential_search to true or false
// Also implements phi(pos), phi(pos, steps), sa(pos)

// === num_bits_type ===
// See above

// === runperm_phi_inv ===
// Need to call rlbwt_to_phi_inv(rlbwt_heads, rlbwt_run_lengths) to get interval encoding
// Otherwise, same as runperm -->
// runperm_phi_inv is just a named specilization which sets store_absolute_positions to true
// also lets set exponential_search to true or false
// Also implements phi_inv(pos), phi_inv(pos, steps), sa(pos)

// === move_phi_inv<> ===
// See above

} // end namespace orbit::rlbwt

#endif /* end of include guard: _PUBLIC_RUNPERM_RLBWT_HPP */