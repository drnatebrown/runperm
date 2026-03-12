#ifndef _MOVE_SPLITTING_HH
#define _MOVE_SPLITTING_HH

#include "internal/common.hpp"
#include <cmath>
#include <cassert>
#include <optional>
#include <iostream>

inline constexpr std::optional<double> DEFAULT_LENGTH_CAPPING = 8.0;
inline constexpr std::optional<ulint> DEFAULT_BALANCING = 16;

struct SplitParams {
    std::optional<double> length_capping;
    std::optional<ulint> balancing;

    SplitParams() : length_capping(DEFAULT_LENGTH_CAPPING), balancing(DEFAULT_BALANCING) {}
    SplitParams(std::optional<double> length_capping, std::optional<ulint> balancing)
    : length_capping(std::move(length_capping)), balancing(std::move(balancing)) {}
};

inline SplitParams NO_SPLITTING = SplitParams(std::nullopt, std::nullopt);
inline SplitParams DEFAULT_SPLITTING = SplitParams(DEFAULT_LENGTH_CAPPING, DEFAULT_BALANCING);
inline SplitParams ONLY_LENGTH_CAPPING = SplitParams(DEFAULT_LENGTH_CAPPING, std::nullopt);
inline SplitParams ONLY_BALANCING = SplitParams(std::nullopt, DEFAULT_BALANCING);

template<class IntVectorType>
struct SplitResult {
    IntVectorType lengths;
    IntVectorType tau_inv;
    ulint max_length;
};

template<class IntVectorType>
inline void split_by_length_capping(
    const IntVectorType& lengths, 
    const IntVectorType& tau_inv, 
    const ulint domain, 
    const double length_capping_factor, 
    SplitResult<IntVectorType>& result
) {
    assert(lengths.size() == tau_inv.size());
    assert(length_capping_factor > 0.0);

    double avg_run_length = static_cast<double>(domain) / static_cast<double>(lengths.size());
    ulint desired_max_allowed_length = static_cast<ulint>(std::ceil(avg_run_length * length_capping_factor));
    uchar bits = bit_width(desired_max_allowed_length);
    ulint max_allowed_length = MAX_VAL(bits);

    size_t new_intervals_upper_bound = std::ceil(static_cast<double>(lengths.size()) / length_capping_factor);

    split_by_max_allowed_length(lengths, tau_inv, domain, max_allowed_length, result, new_intervals_upper_bound);
}

template<class IntVectorType>
inline void split_by_max_allowed_length(
    const IntVectorType& lengths, 
    const IntVectorType& tau_inv,
    const ulint domain,
    const ulint max_allowed_length, 
    SplitResult<IntVectorType>& result,
    const size_t new_intervals_upper_bound = 0
) {
    assert(lengths.size() == tau_inv.size());
    assert(max_allowed_length > 0);
    size_t intervals_after_splitting_upper_bound = 
        (new_intervals_upper_bound == 0) ? domain - 1 : lengths.size() + new_intervals_upper_bound - 1;

    // For each original run, how many new intervals were added up to but not including this run?
    IntVectorType input_splits_exclusive_cumsum(lengths.size(), bit_width(intervals_after_splitting_upper_bound));
    
    // First pass to determine the number of intervals after splitting in input order, 
    // the number of new intervals added up to but not including each run, 
    // and the max length of the new intervals
    size_t num_intervals_after_splitting = 0;
    size_t cumulative_new_intervals = 0;
    result.max_length = 0;
    for (size_t i = 0; i < lengths.size(); ++i) {
        input_splits_exclusive_cumsum[i] = cumulative_new_intervals;

        if (lengths[i] > max_allowed_length) {
            ulint remaining = lengths[i];
            while (remaining > 0) {
                ulint chunk = std::min(remaining, max_allowed_length);
                remaining -= chunk;
                ++num_intervals_after_splitting;
                ++cumulative_new_intervals;
            }
            --cumulative_new_intervals;
            result.max_length = max_allowed_length;
        } else {
            result.max_length = std::max(result.max_length, lengths[i]);
            ++num_intervals_after_splitting;
        }
    }

    result.lengths = IntVectorType(num_intervals_after_splitting, bit_width(result.max_length));
    result.tau_inv = IntVectorType(num_intervals_after_splitting, bit_width(num_intervals_after_splitting - 1));

    // Second pass to fill the lengths and tau_inv arrays, in output order
    size_t curr_tau_inv_idx = 0;
    for (size_t i = 0; i < lengths.size(); ++i) {
        size_t j = tau_inv[i];
        size_t length = lengths[j];
        size_t num_splits = 0;
        if (length > max_allowed_length) {
            ulint remaining = length;
            while (remaining > 0) {
                ulint chunk = std::min(remaining, max_allowed_length);
                remaining -= chunk;
                result.lengths[j + input_splits_exclusive_cumsum[j] + num_splits] = chunk;
                ++num_splits;
            }
            --num_splits;
        } else {
            result.lengths[j + input_splits_exclusive_cumsum[j]] = length;
        }
        
        // Fill the tau_inv array
        size_t curr_tau_inv_val = j + input_splits_exclusive_cumsum[j];
        for (size_t k = 0; k < num_splits + 1; ++k) {
            result.tau_inv[curr_tau_inv_idx] = curr_tau_inv_val + k;
            ++curr_tau_inv_idx;
        }
    }
}

template<class IntVectorType>
inline void split_by_balancing(const IntVectorType& lengths, const IntVectorType& tau_inv, const ulint domain, const ulint balancing_factor, SplitResult<IntVectorType>& result) {
    assert(lengths.size() == tau_inv.size());
    // TODO
    result.lengths = lengths;
    result.tau_inv = tau_inv;
    // result.max_length = 0;
}

#endif