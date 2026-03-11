#ifndef _MOVE_SPLITTING_HH
#define _MOVE_SPLITTING_HH

#include "internal/common.hpp"
#include <cmath>
#include <cassert>
#include <optional>

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

    split_by_max_allowed_length(lengths, tau_inv, domain, max_allowed_length, result);
}

template<class IntVectorType>
inline void split_by_max_allowed_length(
    const IntVectorType& lengths, 
    const IntVectorType& tau_inv,
    const ulint domain,
    const ulint max_allowed_length, 
    SplitResult<IntVectorType>& result
) {
    assert(lengths.size() == tau_inv.size());
    assert(max_allowed_length > 0);

    // First pass to determine the number of intervals after splitting and the max length
    size_t num_intervals_after_splitting = 0;
    // Map from original interval index to new interval index
    IntVectorType old_to_new_interval_idx(lengths.size(), bit_width(domain - 1));
    result.max_length = 0;
    for (size_t i = 0; i < lengths.size(); ++i) {
        old_to_new_interval_idx.set(i, num_intervals_after_splitting);
        if (lengths[i] > max_allowed_length) {
            ulint remaining = lengths[i];
            while (remaining > 0) {
                ulint chunk = std::min(remaining, max_allowed_length);
                remaining -= chunk;
                ++num_intervals_after_splitting;
            }
            result.max_length = max_allowed_length;
        } else {
            result.max_length = std::max(result.max_length, lengths[i]);
            ++num_intervals_after_splitting;
        }
    }

    result.lengths = IntVectorType(num_intervals_after_splitting, bit_width(result.max_length));
    result.tau_inv = IntVectorType(num_intervals_after_splitting, bit_width(num_intervals_after_splitting - 1));

    // Second pass to fill the lengths and tau_inv arrays
    size_t curr_interval_idx = 0;
    for (size_t i = 0; i < lengths.size(); ++i) {
        if (lengths[i] > max_allowed_length) {
            ulint remaining = lengths[i];
            size_t relative_interval_idx = 0;
            while (remaining > 0) {
                ulint chunk = std::min(remaining, max_allowed_length);
                result.lengths[curr_interval_idx] = chunk;
                result.tau_inv[curr_interval_idx] = old_to_new_interval_idx[tau_inv[i]] + relative_interval_idx;
                ++curr_interval_idx;
                ++relative_interval_idx;
                remaining -= chunk;
            }
        } else {
            result.lengths[curr_interval_idx] = lengths[i];
            result.tau_inv[curr_interval_idx] = old_to_new_interval_idx[tau_inv[i]];
            ++curr_interval_idx;
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