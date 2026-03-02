#ifndef _MOVE_SPLITTING_HH
#define _MOVE_SPLITTING_HH

#include "internal/common.hpp"
#include <cmath>

inline constexpr std::optional<double> DEFAULT_LENGTH_CAPPING_FACTOR = 8.0;
inline constexpr std::optional<ulint> DEFAULT_BALANCING_FACTOR = 16;

struct SplitParams {
    std::optional<double> length_capping_factor;
    std::optional<ulint> balancing_factor;

    SplitParams() : length_capping_factor(DEFAULT_LENGTH_CAPPING_FACTOR), balancing_factor(DEFAULT_BALANCING_FACTOR) {}
    SplitParams(std::optional<double> length_capping_factor, std::optional<ulint> balancing_factor)
    : length_capping_factor(std::move(length_capping_factor)), balancing_factor(std::move(balancing_factor)) {}
};

inline SplitParams NO_SPLITTING = SplitParams(std::nullopt, std::nullopt);
inline SplitParams DEFAULT_SPLITTING = SplitParams(DEFAULT_LENGTH_CAPPING_FACTOR, DEFAULT_BALANCING_FACTOR);
inline SplitParams DEFAULT_LENGTH_CAPPING = SplitParams(DEFAULT_LENGTH_CAPPING_FACTOR, std::nullopt);
inline SplitParams DEFAULT_BALANCING = SplitParams(std::nullopt, DEFAULT_BALANCING_FACTOR);

struct SplitResult {
    std::vector<ulint> lengths;
    std::vector<ulint> interval_permutations;
    ulint max_length;
};

// TODO see if double scan is faster
inline void split_by_length_capping(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const ulint domain, const double length_capping_factor, SplitResult& result) {
    assert(lengths.size() == interval_permutation.size());
    assert(length_capping_factor > 0.0);

    double avg_run_length = static_cast<double>(domain) / static_cast<double>(lengths.size());
    ulint max_allowed_length = static_cast<ulint>(std::ceil(avg_run_length * length_capping_factor));

    result.lengths.clear();
    result.interval_permutations.clear();
    result.lengths.reserve(static_cast<size_t>(lengths.size() + lengths.size()/length_capping_factor));
    result.interval_permutations.reserve(static_cast<size_t>(interval_permutation.size() + interval_permutation.size()/length_capping_factor));
    result.max_length = 0;

    for (size_t i = 0; i < lengths.size(); ++i) {
        if (lengths[i] > max_allowed_length) {
            ulint remaining = lengths[i];
            size_t sum_to_curr_chunk = 0;
            while (remaining > 0) {
                ulint chunk = std::min(remaining, max_allowed_length);
                result.lengths.push_back(chunk);
                result.interval_permutations.push_back(interval_permutation[i] + sum_to_curr_chunk);
                remaining -= chunk;
                sum_to_curr_chunk += chunk;
            }
            result.max_length = max_allowed_length;
        } else {
            result.lengths.push_back(lengths[i]);
            result.interval_permutations.push_back(interval_permutation[i]);
            result.max_length = std::max(result.max_length, lengths[i]);
        }
    }
    result.lengths.shrink_to_fit();
    result.interval_permutations.shrink_to_fit();
}

inline void split_by_balancing(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const ulint domain, const ulint balancing_factor, SplitResult& result) {
    assert(lengths.size() == interval_permutation.size());
    // TODO
    result.lengths = lengths;
    result.interval_permutations = interval_permutation;
    // result.max_length = 0;
}

#endif