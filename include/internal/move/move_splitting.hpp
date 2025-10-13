#ifndef _MOVE_SPLITTING_HH
#define _MOVE_SPLITTING_HH

#include "internal/common.hpp"

struct SplitParams {
    std::optional<ulint> max_allowed_length;
    std::optional<ulint> balancing_factor;

    SplitParams() : max_allowed_length(std::nullopt), balancing_factor(std::nullopt) {}
};

struct SplitResult {
    std::vector<ulint> lengths;
    std::vector<ulint> interval_permutations;
    ulint max_length;
};

// TODO see if double scan is faster
void split_by_max_allowed_length(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const ulint max_allowed_length, SplitResult& result) {
    assert(lengths.size() == interval_permutation.size());

    result.lengths.clear();
    result.interval_permutations.clear();
    result.lengths.reserve(static_cast<size_t>(1.5*lengths.size()));
    result.interval_permutations.reserve(static_cast<size_t>(1.5*interval_permutation.size()));
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

void split_by_balancing_factor(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const ulint balancing_factor, SplitResult& result) {
    assert(lengths.size() == interval_permutation.size());
    // TODO
    result.lengths = lengths;
    result.interval_permutations = interval_permutation;
    result.max_length = 0;
}

#endif